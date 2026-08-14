#!/usr/bin/env python3
"""enzyme_characterize.py -- per-reaction heterologous-enzyme variant characterisation.

For ONE reaction (its KO/enzyme) this builds the exact deep-dive the PanRoute webapp shows
when a user clicks a heterologous-expression candidate step:

  1. gather the enzyme's protein variants  : every KEGG gene annotated to the reaction's KO
     (organism attribution comes for free from the KEGG gene id).
  2. length filter                          : drop variants < 50% of the longest variant length.
  3. cluster at 80% identity (mmseqs)       : collapse near-duplicates to representatives; keep the
                                              cluster's member organisms (= who you can clone it from).
  4. kcat + Km (MPEK / MTLKcatKM)           : per representative, from (AA sequence, substrate SMILES,
                                              organism, temperature). MPEK is temperature-aware.
  5. thermostability curve (TemStaPro)      : per representative, P(Tm above 40/45/50/55/60/65 C) +
                                              a melting-range label -> a real thermostability-vs-temp curve.
  6. rank best -> worst                     : by catalytic efficiency (kcat/Km) and thermostability.

Writes docs/data/enzymes/<rid>.json consumed by the client (3D scatter of Km x kcat x thermostability
+ sortable ranked table + per-variant sequence and source organisms).

NOTHING here is invented: sequences are the real KEGG orthologues, kcat/Km are MPEK predictions,
thermostability is TemStaPro. Predictions are model estimates -- labelled as such in the output.

USAGE
  python3 bin/enzyme_characterize.py --rid R06943 --ko K06446 --sub-cid C05662 \
      [--temps 25,37,55] [--max-genes 600] [--max-reps 120] [--stage all|fetch|predict]
"""
import os, sys, json, time, csv, subprocess, urllib.request, urllib.error, argparse, hashlib

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs")
DATA = os.path.join(DOCS, "data")
CACHE = os.path.join(ROOT, "cache", "aaseq")
WORK = os.path.join(ROOT, "cache", "enzyme_work")
OUTDIR = os.path.join(DATA, "enzymes")
MPEK_SH = "/data/metabolic_atlas/HETERO_CANDIDATES/pipeline/predict_kinetics_mpek.sh"
TEMSTAPRO = os.path.join(ROOT, "tools", "TemStaPro", "temstapro")
PROTT5 = "/home/omidard/panGEM_pipeline/mpek/MTLKcatKM/checkpoints/prot_t5_xl_uniref50"
MMSEQS = "/home/omidard/anaconda3/bin/mmseqs"
KEGG = "https://rest.kegg.jp"
MIN_INTERVAL = 1.0 / 3.0
_last = [0.0]
for d in (CACHE, WORK, OUTDIR):
    os.makedirs(d, exist_ok=True)


def _throttle():
    dt = time.time() - _last[0]
    if dt < MIN_INTERVAL:
        time.sleep(MIN_INTERVAL - dt)
    _last[0] = time.time()


def kegg_get(op):
    """GET rest.kegg.jp/<op> with on-disk cache (MISS cached as empty sentinel)."""
    key = hashlib.sha1(op.encode()).hexdigest()[:16]
    cp = os.path.join(CACHE, key)
    if os.path.exists(cp):
        return open(cp).read()
    _throttle()
    try:
        with urllib.request.urlopen(f"{KEGG}/{op}", timeout=40) as r:
            txt = r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        txt = "" if e.code == 404 else None
    except Exception:
        txt = None
    if txt is not None:
        open(cp, "w").write(txt)
    return txt or ""


def ko_genes(ko):
    """org:gene ids annotated to a KO."""
    txt = kegg_get(f"link/genes/ko:{ko}")
    genes = []
    for line in txt.splitlines():
        p = line.split("\t")
        if len(p) == 2:
            genes.append(p[1].strip())
    return genes


def fetch_aaseq(genes):
    """batched /get/<g1>+<g2>.../aaseq -> {gene: sequence}. up to 10 per call."""
    out = {}
    for i in range(0, len(genes), 10):
        batch = genes[i:i + 10]
        rec = kegg_get("get/" + "+".join(batch) + "/aaseq")
        cur, seq = None, []
        for line in rec.splitlines():
            if line.startswith(">"):
                if cur:
                    out[cur] = "".join(seq)
                # header: >org:gene  description
                cur = line[1:].split()[0]
                seq = []
            elif cur:
                seq.append(line.strip())
        if cur:
            out[cur] = "".join(seq)
        if (i // 10) % 20 == 0:
            print(f"  aaseq {min(i+10,len(genes))}/{len(genes)}", flush=True)
    # clean: uppercase, strip non-standard chars
    clean = {}
    for g, s in out.items():
        s = "".join(c for c in s.upper() if c.isalpha())
        if len(s) >= 40:
            clean[g] = s
    return clean


def load_json(p, default):
    try:
        return json.load(open(p))
    except Exception:
        return default


def species_of(gene, tax):
    org = gene.split(":")[0]
    t = tax.get(org)
    return (t[0] if t else org), org


def run_mmseqs(fasta, prefix):
    tmp = os.path.join(WORK, "mm_tmp")
    os.makedirs(tmp, exist_ok=True)
    # 80% identity, coverage 0.5 (bidirectional) -> representative per near-identical cluster
    subprocess.run([MMSEQS, "easy-cluster", fasta, prefix, tmp,
                    "--min-seq-id", "0.8", "-c", "0.5", "-v", "1"], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    clusters = {}   # rep -> [members]
    with open(prefix + "_cluster.tsv") as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")
            clusters.setdefault(rep, []).append(mem)
    return clusters


def run_temstapro(fasta, out_tsv):
    env = dict(os.environ, USE_TF="0", USE_FLAX="0")
    cmd = [sys.executable, TEMSTAPRO, "--input-fasta", fasta,
           "--PT5-model", os.path.dirname(PROTT5),
           "--mean-output", out_tsv, "--temperature-ranges"]
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(out_tsv):
        print("  [temstapro] FAILED:", r.stderr[-500:], flush=True)
        return {}
    therm = {}
    with open(out_tsv) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            cid = row.get("cand_id") or row.get(rd.fieldnames[0])
            curve = {}
            for t in (40, 45, 50, 55, 60, 65):
                v = row.get(f"p{t}")
                if v not in (None, ""):
                    curve[t] = float(v)
            therm[cid] = {"curve": curve, "label": row.get("temstapro_label", "")}
    return therm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rid", required=True)
    ap.add_argument("--ko", required=True, help="comma-separated KO id(s)")
    ap.add_argument("--sub-cid", default="", help="substrate KEGG compound id (for SMILES)")
    ap.add_argument("--sub-smiles", default="")
    ap.add_argument("--sub-name", default="")
    ap.add_argument("--temps", default="37")
    ap.add_argument("--max-genes", type=int, default=600)
    ap.add_argument("--max-reps", type=int, default=120)
    ap.add_argument("--stage", default="all", choices=["all", "fetch", "predict"])
    args = ap.parse_args()

    tax = load_json(os.path.join(DATA, "taxonomy.json"), {})
    smiles_db = load_json(os.path.join(DATA, "smiles.json"), {})
    smiles = args.sub_smiles or smiles_db.get(args.sub_cid, "")
    temps = [int(x) for x in args.temps.split(",") if x.strip()]
    kos = [k.strip() for k in args.ko.split(",") if k.strip()]
    wpre = os.path.join(WORK, args.rid)

    # ---- 1. gather variants (KEGG KO orthologues) ----
    print(f"[{args.rid}] KO {kos} — gathering variants", flush=True)
    genes = []
    for ko in kos:
        genes += ko_genes(ko)
    genes = sorted(set(genes))
    n_raw_genes = len(genes)
    print(f"  {n_raw_genes} KO genes", flush=True)
    if n_raw_genes > args.max_genes:
        # deterministic subsample by hash so runs are reproducible (cluster still dedups)
        genes = sorted(genes, key=lambda g: hashlib.sha1(g.encode()).hexdigest())[:args.max_genes]
        print(f"  capped to {len(genes)} for sequence fetch", flush=True)
    seqs = fetch_aaseq(genes)
    print(f"  {len(seqs)} sequences fetched", flush=True)
    if not seqs:
        print("  no sequences — abort"); return

    # ---- 2. length filter (drop < 50% of the longest) ----
    longest = max(len(s) for s in seqs.values())
    kept = {g: s for g, s in seqs.items() if len(s) >= 0.5 * longest}
    print(f"  longest={longest} aa; {len(kept)}/{len(seqs)} pass the 50%-length filter", flush=True)

    faa = wpre + "_all.faa"
    with open(faa, "w") as fh:
        for g, s in kept.items():
            fh.write(f">{g}\n{s}\n")

    # ---- 3. cluster at 80% identity ----
    clusters = run_mmseqs(faa, wpre)
    print(f"  {len(clusters)} clusters at 80% identity", flush=True)
    # representative = longest member; carry member organisms
    reps = []
    for rep, members in clusters.items():
        best = max(members, key=lambda m: len(kept.get(m, "")))
        orgs = sorted({species_of(m, tax)[0] for m in members})
        reps.append({"rep": best, "seq": kept[best], "members": members,
                     "cluster_size": len(members), "organisms": orgs})
    reps.sort(key=lambda r: -r["cluster_size"])
    if len(reps) > args.max_reps:
        reps = reps[:args.max_reps]
        print(f"  keeping top {len(reps)} clusters by size", flush=True)

    # stash intermediate so 'predict' stage can resume
    interm = wpre + "_reps.json"
    json.dump({"reps": reps, "n_raw_genes": n_raw_genes, "n_seqs": len(seqs),
               "n_after_length": len(kept), "longest": longest}, open(interm, "w"))
    if args.stage == "fetch":
        print("  fetch stage done ->", interm); return

    # ---- 4. kcat + Km (MPEK) ----
    print(f"  MPEK kcat/Km on {len(reps)} variants x temps {temps}", flush=True)
    mp_in = wpre + "_mpek_in.csv"
    with open(mp_in, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["id", "sequence", "smiles", "species", "temp"])
        for i, r in enumerate(reps):
            sp = species_of(r["rep"], tax)[0]
            for t in temps:
                w.writerow([f"v{i}@{t}", r["seq"], smiles, sp, t])
    mp_out = wpre + "_mpek_out.csv"
    kin = {}
    if smiles:
        subprocess.run([MPEK_SH, mp_in, mp_out], check=True)
        with open(mp_out) as fh:
            for row in csv.DictReader(fh):
                kin[row["id"]] = row
    else:
        print("  no substrate SMILES — skipping kinetics", flush=True)

    # ---- 5. thermostability (TemStaPro) ----
    print(f"  TemStaPro thermostability on {len(reps)} variants", flush=True)
    rep_faa = wpre + "_reps.faa"
    with open(rep_faa, "w") as fh:
        for i, r in enumerate(reps):
            fh.write(f">v{i}\n{r['seq']}\n")
    therm = run_temstapro(rep_faa, wpre + "_therm.tsv")

    # ---- 6. assemble + rank ----
    variants = []
    for i, r in enumerate(reps):
        row37 = kin.get(f"v{i}@{temps[0]}", {})
        kcat = _f(row37.get("kcat"))
        km = _f(row37.get("km"))
        eff = (kcat / km) if (kcat and km) else None
        tvals = {t: _f(kin.get(f"v{i}@{t}", {}).get("kcat")) for t in temps}
        th = therm.get(f"v{i}", {})
        variants.append({
            "id": f"v{i}", "rep_gene": r["rep"], "organisms": r["organisms"],
            "cluster_size": r["cluster_size"], "length": len(r["seq"]),
            "kcat": kcat, "km": km, "kcat_km": eff,
            "kcat_by_temp": tvals,
            "thermo": th.get("curve", {}), "thermo_label": th.get("label", ""),
            "sequence": r["seq"],
        })
    # rank best->worst: catalytic efficiency (log kcat/Km) + thermostability (P>55C)
    import math

    def score(v):
        s = 0.0
        if v["kcat_km"]:
            s += math.log10(v["kcat_km"])
        p55 = (v["thermo"] or {}).get("55")
        if p55 is not None:
            s += 1.5 * p55
        return s
    variants.sort(key=score, reverse=True)
    for rank, v in enumerate(variants):
        v["rank"] = rank + 1
        v["score"] = round(score(v), 3)

    out = {
        "rid": args.rid, "ko": kos,
        "substrate": {"cid": args.sub_cid, "name": args.sub_name, "smiles": smiles},
        "temps": temps,
        "n_ko_genes": n_raw_genes, "n_sequences": len(seqs),
        "n_after_length_filter": len(kept), "n_clusters": len(clusters),
        "n_variants": len(variants),
        "method": {"cluster": "mmseqs easy-cluster --min-seq-id 0.8 -c 0.5",
                   "length_filter": ">=50% of longest",
                   "kinetics": "MPEK / MTLKcatKM (kcat 1/s, Km mM)",
                   "thermostability": "TemStaPro ProtT5-XL (P Tm > threshold)"},
        "variants": variants,
    }
    outp = os.path.join(OUTDIR, f"{args.rid}.json")
    json.dump(out, open(outp, "w"))
    print(f"  wrote {outp}: {len(variants)} ranked variants", flush=True)


def _f(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


if __name__ == "__main__":
    main()

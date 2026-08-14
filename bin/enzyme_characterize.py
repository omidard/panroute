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
TSP_DIR = os.path.join(ROOT, "tools", "TemStaPro")
TEMSTAPRO = os.path.join(TSP_DIR, "temstapro")
# ProtT5-XL dir with a local pytorch_model.bin + tokenizer (T5EncoderModel loads the encoder from it,
# so TemStaPro doesn't re-download the half-enc model; same encoder embeddings).
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


_ORGNAMES = {}
def _load_orgnames():
    global _ORGNAMES
    if not _ORGNAMES:
        p = os.path.join(DATA, "kegg_org_names.json")
        if os.path.exists(p):
            _ORGNAMES = json.load(open(p))
    return _ORGNAMES


def species_of(gene, tax):
    org = gene.split(":")[0]
    t = tax.get(org)
    if t:
        return t[0], org
    # orthologues can live in genomes outside our tracked set -> resolve the KEGG code to a species
    # name from the bundled kegg_org_names map (so the UI never shows a bare code like "sgra").
    return _load_orgnames().get(org, org), org


def run_mmseqs(fasta, prefix):
    tmp = prefix + "_mmtmp"                       # per-prefix tmp so batch (repeated) calls don't collide
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
    cmd = [sys.executable, TEMSTAPRO, "-f", fasta, "-d", PROTT5, "-t", TSP_DIR,
           "--more-thresholds", "--mean-output", out_tsv]
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(out_tsv):
        print("  [temstapro] FAILED:", r.stderr[-500:], flush=True)
        return {}
    therm = {}
    with open(out_tsv) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            cid = row.get("protein_id") or row.get(rd.fieldnames[0])
            curve = {}
            for t in (40, 45, 50, 55, 60, 65, 70, 75, 80):   # P(Tm above threshold), raw probability
                v = row.get(f"t{t}_raw")
                if v not in (None, ""):
                    curve[t] = float(v)
            therm[cid] = {"curve": curve, "label": row.get("left_hand_label") or row.get("right_hand_label", ""),
                          "thermophilicity": row.get("thermophilicity", "")}
    return therm


def _f(x):
    try:
        v = float(x); return v if v == v else None
    except (TypeError, ValueError):
        return None


METHOD = {"cluster": "mmseqs easy-cluster --min-seq-id 0.8 -c 0.5",
          "length_filter": ">=50% of longest",
          "kinetics": "MPEK / MTLKcatKM (kcat 1/s, Km mM)",
          "thermostability": "TemStaPro ProtT5-XL (P Tm > threshold)"}


def gather_reps(kos, max_genes, max_reps, tax, tag):
    """KEGG orthologues -> length-filter (>=50% longest) -> 80%-identity clusters -> representatives."""
    genes = []
    for ko in kos:
        genes += ko_genes(ko)
    genes = sorted(set(genes))
    n_raw = len(genes)
    print(f"  [{tag}] {n_raw} KO genes", flush=True)
    if n_raw > max_genes:                              # deterministic subsample; clustering still dedups
        genes = sorted(genes, key=lambda g: hashlib.sha1(g.encode()).hexdigest())[:max_genes]
    seqs = fetch_aaseq(genes)
    print(f"  [{tag}] {len(seqs)} sequences fetched", flush=True)
    if not seqs:
        return [], {"n_ko_genes": n_raw, "n_sequences": 0, "n_after_length": 0, "n_clusters": 0}
    longest = max(len(s) for s in seqs.values())
    kept = {g: s for g, s in seqs.items() if len(s) >= 0.5 * longest}
    faa = os.path.join(WORK, tag + "_all.faa")
    with open(faa, "w") as fh:
        for g, s in kept.items():
            fh.write(f">{g}\n{s}\n")
    clusters = run_mmseqs(faa, os.path.join(WORK, tag))
    reps = []
    for rep, members in clusters.items():
        best = max(members, key=lambda m: len(kept.get(m, "")))
        orgs = sorted({species_of(m, tax)[0] for m in members})
        reps.append({"rep": best, "seq": kept[best], "cluster_size": len(members), "organisms": orgs})
    reps.sort(key=lambda r: -r["cluster_size"])
    reps = reps[:max_reps]
    print(f"  [{tag}] {len(clusters)} clusters at 80% id -> {len(reps)} representatives", flush=True)
    return reps, {"n_ko_genes": n_raw, "n_sequences": len(seqs),
                  "n_after_length": len(kept), "n_clusters": len(clusters)}


def run_mpek(rows, workprefix):
    """rows=[[id,seq,smiles,species,temp],...] -> {id: mpek_row}. One MPEK invocation (model loads once)."""
    if not rows:
        return {}
    mp_in, mp_out = workprefix + "_mpek_in.csv", workprefix + "_mpek_out.csv"
    with open(mp_in, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["id", "sequence", "smiles", "species", "temp"]); w.writerows(rows)
    subprocess.run([MPEK_SH, mp_in, mp_out], check=True)
    kin = {}
    with open(mp_out) as fh:
        for row in csv.DictReader(fh):
            kin[row["id"]] = row
    return kin


def assemble_variants(reps, kin, therm, temps, idp):
    """reps + MPEK rows (id `<idp>_v<i>_<t>`) + TemStaPro (id `<idp>_v<i>`) -> ranked variant list."""
    import math
    variants = []
    for i, r in enumerate(reps):
        base = f"{idp}_v{i}"
        row0 = kin.get(f"{base}_{temps[0]}", {})
        kcat, km = _f(row0.get("kcat")), _f(row0.get("km"))
        eff = (kcat / km) if (kcat and km) else None
        tvals = {t: _f(kin.get(f"{base}_{t}", {}).get("kcat")) for t in temps}
        th = therm.get(base, {})
        variants.append({"id": f"v{i}", "rep_gene": r["rep"], "organisms": r["organisms"],
                         "cluster_size": r["cluster_size"], "length": len(r["seq"]),
                         "kcat": kcat, "km": km, "kcat_km": eff, "kcat_by_temp": tvals,
                         "thermo": th.get("curve", {}), "thermo_label": th.get("label", ""),
                         "thermophilicity": th.get("thermophilicity", ""),
                         "sequence": r["seq"]})

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
        v["rank"] = rank + 1; v["score"] = round(score(v), 3)
    return variants


def write_bundle(rid, kos, sub, temps, stats, variants):
    out = {"rid": rid, "ko": kos, "substrate": sub, "temps": temps,
           "n_ko_genes": stats["n_ko_genes"], "n_sequences": stats["n_sequences"],
           "n_after_length_filter": stats["n_after_length"], "n_clusters": stats["n_clusters"],
           "n_variants": len(variants), "method": METHOD, "variants": variants}
    outp = os.path.join(OUTDIR, f"{rid}.json")
    json.dump(out, open(outp, "w"))
    print(f"  wrote {outp}: {len(variants)} ranked variants", flush=True)
    return outp


def _smiles(sub_cid, sub_smiles, smiles_db):
    return sub_smiles or smiles_db.get(sub_cid, "")


def characterize_single(rid, kos, sub_cid, sub_name, sub_smiles, temps, max_genes, max_reps, tax, smiles_db):
    smiles = _smiles(sub_cid, sub_smiles, smiles_db)
    print(f"[{rid}] KO {kos} substrate {sub_name or sub_cid}", flush=True)
    reps, stats = gather_reps(kos, max_genes, max_reps, tax, rid)
    if not reps:
        print("  no sequences — abort"); return None
    rows = []
    for i, r in enumerate(reps):
        sp = species_of(r["rep"], tax)[0]
        for t in temps:
            rows.append([f"{rid}_v{i}_{t}", r["seq"], smiles, sp, t])
    kin = run_mpek(rows, os.path.join(WORK, rid)) if smiles else {}
    faa = os.path.join(WORK, rid + "_reps.faa")
    with open(faa, "w") as fh:
        for i, r in enumerate(reps):
            fh.write(f">{rid}_v{i}\n{r['seq']}\n")
    therm = run_temstapro(faa, os.path.join(WORK, rid + "_therm.tsv"))
    variants = assemble_variants(reps, kin, therm, temps, rid)
    return write_bundle(rid, kos, {"cid": sub_cid, "name": sub_name, "smiles": smiles}, temps, stats, variants)


def characterize_batch(steps, temps, max_genes, max_reps, tax, smiles_db):
    """Whole-pathway: gather+cluster every step, then ONE MPEK + ONE TemStaPro over ALL variants (each
    gigabyte-scale model loads once instead of once-per-reaction), then split back to per-reaction bundles.
    This is why the UI runs the whole pathway from a single button rather than per enzyme."""
    print(f"[pathway] {len(steps)} reactions — one MPEK + one TemStaPro pass for the whole pathway", flush=True)
    prepared, all_rows, faa_lines = [], [], []
    for st in steps:
        rid = st.get("rid", "")
        kos = [k for k in (st.get("ko") or "").split(",") if k]
        if not rid or not kos:
            print(f"  [{rid or '?'}] no KO — skip"); continue
        if os.path.exists(os.path.join(OUTDIR, f"{rid}.json")):
            print(f"  [{rid}] already computed — skip"); continue
        smiles = _smiles(st.get("sub_cid", ""), "", smiles_db)
        reps, stats = gather_reps(kos, max_genes, max_reps, tax, rid)
        if not reps:
            print(f"  [{rid}] no sequences — skip"); continue
        for i, r in enumerate(reps):
            sp = species_of(r["rep"], tax)[0]
            faa_lines.append(f">{rid}_v{i}\n{r['seq']}\n")
            if smiles:
                for t in temps:
                    all_rows.append([f"{rid}_v{i}_{t}", r["seq"], smiles, sp, t])
        prepared.append({"rid": rid, "kos": kos, "reps": reps, "stats": stats,
                         "sub": {"cid": st.get("sub_cid", ""), "name": st.get("sub_name", ""), "smiles": smiles}})
    if not prepared:
        print("[pathway] nothing to characterize (all cached or no KOs)"); return []
    print(f"[pathway] MPEK on {len(all_rows)} rows across {len(prepared)} reactions (model loaded once)", flush=True)
    kin = run_mpek(all_rows, os.path.join(WORK, "pathway"))
    nvar = sum(len(p["reps"]) for p in prepared)
    print(f"[pathway] TemStaPro on {nvar} variants across {len(prepared)} reactions (model loaded once)", flush=True)
    faa = os.path.join(WORK, "pathway_reps.faa")
    open(faa, "w").write("".join(faa_lines))
    therm = run_temstapro(faa, os.path.join(WORK, "pathway_therm.tsv"))
    written = []
    for p in prepared:
        variants = assemble_variants(p["reps"], kin, therm, temps, p["rid"])
        written.append(write_bundle(p["rid"], p["kos"], p["sub"], temps, p["stats"], variants))
    print(f"[pathway] done — {len(written)} reaction bundles written", flush=True)
    return written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rid"); ap.add_argument("--ko", default="")
    ap.add_argument("--sub-cid", default=""); ap.add_argument("--sub-smiles", default=""); ap.add_argument("--sub-name", default="")
    ap.add_argument("--steps-json", default="", help="JSON file [{rid,ko,sub_cid,sub_name}] -> whole-pathway batch")
    ap.add_argument("--temps", default="37")
    ap.add_argument("--max-genes", type=int, default=600)
    ap.add_argument("--max-reps", type=int, default=120)
    ap.add_argument("--stage", default="all")          # accepted for back-compat, ignored
    args = ap.parse_args()
    tax = load_json(os.path.join(DATA, "taxonomy.json"), {})
    smiles_db = load_json(os.path.join(DATA, "smiles.json"), {})
    temps = [int(x) for x in args.temps.split(",") if x.strip()]
    if args.steps_json:
        characterize_batch(load_json(args.steps_json, []), temps, args.max_genes, args.max_reps, tax, smiles_db)
    else:
        kos = [k.strip() for k in args.ko.split(",") if k.strip()]
        characterize_single(args.rid, kos, args.sub_cid, args.sub_name, args.sub_smiles,
                            temps, args.max_genes, args.max_reps, tax, smiles_db)


if __name__ == "__main__":
    main()

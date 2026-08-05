#!/usr/bin/env python3
"""Precompute the static KEGG data bundle for the CLIENT-SIDE (browser) engine, so
github.io can run ANY start->end bioconversion in realtime with no backend.

Outputs (to docs/data/):
  network.json    carbon-skeleton graph: compounds + directed edges + reaction->KO/EC/dir
  compounds.json  {cid: name}                 (autocomplete + labels)
  taxonomy.json   {orgcode: [species,gram,domain]}   (genome gating -> species)
  ko/<KO>.json    {orgs:[...]}                 per-KO genome list (fetched on demand)
Thermo ΔG is layered in separately (bin/precompute_thermo.py) — until then the client uses
KEGG arrow directionality (already in network.json).

Resumable via the KeggClient disk cache. Run in the background; ~20-40 min first time.
"""
import sys, os, json, time
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from panroute.keggfetch import KeggClient, parse_reaction, parse_compound, parse_equation, carbon_count
from panroute.network import load_currency
from panroute.taxonomy import load_taxonomy, gram, species

OUT = os.path.join(ROOT, "docs", "data")
os.makedirs(OUT, exist_ok=True)
os.makedirs(os.path.join(OUT, "ko"), exist_ok=True)
cl = KeggClient(os.path.join(ROOT, "cache"))
currency = load_currency(os.path.join(ROOT, "assets", "currency_metabolites.tsv"))


def build_network_bundle():
    rxn_ids = [r.replace("rn:", "") for r, _ in cl.list_db("reaction")]
    sys.stderr.write(f"[net] fetching {len(rxn_ids)} reactions...\n")
    parsed = {}
    for i in range(0, len(rxn_ids), 10):
        recs = cl.get_entries([f"rn:{r}" for r in rxn_ids[i:i+10]])
        for rid in rxn_ids[i:i+10]:
            rec = recs.get(f"rn:{rid}")
            if rec:
                p = parse_reaction(rec)
                if p["id"]:
                    parsed[p["id"]] = p
        if i % 1000 == 0:
            sys.stderr.write(f"  {i}/{len(rxn_ids)}\n")
    # compounds referenced
    cids = set()
    for p in parsed.values():
        s, pr, _ = parse_equation(p["equation"])
        for _, c in s + pr:
            cids.add(c)
        for _, a, b in p["rclass"]:
            cids.add(a); cids.add(b)
    sys.stderr.write(f"[cpd] fetching {len(cids)} compounds...\n")
    comp = {}
    cids = [c for c in cids if c.startswith("C")]
    for i in range(0, len(cids), 10):
        recs = cl.get_entries([f"cpd:{c}" for c in cids[i:i+10]])
        for c in cids[i:i+10]:
            rec = recs.get(f"cpd:{c}")
            if rec:
                cp = parse_compound(rec)
                comp[c] = {"n": cp["name"], "c": carbon_count(cp["formula"])}
            else:
                comp[c] = {"n": c, "c": 0}

    # skeleton edges (RCLASS pairs, currency-excluded, carbon on both sides)
    def is_cur(c): return c in currency
    edges = []          # [src, dst, rid]
    rxn = {}            # rid -> {k:[kos], e:ec, d:dir}
    for rid, p in parsed.items():
        subs, prods, rev = parse_equation(p["equation"])
        sub, prod = {c for _, c in subs}, {c for _, c in prods}
        direction = "b" if rev else "f"
        rxn[rid] = {"k": list(dict.fromkeys(p["kos"])), "e": (p["ec"][0] if p["ec"] else ""), "d": direction}
        seen = set()
        for _rc, a, b in p["rclass"]:
            if a in sub and b in prod: s, d = a, b
            elif b in sub and a in prod: s, d = b, a
            else: continue
            if is_cur(s) or is_cur(d): continue
            if comp.get(s, {}).get("c", 0) == 0 or comp.get(d, {}).get("c", 0) == 0: continue
            if (rid, s, d) in seen: continue
            seen.add((rid, s, d))
            if direction in ("f", "b"): edges.append([s, d, rid])
            if direction == "b": edges.append([d, s, rid])
    json.dump({"compounds": comp, "edges": edges, "rxn": rxn,
               "kegg_release": cl.release()},
              open(f"{OUT}/network.json", "w"))
    json.dump({c: v["n"] for c, v in comp.items()}, open(f"{OUT}/compounds.json", "w"))
    sys.stderr.write(f"[net] {len(comp)} compounds, {len(edges)} edges, {len(rxn)} reactions -> network.json "
                     f"({os.path.getsize(f'{OUT}/network.json')//1024//1024} MB)\n")


def build_taxonomy():
    on, od, ol = load_taxonomy("/data/bioconversion/data")
    tax = {}
    for code, name in on.items():
        if od.get(code) in ("Bacteria", "Archaea"):
            tax[code] = [species(name), gram(ol.get(code, []), od.get(code)), od.get(code)]
    json.dump(tax, open(f"{OUT}/taxonomy.json", "w"))
    sys.stderr.write(f"[tax] {len(tax)} prokaryote genomes -> taxonomy.json\n")


if __name__ == "__main__":
    t = time.time()
    build_network_bundle()
    build_taxonomy()
    sys.stderr.write(f"[done] client bundle in {int(time.time()-t)}s (per-KO genome lists: run precompute_kogenomes.py)\n")

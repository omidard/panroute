#!/usr/bin/env python3
"""Precompute per-KO genome lists for the client-side engine: for every KO used by a
metabolic reaction in the bundled network, fetch the set of prokaryote genomes carrying it,
and write docs/data/ko/<KO>.json = {"orgs": [codes...]}. The browser fetches only the KOs a
query's routes need (a few hundred), so the genome-gating runs client-side on demand.

Reads docs/data/network.json (reaction KOs) + docs/data/taxonomy.json (prokaryote filter).
Resumable: skips KOs already written. ~30-60 min first time (rate-limited)."""
import sys, os, json
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from panroute.keggfetch import KeggClient

OUT = os.path.join(ROOT, "docs", "data", "ko")
os.makedirs(OUT, exist_ok=True)
cl = KeggClient(os.path.join(ROOT, "cache"))

net = json.load(open(os.path.join(ROOT, "docs", "data", "network.json")))
tax = json.load(open(os.path.join(ROOT, "docs", "data", "taxonomy.json")))
prok = set(tax)                                   # prokaryote org codes

# all KOs used by reactions + feedstock rules
kos = set()
for r in net["rxn"].values():
    kos.update(r.get("k", []))
import re
fr = json.load(open(os.path.join(ROOT, "assets", "feedstock_rules.json")))["feedstocks"]
for f in fr.values():
    for key in ("uptake", "overflow_ambiguous"):
        if f.get(key):
            kos.update(re.findall(r"K\d{5}", f[key]["definition"]))

kos = sorted(kos)
sys.stderr.write(f"[ko] {len(kos)} KOs to export\n")
done = 0
for i, ko in enumerate(kos):
    fp = os.path.join(OUT, f"{ko}.json")
    if os.path.exists(fp):
        continue
    orgs = set()
    for _src, gene in cl.link("genes", f"ko:{ko}"):
        code = gene.split(":")[0]
        if code in prok:
            orgs.add(code)
    json.dump({"orgs": sorted(orgs)}, open(fp, "w"))
    done += 1
    if i % 200 == 0:
        sys.stderr.write(f"  {i}/{len(kos)} ({done} fetched)\n")
sys.stderr.write(f"[ko] done — {len(kos)} KO files in docs/data/ko/\n")

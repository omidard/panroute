#!/usr/bin/env python3
"""Per-reaction genome prevalence = number of KEGG genomes that encode the reaction
(union of its KOs' genome lists, OR over isozymes). This is the biological-plausibility
weight the route search needs: atom-conservation (RCLASS) admits carbon-skeleton edges that
no organism actually runs (e.g. R00327 glucose->6-acetyl-glucose is in 0 genomes), which let
hop-count search find absurd shortcuts. Weighting edges by -log(prevalence) makes a real
10-step pathway (all high-prevalence reactions) cheaper than a 3-step path through a rare one.

Output docs/data/rxnprev.json = {rid: genome_count} for every reaction on a network edge."""
import json, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "data")
KO = os.path.join(OUT, "ko")
net = json.load(open(f"{OUT}/network.json"))
edge_rxns = sorted({e[2] for e in net["edges"]})

# cache KO -> genome count-set (load each ko file once)
ko_orgs = {}
def orgs(k):
    if k not in ko_orgs:
        p = os.path.join(KO, f"{k}.json")
        ko_orgs[k] = set(json.load(open(p))["orgs"]) if os.path.exists(p) else set()
    return ko_orgs[k]

prev = {}
for i, rid in enumerate(edge_rxns):
    kos = net["rxn"].get(rid, {}).get("k", [])
    S = set()
    for k in kos:
        S |= orgs(k)
    prev[rid] = len(S)
    if i % 2000 == 0:
        print(f"  {i}/{len(edge_rxns)}")
json.dump(prev, open(f"{OUT}/rxnprev.json", "w"))
nz = sum(1 for v in prev.values() if v > 0)
print(f"[rxnprev] {len(prev)} reactions · {nz} with >=1 genome · {len(prev)-nz} with ZERO genomes "
      f"(carbon-skeleton edges no organism runs) -> rxnprev.json ({os.path.getsize(f'{OUT}/rxnprev.json')//1024} KB)")

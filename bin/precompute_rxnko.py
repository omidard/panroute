#!/usr/bin/env python3
"""Per-reaction KO requirement grouped by EC component, so multi-subunit / multi-enzyme
complexes are gated correctly (AND across catalytic components, OR within = isozymes) instead
of the naive OR that over-credits any genome carrying a single shared subunit.

KEGG reaction ORTHOLOGY lists each KO with its [EC:...]; a genuine complex (e.g. pyruvate
dehydrogenase) spans several ECs (E1 1.2.4.1 · E2 2.3.1.12 · E3 1.8.1.4). We emit, for
reactions spanning >=2 EC groups, rxnko[rid] = [[kos of EC group 1], [kos of EC group 2], ...]
(AND across the groups). Single-EC reactions are omitted — the engine ORs their KOs as before.

Output docs/data/rxnko.json."""
import json, os, re, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from panroute.keggfetch import KeggClient
OUT = os.path.join(ROOT, "docs", "data")
net = json.load(open(f"{OUT}/network.json"))
cl = KeggClient(os.path.join(ROOT, "cache"), offline=True)

rxnko = {}
n_multi = 0
for rid in net["rxn"]:
    if len(net["rxn"][rid].get("k", [])) < 2:
        continue
    rec = cl.get_entries([f"rn:{rid}"]).get(f"rn:{rid}", "")
    groups = {}          # ec -> [kos]
    field = None
    for line in rec.splitlines():
        key = line[:12].strip()
        if key:
            field = key
        if field != "ORTHOLOGY":
            continue
        m = re.search(r"(K\d{5})", line)
        if not m:
            continue
        ko = m.group(1)
        ecm = re.search(r"\[EC:([^\]]+)\]", line)
        ec = ecm.group(1).split()[0] if ecm else "noec"
        groups.setdefault(ec, [])
        if ko not in groups[ec]:
            groups[ec].append(ko)
    # only store when there are >=2 genuine EC groups (a complex / multienzyme system)
    ecgroups = [v for k, v in groups.items() if k != "noec"]
    if len(ecgroups) >= 2:
        rxnko[rid] = [sorted(g) for g in ecgroups]
        n_multi += 1
json.dump(rxnko, open(f"{OUT}/rxnko.json", "w"))
print(f"[rxnko] {n_multi} multi-EC-component reactions -> rxnko.json "
      f"({os.path.getsize(f'{OUT}/rxnko.json')//1024} KB)")
print("  PDH R00209:", rxnko.get("R00209"))

#!/usr/bin/env python3
"""Build docs/data/smiles.json = {kegg_cid: SMILES} for the bundled compounds, so the route
report page can draw chemical structures client-side (SmilesDrawer). KEGG C -> MNX (chem_xref)
-> SMILES (chem_prop), from MetaNetX (open, redistributable)."""
import json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "data")
MNX = "/data/mnx_tmp"

net = json.load(open(f"{OUT}/network.json"))
want = set(net["compounds"].keys())
sys.stderr.write(f"[smiles] {len(want)} compounds wanted\n")

# KEGG C -> MNX
kegg2mnx = {}
with open(f"{MNX}/chem_xref.tsv") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 2 or not p[0].startswith("kegg.compound:"):
            continue
        cid = p[0].split(":", 1)[1]
        if cid in want:
            kegg2mnx.setdefault(cid, p[1])
sys.stderr.write(f"[smiles] {len(kegg2mnx)} mapped KEGG->MNX\n")

# MNX -> SMILES (col 9, 0-indexed 8)
need_mnx = set(kegg2mnx.values())
mnx_smiles = {}
with open(f"{MNX}/chem_prop.tsv") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if p[0] in need_mnx and len(p) >= 9 and p[8]:
            mnx_smiles[p[0]] = p[8]

smiles = {}
for cid, mnx in kegg2mnx.items():
    s = mnx_smiles.get(mnx)
    if s and "*" not in s and len(s) < 400:      # skip polymers/generic (*) and huge
        smiles[cid] = s
json.dump(smiles, open(f"{OUT}/smiles.json", "w"))
sys.stderr.write(f"[smiles] wrote {len(smiles)} structures -> smiles.json "
                 f"({os.path.getsize(f'{OUT}/smiles.json')//1024} KB)\n")

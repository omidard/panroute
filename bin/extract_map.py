#!/usr/bin/env python3
"""Extract the KEGG global metabolic map (map01100) layout for the live web UI:
compound node positions + REAL reaction polylines (the lines KEGG actually drew, which
bend through the network), plus EC->reaction and compound-pair indices so a route's
reactions can be resolved to on-map geometry even when the exact R-number is a variant
not drawn on the map (map01100 draws a curated subset; e.g. acetolactate synthase is
R00226 on the map, not R04672).

Input : reference KGML rn01100 (all reactions) + map01100.png (background image).
Output: assets/map01100/layout.json  +  copy of map01100.png
"""
import xml.etree.ElementTree as ET, json, os, sys, shutil

KGML = sys.argv[1] if len(sys.argv) > 1 else "/tmp/rn01100.kgml"
PNG  = sys.argv[2] if len(sys.argv) > 2 else "/home/omidard/repos/panGPR/data/maps/kegg/png/map01100.png"
OUT  = "/data/bioconversion/panroute/assets/map01100"
os.makedirs(OUT, exist_ok=True)

root = ET.parse(KGML).getroot()

compounds = {}          # cid -> [x, y]
reactions = {}          # rid -> [[x,y], ...] (first/main polyline)
rxn_ec = {}             # rid -> [EC, ...]  (from associated <reaction>/entry ec)
ec_rxn = {}             # EC -> [rid, ...] present on the map
pair_rxn = {}           # "Ca|Cb" (sorted) -> [rid, ...] present on the map

# compound positions (circles)
for e in root.findall("entry"):
    if e.get("type") == "compound":
        cid = e.get("name", "").replace("cpd:", "").split()[0]
        for g in e.findall("graphics"):
            if g.get("x") and g.get("y"):
                compounds[cid] = [float(g.get("x")), float(g.get("y"))]

# reaction polylines (entries carrying a reaction attr + line graphics)
for e in root.findall("entry"):
    ra = e.get("reaction")
    if not ra:
        continue
    line = None
    for g in e.findall("graphics"):
        if g.get("type") == "line" and g.get("coords"):
            nums = [float(x) for x in g.get("coords").split(",")]
            line = list(zip(nums[0::2], nums[1::2]))
            line = [[x, y] for x, y in line]
            break
    if line is None:
        continue
    for r in ra.split():
        r = r.replace("rn:", "")
        if r not in reactions:
            reactions[r] = line

# reaction -> substrates/products (from <reaction> elements) for pair index + EC
for rn in root.findall("reaction"):
    for r in rn.get("name", "").split():
        r = r.replace("rn:", "")
        subs = [s.get("name", "").replace("cpd:", "") for s in rn.findall("substrate")]
        prods = [p.get("name", "").replace("cpd:", "") for p in rn.findall("product")]
        for a in subs:
            for b in prods:
                key = "|".join(sorted([a, b]))
                pair_rxn.setdefault(key, [])
                if r not in pair_rxn[key]:
                    pair_rxn[key].append(r)

layout = {
    "image": {"file": "map01100.png", "width": 4961, "height": 3199},
    "compounds": compounds,
    "reactions": reactions,
    "pair_rxn": pair_rxn,
}
json.dump(layout, open(f"{OUT}/layout.json", "w"))
if os.path.exists(PNG):
    shutil.copy(PNG, f"{OUT}/map01100.png")

print(f"compounds: {len(compounds)} | reaction polylines: {len(reactions)} | "
      f"compound-pairs: {len(pair_rxn)}")
print(f"wrote {OUT}/layout.json ({os.path.getsize(f'{OUT}/layout.json')//1024} KB) + map01100.png")

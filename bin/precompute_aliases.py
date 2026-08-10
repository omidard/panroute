#!/usr/bin/env python3
"""Compound synonym groups so a search targets all forms of ONE metabolite (generic +
stereo/anomeric/protonation variants), WITHOUT merging chemically distinct metabolites.

Method: group by (base name with stereo/config descriptors stripped) + (carbon count),
then require identical molecular FORMULA within a group (from MetaNetX, fallback KEGG cache).
The NAME is the correct discriminator: (R)-/(S)-/generic acetoin share the base name "acetoin"
and merge; D-glucose vs D-galactose have DIFFERENT base names and stay separate (the InChIKey
connectivity block is stereo-blind and wrongly merged them). Generic class terms (D-Hexose,
D-Aldose) keep their own base name and do not merge into a specific sugar.

Output docs/data/aliases.json = {cid: [all cids in its group]}."""
import json, os, re
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "data")
MNX = "/data/mnx_tmp"
net = json.load(open(f"{OUT}/network.json"))
C = net["compounds"]

# molecular formula per KEGG compound (MetaNetX chem_prop col 4 via chem_xref)
kegg2mnx = {}
with open(f"{MNX}/chem_xref.tsv") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) >= 2 and p[0].startswith("kegg.compound:"):
            cid = p[0].split(":", 1)[1]
            if cid in C:
                kegg2mnx.setdefault(cid, p[1])
mnx_formula = {}
need = set(kegg2mnx.values())
with open(f"{MNX}/chem_prop.tsv") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if p[0] in need and len(p) >= 4 and p[3]:
            mnx_formula[p[0]] = p[3]
def formula(cid):
    return mnx_formula.get(kegg2mnx.get(cid, ""), "")

DESC = re.compile(
    r"^\(\s*[0-9RSEZ,\+\-'a-z ]+\)\s*-?\s*"                          # (R)- (2R,3S)- (+)-
    r"|^(meso|cis|trans|syn|anti|DL|D|L|alpha|beta|gamma|delta|epsilon|n|sn|o|m|p|"
    r"threo|erythro|allo|xylo|arabino|lyxo|ribo|gluco|galacto|manno)\s*-\s*",
    re.IGNORECASE)
def base(name):
    n = name.split(";")[0].strip(); prev = None
    while prev != n:
        prev = n
        n = DESC.sub("", n)
    return re.sub(r"\s+", " ", n).lower().strip()

groups = {}
for cid, v in C.items():
    if not v.get("n"):
        continue
    groups.setdefault((base(v["n"]), v.get("c", 0)), []).append(cid)

aliases = {}
n_multi = 0
for (bn, carb), cids in groups.items():
    if len(cids) < 2 or carb <= 0 or not bn:
        continue
    # split by molecular formula so only true same-compound variants merge
    byf = {}
    for c in cids:
        byf.setdefault(formula(c) or f"c{carb}", []).append(c)
    for f, members in byf.items():
        if len(members) > 1:
            n_multi += 1
            for c in members:
                aliases[c] = sorted(members)
json.dump(aliases, open(f"{OUT}/aliases.json", "w"))
print(f"[aliases] {len(aliases)} compounds in {n_multi} groups (base-name + formula) -> aliases.json "
      f"({os.path.getsize(f'{OUT}/aliases.json')//1024} KB)")
comp = {c: v["n"] for c, v in C.items()}
for probe in ["C00031", "C00466", "C03044", "C00186"]:
    print(f"  {probe} ({comp.get(probe)}):", [(c, comp.get(c)) for c in aliases.get(probe, [probe])])

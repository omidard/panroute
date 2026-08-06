#!/usr/bin/env python3
"""Compound identity groups so a search targets ALL forms of a metabolite (stereoisomers,
protonation/tautomer variants). Done by CHEMISTRY, not id-strings: KEGG compounds are grouped
by their InChIKey connectivity block (first 14 chars, stereo/charge-independent) via MetaNetX,
with name-base grouping as a fallback for compounds lacking an InChIKey.

Why it matters: KEGG splits e.g. acetoin into C00466 (generic, NO KO-annotated reactions),
C00810 ((R)-, the canonical alsD product WITH KOs) and C01769 ((S)-). Searching the generic
form alone returns 0 genomes. Merging by identity fixes this for every such metabolite.

Output docs/data/aliases.json = {cid: [all cids in its identity group]}."""
import json, os, re
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "data")
MNX = "/data/mnx_tmp"
net = json.load(open(f"{OUT}/network.json"))
C = net["compounds"]
want = set(C)

# KEGG C -> MNX
kegg2mnx = {}
with open(f"{MNX}/chem_xref.tsv") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) >= 2 and p[0].startswith("kegg.compound:"):
            cid = p[0].split(":", 1)[1]
            if cid in want:
                kegg2mnx.setdefault(cid, p[1])
# MNX -> InChIKey connectivity block
need = set(kegg2mnx.values())
mnx_ikey = {}
with open(f"{MNX}/chem_prop.tsv") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if p[0] in need and len(p) >= 8 and p[7].startswith("InChIKey="):
            mnx_ikey[p[0]] = p[7].split("=", 1)[1].split("-")[0]     # 14-char connectivity block

# group by InChIKey block
by_ikey = {}
for cid, mnx in kegg2mnx.items():
    ik = mnx_ikey.get(mnx)
    if ik:
        by_ikey.setdefault(ik, set()).add(cid)

# name-base fallback (for compounds without an InChIKey)
def base(name):
    n = name.split(";")[0].strip(); prev = None
    while prev != n:
        prev = n
        n = re.sub(r"^\(\s*[0-9RSEZ,\+\-'a-z ]+\)\s*-?\s*", "", n)
        n = re.sub(r"^(meso|cis|trans|DL|D|L|alpha|beta|n|sn|threo|erythro)\s*-\s*", "", n, flags=re.I)
    return re.sub(r"\s+", " ", n).lower().strip()
by_name = {}
for cid, v in C.items():
    if cid in mnx_ikey.get(kegg2mnx.get(cid, ""), ""):   # already has ikey group -> skip name
        pass
    if v.get("n"):
        by_name.setdefault((base(v["n"]), v.get("c", 0)), set()).add(cid)

# union: each compound's group = its ikey group ∪ its name group
group_of = {}
def union_groups(cids):
    for c in cids:
        group_of.setdefault(c, set()).update(cids)
for s in by_ikey.values():
    if len(s) > 1: union_groups(s)
for (bn, c), s in by_name.items():
    if len(s) > 1 and c > 0: union_groups(s)
# transitive closure (name+ikey overlaps)
changed = True
while changed:
    changed = False
    for c in list(group_of):
        g = group_of[c]
        merged = set(g)
        for m in g:
            merged |= group_of.get(m, {m})
        if merged != g:
            for m in merged: group_of[m] = merged
            changed = True

aliases = {c: sorted(g) for c, g in group_of.items() if len(g) > 1}
json.dump(aliases, open(f"{OUT}/aliases.json", "w"))
n_groups = len({frozenset(g) for g in aliases.values()})
print(f"[aliases] {len(aliases)} compounds in {n_groups} identity groups "
      f"(InChIKey + name) -> aliases.json ({os.path.getsize(f'{OUT}/aliases.json')//1024} KB)")

# --- SCAN: how many 'false-zero' compounds does this rescue? ---
ko_producer = set()   # compounds with >=1 KO-annotated producing reaction
for s, d, rid in net["edges"]:
    if net["rxn"].get(rid, {}).get("k"):
        ko_producer.add(d)
rescued = 0
examples = []
for c, g in aliases.items():
    if c not in ko_producer and any(m in ko_producer for m in g):
        rescued += 1
        if len(examples) < 10:
            examples.append(f"{c}({C[c]['n'][:20]})")
print(f"[scan] {rescued} compounds had NO KO-producing reaction of their own but a merged "
      f"identity-partner does — these would return a FALSE 0 if searched alone, now fixed.")
print("  examples:", ", ".join(examples))
print("  acetoin group:", aliases.get("C00466"), "| butanediol:", aliases.get("C03044"),
      "| L-lactate:", aliases.get("C00186"))

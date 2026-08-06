#!/usr/bin/env python3
"""Build compound stereoisomer/synonym groups so a search for e.g. 'acetoin' (C00466) also
targets (R)-acetoin (C00810) and (S)-acetoin (C01769) — otherwise the canonical route (which
often produces a specific stereoisomer, and whose generic form has no KO-annotated reactions)
is missed and returns 0 genomes.

Groups = compounds sharing the same base name (stereo/optical descriptors stripped) AND the
same carbon count. Output docs/data/aliases.json = {cid: [all cids in its group]}."""
import json, os, re
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "data")
net = json.load(open(f"{OUT}/network.json"))
C = net["compounds"]

def base(name):
    n = name.split(";")[0].strip()
    prev = None
    while prev != n:
        prev = n
        n = re.sub(r"^\(\s*[0-9RSEZ,\+\-'a-z ]+\)\s*-?\s*", "", n)            # (R)-, (2R,3S)-, (+)-
        n = re.sub(r"^(meso|cis|trans|DL|D|L|alpha|beta|n|sn|o|m|p|R|S|E|Z|threo|erythro)\s*-\s*", "", n, flags=re.I)
    return re.sub(r"\s+", " ", n).lower().strip()

groups = {}
for cid, v in C.items():
    if not v.get("n"):
        continue
    key = (base(v["n"]), v.get("c", 0))
    groups.setdefault(key, []).append(cid)

aliases = {}
n_multi = 0
for key, cids in groups.items():
    if len(cids) > 1 and key[1] > 0:          # only real multi-member carbon groups
        n_multi += 1
        for c in cids:
            aliases[c] = sorted(cids)
json.dump(aliases, open(f"{OUT}/aliases.json", "w"))
print(f"[aliases] {len(aliases)} compounds in {n_multi} multi-member groups -> aliases.json "
      f"({os.path.getsize(f'{OUT}/aliases.json')//1024} KB)")
print("acetoin group:", aliases.get("C00466"))
print("2,3-butanediol group:", aliases.get("C03044"))

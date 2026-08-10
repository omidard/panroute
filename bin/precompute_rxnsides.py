#!/usr/bin/env python3
"""Reactant/product sides per reaction, so route feasibility can apply ΔrG'° with the
correct SIGN for the direction each step actually traverses (audit A5).

thermo.json stores the as-written ΔG (left => right). A route step may cross a reaction in
reverse (a right-side compound -> a left-side compound); its effective ΔG is then -dg, not
+dg. Without sides the engine applied the as-written sign to every step, producing thousands
of false-feasible / false-infeasible calls on reverse traversals.

Output docs/data/rxnsides.json = {rid: {"L": [cids...], "R": [cids...]}} using raw KEGG
compound ids (the engine canonicalises both sides and the step endpoints with aliases.json)."""
import sys, os, json, glob, re
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from panroute.keggfetch import parse_reaction

OUT = os.path.join(ROOT, "docs", "data")
net = json.load(open(f"{OUT}/network.json"))
edge_rxns = {e[2] for e in net["edges"]}

CID = re.compile(r"C\d{5}")
def sides(eq):
    if "<=>" in eq: L, R = eq.split("<=>")
    elif "=>" in eq: L, R = eq.split("=>")
    elif "<=" in eq: L, R = eq.split("<=")
    else: return None
    return CID.findall(L), CID.findall(R)

out = {}
for f in glob.glob(os.path.join(ROOT, "cache", "get__rn_*.txt")):
    for block in open(f).read().split("///"):
        p = parse_reaction(block)
        if not (p and p.get("id") and p["id"] in edge_rxns): continue
        s = sides(p.get("equation", ""))
        if s and s[0] and s[1]:
            out[p["id"]] = {"L": sorted(set(s[0])), "R": sorted(set(s[1]))}
json.dump(out, open(f"{OUT}/rxnsides.json", "w"))
miss = len(edge_rxns) - len(out)
print(f"[rxnsides] {len(out)}/{len(edge_rxns)} edge reactions have parsed sides "
      f"({miss} without a two-sided equation) -> rxnsides.json "
      f"({os.path.getsize(f'{OUT}/rxnsides.json')//1024} KB)")
for r in ["R00754", "R00200", "R00703"]:
    print(f"  {r}: {out.get(r)}")

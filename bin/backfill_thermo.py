#!/usr/bin/env python3
"""Backfill ΔrG'° for routable reactions that the first bundle build dropped.

Root cause (audit A/H): Thermo.__call__ caches the arrow fallback tuple ('both', None, ...)
whenever _from_cc fails transiently (equilibrator/compound-cache not warm on the first pass),
and its `if rid in self._dg: return` short-circuit then blocks any retry forever. Central,
perfectly balanced reactions (pyruvate kinase R00200, LDH R00703, ADH R00754, aldehyde-DH
R00228, ...) were left with no numeric ΔG and dropped from docs/data/thermo.json.

This script force-retries _from_cc for every edge reaction whose cached ΔG is None (or is
absent), keeps numeric results, flags |ΔG|>200 kJ/mol as suspect (usually an unbalanced /
generic-R-group equation) without discarding it, and rewrites docs/data/thermo.json + the
thermo_dg.json cache. Idempotent: re-running only touches still-missing reactions."""
import sys, os, json, time, glob
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from panroute.keggfetch import parse_reaction
from panroute.thermo import Thermo

OUT = os.path.join(ROOT, "docs", "data")
net = json.load(open(f"{OUT}/network.json"))
edge_rxns = sorted({e[2] for e in net["edges"]})

# parse every reaction we can from the bundled KEGG cache (get__rn_*.txt)
parsed = {}
for f in glob.glob(os.path.join(ROOT, "cache", "get__rn_*.txt")):
    for block in open(f).read().split("///"):
        p = parse_reaction(block)
        if p and p.get("id"):
            parsed[p["id"]] = p
sys.stderr.write(f"[backfill] parsed {len(parsed)} reactions; {len(edge_rxns)} routable edges\n")

th = Thermo(parsed, consensus_path="/data/bioconversion/thermo/directionality_consensus.json",
            cache_path=os.path.join(ROOT, "cache", "thermo_dg.json"), use_equilibrator=True)

# which edge reactions currently have NO numeric ΔG (missing or None in cache)?
def numeric(rid):
    v = th._dg.get(rid)
    return isinstance(v, (list, tuple)) and isinstance(v[1], (int, float))
todo = [r for r in edge_rxns if r in parsed and not numeric(r)]
sys.stderr.write(f"[backfill] {len(todo)} edge reactions missing numeric ΔG -> retrying equilibrator\n")

t0 = time.time(); fixed = 0
for i, rid in enumerate(todo):
    res = th._from_cc(rid)            # force a fresh equilibrator attempt
    if res and isinstance(res[1], (int, float)):
        th._dg[rid] = res            # overwrite the stale None/arrow tuple
        fixed += 1
    if i % 100 == 0 and i:
        th.save()
        sys.stderr.write(f"  {i}/{len(todo)} · +{fixed} recovered · {int(time.time()-t0)}s\n")
th.save()

# rewrite the flat docs bundle; tag |ΔG|>200 suspects in a sidecar (not discarded)
flat = {r: round(v[1], 3) for r, v in th._dg.items() if isinstance(v[1], (int, float))}
json.dump(flat, open(f"{OUT}/thermo.json", "w"))
suspect = sorted(r for r, dg in flat.items() if abs(dg) > 200)
json.dump(suspect, open(f"{OUT}/thermo_suspect.json", "w"))
sys.stderr.write(f"[backfill] done · recovered {fixed} · thermo.json now {len(flat)} numeric "
                 f"({len(suspect)} |ΔG|>200 flagged suspect) · {int(time.time()-t0)}s\n")
for r in ["R00200", "R00703", "R00754", "R00228"]:
    sys.stderr.write(f"  {r}: {flat.get(r)}\n")

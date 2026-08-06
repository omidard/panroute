#!/usr/bin/env python3
"""Compute our own eQuilibrator component-contribution ΔrG′° for the routable reactions
(those in the carbon-skeleton network edges) and refresh docs/data/thermo.json + the ΔG in
rxninfo.json. Resumable via cache/thermo_dg.json. Slow (eQuilibrator per reaction); run in
the background — it improves ΔG coverage incrementally."""
import sys, os, json, time
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from panroute.keggfetch import KeggClient, parse_reaction
from panroute.thermo import Thermo

OUT = os.path.join(ROOT, "docs", "data")
cl = KeggClient(os.path.join(ROOT, "cache"), offline=True)   # reactions already cached
net = json.load(open(f"{OUT}/network.json"))

# routable reactions = those used by carbon-skeleton edges (what routes can traverse)
edge_rxns = sorted({e[2] for e in net["edges"]})
sys.stderr.write(f"[thermo] {len(edge_rxns)} routable reactions to price\n")

# parsed reactions (equations) from the KEGG cache
parsed = {}
for rid in edge_rxns:
    rec = cl.get_entries([f"rn:{rid}"]).get(f"rn:{rid}")
    if rec:
        p = parse_reaction(rec)
        if p["id"]:
            parsed[p["id"]] = p

thermo = Thermo(parsed, consensus_path="/data/bioconversion/thermo/directionality_consensus.json",
                cache_path=os.path.join(ROOT, "cache", "thermo_dg.json"), use_equilibrator=True)

t0 = time.time()
for i, rid in enumerate(edge_rxns):
    thermo(rid)                       # computes + caches
    if i % 100 == 0:
        thermo.save()
        num = sum(1 for v in thermo._dg.values() if isinstance(v[1], (int, float)))
        json.dump({r: v[1] for r, v in thermo._dg.items() if isinstance(v[1], (int, float))},
                  open(f"{OUT}/thermo.json", "w"))
        sys.stderr.write(f"  {i}/{len(edge_rxns)} · {num} numeric ΔG · {int(time.time()-t0)}s\n")
thermo.save()
json.dump({r: v[1] for r, v in thermo._dg.items() if isinstance(v[1], (int, float))},
          open(f"{OUT}/thermo.json", "w"))
sys.stderr.write(f"[thermo] done · {int(time.time()-t0)}s\n")

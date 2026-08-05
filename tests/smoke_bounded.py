#!/usr/bin/env python3
"""Fast bounded validation of the core (network + retro) on an explicit reaction set —
no slow network expansion. Verifies: (a) no currency leak, (b) acetone reachable from
acetyl-CoA, (c) the canonical 3-step route is recovered."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from panroute.keggfetch import KeggClient
from panroute.network import build_network, load_currency
from panroute.retro import reachable, enumerate_routes

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cl = KeggClient(os.path.join(ROOT, "cache"))
currency = load_currency(os.path.join(ROOT, "assets", "currency_metabolites.tsv"))

# explicit clostridial acetone neighbourhood + a couple of alternatives
RXN = ["R00238", "R00463", "R01357", "R01359", "R01366", "R10707", "R11026"]
C_ACCOA, C_ACETONE, C_AACOA, C_AA = "C00024", "C00207", "C00332", "C00164"

net = build_network(RXN, cl, currency, keep_endpoints={C_ACCOA, C_ACETONE}, thermo=None)
print("stats:", net.stats())
print("\nreactions resolved + skeleton pairs:")
for rid, p in net.reactions.items():
    pairs = [f"{a}->{b}" for _rc, a, b in p["rclass"]]
    print(f"  {rid}: {p['equation'][:60]:60s}  rclass={pairs}")

print("\ncurrency-leak check (NAD/CoA/ATP degree must be 0):")
for cur in ("C00003", "C00010", "C00002"):
    deg = len(net.edges_out.get(cur, [])) + len(net.edges_in.get(cur, []))
    print(f"  {cur}: {deg}")
    assert deg == 0, f"CURRENCY LEAK {cur}"

d = reachable(net, C_ACCOA, C_ACETONE)
print(f"\nmin length acetyl-CoA->acetone: {d}")
routes = enumerate_routes(net, C_ACCOA, C_ACETONE, max_len=6, max_routes=20)
print(f"routes: {len(routes)}")
for i, r in enumerate(routes, 1):
    chain = " -> ".join(net.compounds.get(c, {}).get("name", c).split(";")[0][:20] for c in r.path)
    print(f"  [{i}] len={r.length} via {[ '/'.join(x['rid'] for x in st['reactions']) for st in r.steps ]}")
    print(f"       {chain}")

assert d is not None, "acetone unreachable"
assert any(C_AA in r.path and C_AACOA in r.path for r in routes), "canonical route missing"
print("\nBOUNDED SMOKE TEST PASSED")

#!/usr/bin/env python3
"""Smoke test: does the carbon-skeleton retro-search recover the real acetone route
(acetyl-CoA -> acetoacetyl-CoA -> acetoacetate -> acetone) from live KEGG data, and does
it AVOID currency-metabolite shortcuts?  Runs on a bounded subnetwork around acetone."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from panroute.keggfetch import KeggClient
from panroute.network import build_network, load_currency
from panroute.retro import expand_subnetwork, reachable, enumerate_routes

CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cache")
ASSETS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets")

ACETYL_COA = "C00024"
ACETONE = "C00207"

def main():
    cl = KeggClient(CACHE)
    print("KEGG release:", cl.release())
    currency = load_currency(os.path.join(ASSETS, "currency_metabolites.tsv"))
    print(f"currency metabolites: {len(currency)}")

    print("expanding bounded subnetwork around acetone (depth 3)...")
    rxns = expand_subnetwork(cl, [ACETONE], max_depth=3, currency=currency)
    print(f"  reactions in subnetwork: {len(rxns)}")

    net = build_network(rxns, cl, currency,
                        keep_endpoints={ACETYL_COA, ACETONE}, thermo=None)
    print("network stats:", net.stats())

    # sanity: currency compounds must NOT be nodes with through-traffic
    for cur in ("C00003", "C00010", "C00002"):   # NAD, CoA, ATP
        deg = len(net.edges_out.get(cur, [])) + len(net.edges_in.get(cur, []))
        print(f"  currency {cur} degree in skeleton graph: {deg}  (expect 0)")
        assert deg == 0, f"CURRENCY LEAK: {cur} is a routing node!"

    d = reachable(net, ACETYL_COA, ACETONE)
    print(f"min route length acetyl-CoA -> acetone: {d}")
    assert d is not None, "acetone not reachable from acetyl-CoA!"

    routes = enumerate_routes(net, ACETYL_COA, ACETONE, max_len=6, max_routes=25)
    print(f"routes found: {len(routes)}")
    for i, r in enumerate(routes[:8], 1):
        chain = " -> ".join(net.compounds.get(c, {}).get("name", c).split(";")[0][:22]
                            for c in r.path)
        rxns_per_step = [ "/".join(x["rid"] for x in st["reactions"]) for st in r.steps ]
        print(f"  [{i}] len={r.length}  {chain}")
        print(f"       steps: {rxns_per_step}")
    # the canonical route must appear (3 steps, passes through acetoacetate C00164)
    got_canonical = any(ACETONE in r.path and "C00164" in r.path and r.length <= 4
                        for r in routes)
    print("canonical acetoacetate route present:", got_canonical)
    assert got_canonical, "did not recover the acetoacetate->acetone route"
    print("\nSMOKE TEST PASSED")

if __name__ == "__main__":
    main()

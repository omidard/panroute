#!/usr/bin/env python3
"""Streaming engine for the live web UI. run_query() is a GENERATOR that runs the real
PanRoute pipeline and yields events as it goes, so the frontend can animate true progress
(retro-search product->substrate on the KEGG map, organism discovery, per-route
thermodynamic feasibility). Every event carries real data — no mocks."""
from __future__ import annotations
import os
from collections import Counter

from .keggfetch import KeggClient
from .network import build_network, load_currency, fetch_reactions
from .thermo import Thermo
from .retro import expand_subnetwork, reachable, enumerate_routes
from .taxonomy import load_taxonomy
from . import genomes as G
from .feedstock import FeedstockGate, feedstock_kos
from . import mapviz

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "cache")
ASSETS = os.path.join(ROOT, "assets")
DATA = "/data/bioconversion/data"


def _cpd_name(net, cid):
    return net.compounds.get(cid, {}).get("name", cid).split(";")[0]


def run_query(start, end, feedstock=None, expand_depth=3, max_len=5, max_routes=60,
              thermo_routes=24, offline=False):
    """Yield (event_type, payload) tuples. Terminal event is ('done', summary)."""
    cl = KeggClient(CACHE, offline=offline)
    lay = mapviz.load_layout()
    currency = load_currency(os.path.join(ASSETS, "currency_metabolites.tsv"))
    keep = {start, end}

    yield "phase", {"msg": "expanding carbon-skeleton network from the product", "pct": 5}
    rxns = expand_subnetwork(cl, [end], max_depth=expand_depth, currency=currency)
    rxns |= expand_subnetwork(cl, [start], max_depth=1, currency=currency)

    parsed = fetch_reactions(cl, rxns)
    net = build_network(rxns, cl, currency, keep_endpoints=keep, thermo=None)

    yield "endpoints", {
        "start": {"cid": start, "name": _cpd_name(net, start), "xy": mapviz.endpoint_xy(start, lay)},
        "end":   {"cid": end,   "name": _cpd_name(net, end),   "xy": mapviz.endpoint_xy(end, lay)},
        "map": lay["image"]}

    yield "phase", {"msg": "searching routes product → feedstock", "pct": 20}
    d = reachable(net, start, end)
    routes = enumerate_routes(net, start, end, max_len=max_len, max_routes=max_routes)
    if not routes:
        yield "done", {"error": f"no route found from {start} to {end} within {max_len} steps"}
        return

    # attach names + resolve map geometry
    rroutes = []
    for i, r in enumerate(routes):
        rj = {"id": i, "length": r.length,
              "path": [{"cid": c, "name": _cpd_name(net, c)} for c in r.path],
              "steps": [{"from": s["from"], "to": s["to"],
                         "reactions": s["reactions"],
                         "enzymes": "/".join(sorted({k for x in s["reactions"] for k in x["ec"]})[:2])}
                        for s in r.steps]}
        rj["map"] = mapviz.resolve_route(rj, lay)
        rroutes.append(rj)

    # animate the shortest route on the map, retro order (product -> substrate).
    # enrich each step with compound names + on-map positions so the frontend can draw
    # on-map polylines AND off-map peripheral chips.
    shortest = min(rroutes, key=lambda x: x["length"])
    nm = {c: _cpd_name(net, c) for r in rroutes for c in [p["cid"] for p in r["path"]]}
    steps = shortest["map"]["steps"]
    for k, st in enumerate(steps):
        st = dict(st)
        st["from_name"] = nm.get(st["from"], st["from"])
        st["to_name"] = nm.get(st["to"], st["to"])
        st["from_xy"] = mapviz.endpoint_xy(st["from"], lay)
        st["to_xy"] = mapviz.endpoint_xy(st["to"], lay)
        yield "explore", {"step": st, "index": k, "total": len(steps),
                          "pct": 20 + int(25 * (k + 1) / max(1, len(steps)))}

    yield "routes", {"routes": rroutes, "shortest_len": d, "n_routes": len(rroutes)}

    # per-route thermodynamic feasibility (real eQuilibrator ΔG), streamed live
    yield "phase", {"msg": "computing thermodynamic feasibility of each route", "pct": 48}
    thermo = Thermo(parsed, consensus_path="/data/bioconversion/thermo/directionality_consensus.json",
                    cache_path=os.path.join(CACHE, "thermo_dg.json"), use_equilibrator=True)
    feas = {}
    for r in sorted(rroutes, key=lambda x: x["length"])[:thermo_routes]:
        ok = True; dgs = []
        for st in r["steps"]:
            rid = st["reactions"][0]["rid"]
            direction, dg, src = thermo(rid)
            dgs.append(dg)
            if direction == "r":            # can only run opposite the needed direction
                ok = False
        feas[r["id"]] = {"feasible": ok, "dG": [d for d in dgs if isinstance(d, (int, float))]}
        yield "thermo", {"route_id": r["id"], "feasible": ok, "length": r["length"],
                         "dG_sum": round(sum(x for x in dgs if isinstance(x, (int, float))), 1)}
    thermo.save()

    yield "phase", {"msg": "gating all KEGG genomes on the routes", "pct": 62}
    org_name, org_domain, org_lineage = load_taxonomy(DATA)
    route_kos = G.collect_route_kos(routes)
    fkos = feedstock_kos(os.path.join(ASSETS, "feedstock_rules.json"), feedstock) if feedstock else set()
    ko_orgs = G.fetch_ko_orgs(cl, route_kos | fkos, org_domain)
    org_kos = G.build_org_kos(ko_orgs)
    gate_species, per_genome = G.gate_all(routes, org_kos, org_name, org_domain, org_lineage)

    fg = FeedstockGate(os.path.join(ASSETS, "feedstock_rules.json")) if feedstock else None
    # stream organisms (encoding a route), richest first
    orgs = [row for row in per_genome if row["encodes_route"]]
    # collapse to species (best strain)
    by_sp = {}
    for row in orgs:
        s = row["species"]
        if s not in by_sp or row["n_routes"] > by_sp[s]["n_routes"]:
            by_sp[s] = row
    sp_rows = sorted(by_sp.values(), key=lambda r: -r["n_routes"])
    yield "phase", {"msg": f"found {len(sp_rows)} species with a native route", "pct": 80}
    for row in sp_rows[:400]:
        code = row["code"]
        feed = fg.status(feedstock, org_kos.get(code, set())) if fg else "n/a"
        # does this genome encode a thermodynamically-feasible route?
        tfeas = any(feas.get(i, {}).get("feasible", True) for i in row["route_idx"])
        yield "organism", {"species": row["species"], "domain": row["domain"], "gram": row["gram"],
                           "n_routes": row["n_routes"], "route_idx": row["route_idx"][:12],
                           "feedstock": feed, "thermo_feasible": tfeas}

    g = gate_species["T2"]
    gc = Counter(v["gram"] for v in g.values())
    up = ov = 0
    if fg:
        for s, v in g.items():
            stt = fg.status(feedstock, org_kos.get(v["example_code"], set()))
            up += stt == "uptake"; ov += stt == "overflow_capable"
    yield "done", {
        "n_routes": len(rroutes), "shortest": d,
        "T0": len(gate_species["T0"]), "T2": len(g),
        "T3": up if fg else None, "overflow_excluded": ov if fg else None,
        "gram": {"Gpos": gc["Gpos"], "Gneg": gc["Gneg"], "Arch": gc["Arch"], "Other": gc["Other"]},
        "kegg_release": net.kegg_release,
    }

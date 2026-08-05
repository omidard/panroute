#!/usr/bin/env python3
"""PanRoute end-to-end CLI: start metabolite + end metabolite -> feasible routes across
all KEGG genomes, with honest tiered reporting.

    python -m panroute.cli --start C00024 --end C00207 \
        --feedstock C00033 --out results/ --expand-depth 3

Sub-commands mirror the Nextflow processes so each stage is independently runnable and
cacheable (build-network / thermo / search / gate / report), or `run` does all of them.
"""
from __future__ import annotations
import argparse, os, sys, json

from .keggfetch import KeggClient
from .network import (build_network, load_currency, fetch_reactions,
                      save_network, load_network)
from .thermo import Thermo
from .retro import expand_subnetwork, reachable, enumerate_routes
from .taxonomy import load_taxonomy
from . import genomes as G
from .feedstock import FeedstockGate, feedstock_kos
from . import report as R
from .validate import validate

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def _paths(a):
    return {"cache": a.cache or os.path.join(ROOT, "cache"),
            "assets": a.assets or os.path.join(ROOT, "assets"),
            "data": a.data or "/data/bioconversion/data",
            "out": a.out or os.path.join(ROOT, "results")}


def cmd_run(a):
    P = _paths(a)
    os.makedirs(P["out"], exist_ok=True)
    cl = KeggClient(P["cache"], offline=a.offline)
    currency = load_currency(os.path.join(P["assets"], "currency_metabolites.tsv"))
    keep = {a.start, a.end}

    # 1. reaction set: explicit file, or bounded reverse-expansion from the product
    if a.reactions:
        rxns = set(l.strip() for l in open(a.reactions) if l.strip().startswith("R"))
    else:
        sys.stderr.write(f"[build] expanding subnetwork from {a.end} (depth {a.expand_depth})...\n")
        rxns = expand_subnetwork(cl, [a.end], max_depth=a.expand_depth, currency=currency)
        # also seed from start so short forward links are captured
        rxns |= expand_subnetwork(cl, [a.start], max_depth=1, currency=currency)
    sys.stderr.write(f"[build] {len(rxns)} reactions\n")

    # 2. thermo over the reaction set
    parsed = fetch_reactions(cl, rxns)
    thermo = Thermo(parsed,
                    consensus_path=a.consensus,
                    cache_path=os.path.join(P["cache"], "thermo_dg.json"),
                    use_equilibrator=not a.no_thermo)

    # 3. build carbon-skeleton network (thermo-gated edges)
    net = build_network(rxns, cl, currency, keep_endpoints=keep,
                        min_shared_c=a.min_shared_c, thermo=thermo)
    thermo.save()
    save_network(net, os.path.join(P["out"], f"network_{a.start}_{a.end}.pkl"))
    sys.stderr.write(f"[net] {net.stats()}\n")

    # 4. retro search
    d = reachable(net, a.start, a.end)
    if d is None:
        sys.stderr.write(f"[search] NO route {a.start} -> {a.end} in the skeleton network\n")
    routes = enumerate_routes(net, a.start, a.end, max_len=a.max_len, max_routes=a.max_routes)
    sys.stderr.write(f"[search] {len(routes)} routes (shortest {d})\n")

    # 5. genome gating (fetch only the KOs the routes use, + feedstock KOs)
    org_name, org_domain, org_lineage = load_taxonomy(P["data"])
    route_kos = G.collect_route_kos(routes)
    feed_kos = set()
    if a.feedstock:
        feed_kos = feedstock_kos(os.path.join(P["assets"], "feedstock_rules.json"), a.feedstock)
    sys.stderr.write(f"[gate] fetching {len(route_kos | feed_kos)} KO->genome maps...\n")
    ko_orgs = G.fetch_ko_orgs(cl, route_kos | feed_kos, org_domain)
    org_kos = G.build_org_kos(ko_orgs)
    gate_species, per_genome = G.gate_all(routes, org_kos, org_name, org_domain, org_lineage)
    sys.stderr.write(f"[gate] T2 species encoding a full route: {len(gate_species['T2'])}\n")

    # 6. feedstock direction-aware status per T2 species
    feed_status = None
    if a.feedstock:
        fg = FeedstockGate(os.path.join(P["assets"], "feedstock_rules.json"))
        feed_status = {}
        # map species -> representative org's KO set
        sp_code = {v["example_code"]: sp for sp, v in gate_species["T2"].items()}
        for code, sp in sp_code.items():
            feed_status[sp] = fg.status(a.feedstock, org_kos.get(code, set()))

    # 7. validation
    val = None
    truth = os.path.join(P["assets"], "validation_truth.tsv")
    if os.path.exists(truth):
        val = validate(a.end, gate_species["T2"], truth)

    # 8. report
    rep = R.funnel(a.start, a.end, routes, gate_species, per_genome,
                   thermo.coverage(), feedstock_status=feed_status,
                   validation=val, kegg_release=net.kegg_release)
    path = R.emit(P["out"], a.start, a.end, net, routes, gate_species, per_genome, rep)
    # per-genome table
    import csv
    with open(os.path.join(P["out"], f"per_genome_{a.start}_{a.end}.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["code", "name", "species", "domain", "gram",
                                           "encodes_route", "n_routes", "terminal_only"])
        w.writeheader()
        for row in per_genome:
            w.writerow({k: row[k] for k in w.fieldnames})
    print(json.dumps(rep, indent=2))
    sys.stderr.write(f"[done] {path}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="panroute")
    ap.add_argument("--start", required=True, help="start/feedstock KEGG compound id, e.g. C00024")
    ap.add_argument("--end", required=True, help="end/product KEGG compound id, e.g. C00207")
    ap.add_argument("--feedstock", help="feedstock cid for direction-aware uptake gating (e.g. C00033 acetate)")
    ap.add_argument("--reactions", help="file of reaction ids (one per line); else bounded expansion")
    ap.add_argument("--expand-depth", type=int, default=3)
    ap.add_argument("--max-len", type=int, default=12)
    ap.add_argument("--max-routes", type=int, default=200)
    ap.add_argument("--min-shared-c", type=int, default=1)
    ap.add_argument("--no-thermo", action="store_true", help="skip eQuilibrator (KEGG-arrow direction only)")
    ap.add_argument("--consensus", default="/data/bioconversion/thermo/directionality_consensus.json")
    ap.add_argument("--offline", action="store_true", help="use cache only; never hit KEGG")
    ap.add_argument("--cache"); ap.add_argument("--assets"); ap.add_argument("--data"); ap.add_argument("--out")
    a = ap.parse_args(argv)
    cmd_run(a)


if __name__ == "__main__":
    main()

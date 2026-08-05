#!/usr/bin/env python3
"""Honest tiered reporting (fix #4) + webapp-ready JSON emit.

Never presents a bare 'can produce' number. Every artifact carries: the funnel
(terminal-potential -> route-exists -> genome-encodes -> +feedstock -> +validated), the
thermo coverage, the KEGG release, and the standing caveats (gene != flux; KEGG sampling
bias; taxonomy heuristic; KO-OR simplification).
"""
from __future__ import annotations
import json, os

CAVEATS = [
    "Genome POTENTIAL only: a KO hit means the enzyme is ENCODED, not expressed, regulated, "
    "kinetically competent, or carrying flux. Not proof of production, titre or rate.",
    "KEGG genomes are culture/interest-biased, not a census of prokaryotic diversity; counts "
    "are lower bounds skewed to well-studied clades and drift with KEGG releases.",
    "Gram/taxonomy is a heuristic phylum mapping, not curated.",
    "Reaction KO requirements use OR-over-isozymes (v1); lone subunits of multi-subunit "
    "complexes may be over-credited unless a MODULE subunit-AND refinement was applied.",
    "Thermodynamics is standard-state (pH 7, I=0.25 M); in vivo concentrations can shift "
    "directionality. Routes flagged direction_source=kegg_arrow lack a computed ΔG.",
]


def funnel(start, end, routes, gate_species, per_genome, thermo_cov,
           feedstock_status=None, validation=None, kegg_release="unknown"):
    from collections import Counter
    t2 = gate_species["T2"]
    t0 = gate_species["T0"]

    def gc(m):
        c = Counter(v["gram"] for v in m.values())
        return {"species": len(m), "Gpos": c["Gpos"], "Gneg": c["Gneg"],
                "Arch": c["Arch"], "Other": c["Other"]}

    tiers = {
        "T0_terminal_enzyme_potential": gc(t0),
        "T1_route_exists": bool(routes),
        "T2_genome_encodes_full_route": gc(t2),
    }
    if feedstock_status is not None:
        up = {s for s, st in feedstock_status.items() if st == "uptake" and s in t2}
        ov = {s for s, st in feedstock_status.items() if st == "overflow_capable" and s in t2}
        tiers["T3_plus_feedstock_uptake"] = {
            "species": len(up),
            "species_overflow_only_excluded": len(ov)}
    if validation is not None:
        tiers["T4_validation"] = validation
    return {
        "query": {"start": start, "end": end},
        "kegg_release": kegg_release,
        "n_routes": len(routes),
        "shortest_route_len": min((r.length for r in routes), default=None),
        "thermo_coverage": thermo_cov,
        "tiers": tiers,
        "caveats": CAVEATS,
    }


def emit(outdir, start, end, net, routes, gate_species, per_genome, report_obj):
    os.makedirs(outdir, exist_ok=True)
    # 1. funnel/report
    json.dump(report_obj, open(f"{outdir}/report_{start}_{end}.json", "w"), indent=2)

    # 2. routes JSON (webapp: route DAG)
    routes_out = []
    for i, r in enumerate(routes):
        routes_out.append({
            "id": i, "length": r.length, "dg_total": r.dg_total,
            "dg_known_frac": round(r.dg_known_frac, 3),
            "path": [{"cid": c, "name": net.compounds.get(c, {}).get("name", c).split(";")[0]}
                     for c in r.path],
            "steps": [{"from": st["from"], "to": st["to"],
                       "reactions": [{"rid": x["rid"], "kos": x["kos"], "ec": x["ec"],
                                      "dg": x["dg"], "dg_source": x["dg_source"]}
                                     for x in st["reactions"]]}
                      for st in r.steps],
            "modules": sorted(r.modules)})
    json.dump({"query": {"start": start, "end": end}, "routes": routes_out},
              open(f"{outdir}/routes_{start}_{end}.json", "w"), indent=2)

    # 3. species table (T2, csv)
    import csv
    with open(f"{outdir}/species_{start}_{end}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "domain", "gram", "n_routes_encoded", "example_kegg_code", "tier"])
        for sp, v in sorted(gate_species["T2"].items()):
            w.writerow([sp, v["domain"], v["gram"], v["n_routes"], v["example_code"], "T2"])

    # 4. webapp payload (single self-describing bundle)
    json.dump({
        "query": {"start": start, "end": end,
                  "start_name": net.compounds.get(start, {}).get("name", start),
                  "end_name": net.compounds.get(end, {}).get("name", end)},
        "report": report_obj,
        "routes": routes_out,
        "network_stats": net.stats(),
    }, open(f"{outdir}/webapp_{start}_{end}.json", "w"), indent=2)
    return f"{outdir}/report_{start}_{end}.json"

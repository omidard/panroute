#!/usr/bin/env python3
"""Aggregate per-query PanRoute reports into one summary table (stdout JSON)."""
import sys, json

rows = []
for path in sys.argv[1:]:
    try:
        d = json.load(open(path))
    except Exception as e:
        sys.stderr.write(f"skip {path}: {e}\n"); continue
    t = d.get("tiers", {})
    rows.append({
        "start": d["query"]["start"], "end": d["query"]["end"],
        "kegg_release": d.get("kegg_release"),
        "n_routes": d.get("n_routes"),
        "shortest_route_len": d.get("shortest_route_len"),
        "T0_terminal_species": t.get("T0_terminal_enzyme_potential", {}).get("species"),
        "T2_route_species": t.get("T2_genome_encodes_full_route", {}).get("species"),
        "T3_feedstock_species": t.get("T3_plus_feedstock_uptake", {}).get("species"),
        "overflow_excluded": t.get("T3_plus_feedstock_uptake", {}).get("species_overflow_only_excluded"),
        "validation": t.get("T4_validation"),
        "thermo_fraction_real_dg": d.get("thermo_coverage", {}).get("fraction_real_dg"),
    })
print(json.dumps({"queries": rows, "n": len(rows)}, indent=2))

#!/usr/bin/env python3
"""Genome gating: which prokaryotic genomes ENCODE a complete enumerated route.

Efficient by design — we fetch KO->genome membership only for the KOs that actually
appear in the enumerated routes (+ feedstock), not all KOs for all genomes. For one
product that is ~10^2 KOs, a few hundred cached KEGG calls.

Route KO-satisfaction (v1):
  genome encodes route  <=>  for EVERY step, the genome has >=1 KO among the KOs of at
  least one reaction realising that step   (OR over isozymes and over parallel reactions).
  NOTE (documented simplification): reaction KO sets are treated as OR. Multi-subunit
  complexes (subunit-AND) are refined via KEGG MODULE grammar in genomes_modules.py when a
  reaction maps to a module; otherwise OR may over-credit a lone subunit. Reported in the
  quality notes, never hidden.

Tiers (the honest funnel, fix #4):
  T1  route exists at all (organism-independent)         -> boolean
  T2  genome encodes >=1 full route                      -> headline honest count
  T0  terminal-enzyme only (last step KOs)               -> the OLD inflated metric, for contrast
"""
from __future__ import annotations
from collections import defaultdict, Counter
from .taxonomy import gram, species, PROK


def collect_route_kos(routes):
    kos = set()
    for r in routes:
        for st in r.steps:
            for rx in st["reactions"]:
                kos.update(rx["kos"])
    return kos


def fetch_ko_orgs(client, kos, org_domain):
    """KO -> set(prokaryote org codes) via link/genes/ko:*  (cached)."""
    ko_orgs = {}
    for ko in sorted(kos):
        orgs = set()
        for _src, gene in client.link("genes", f"ko:{ko}"):
            code = gene.split(":")[0]
            if org_domain.get(code) in PROK:
                orgs.add(code)
        ko_orgs[ko] = orgs
    return ko_orgs


def build_org_kos(ko_orgs):
    org_kos = defaultdict(set)
    for ko, orgs in ko_orgs.items():
        for o in orgs:
            org_kos[o].add(ko)
    return org_kos


def step_ko_set(step):
    """Union of KOs across all reactions that realise a step (OR semantics)."""
    s = set()
    for rx in step["reactions"]:
        s.update(rx["kos"])
    return s


def genome_encodes_route(org_ko: set, route) -> bool:
    for st in route.steps:
        ks = step_ko_set(st)
        if not ks:                       # step has no KO annotation at all -> cannot credit
            return False
        if org_ko.isdisjoint(ks):
            return False
    return True


def gate_all(routes, org_kos, org_name, org_domain, org_lineage):
    """Return {tier: {species: {gram,domain,example_code,routes:[idx...]}}} and per-genome rows."""
    # candidate genomes = any that have >=1 route KO
    all_route_kos = collect_route_kos(routes)
    cand = {o for o, ks in org_kos.items() if ks & all_route_kos}

    # T0 terminal enzyme = last-step KO set of any route
    terminal_kos = set()
    for r in routes:
        if r.steps:
            terminal_kos |= step_ko_set(r.steps[-1])

    per_genome = []
    t2_species, t0_species = {}, {}
    for o in cand:
        oko = org_kos.get(o, set())
        enc = [i for i, r in enumerate(routes) if genome_encodes_route(oko, r)]
        term = bool(oko & terminal_kos)
        sp = species(org_name.get(o, o))
        g = gram(org_lineage.get(o, []), org_domain.get(o))
        dom = org_domain.get(o)
        per_genome.append({"code": o, "name": org_name.get(o, o), "species": sp,
                           "domain": dom, "gram": g,
                           "encodes_route": bool(enc), "n_routes": len(enc),
                           "route_idx": enc, "terminal_only": term and not enc})
        if enc and sp not in t2_species:
            t2_species[sp] = {"gram": g, "domain": dom, "example_code": o, "n_routes": len(enc)}
        if term and sp not in t0_species:
            t0_species[sp] = {"gram": g, "domain": dom, "example_code": o}
    return {"T0": t0_species, "T2": t2_species}, per_genome


def species_gram_counts(species_map):
    c = Counter(v["gram"] for v in species_map.values())
    return {"n_species": len(species_map), "Gpos": c["Gpos"], "Gneg": c["Gneg"],
            "Arch": c["Arch"], "Other": c["Other"]}

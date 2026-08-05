#!/usr/bin/env python3
"""Build the carbon-skeleton reaction network: a DIRECTED multigraph over compounds
where an edge X -> Y means "some reaction can, in a thermodynamically-allowed direction,
convert carbon-skeleton X into carbon-skeleton Y", carrying the reaction id, the KO
requirement, EC and ΔrG′°.

Currency metabolites (cofactors, Pi, CO2 hub, …) are excluded as intermediates so the
search cannot shortcut through the cofactor pool — EXCEPT when a currency metabolite is
the user's explicit start/end (e.g. methanol or CO2 as feedstock). Edges are derived
ONLY from KEGG RCLASS atom-conserved substrate<->product pairs, never from arbitrary
substrate x product combinations. See docs/ARCHITECTURE.md §1.
"""
from __future__ import annotations
import os, pickle
from collections import defaultdict
from dataclasses import dataclass, field

from .keggfetch import (KeggClient, parse_reaction, parse_compound,
                        parse_equation, carbon_count)


@dataclass
class Edge:
    src: str            # compound consumed
    dst: str            # compound produced
    reaction: str       # Rxxxxx
    rclass: str
    direction: str      # 'f' (as-written) or 'r' (reverse-of-written) used for this edge
    kos: tuple          # KO alternatives catalysing the reaction (OR semantics, v1)
    ec: tuple
    dg: float | None = None            # ΔrG′° in the direction of this edge (kJ/mol)
    dg_source: str = "none"


@dataclass
class Network:
    edges_out: dict = field(default_factory=lambda: defaultdict(list))  # cpd -> [Edge]
    edges_in: dict = field(default_factory=lambda: defaultdict(list))   # cpd -> [Edge]
    compounds: dict = field(default_factory=dict)   # cid -> {name,formula,carbons}
    reactions: dict = field(default_factory=dict)   # rid -> parsed reaction
    currency: set = field(default_factory=set)
    kegg_release: str = "unknown"

    def add_edge(self, e: Edge):
        self.edges_out[e.src].append(e)
        self.edges_in[e.dst].append(e)

    def stats(self) -> dict:
        n_edges = sum(len(v) for v in self.edges_out.values())
        return {"compounds": len(self.compounds), "reactions": len(self.reactions),
                "directed_edges": n_edges, "currency_excluded": len(self.currency),
                "kegg_release": self.kegg_release}


def fetch_reactions(client: KeggClient, reaction_ids) -> dict:
    """Fetch + parse a set of reaction ids -> {rid: parsed}. Cached; used to seed Thermo
    before building the network (build_network re-reads from cache, so no double fetch)."""
    recs = client.get_entries([f"rn:{r}" for r in reaction_ids])
    out = {}
    for rid in reaction_ids:
        rec = recs.get(f"rn:{rid}")
        if rec:
            p = parse_reaction(rec)
            if p["id"]:
                out[p["id"]] = p
    return out


def load_currency(path: str) -> set:
    cur = set()
    with open(path) as fh:
        next(fh)                                    # header
        for line in fh:
            cid = line.split("\t", 1)[0].strip()
            # C00024 acetyl-CoA is explicitly KEPT (real skeleton node)
            if cid and cid != "C00024":
                cur.add(cid)
    return cur


def build_network(reaction_ids, client: KeggClient, currency: set,
                  keep_endpoints: set | None = None,
                  min_shared_c: int = 1,
                  thermo=None) -> Network:
    """Assemble the network from a set of reaction ids.

    keep_endpoints: currency ids that must remain routable (user start/end).
    min_shared_c:   require min(C(src),C(dst)) >= this (RCLASS already guarantees atom
                    conservation; this is a backstop, default 1 = permissive for C1).
    thermo:         optional callable rid -> (direction, dg_f, source) where direction in
                    {'f','r','both','unknown'} and dg_f is ΔrG′° of the as-written forward.
    """
    keep_endpoints = keep_endpoints or set()
    net = Network(currency=set(currency), kegg_release=client.release())

    # 1. fetch + parse reactions
    recs = client.get_entries([f"rn:{r}" for r in reaction_ids])
    parsed = {}
    for rid in reaction_ids:
        rec = recs.get(f"rn:{rid}")
        if not rec:
            continue
        p = parse_reaction(rec)
        if p["id"]:
            parsed[p["id"]] = p
    net.reactions = parsed

    # 2. fetch compound formulas for everything referenced
    cids = set()
    for p in parsed.values():
        subs, prods, _ = parse_equation(p["equation"])
        for _, c in subs + prods:
            cids.add(c)
        for _, a, b in p["rclass"]:
            cids.add(a); cids.add(b)
    crecs = client.get_entries([f"cpd:{c}" for c in cids if c.startswith("C")])
    for c in cids:
        rec = crecs.get(f"cpd:{c}")
        if rec:
            cp = parse_compound(rec)
            net.compounds[c] = {"name": cp["name"], "formula": cp["formula"],
                                "carbons": carbon_count(cp["formula"])}
        else:
            net.compounds[c] = {"name": c, "formula": "", "carbons": 0}

    def is_currency(c):
        return c in currency and c not in keep_endpoints

    # 3. build edges from RCLASS atom-conserved pairs
    for rid, p in parsed.items():
        subs, prods, rev = parse_equation(p["equation"])
        sub_ids = {c for _, c in subs}
        prod_ids = {c for _, c in prods}
        # thermodynamic direction for this reaction
        if thermo:
            direction, dg_f, dgsrc = thermo(rid)
        else:
            direction, dg_f, dgsrc = ("both" if rev else "f"), None, "kegg_arrow"
        kos = tuple(dict.fromkeys(p["kos"]))
        ec = tuple(p["ec"])
        seen_pairs = set()
        for rc, a, b in p["rclass"]:
            # orient the pair to (substrate_side, product_side) as written
            if a in sub_ids and b in prod_ids:
                s, d = a, b
            elif b in sub_ids and a in prod_ids:
                s, d = b, a
            else:
                # pair not cleanly on opposite sides (e.g. same-side or generic) -> skip
                continue
            if is_currency(s) or is_currency(d):
                continue
            cs, cd = net.compounds.get(s, {}).get("carbons", 0), net.compounds.get(d, {}).get("carbons", 0)
            if cs == 0 or cd == 0:
                continue                     # need carbon on both sides of a skeleton edge
            if min(cs, cd) < min_shared_c:
                continue
            key = (rid, s, d)
            if key in seen_pairs:
                continue
            seen_pairs.add(key)
            # emit directed edges consistent with allowed thermodynamic direction
            #   as-written forward  = s -> d  (dg = dg_f)
            #   reverse             = d -> s  (dg = -dg_f)
            if direction in ("f", "both", "unknown"):
                net.add_edge(Edge(s, d, rid, rc, "f", kos, ec, dg_f, dgsrc))
            if direction in ("r", "both", "unknown"):
                dg_r = (-dg_f) if isinstance(dg_f, (int, float)) else None
                net.add_edge(Edge(d, s, rid, rc, "r", kos, ec, dg_r, dgsrc))
    return net


def save_network(net: Network, path: str):
    with open(path, "wb") as fh:
        pickle.dump(net, fh)


def load_network(path: str) -> Network:
    with open(path, "rb") as fh:
        return pickle.load(fh)

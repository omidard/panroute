#!/usr/bin/env python3
"""Retrosynthetic route search over the carbon-skeleton network.

A *route* is a simple carbon-skeleton path  start -> c1 -> ... -> end.  Each step
(u->v) is realised by one or more reactions (alternative enzymes); a genome later
"encodes" the route iff for EVERY step it has >=1 realising reaction whose KO
requirement it satisfies (see genomes.py).

Search is anchored at the END product and pruned by reverse-BFS distance so we only
expand compounds that can still reach the start within the length budget — this is the
'retro' framing the user asked for (start from the product, work back to the feedstock)
and keeps enumeration tractable on the full KEGG network.
"""
from __future__ import annotations
from collections import deque, defaultdict
from dataclasses import dataclass, field


@dataclass
class Route:
    path: list                      # [start, c1, ..., end]
    steps: list                     # per step: {"from","to","reactions":[{rid,kos,ec,dg,dir}]}
    length: int = 0
    dg_total: float | None = None
    dg_known_frac: float = 0.0
    modules: set = field(default_factory=set)

    def reaction_ids(self):
        s = set()
        for st in self.steps:
            for r in st["reactions"]:
                s.add(r["rid"])
        return s


def reverse_bfs_dist(net, end: str) -> dict:
    """Min #edges from any compound to `end`, following directed edges forward
    (computed on the reversed graph starting at end)."""
    dist = {end: 0}
    q = deque([end])
    while q:
        v = q.popleft()
        for e in net.edges_in.get(v, []):      # edges u->v ; predecessor u
            u = e.src
            if u not in dist:
                dist[u] = dist[v] + 1
                q.append(u)
    return dist


def reachable(net, start: str, end: str) -> int | None:
    """Minimum route length (in reactions) from start to end, or None."""
    dist = reverse_bfs_dist(net, end)
    return dist.get(start)


def _collapse_steps(net, path):
    """For a compound path, gather all reactions realising each consecutive step."""
    steps = []
    modules = set()
    dg_vals = []
    for u, v in zip(path, path[1:]):
        rxns = []
        seen = set()
        for e in net.edges_out.get(u, []):
            if e.dst == v and e.reaction not in seen:
                seen.add(e.reaction)
                rxns.append({"rid": e.reaction, "kos": list(e.kos), "ec": list(e.ec),
                             "dg": e.dg, "dir": e.direction, "dg_source": e.dg_source})
                modules.update(net.reactions.get(e.reaction, {}).get("modules", []))
        # best (most negative known) ΔG for the step, for scoring
        kn = [r["dg"] for r in rxns if isinstance(r["dg"], (int, float))]
        dg_vals.append(min(kn) if kn else None)
        steps.append({"from": u, "to": v, "reactions": rxns})
    known = [d for d in dg_vals if d is not None]
    dg_total = sum(known) if known else None
    dg_frac = (len(known) / len(dg_vals)) if dg_vals else 0.0
    return steps, modules, dg_total, dg_frac


def enumerate_routes(net, start: str, end: str, max_len: int = 12,
                     max_routes: int = 200, max_expansions: int = 500000) -> list[Route]:
    """Enumerate distinct simple carbon-skeleton routes start->end, length <= max_len.

    Uses reverse-BFS distance-to-end to prune: from a node at path-length L, only expand
    successors whose dist_to_end <= (max_len - L - 1). Returns routes ranked by
    (length, -dg_known_frac, dg_total)."""
    dist = reverse_bfs_dist(net, end)
    if start not in dist or dist[start] > max_len:
        return []
    routes: list[Route] = []
    expansions = 0

    # DFS with on-path set to keep routes simple (no repeated compound)
    stack = [(start, [start], {start})]
    while stack and len(routes) < max_routes and expansions < max_expansions:
        node, path, onpath = stack.pop()
        if node == end and len(path) > 1:
            steps, modules, dg_total, dg_frac = _collapse_steps(net, path)
            routes.append(Route(path=path, steps=steps, length=len(path) - 1,
                                dg_total=dg_total, dg_known_frac=dg_frac, modules=modules))
            continue
        remaining = max_len - (len(path) - 1)
        if remaining <= 0:
            continue
        # unique successors reachable within budget
        nxt = {}
        for e in net.edges_out.get(node, []):
            v = e.dst
            if v in onpath:
                continue
            d = dist.get(v)
            if d is None or d > remaining - 1:
                continue
            nxt[v] = True
        expansions += 1
        # push successors (closer-to-end first so best routes surface early)
        for v in sorted(nxt, key=lambda x: dist.get(x, 1e9), reverse=True):
            stack.append((v, path + [v], onpath | {v}))

    routes.sort(key=lambda r: (r.length, -r.dg_known_frac,
                               r.dg_total if r.dg_total is not None else 1e9))
    return routes[:max_routes]


def expand_subnetwork(client, seed_cpds, max_depth: int, currency: set,
                      max_reactions: int = 20000) -> set:
    """Reverse-expand a bounded reaction set around seed compounds by walking KEGG
    link/reaction/cpd out to `max_depth` carbon-skeleton hops. Used for on-demand /
    bounded queries (and the smoke test) without fetching all of KEGG.

    Returns a set of reaction ids. Skeleton neighbours are read from each reaction's
    RCLASS pairs, skipping currency compounds, so expansion follows carbon, not cofactors.
    """
    from .keggfetch import parse_reaction
    rxn_ids = set()
    frontier = set(seed_cpds)
    visited_cpd = set()
    for _ in range(max_depth):
        if not frontier or len(rxn_ids) >= max_reactions:
            break
        # reactions touching any frontier compound
        new_rxns = set()
        for c in frontier:
            for _src, tgt in client.link("reaction", f"cpd:{c}"):
                rid = tgt.split(":")[-1]
                if rid.startswith("R"):
                    new_rxns.add(rid)
        new_rxns -= rxn_ids
        rxn_ids |= new_rxns
        visited_cpd |= frontier
        # next frontier = carbon-skeleton neighbours of these reactions
        recs = client.get_entries([f"rn:{r}" for r in new_rxns])
        nxt = set()
        for r in new_rxns:
            rec = recs.get(f"rn:{r}")
            if not rec:
                continue
            for _rc, a, b in parse_reaction(rec)["rclass"]:
                for x in (a, b):
                    if x.startswith("C") and x not in currency and x not in visited_cpd:
                        nxt.add(x)
        frontier = nxt
    return rxn_ids

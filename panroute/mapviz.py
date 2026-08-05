#!/usr/bin/env python3
"""Resolve a PanRoute route onto the KEGG global map (map01100) geometry for the live UI.

For each step we return REAL geometry, in this order of preference (honest — never a fake
hub line):
  1. the reaction's own drawn polyline (R-number present on the map), OR
  2. an on-map reaction that performs the same compound-pair transformation (map01100 draws
     a curated subset, so R04672 acetolactate synthase -> use R00226's polyline), OR
  3. a direct connector between the two compounds' real node positions (both on map), OR
  4. 'offmap': at least one endpoint is a specialised metabolite not on the core map ->
     the frontend renders it as a peripheral chip with a dashed connector.
"""
from __future__ import annotations
import json, os, math

_LAYOUT = None


def load_layout(path=None):
    global _LAYOUT
    if _LAYOUT is None:
        path = path or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                    "assets", "map01100", "layout.json")
        _LAYOUT = json.load(open(path))
    return _LAYOUT


def _nearest_onmap(cid, lay):
    return lay["compounds"].get(cid)


def resolve_step(from_cid, to_cid, reaction, lay):
    rxns = lay["reactions"]
    comp = lay["compounds"]
    # 1. exact reaction polyline
    if reaction in rxns:
        return {"kind": "polyline", "coords": rxns[reaction], "reaction": reaction}
    # 2. same-transformation reaction drawn on the map
    key = "|".join(sorted([from_cid, to_cid]))
    for r in lay["pair_rxn"].get(key, []):
        if r in rxns:
            return {"kind": "polyline", "coords": rxns[r], "reaction": r}
    # 3. both compounds on map -> connector between real node positions
    a, b = comp.get(from_cid), comp.get(to_cid)
    if a and b:
        return {"kind": "connector", "coords": [a, b], "reaction": reaction}
    # 4. off map
    return {"kind": "offmap", "coords": None, "reaction": reaction}


def resolve_route(route, lay=None):
    """route: dict with 'path' [{cid,name}] and 'steps' [{from,to,reactions:[{rid,...}]}].
    Returns per-node map positions (or None if off-map) and per-step resolved geometry,
    ordered from END back to START (retro direction) for the animation."""
    lay = lay or load_layout()
    comp = lay["compounds"]
    nodes = []
    for p in route["path"]:
        xy = comp.get(p["cid"])
        nodes.append({"cid": p["cid"], "name": p["name"], "xy": xy, "onmap": xy is not None})
    steps = []
    for st in route["steps"]:
        rid = st["reactions"][0]["rid"] if st.get("reactions") else st.get("reaction", "")
        geo = resolve_step(st["from"], st["to"], rid, lay)
        steps.append({"from": st["from"], "to": st["to"], **geo,
                      "enzymes": st.get("enzymes", "")})
    # retro order = reverse (product first)
    return {"nodes": list(reversed(nodes)), "steps": list(reversed(steps))}


def endpoint_xy(cid, lay=None):
    lay = lay or load_layout()
    return lay["compounds"].get(cid)


def onmap_fraction(route, lay=None):
    lay = lay or load_layout()
    tot = len(route["steps"]) or 1
    on = sum(1 for st in route["steps"]
             if resolve_step(st["from"], st["to"],
                             st["reactions"][0]["rid"] if st.get("reactions") else "",
                             lay)["kind"] in ("polyline", "connector"))
    return on / tot

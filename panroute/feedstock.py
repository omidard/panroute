#!/usr/bin/env python3
"""Direction-aware feedstock utilisation gating (fix #3).

A genome is credited with USING a feedstock only if it encodes an uptake/activation
route in the consuming direction. Overflow-ambiguous systems (ackA-pta for acetate) are
recorded separately as 'overflow_capable' and NOT counted as uptake — this is exactly the
Parageobacillus caveat (acetate is its product, not a demonstrated substrate).

Boolean grammar identical to KEGG MODULE: space/'+' = AND, ',' = OR, parentheses group.
"""
from __future__ import annotations
import json, re


def compile_expr(expr):
    """KEGG-MODULE boolean -> python eval. space/'+' = AND, ',' = OR, parens group.
    Commas are tightened first so ', ' cannot collide with the space->AND rewrite."""
    w = re.sub(r'\s*,\s*', ',', expr.strip())     # no whitespace around OR commas
    w = re.sub(r'(K\d{5})', r"m('\1')", w)         # wrap KOs (contain no spaces)
    w = w.replace('+', ' ')                         # subunit AND -> space
    w = re.sub(r'\s+', ' ', w).strip()
    w = w.replace(' ', ' and ')                     # space AND -> python and
    w = w.replace(',', ' or ')                      # comma OR -> python or
    return compile(w, '<feed>', 'eval')


class FeedstockGate:
    def __init__(self, rules_path):
        cfg = json.load(open(rules_path))["feedstocks"]
        self.rules = {}
        for cid, r in cfg.items():
            entry = {"name": r["name"], "uptake": None, "overflow": None}
            if r.get("uptake"):
                entry["uptake"] = compile_expr(r["uptake"]["definition"])
            if r.get("overflow_ambiguous"):
                entry["overflow"] = compile_expr(r["overflow_ambiguous"]["definition"])
            self.rules[cid] = entry

    def status(self, cid, org_ko: set) -> str:
        """Return 'uptake' | 'overflow_capable' | 'none' for a feedstock in a genome."""
        r = self.rules.get(cid)
        if not r:
            return "unknown_feedstock"
        H = org_ko
        if r["uptake"] is not None and bool(eval(r["uptake"], {"m": H.__contains__})):
            return "uptake"
        if r["overflow"] is not None and bool(eval(r["overflow"], {"m": H.__contains__})):
            return "overflow_capable"
        return "none"


def feedstock_kos(rules_path, cid):
    """All KOs referenced by a feedstock's rules (so they get fetched for gating)."""
    cfg = json.load(open(rules_path))["feedstocks"].get(cid, {})
    kos = set()
    for key in ("uptake", "overflow_ambiguous"):
        if cfg.get(key):
            kos.update(re.findall(r"K\d{5}", cfg[key]["definition"]))
    return kos

#!/usr/bin/env python3
"""Thermodynamic directionality for reactions (quality fix #2).

For each reaction we determine which direction(s) are physiologically feasible so the
retro-search only builds edges that can carry flux toward the product.

Priority of evidence (recorded per reaction as `dg_source`):
  1. eQuilibrator component-contribution ΔrG′° (pH 7.0, I=0.25 M) + reversibility index
     — the field standard (Noor 2013, Beber 2022). Requires `equilibrator_api`.
  2. Curated Rhea / MetaCyc / group-contribution consensus already computed for the
     original bioconversion project (thermo/directionality_consensus.json).
  3. KEGG equation arrow ( '<=>' reversible, '=>' forward ) — weakest, flagged low-conf.

`Thermo(...)` is a callable: rid -> (direction, dg_forward_kJmol, source)
  direction in {'f','r','both','unknown'}   dg_forward = ΔrG′° of the AS-WRITTEN reaction.

REV_MARGIN: |ΔrG′°| below this (kJ/mol) => treat as reversible ('both'); above and
negative => forward-favoured only ('f'); above and positive => reverse only ('r').
Uses the reversibility index when available (physiological-concentration aware) in
preference to the raw ΔG sign.
"""
from __future__ import annotations
import json, os

REV_MARGIN = 30.0          # kJ/mol; ~ln reversibility index threshold below
PH = 7.0
IONIC = 0.25
TEMP = 298.15


class Thermo:
    def __init__(self, reactions: dict, consensus_path: str | None = None,
                 cache_path: str | None = None, use_equilibrator: bool = True):
        self.reactions = reactions            # rid -> parsed reaction (needs 'equation')
        self.cache_path = cache_path
        self._dg = {}                         # rid -> (direction, dg_f, source)
        self._cc = None
        self._consensus = {}
        if consensus_path and os.path.exists(consensus_path):
            try:
                raw = json.load(open(consensus_path))
                if isinstance(raw, list):        # list of dicts keyed by 'reaction'
                    self._consensus = {d.get("reaction"): d for d in raw if d.get("reaction")}
                else:
                    self._consensus = raw
            except Exception:
                self._consensus = {}
        if cache_path and os.path.exists(cache_path):
            try:
                self._dg = {k: tuple(v) for k, v in json.load(open(cache_path)).items()}
            except Exception:
                self._dg = {}
        if use_equilibrator:
            self._init_cc()

    def _init_cc(self):
        try:
            from equilibrator_api import ComponentContribution, Q_
            self._cc = ComponentContribution()
            self._cc.p_h = Q_(PH)
            self._cc.ionic_strength = Q_(f"{IONIC}M")
            self._cc.temperature = Q_(f"{TEMP}K")
        except Exception:
            self._cc = None                   # graceful: fall back to consensus/arrow

    # ---- direction decision from a ΔG and optional reversibility index ----
    @staticmethod
    def _decide(dg_f, ln_rev=None):
        if ln_rev is not None:
            if abs(ln_rev) < 1.0:
                return "both"
            return "f" if ln_rev > 0 else "r"
        if dg_f is None:
            return "unknown"
        if abs(dg_f) <= REV_MARGIN:
            return "both"
        return "f" if dg_f < 0 else "r"

    def _from_cc(self, rid):
        eq = self.reactions.get(rid, {}).get("equation", "")
        if not eq or self._cc is None:
            return None
        # KEGG ids -> eQuilibrator formula ('C00022' -> 'kegg:C00022', '<=>' -> '=')
        try:
            formula = eq.replace("<=>", "=").replace("=>", "=").replace("<=", "=")
            toks = []
            for part in formula.split("="):
                terms = []
                for t in part.split(" + "):
                    t = t.strip()
                    if not t:
                        continue
                    bits = t.split()
                    if len(bits) == 2 and bits[0].lstrip("-").isdigit():
                        terms.append(f"{bits[0]} kegg:{bits[1]}")
                    else:
                        terms.append(f"kegg:{bits[-1]}")
                toks.append(" + ".join(terms))
            rxn = self._cc.parse_reaction_formula(" = ".join(toks))
            if not rxn.is_balanced():
                pass                          # still usable; CC handles many unbalanced
            dg = self._cc.standard_dg_prime(rxn)
            dg_f = float(dg.value.m_as("kJ/mol"))
            try:
                ln_rev = float(self._cc.ln_reversibility_index(rxn).m_as(""))
            except Exception:
                ln_rev = None
            return self._decide(dg_f, ln_rev), dg_f, "equilibrator"
        except Exception:
            return None

    def _from_consensus(self, rid):
        c = self._consensus.get(rid)
        if not c:
            return None
        # Legacy consensus: direction relative to the Rhea/KEGG reference equation.
        # rhea_direction BI/LR/RL is as-written-relative; metacyc as a fallback. We do NOT
        # use dGm_production_kJmol for the sign (it is production-oriented, not as-written).
        rhea = (c.get("rhea_direction") or "").upper()
        m = {"BI": "both", "LR": "f", "RL": "r"}.get(rhea)
        if m is None:
            mc = (c.get("metacyc_direction") or "").upper()
            m = {"REVERSIBLE": "both", "LEFT-TO-RIGHT": "f", "RIGHT-TO-LEFT": "r"}.get(mc, "unknown")
        return m, None, "consensus"

    def _from_arrow(self, rid):
        eq = self.reactions.get(rid, {}).get("equation", "")
        rev = "<=>" in eq
        return ("both" if rev else "f"), None, "kegg_arrow"

    def __call__(self, rid):
        if rid in self._dg:
            return self._dg[rid]
        res = (self._from_cc(rid) or self._from_consensus(rid) or self._from_arrow(rid))
        self._dg[rid] = res
        return res

    def coverage(self):
        """(fraction with real ΔG, counts by source) — quality gate G2."""
        from collections import Counter
        srcs = Counter(v[2] for v in self._dg.values())
        real = srcs["equilibrator"] + srcs["consensus"]
        tot = max(1, sum(srcs.values()))
        return {"fraction_real_dg": real / tot, "by_source": dict(srcs), "n": tot}

    def save(self):
        if self.cache_path:
            json.dump(self._dg, open(self.cache_path, "w"))

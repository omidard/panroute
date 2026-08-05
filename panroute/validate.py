#!/usr/bin/env python3
"""Validation harness (fix #5): benchmark the T2 'genome encodes a route' gate against a
curated truth table of known producers / non-producers. Reports precision, recall, and the
specific misses so a run that disagrees with the literature is visible (quality gate G4).

Matching is by species-name substring (KEGG species names vary), reported transparently.
"""
from __future__ import annotations
import csv, os


def load_truth(path, end_cid):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["end_cid"] == end_cid:
                rows.append(r)
    return rows


def _match(truth_species, t2_species_set):
    ts = truth_species.lower()
    for sp in t2_species_set:
        s = sp.lower()
        if s in ts or ts in s or (ts.split()[0] == s.split()[0] and ts.split()[0] not in
                                  ("clostridium", "bacillus", "escherichia", "pseudomonas")):
            return True
    # genus+species token overlap
    tt = set(ts.split()[:2])
    for sp in t2_species_set:
        if tt <= set(sp.lower().split()):
            return True
    return False


def validate(end_cid, t2_species_map, truth_path):
    truth = load_truth(truth_path, end_cid)
    if not truth:
        return {"n_truth": 0, "note": "no truth entries for this product"}
    t2 = set(t2_species_map.keys())
    tp = fp = tn = fn = 0
    misses = {"false_negative": [], "false_positive": []}
    for r in truth:
        hit = _match(r["species"], t2)
        if r["label"] == "producer":
            if hit:
                tp += 1
            else:
                fn += 1; misses["false_negative"].append(r["species"])
        else:  # non_producer
            if hit:
                fp += 1; misses["false_positive"].append(r["species"])
            else:
                tn += 1
    prec = tp / (tp + fp) if (tp + fp) else None
    rec = tp / (tp + fn) if (tp + fn) else None
    return {"n_truth": len(truth), "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "precision": prec, "recall": rec, "misses": misses}

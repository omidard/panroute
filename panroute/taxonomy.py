#!/usr/bin/env python3
"""KEGG organism taxonomy: code -> (domain, lineage, name) from the genome list and the
br08610 organism brite. Gram is a HEURISTIC phylum mapping and is labelled as such in all
outputs (a known caveat, not ground truth)."""
from __future__ import annotations
import os, re

PROK = {"Bacteria", "Archaea"}


def load_taxonomy(data_dir: str):
    org_name = {}
    gen = os.path.join(data_dir, "kegg_genome.tsv")
    for line in open(gen):
        p = line.rstrip("\n").split("\t")
        if len(p) < 2 or ";" not in p[1]:
            continue
        code, name = p[1].split(";", 1)
        org_name[code.strip()] = name.strip()
    valid = set(org_name)
    org_domain, org_lineage = {}, {}
    stack, domain = {}, None
    for line in open(os.path.join(data_dir, "br08610.keg")):
        if not line or line[0] not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
            continue
        lvl, rest = line[0], line[1:].rstrip("\n")
        m = re.match(r"\s*([a-z][a-z0-9]{1,4})  (\S.*)$", rest)
        if m and m.group(1) in valid:
            c = m.group(1)
            org_domain[c] = domain
            org_lineage[c] = [stack[k] for k in "BCDEFGHIJKLMNOPQRS"
                              if k in stack and ord(k) < ord(lvl)]
        else:
            nm = rest.strip()
            if lvl == "A":
                domain = nm
            stack[lvl] = nm
            for k in [c for c in stack if ord(c) > ord(lvl)]:
                del stack[k]
    return org_name, org_domain, org_lineage


def gram(lineage, domain):
    s = " ".join(lineage)
    if domain == "Archaea":
        return "Arch"
    if any(p in s for p in ["Mycoplasm", "Mollicutes", "Tenericutes"]):
        return "Other"
    if any(p in s for p in ["Chloroflex", "Deinococc", "Thermus"]):
        return "Gneg"
    if any(p in s for p in ["Bacillati", "Bacillota", "Firmicutes", "Actinomycet", "Actinobacteri"]):
        return "Gpos"
    if any(p in s for p in ["Pseudomonadati", "Pseudomonadota", "Proteobacteria", "Bacteroid",
                            "Campylobacter", "Cyanobacteri", "Spirochaet", "Chlorobi", "Thermodesulfo",
                            "Verrucomicrobia", "Planctomycet", "Aquific", "Deferribacter", "Fusobacteri",
                            "Chlamydi", "Nitrospir", "Acidobacteri", "Thermotog", "Synergist",
                            "Gemmatimonad", "Bdellovibrio", "Myxococc"]):
        return "Gneg"
    return "Other"


def species(name):
    t = name.replace("'", "").split()
    if not t:
        return name
    if t[0] == "Candidatus" and len(t) >= 3:
        return " ".join(t[:3])
    if len(t) >= 2:
        return " ".join(t[:2])
    return t[0]

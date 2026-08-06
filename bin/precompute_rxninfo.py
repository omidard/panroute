#!/usr/bin/env python3
"""Precompute per-reaction cross-database info for the route report pages:
KEGG <-> BiGG <-> Rhea <-> MetaCyc <-> SEED xrefs (via MetaNetX reac_xref) + directionality
from each source (KEGG arrow, Rhea LR/RL/BI, MetaCyc) + OUR eQuilibrator ΔG.

Outputs docs/data/rxninfo.json = {rid: {ec, eq, mnx, bigg[], rhea, rhea_dir, seed[],
metacyc, metacyc_dir, biocyc[], kegg_dir, our_dg, our_dir, our_src}}.
Only reactions present in the bundled network are emitted (routes can only use those).
Fast (no network calls). Sources are KEGG-independent open DBs where possible."""
import json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BC = "/data/bioconversion"
OUT = os.path.join(ROOT, "docs", "data")

net = json.load(open(os.path.join(OUT, "network.json")))
rxn_ids = set(net["rxn"].keys())
sys.stderr.write(f"[rxninfo] {len(rxn_ids)} network reactions\n")

# --- MetaNetX reac_xref: build MNXR -> {db: ids} and kegg R -> MNXR (+ equation) ---
DBMAP = {"bigg.reaction": "bigg", "biggR": "bigg", "rhea": "rhea", "rheaR": "rhea",
         "seed.reaction": "seed", "seedR": "seed", "metacyc.reaction": "metacyc",
         "metacycR": "metacyc", "biocyc": "biocyc", "sabiork.reaction": "sabiork"}
mnx2ids = {}           # MNXR -> {db: set(id)}
kegg2mnx = {}          # R -> MNXR
kegg_eq = {}           # R -> equation string
with open(f"{BC}/thermo/reac_xref.tsv") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 2:
            continue
        xref, mnx = p[0], p[1]
        if ":" not in xref:
            continue
        pre, xid = xref.split(":", 1)
        if pre in ("kegg.reaction", "keggR"):
            if xid.startswith("R"):
                kegg2mnx[xid] = mnx
                if len(p) >= 3 and "||" in p[2]:
                    parts = p[2].split("||")
                    if len(parts) >= 2 and "<=>" in parts[1] or (len(parts) >= 2 and "=" in parts[1]):
                        kegg_eq[xid] = parts[1]
        db = DBMAP.get(pre)
        if db:
            mnx2ids.setdefault(mnx, {}).setdefault(db, set()).add(xid)
sys.stderr.write(f"[rxninfo] MNX index: {len(mnx2ids)} MNXR, {len(kegg2mnx)} KEGG mapped\n")

# --- Rhea direction for a KEGG reaction ---
rhea_dir = {}          # R -> (rhea_id, direction)
with open(f"{BC}/thermo/rhea2kegg.tsv") as fh:
    next(fh)
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if len(c) >= 4 and c[3].startswith("R"):
            rhea_dir.setdefault(c[3], (c[2], c[1]))   # (master_id, BI/LR/RL)

# --- MetaCyc direction (curated subset) ---
metacyc = {}
try:
    for r in json.load(open(f"{BC}/thermo/metacyc_directions.json")):
        pass
except Exception:
    pass
mc = json.load(open(f"{BC}/thermo/metacyc_directions.json"))
if isinstance(mc, list):
    metacyc = {d["reaction"]: d for d in mc if d.get("reaction")}
elif isinstance(mc, dict):
    metacyc = mc

# --- our eQuilibrator ΔG (numeric, real) ---
our = {}
dgc = json.load(open(os.path.join(ROOT, "cache", "thermo_dg.json")))
for rid, v in dgc.items():
    direction, dg, src = v
    if src in ("equilibrator", "consensus") and isinstance(dg, (int, float)):
        our[rid] = {"dg": round(dg, 1), "dir": direction, "src": src}
    elif src == "equilibrator":
        our[rid] = {"dg": None, "dir": direction, "src": src}

RHEA_D = {"BI": "reversible", "LR": "left→right", "RL": "right→left"}
info = {}
for rid in rxn_ids:
    mnx = kegg2mnx.get(rid)
    ids = mnx2ids.get(mnx, {}) if mnx else {}
    rd = rhea_dir.get(rid)
    md = metacyc.get(rid)
    o = our.get(rid)
    info[rid] = {
        "ec": net["rxn"][rid].get("e", ""),
        "eq": kegg_eq.get(rid, ""),
        "kegg_dir": "reversible" if net["rxn"][rid].get("d") == "b" else "left→right",
        "mnx": mnx,
        "bigg": sorted(ids.get("bigg", []))[:6],
        "seed": sorted(ids.get("seed", []))[:4],
        "biocyc": sorted(ids.get("biocyc", []))[:4],
        "rhea": rd[0] if rd else None,
        "rhea_dir": RHEA_D.get(rd[1]) if rd else None,
        "metacyc": md.get("metacyc_id") if md else None,
        "metacyc_dir": (md.get("direction") or "").lower().replace("physiol-", "") if md else None,
        "our_dg": o["dg"] if o else None,
        "our_dir": o["dir"] if o else None,
    }
json.dump(info, open(f"{OUT}/rxninfo.json", "w"))
n_bigg = sum(1 for v in info.values() if v["bigg"])
n_rhea = sum(1 for v in info.values() if v["rhea"])
n_dg = sum(1 for v in info.values() if v["our_dg"] is not None)
sys.stderr.write(f"[rxninfo] wrote {len(info)} reactions -> rxninfo.json "
                 f"({os.path.getsize(f'{OUT}/rxninfo.json')//1024} KB) | bigg:{n_bigg} rhea:{n_rhea} ourΔG:{n_dg}\n")

# thermo.json for the engine (numeric ΔG only)
thermo = {rid: o["dg"] for rid, o in our.items() if isinstance(o.get("dg"), (int, float))}
json.dump(thermo, open(f"{OUT}/thermo.json", "w"))
sys.stderr.write(f"[rxninfo] thermo.json: {len(thermo)} reactions with numeric ΔG\n")

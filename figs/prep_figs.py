#!/usr/bin/env python3
"""Prepare figure data for a PanRoute query: (1) funnel/overview tables, (2) the route
DAG (nodes + enzyme-labelled edges with biological gene symbols). Writes TSVs to figs/data/."""
import sys, os, json, csv, re
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from panroute.keggfetch import KeggClient

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "cache")
START, END = sys.argv[1], sys.argv[2]
RESDIR = sys.argv[3]                      # e.g. results/butanediol
OUT = os.path.join(ROOT, "figs", "data"); os.makedirs(OUT, exist_ok=True)

cl = KeggClient(CACHE)
report = json.load(open(f"{RESDIR}/report_{START}_{END}.json"))
routes = json.load(open(f"{RESDIR}/routes_{START}_{END}.json"))["routes"]

# ---- metabolite stereoisomer collapse (one conceptual node per metabolite) ----
CANON = {
    "C00900": ("C_ala", "2-acetolactate"), "C06010": ("C_ala", "2-acetolactate"),
    "C00466": ("C_ain", "acetoin"),        "C00810": ("C_ain", "acetoin"),
    "C03044": ("C_bdo", "2,3-butanediol"), "C03046": ("C_bdo", "2,3-butanediol"),
    "C20657": ("C_bdo", "2,3-butanediol"),
}
def canon(cid, name):
    if cid in CANON:
        return CANON[cid]
    return cid, name.split(";")[0]

# ---- EC -> biological enzyme gene symbol (house standard: bio names, not DB ids) ----
EC2SYM = {
    "2.2.1.6": "alsS", "4.1.1.5": "alsD", "1.1.1.4": "budC", "1.1.1.76": "budC",
    "1.1.1.303": "budC", "1.1.1.304": "budC", "2.3.1.190": "acetoin synth.",
    "4.1.1.1": "pdc", "1.2.7.1": "por", "1.2.8.1": "por", "1.2.1.10": "mhpF/adhE",
    "2.8.3.8": "ctf", "2.3.1.8": "pta", "1.2.5.1": "poxB", "1.2.3.3": "poxB",
    "4.1.3.25": "citE", "4.1.3.46": "mcl", "2.3.3.9": "aceB", "1.1.1.40": "maeB",
    "4.1.3.39": "hoa", "4.1.2.36": "hoa", "1.1.3.2": "lox", "1.2.4.1": "pdhA",
    "4.1.1.5 ": "alsD",
}
def enzyme_label(ec, kos):
    if ec in EC2SYM:
        return EC2SYM[ec]
    # fall back to a cleaned KO symbol
    if kos:
        rec = cl.get_entries([f"ko:{kos[0]}"]).get(f"ko:{kos[0]}", "")
        for line in rec.splitlines():
            if line.startswith("NAME"):
                nm = line[12:].strip()
                head = nm.split(";")[0]
                m = re.match(r"^([a-z][a-zA-Z0-9]{1,5})(,|$)", head)
                if m:
                    return m.group(1)
    return ec or "?"

# ---- 1. funnel table ----
t = report["tiers"]
def sp(k): return t[k]["species"] if isinstance(t.get(k), dict) and "species" in t[k] else None
with open(f"{OUT}/{END}_funnel.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["tier", "label", "species"])
    w.writerow(["T0", "has terminal enzyme", sp("T0_terminal_enzyme_potential")])
    w.writerow(["T2", "encodes full route", sp("T2_genome_encodes_full_route")])
    w.writerow(["T3", "+ acetate uptake", t["T3_plus_feedstock_uptake"]["species"]])

# ---- 2. gram composition of T2 ----
g = t["T2_genome_encodes_full_route"]
with open(f"{OUT}/{END}_gram.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["group", "species"])
    for k, lab in [("Gpos", "Gram-positive"), ("Gneg", "Gram-negative"), ("Arch", "Archaea"), ("Other", "Other")]:
        if g.get(k, 0): w.writerow([lab, g[k]])

# ---- 3. meta (validation, counts, feedstock, overflow excluded) ----
val = t.get("T4_validation", {})
meta = {"start": START, "end": END,
        "n_routes": report["n_routes"], "shortest": report["shortest_route_len"],
        "T0": sp("T0_terminal_enzyme_potential"), "T2": sp("T2_genome_encodes_full_route"),
        "T3": t["T3_plus_feedstock_uptake"]["species"],
        "overflow_excluded": t["T3_plus_feedstock_uptake"]["species_overflow_only_excluded"],
        "precision": val.get("precision"), "recall": val.get("recall"),
        "kegg_release": report.get("kegg_release")}
json.dump(meta, open(f"{OUT}/{END}_meta.json", "w"), indent=2)

# ---- 4. notable species (validated producers first, then a spread) ----
truth = {}
tp = os.path.join(ROOT, "assets", "validation_truth.tsv")
for r in csv.DictReader(open(tp), delimiter="\t"):
    if r["end_cid"] == END:
        truth[r["species"].lower()] = r["label"]
sp_rows = list(csv.DictReader(open(f"{RESDIR}/species_{START}_{END}.csv")))
notable = []
for row in sp_rows:
    lab = None
    for ts, l in truth.items():
        if ts in row["species"].lower() or row["species"].lower() in ts:
            lab = l; break
    if lab == "producer":
        notable.append((row["species"], row["domain"], row["gram"], "validated"))
seen = {n[0] for n in notable}
for row in sp_rows:
    if len(notable) >= 12: break
    if row["species"] not in seen:
        notable.append((row["species"], row["domain"], row["gram"], "")); seen.add(row["species"])
with open(f"{OUT}/{END}_notable.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["species", "domain", "gram", "flag"])
    for n in notable: w.writerow(n)

# ---- 5. route DAG: nodes + enzyme-labelled edges (stereoisomers collapsed) ----
START_C, _ = canon(START, "")
END_C, _ = canon(END, "")
nodes = {}
edges = {}
for rt in routes:
    L = rt["length"]
    for p in rt["path"]:
        cid, nm = canon(p["cid"], p["name"])
        nodes.setdefault(cid, {"cid": cid, "name": nm, "role": "intermediate", "min_len": L})
        nodes[cid]["min_len"] = min(nodes[cid]["min_len"], L)
    for st in rt["steps"]:
        fc, _ = canon(st["from"], ""); tc, _ = canon(st["to"], "")
        if fc == tc:                       # self-edge from stereoisomer collapse -> skip
            continue
        key = (fc, tc)
        rx = st["reactions"][0]
        ec = rx["ec"][0] if rx["ec"] else ""
        lab = enzyme_label(ec, rx["kos"])
        e = edges.setdefault(key, {"from": fc, "to": tc, "reaction": rx["rid"], "ec": ec,
                                   "enzymes": set(), "on_shortest": False, "n_reactions": 0})
        e["enzymes"].add(lab)
        e["n_reactions"] = len({r["rid"] for r in st["reactions"]})
        if L == report["shortest_route_len"]:
            e["on_shortest"] = True
nodes[START_C]["role"] = "start"; nodes[END_C]["role"] = "end"
with open(f"{OUT}/{END}_nodes.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["cid", "name", "role", "min_len"])
    for n in nodes.values(): w.writerow([n["cid"], n["name"], n["role"], n["min_len"]])
with open(f"{OUT}/{END}_edges.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["from", "to", "reaction", "ec", "enzymes", "n_reactions", "on_shortest"])
    for e in edges.values():
        w.writerow([e["from"], e["to"], e["reaction"], e["ec"],
                    "/".join(sorted(e["enzymes"])[:2]) or "?", e["n_reactions"],
                    int(e["on_shortest"])])
print(f"prep done: {len(nodes)} nodes, {len(edges)} edges, {len(notable)} notable species")
print("nodes:", ", ".join(sorted(n["name"][:16] for n in nodes.values())))

#!/usr/bin/env python3
"""Per-reaction KO requirement grouped by EC component, so genuine multi-subunit complexes
are gated with AND across catalytic components (OR within = isozymes) instead of the naive OR
that over-credits any genome carrying one shared subunit.

KEGG reaction ORTHOLOGY lists each KO with its [EC:...]; a complex (e.g. pyruvate
dehydrogenase) spans several ECs. BUT different EC numbers do NOT imply a complex — most
often they are ISOZYMES / alternative enzymes (hexokinase 2.7.1.1 vs glucokinase 2.7.1.2 on
R00299), each catalysing the WHOLE reaction. Blindly ANDing them makes central metabolism
(glycolysis, TCA) gate to ~0 genomes (a cell with only glucokinase fails hexose kinase).

Discriminator (DATA-DRIVEN co-occurrence): true complex subunits are co-inherited, so the
fraction of genomes carrying all EC groups (AND) over genomes carrying any (OR) is HIGH;
isozymes are substitutable, so that ratio is LOW. Calibration: PDH R00209 = 0.88; hexokinase
0.01, PGI 0.00, phosphoglycerate-mutase 0.00, GAPDH 0.10. We keep the AND split ONLY when
AND/OR >= COOC_MIN, erring toward OR (mild over-credit) rather than under-crediting a real
pathway to zero.

Output docs/data/rxnko.json = {rid: [[EC-group-1 KOs], [EC-group-2 KOs], ...]} (AND of groups)."""
import json, os, re, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from panroute.keggfetch import KeggClient
OUT = os.path.join(ROOT, "docs", "data")
KO = os.path.join(OUT, "ko")
net = json.load(open(f"{OUT}/network.json"))
cl = KeggClient(os.path.join(ROOT, "cache"), offline=True)
COOC_MIN = 0.5    # AND/OR ratio above which EC groups are treated as an obligate complex

_ko = {}
def orgs(k):
    if k not in _ko:
        p = os.path.join(KO, f"{k}.json")
        _ko[k] = set(json.load(open(p))["orgs"]) if os.path.exists(p) else set()
    return _ko[k]

rxnko = {}
n_multi = n_iso = 0
for rid in net["rxn"]:
    if len(net["rxn"][rid].get("k", [])) < 2:
        continue
    rec = cl.get_entries([f"rn:{rid}"]).get(f"rn:{rid}", "")
    groups = {}          # ec -> [kos]
    field = None
    for line in rec.splitlines():
        key = line[:12].strip()
        if key:
            field = key
        if field != "ORTHOLOGY":
            continue
        m = re.search(r"(K\d{5})", line)
        if not m:
            continue
        ko = m.group(1)
        ecm = re.search(r"\[EC:([^\]]+)\]", line)
        ec = ecm.group(1).split()[0] if ecm else "noec"
        groups.setdefault(ec, [])
        if ko not in groups[ec]:
            groups[ec].append(ko)
    ecgroups = [sorted(v) for k, v in groups.items() if k != "noec"]
    if len(ecgroups) < 2:
        continue
    # co-occurrence gate: keep AND only if the groups are co-inherited (a real complex)
    per = [set().union(*[orgs(k) for k in g]) if g else set() for g in ecgroups]
    OR = set().union(*per) if per else set()
    AND = set.intersection(*per) if per else set()
    ratio = (len(AND) / len(OR)) if OR else 0
    if ratio >= COOC_MIN:
        rxnko[rid] = ecgroups
        n_multi += 1
    else:
        n_iso += 1      # substitutable isozymes -> omit, engine ORs the KOs
json.dump(rxnko, open(f"{OUT}/rxnko.json", "w"))
print(f"[rxnko] {n_multi} obligate-complex reactions kept as AND (ratio>={COOC_MIN}); "
      f"{n_iso} multi-EC reactions were ISOZYMES -> collapsed to OR. -> rxnko.json "
      f"({os.path.getsize(f'{OUT}/rxnko.json')//1024} KB)")
print("  PDH R00209:", rxnko.get("R00209"), "| hexokinase R00299:", rxnko.get("R00299"),
      "| PGI R00771:", rxnko.get("R00771"))

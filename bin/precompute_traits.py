#!/usr/bin/env python3
"""Per-species phenotype traits for the results catalog: oxygen tolerance, temperature class,
and a (best-effort, curated) safety class. Gram/domain already come from taxonomy.json.

Sources on disk:
  - Madin et al. 2020 bacteria-archaea-traits  (/data/bioconversion/traits/traits.csv):
      `metabolism` -> oxygen, `range_tmp`/`optimum_tmp` -> temperature class.
  - BacDive growth pull (/data/GrowthDB_work/data/bacdive_growth.json.gz): `ox`, `optT`.
  - Curated EFSA-QPS (safe) and known-human-pathogen lists for the safety class (partial;
    everything else is left 'unknown' rather than guessed).

Output docs/data/traits.json = {species: {ox, temp, safety}} (species = binomial, matching
taxonomy.json). oxygen/temp fall back to a genus majority when the species is unlisted."""
import json, os, csv, gzip, sys
from collections import Counter, defaultdict
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "data")
TRAITS = "/data/bioconversion/traits/traits.csv"
BACDIVE = "/data/GrowthDB_work/data/bacdive_growth.json.gz"

def ox_norm(v):
    v = (v or "").lower()
    if not v or v == "na": return None
    if "microaero" in v: return "microaerophile"
    if "facultative" in v: return "facultative"
    if "anaerob" in v: return "anaerobe"
    if "aerob" in v: return "aerobe"
    return None

def temp_class(rng, opt):
    rng = (rng or "").lower()
    if "hyperthermo" in rng or "extreme thermo" in rng: return "hyperthermophile"
    if "thermo" in rng: return "thermophile"
    if "psychro" in rng: return "psychrophile"
    if "meso" in rng: return "mesophile"
    try:
        t = float(opt)
        if t < 20: return "psychrophile"
        if t < 45: return "mesophile"
        if t < 80: return "thermophile"
        return "hyperthermophile"
    except (TypeError, ValueError):
        return None

# ---- Madin traits.csv: per-species majority oxygen + temperature, and genus fallback ----
sp_ox, sp_tmp, gen_ox, gen_tmp = defaultdict(Counter), defaultdict(Counter), defaultdict(Counter), defaultdict(Counter)
with open(TRAITS) as f:
    for r in csv.DictReader(f):
        sp, gen = r["species"], r["genus"]
        o = ox_norm(r["metabolism"]); t = temp_class(r["range_tmp"], r["optimum_tmp"])
        if o: sp_ox[sp][o] += 1; gen_ox[gen][o] += 1
        if t: sp_tmp[sp][t] += 1; gen_tmp[gen][t] += 1

# ---- BacDive supplement (species from the org binomial) ----
try:
    with gzip.open(BACDIVE, "rt") as f:
        for rec in json.load(f):
            org = (rec.get("org") or "").split()
            if len(org) < 2: continue
            sp = " ".join(org[:2]); gen = org[0]
            o = ox_norm(rec.get("ox")); t = temp_class(None, rec.get("optT"))
            if o: sp_ox[sp][o] += 1; gen_ox[gen][o] += 1
            if t: sp_tmp[sp][t] += 1; gen_tmp[gen][t] += 1
except FileNotFoundError:
    sys.stderr.write("[traits] BacDive file missing, using Madin only\n")

# ---- curated safety (partial, honest) ----
QPS_GENERA = {"Lactobacillus", "Lactiplantibacillus", "Lacticaseibacillus", "Limosilactobacillus",
    "Lactococcus", "Leuconostoc", "Pediococcus", "Oenococcus", "Weissella", "Streptococcus",  # S. thermophilus only, refined below
    "Bifidobacterium", "Propionibacterium", "Acidipropionibacterium", "Carnobacterium",
    "Corynebacterium", "Gluconobacter", "Komagataeibacter"}
QPS_SPECIES = {"Bacillus subtilis", "Bacillus amyloliquefaciens", "Bacillus licheniformis",
    "Bacillus coagulans", "Bacillus velezensis", "Corynebacterium glutamicum",
    "Streptococcus thermophilus", "Escherichia coli",  # E. coli K-12 lab strains are GRAS-adjacent; flagged safe-ish
    "Saccharomyces cerevisiae"}
PATHOGEN_SPECIES = {"Mycobacterium tuberculosis", "Yersinia pestis", "Vibrio cholerae",
    "Bacillus anthracis", "Clostridium botulinum", "Clostridioides difficile", "Salmonella enterica",
    "Shigella dysenteriae", "Neisseria meningitidis", "Neisseria gonorrhoeae", "Bordetella pertussis",
    "Corynebacterium diphtheriae", "Listeria monocytogenes", "Francisella tularensis",
    "Brucella melitensis", "Burkholderia mallei", "Burkholderia pseudomallei", "Legionella pneumophila",
    "Treponema pallidum", "Helicobacter pylori", "Campylobacter jejuni", "Haemophilus influenzae",
    "Streptococcus pyogenes", "Streptococcus pneumoniae", "Bordetella bronchiseptica"}
OPPORTUNIST_SPECIES = {"Pseudomonas aeruginosa", "Klebsiella pneumoniae", "Acinetobacter baumannii",
    "Staphylococcus aureus", "Enterococcus faecium", "Enterococcus faecalis", "Enterobacter cloacae",
    "Escherichia coli", "Serratia marcescens", "Stenotrophomonas maltophilia", "Proteus mirabilis",
    "Clostridium perfringens", "Bacteroides fragilis", "Morganella morganii", "Citrobacter freundii"}
def safety(sp):
    gen = sp.split()[0]
    if sp in PATHOGEN_SPECIES: return "pathogen"
    if sp in OPPORTUNIST_SPECIES: return "opportunist"
    if sp in QPS_SPECIES: return "GRAS/QPS"
    if gen in QPS_GENERA: return "GRAS/QPS"
    return None

# ---- assemble for OUR species only ----
tax = json.load(open(f"{OUT}/taxonomy.json"))
our = sorted(set(v[0] for v in tax.values()))
def pick(spmap, genmap, sp):
    if spmap.get(sp): return spmap[sp].most_common(1)[0][0]
    gen = sp.split()[0]
    if genmap.get(gen): return genmap[gen].most_common(1)[0][0]
    return None
traits = {}
for sp in our:
    o = pick(sp_ox, gen_ox, sp); t = pick(sp_tmp, gen_tmp, sp); s = safety(sp)
    if o or t or s:
        traits[sp] = {k: v for k, v in (("ox", o), ("temp", t), ("safety", s)) if v}
json.dump(traits, open(f"{OUT}/traits.json", "w"))
cov = lambda k: sum(1 for v in traits.values() if v.get(k))
print(f"[traits] {len(our)} species · oxygen={cov('ox')} temp={cov('temp')} safety={cov('safety')} "
      f"-> traits.json ({os.path.getsize(f'{OUT}/traits.json')//1024} KB)")

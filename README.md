# PanRoute

**Retrosynthetic route search across all KEGG genomes.** Given a *start* metabolite
(feedstock) and an *end* metabolite (target), PanRoute finds every thermodynamically
feasible native route through the KEGG reaction network, then reports which prokaryotic
genomes in KEGG *encode* a complete route — with an honest, tiered confidence funnel.

It replaces the hand-written-route survey in `/data/bioconversion`. See
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for the full design and the five flaws it
fixes.

---

## Why it's different from a naive KEGG survey

- **Real retrosynthetic search** (not hand-written pathways): reverse traversal from the
  product back to the feedstock over the whole KEGG reaction network.
- **Carbon-skeleton graph**: edges come only from KEGG **RCLASS** atom-conserved pairs,
  with currency metabolites (ATP/NAD(P)(H)/CoA/CO₂/…) excluded — so the search cannot
  shortcut through the cofactor pool. C1 metabolism (methanol/formate) is handled.
- **Thermodynamically gated** every edge (eQuilibrator ΔrG′°; falls back to a Rhea/MetaCyc
  consensus, then KEGG arrows — coverage is reported, never hidden).
- **Direction-aware feedstock gating**: `ackA-pta`-only genomes are flagged
  *acetate-overflow*, not counted as uptake (the Parageobacillus caveat).
- **Honest tiered funnel**: terminal-enzyme potential → route exists → genome encodes route
  → + feedstock uptake → + validated. Counts are labelled "genome *encodes* a route,"
  never "can produce."
- **Validated** against a curated producer/non-producer truth table.

## Install

On this machine the **base conda env already has the dependencies** (`networkx` +
`equilibrator_api`), so nothing to install — `-profile standard` uses base Python directly.

For a **fresh machine**, build the env from the spec and use `-profile conda`:
```bash
conda env create -f conf/environment.yml    # networkx + equilibrator-api
```
`equilibrator-api` is optional; without it PanRoute falls back to the Rhea/MetaCyc
consensus + KEGG arrows and reports reduced thermodynamic coverage (never fails).

## Run — Nextflow (recommended, self-contained, resumable)

```bash
# a batch of bioconversions
nextflow run main.nf --queries conf/queries.csv -profile standard

# a single query
nextflow run main.nf --start C00024 --end C00207 --feedstock C00033 -profile standard

# fast smoke test (cached, no thermo)
nextflow run main.nf -profile test
```
KEGG fetches are cached on disk under `cache/` and reused across runs and queries; use
`-resume` to skip completed queries. `-profile slurm` for the cluster.

## Run — CLI (one query, no Nextflow)

```bash
python -m panroute.cli --start C00024 --end C00207 --feedstock C00033 --out results/
```

Key options: `--expand-depth` (bounded network expansion depth, default 3),
`--max-len`/`--max-routes` (route enumeration bounds), `--no-thermo` (skip eQuilibrator),
`--reactions FILE` (use an explicit reaction set instead of expansion), `--offline`
(cache only).

## Inputs

Metabolites are **KEGG compound ids** (`C#####`). Look them up at
<https://rest.kegg.jp/find/compound/acetone>. A query is `name,start,end,feedstock`
(feedstock optional; enables direction-aware uptake gating).

## Outputs (per query, under `results/<name>/`)

| file | contents |
|------|----------|
| `report_<s>_<e>.json`   | the tiered funnel + thermo coverage + KEGG release + caveats |
| `routes_<s>_<e>.json`   | every enumerated route (compound path + per-step reactions/KOs/EC/ΔG) |
| `species_<s>_<e>.csv`   | T2 species that encode ≥1 full route (domain, Gram, #routes) |
| `per_genome_<s>_<e>.csv`| per-genome gating detail |
| `webapp_<s>_<e>.json`   | self-describing bundle for the web UI |

Aggregate: `results/panroute_summary.json` (funnel table across queries) and
`results/panroute_index.json` (webapp index).

## Interpreting the funnel (read this)

`T0` (terminal enzyme) overcounts — it only checks the last step. `T2` (encodes a full
route) is the honest headline. `T3` adds feedstock uptake. Example (acetone from acetate):
`T0 455 → T2 265 → T3 245` species. **This is genome potential, not proof of production**:
a KO hit means the enzyme is *encoded*, not expressed, regulated, or carrying flux. KEGG
genomes are culture-biased. Confirm with GEM/FBA and experiment.

## Tests / quality gates

```bash
python tests/smoke_bounded.py     # core: no currency leak, recovers the acetone route
python tests/smoke_acetone.py     # end-to-end on a live bounded subnetwork
```
Gates G1–G5 (currency-leak, thermo coverage, no-silent-drops, validation, reproducibility)
are described in `docs/ARCHITECTURE.md §7`.

## Web app

`webapp/` serves the precomputed `webapp_*.json` bundles (route DAG + species funnel).
See `webapp/README.md`. **Licensing:** public serving of KEGG-derived data is restricted;
the data layer is abstracted so KEGG can be swapped for MetaCyc/BiGG/ModelSEED for a public
deployment.

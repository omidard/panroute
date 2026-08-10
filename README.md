# PanRoute

**Retrosynthetic route search across all KEGG genomes.** Given a *start* metabolite
(feedstock) and an *end* metabolite (target), PanRoute finds every thermodynamically
feasible native route through the KEGG reaction network, then reports which prokaryotic
genomes in KEGG *encode* a complete route.



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

The shipped app is the **client-side JS engine** at <https://omidard.github.io/panroute/>
(`docs/engine.js`, node-testable via an injectable `loader`). It mirrors the Python engine
but is the canonical runtime: it runs *any* start→end pair live in the browser over the
bundled static KEGG data (`docs/data/`), with the same aliasing, multi-EC subunit-AND gating,
and signed-ΔG feasibility as the pipeline. The Python CLI/Nextflow path is for batch/offline
reproduction of the data bundle. **Licensing:** public serving of KEGG-derived data is
restricted; the data layer is abstracted so KEGG can be swapped for MetaCyc/BiGG/ModelSEED.

### Reproducing the data bundle

`docs/data/` is rebuilt by the `bin/precompute_*.py` scripts (each is idempotent and cached):
`precompute_clientdata` (network/taxonomy), `precompute_kogenomes` (per-KO genome lists),
`precompute_aliases` (compound-identity groups), `precompute_thermo` + `backfill_thermo`
(eQuilibrator ΔrG′°), `precompute_rxnsides` (L/R sides for signed ΔG), `precompute_rxnko`
(complex subunit groups), `precompute_rxninfo`/`precompute_smiles`. Every correction is a
script here, not an ad-hoc data patch, so a rebuild reproduces the corrected bundle.

### Known limitations (honest caveats)

- **Distant pairs are FBA territory.** Bounded retrosynthetic search (default ≤6 steps) will
  not find a native long pathway (e.g. glucose→pyruvate is ~10-step glycolysis). When no
  genome encodes a route *within the horizon*, the UI says so explicitly rather than implying
  the conversion is impossible or "engineered only."
- **~2,000 reactions are uncreditable to genomes.** KEGG annotates some reactions only with
  eukaryote-specific KOs (e.g. K00001/K00006 list only insects/human); bacteria use different
  KOs for the same EC. Where a reaction carries *only* such KOs, no prokaryote is credited.
  Route-critical central metabolism is unaffected (verified populated).
- **Same-EC obligate heterodimers may be over-credited.** Multi-EC complexes get subunit-AND
  gating (`rxnko.json`), but two subunits sharing one EC (e.g. AHAS large+regulatory, both
  EC 2.2.1.6) cannot be split by EC alone, so an OR over them can over-credit a genome that
  has only one subunit. Bounded (~25% on the few affected reactions).
- **L-/D-stereoisomer merges** (e.g. lactate) are intentional at the routing layer; strain
  stereospecificity is not resolved.
- **ΔG coverage ~92%** of routable reactions; `|ΔG|>200 kJ/mol` (unbalanced/generic-R-group
  equations) is flagged suspect (`thermo_suspect.json`) and treated as *unknown*, never gated.
  A genome KO hit means the enzyme is *encoded* — not expressed, regulated, or carrying flux.

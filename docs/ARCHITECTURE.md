# PanRoute — retrosynthetic route search across all KEGG genomes

**One line:** given a *start* metabolite (feedstock) and an *end* metabolite (target product),
find every thermodynamically-feasible native route through the KEGG reaction network, then
determine which prokaryotic genomes in KEGG encode a complete such route — with honest,
tiered confidence.

This replaces the hand-written-route survey in `/data/bioconversion` (see that project's
`REPORT.md`). It fixes the flaws catalogued there:

| # | Old flaw | Fix here |
|---|----------|----------|
| 1 | Routes hand-curated & sparse (acetone = 1 route) | **Real retro-search** over the whole KEGG reaction network |
| 2 | No thermodynamics on the full-pathway layer | **eQuilibrator ΔrG′° gates every edge**, inside the search |
| 3 | Feedstock gating = KO presence (ackA-pta counted as uptake) | **Direction-aware** feedstock module (overflow vs uptake) |
| 4 | Counts presented as "can produce" | **Tiered funnel**, relabelled "genome *encodes* a route" |
| 5 | No experimental grounding | **Validation harness** vs curated producers/non-producers |

Plus the deeper flaws (gene≠flux, KEGG sampling bias, taxonomy heuristics) are surfaced as
explicit caveats in every output, never hidden.

---

## 1. The core problem and why naive search fails

A KEGG reaction network is a bipartite graph (compounds ↔ reactions). A naive
compound→compound BFS treats **currency metabolites** (ATP, NAD(P)(H), CoA, H₂O, CO₂, Pi,
PPi, H⁺, NH₃, O₂, glutamate/2-oxoglutarate as amino donors, …) as ordinary nodes, so
*everything* becomes reachable from *everything* through the cofactor pool. Result: garbage
"routes" that share no carbon.

**Rigorous fix — carbon-skeleton graph.** An edge S→P is added **only** when S and P are an
**atom-conserved pair** that transfers a real carbon skeleton. We determine this three ways,
combined (belt-and-suspenders):

1. **KEGG RCLASS pairs.** Each reaction record lists `RCLASS  RCxxxxx  Ca_Cb` — the
   atom-mapped substrate↔product pairs. We only build edges from these pairs (never from
   arbitrary substrate×product combinations).
2. **Currency blocklist.** A curated set (`assets/currency_metabolites.tsv`) removes cofactor
   pairs (NAD⁺/NADH, CoA/acyl-CoA carrier side, ATP/ADP, …). CoA *carrier* pairs are dropped;
   the acyl carbon skeleton pair is kept.
3. **Carbon-conservation threshold.** Using compound formulas, an edge requires
   `shared_carbons ≥ MIN_C_SHARED` (default 2) and both endpoints carbon-containing. This
   catches currency pairs RCLASS/blocklist miss and any decarboxylation bookkeeping.

An edge carries: reaction id, direction(s) allowed (thermo), the KO requirement (isozyme OR /
subunit AND), EC, ΔrG′°, and provenance.

## 2. Thermodynamics (fix #2)

`panroute/thermo.py`. Primary engine = **eQuilibrator** (`equilibrator_api`,
component-contribution) keyed by RHEA id (from reaction `DBLINKS`) or by the KEGG equation.
For each reaction we compute ΔrG′° at pH 7.0, ionic strength 0.25 M, and a **reversibility
index** (Noor et al. 2012). An edge is traversable:

- **forward only** if ΔrG′° ≤ −RT·ln(Γ) margin (strongly favourable),
- **reverse only** if strongly unfavourable,
- **both** if within the reversible window.

Retro-search walks reactions *backwards* (product→substrate), so a route is feasible only if
each reaction can run in the *production* (forward-to-target) direction. **Fallback** when
eQuilibrator lacks a compound: the existing Rhea/MetaCyc/group-contribution consensus in
`/data/bioconversion/thermo/directionality_consensus.json`, then KEGG `<=>`/`=>` arrow as last
resort (flagged `direction_source=kegg_arrow`, lower confidence).

## 3. Retro-search (fix #1)

`panroute/retro.py`. Reverse traversal from `end` compound toward `start` compound on the
carbon-skeleton graph:

- **Reachability**: reverse BFS to confirm a path exists and get minimum length.
- **Route enumeration**: bounded DFS / Yen's k-shortest-paths returning up to `--max-routes`
  distinct reaction sequences with length ≤ `--max-len` (default 12). Each route = ordered list
  of (reaction, direction, KO-requirement).
- **Ranking**: routes scored by (length, total |ΔrG′°| headroom, number of currency-independent
  steps, pathway-coherence bonus if steps share a KEGG MODULE/map).
- **De-duplication**: routes equal under reaction-set are merged; isozyme choices collapse.

Output per (start,end): `routes.json` (the DAG of feasible routes + per-route metadata).

## 4. Genome gating (fix #4 tiering)

`panroute/genomes.py`. For each enumerated route, the **KO requirement** is a boolean over the
genome's KO set: `AND over steps( OR over isozymes( AND over subunits ) )`. A genome *encodes*
the route iff the boolean is satisfiable. We report, as a **funnel**:

- **T0 terminal-enzyme potential** — genome has the last step only (the old, inflated metric; kept for comparison).
- **T1 route exists** — a thermo-feasible network route start→end exists at all (organism-independent).
- **T2 genome encodes ≥1 full route** — the headline, honest count.
- **T3 + feedstock uptake** — T2 ∧ direction-aware feedstock utilisation (fix #3).
- **T4 + validated** — T3 ∧ passes the validation benchmark class (fix #5).

Every count is a **distinct species** (strain-collapsed) with domain/Gram, but taxonomy is
labelled *heuristic* and the KEGG-sampling caveat is printed on every artifact.

## 5. Feedstock direction-awareness (fix #3)

`panroute/feedstock.py`. A feedstock is *usable* only if the genome has an uptake/activation
route in the **consuming** direction:
- **Acetate**: ACS (`K01895`, one-step, irreversible uptake) → uptake-competent. ackA+pta
  (`K00925`+`K00625`) alone → flagged `acetate:overflow_capable` (runs both ways; not proof of
  growth on acetate). Configurable strictness.
- **Methanol**: require a *methanol dehydrogenase* AND an assimilation module (RuMP or serine or
  the WLP methyl branch), not mdh alone.
- General feedstocks: uptake defined in `assets/feedstock_rules.json`, extensible.

## 6. Validation (fix #5)

`panroute/validate.py` + `assets/validation_truth.tsv`: curated known **producers** and
**non-producers** per product (with DOIs). Reports precision/recall of the T2/T3 gate against
truth, per product. A run that scores poorly on validation is flagged in the report and (in CI)
fails the quality gate.

## 7. Quality gates (the "high quality" mandate)

Reviewer ≠ author (institute principle): the search author-code and the validator are separate
modules with separate truth inputs. Hard gates in `tests/` and in the Nextflow `VALIDATE`
process:

- **G1 network sanity** — no route may pass through a currency metabolite; unit-tested on the
  glycolysis/TCA backbone (known atom-mapped answers).
- **G2 thermo coverage** — ≥X% of edges have a real ΔG (not just KEGG arrow); report the %.
- **G3 no silent drops** — any KO/reaction/compound that failed to fetch is logged and counted;
  the run refuses to present a clean number while >Y% of required data is missing.
- **G4 validation** — precision/recall vs truth ≥ threshold per product.
- **G5 reproducibility** — every KEGG object cached with its fetch date + KEGG release; a run is
  reproducible from cache; counts are stamped with the KEGG release.

## 8. Pipeline (Nextflow DSL2)

`main.nf`, fully parametrized — runs *any* start→end pair:

```
FETCH_TARGETS      one-off: resolve start/end names→KEGG cpd ids
BUILD_NETWORK      fetch/refresh reactions+rclass+compounds (cached, resumable)  → network.pkl
THERMO_ANNOTATE    eQuilibrator ΔG per reaction                                  → thermo.parquet
RETRO_SEARCH       enumerate feasible routes start→end                           → routes.json
FETCH_GENOME_KOS   per-genome KO sets (cached; the heavy, resumable step)        → genome_kos.parquet
GATE_GENOMES       route boolean over every genome                              → per_genome.parquet
FEEDSTOCK_GATE     direction-aware feedstock utilisation                         → feedstock.parquet
VALIDATE           benchmark vs truth (hard gate)                                → validation.json
REPORT             tiered funnel + figure + species tables + webapp JSON         → results/
```

`nextflow.config` profiles: `standard` (local), `slurm` (compute), `test` (tiny stub set).
Params in `conf/params.yaml`. Resumable: `-resume` reuses cached KEGG objects.

## 9. Webapp (phase 2)

`webapp/` — user enters start + end metabolite → the precomputed route-graph + species table is
served (static JSON emitted by REPORT for common pairs; on-demand for novel pairs via a thin
API that runs RETRO_SEARCH on the cached network). Route DAG rendered with Cytoscape.js; species
table with the funnel + caveats. **Licensing note:** public serving of KEGG-derived content is
constrained; the data layer is abstracted (`panroute/backends/`) so KEGG can be swapped for
MetaCyc/BiGG/ModelSEED for a public deployment.

## 10. What this is and isn't (honesty, always printed)

**Is:** a genome-*potential* screen — "which sequenced prokaryotes *encode* a thermodynamically
feasible native route." A ranked hypothesis generator.
**Isn't:** proof of production. Gene presence ≠ expression ≠ flux ≠ titre. No regulation, no in
vivo concentrations, no rates. KEGG genomes are culture-biased (not a census). Downstream GEM/FBA
(the metabolic_atlas layer) and experiment are required to confirm.

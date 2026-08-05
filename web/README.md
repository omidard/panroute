# PanRoute — live web UI

A local, futuristic web app for running a bioconversion retro-search and watching it happen
**live on the real KEGG global metabolic map** (map01100), Google-Maps style.

> **Internal / private.** This app serves KEGG-derived data. Keep the repo private and run it
> only on internal machines (that is why licensing is safe). Do **not** deploy it publicly.

## Run

```bash
# from the repo root (base conda env has the deps: networkx + equilibrator_api)
python -m server.app
# open http://localhost:8000
```
Set `PANROUTE_PORT` to change the port. The app needs the KEGG genome/taxonomy data at
`/data/bioconversion/data/` (kegg_genome.tsv, br08610.keg) and the map asset
`assets/map01100/` (shipped in the repo).

## What you get

- **Type two metabolites** (product + start/feedstock) — autocomplete resolves them to KEGG
  compound ids.
- **Live KEGG-map trace** — the route is drawn from the **product back toward the feedstock**
  along the *real* reaction polylines KEGG drew (bending through the network). Specialised
  products that are off the core map appear as peripheral chips with dashed connectors.
- **Left panel, live** — organisms discovered as the genomes are gated, each with its Gram/
  domain, how many routes it can use, and a **thermodynamic-feasibility** badge (real
  eQuilibrator ΔG).
- **Results** — the honest genome funnel (terminal-enzyme → encodes-route → +feedstock),
  composition donut, and a searchable species list. **Everything is clickable**: a species
  opens its routes (metabolite chain + enzymes + reactions + ΔG); a funnel tier explains what
  it means.

## How the live map is true, not mocked

`panroute/mapviz.py` resolves each route step to real geometry: the reaction's own drawn
polyline if present; else an on-map reaction performing the same compound-pair transformation
(map01100 draws a curated subset — e.g. acetolactate synthase is R00226 on the map, not the
variant R04672); else a connector between the two compounds' real node positions; else a
peripheral chip. The blue highlight is animated along these real coordinates.

## Architecture

- **Backend** `server/app.py` (stdlib only): serves the frontend + streams true progress from
  `panroute/engine.py:run_query()` (a generator) over Server-Sent Events.
- **Frontend** `web/`: `map.js` (KEGG-map SVG + draw-on animation), `app.js` (SSE, panels,
  results, drawer), `style.css` (dark "night-map" theme).
- **Map asset** `assets/map01100/` built by `bin/extract_map.py` from the reference KGML.

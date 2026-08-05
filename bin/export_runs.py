#!/usr/bin/env python3
"""Export precomputed bioconversion runs as static JSON bundles for the public github.io
site (client-side replay — no backend). Each bundle is the full ordered event stream from
panroute.engine.run_query, so the frontend replays the exact live experience statically.

Usage: export_runs.py  (edits the RUNS list below), writes docs/runs/<slug>.json + index.json
"""
import sys, os, json, re
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from panroute.engine import run_query

OUT = os.path.join(ROOT, "docs", "runs")
os.makedirs(OUT, exist_ok=True)

# (slug, title, start_cid, end_cid, feedstock_cid, max_len)
RUNS = [
    ("butanediol-from-acetate", "2,3-Butanediol from acetate", "C00022", "C03044", "C00033", 4),
    ("acetone-from-acetate",    "Acetone from acetate",        "C00024", "C00207", "C00033", 4),
]


def slugify(s):
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


def export(slug, title, start, end, feed, max_len):
    events = []
    meta = {"slug": slug, "title": title}
    for ev, payload in run_query(start, end, feedstock=feed, max_len=max_len, thermo_routes=24):
        events.append({"event": ev, "data": payload})
        if ev == "endpoints":
            meta["start_name"] = payload["start"]["name"]
            meta["end_name"] = payload["end"]["name"]
        if ev == "done":
            meta["T2"] = payload.get("T2")
            meta["n_routes"] = payload.get("n_routes")
    json.dump({"meta": meta, "events": events}, open(os.path.join(OUT, f"{slug}.json"), "w"))
    return meta


def main():
    index = []
    for slug, title, start, end, feed, ml in RUNS:
        print(f"exporting {slug} ...", flush=True)
        m = export(slug, title, start, end, feed, ml)
        index.append(m)
        print(f"  {m.get('end_name')} <- {m.get('start_name')} | "
              f"{m.get('T2')} species | {m.get('n_routes')} routes")
    json.dump({"runs": index}, open(os.path.join(OUT, "index.json"), "w"), indent=2)
    print(f"wrote {len(index)} runs + index.json -> {OUT}")


if __name__ == "__main__":
    main()

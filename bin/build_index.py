#!/usr/bin/env python3
"""Build a webapp index of all available precomputed queries (webapp_*.json bundles)."""
import sys, json, os, glob

root = sys.argv[1] if len(sys.argv) > 1 else "."
index = []
for f in glob.glob(os.path.join(root, "**", "webapp_*.json"), recursive=True):
    try:
        d = json.load(open(f))
    except Exception:
        continue
    q = d.get("query", {})
    index.append({
        "start": q.get("start"), "end": q.get("end"),
        "start_name": q.get("start_name"), "end_name": q.get("end_name"),
        "n_routes": len(d.get("routes", [])),
        "file": os.path.relpath(f, root),
    })
print(json.dumps({"queries": index}, indent=2))

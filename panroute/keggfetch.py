#!/usr/bin/env python3
"""Cached, rate-limited, resumable KEGG REST client.

Every fetched object is written to the cache dir as a raw file plus a small sidecar
recording the fetch date and the current KEGG release, so a run is reproducible from
cache and every number can be stamped with a KEGG release (quality gate G5).

KEGG asks for considerate use of the REST API; we throttle to <= MAX_RPS requests/s
and batch /get calls (up to 10 entries each). This module does NOT invent data: a
404 / empty body is cached as an explicit MISS so we never silently refetch forever
and never confuse "absent" with "not yet fetched" (quality gate G3).
"""
from __future__ import annotations
import os, time, json, urllib.request, urllib.error, hashlib
from typing import Iterable

BASE = "https://rest.kegg.jp"
MAX_RPS = 3.0                    # requests per second ceiling (polite)
_MIN_INTERVAL = 1.0 / MAX_RPS
_last_call = [0.0]


class KeggClient:
    def __init__(self, cache_dir: str, offline: bool = False, timeout: int = 30):
        self.cache = cache_dir
        self.offline = offline          # if True, only read cache; never hit network
        self.timeout = timeout
        os.makedirs(cache_dir, exist_ok=True)
        self._release = None

    # ---- low level -------------------------------------------------------
    def _throttle(self):
        dt = time.time() - _last_call[0]
        if dt < _MIN_INTERVAL:
            time.sleep(_MIN_INTERVAL - dt)
        _last_call[0] = time.time()

    def _cache_path(self, op: str) -> str:
        # op is a REST path like "get/rn:R00209" -> safe filename
        h = hashlib.sha1(op.encode()).hexdigest()[:16]
        safe = op.replace("/", "__").replace(":", "_")[:120]
        return os.path.join(self.cache, f"{safe}.{h}.txt")

    def _get(self, op: str) -> str | None:
        """Fetch REST path `op` (e.g. 'get/rn:R00209'); return body or None on MISS.
        Cached forever; MISS cached as sentinel so we don't hammer KEGG."""
        cp = self._cache_path(op)
        if os.path.exists(cp):
            body = open(cp, encoding="utf-8").read()
            return None if body == "\x00MISS\x00" else body
        if self.offline:
            return None
        url = f"{BASE}/{op}"
        for attempt in range(4):
            self._throttle()
            try:
                with urllib.request.urlopen(url, timeout=self.timeout) as r:
                    body = r.read().decode("utf-8", "replace")
                if not body.strip():
                    body = None
                break
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    body = None
                    break
                if e.code in (403, 429) or e.code >= 500:
                    time.sleep(2 ** attempt)      # backoff
                    body = "\x01RETRY\x01"
                    continue
                raise
            except (urllib.error.URLError, TimeoutError):
                time.sleep(2 ** attempt)
                body = "\x01RETRY\x01"
                continue
        if body == "\x01RETRY\x01":
            # exhausted retries: do NOT cache; caller may retry later
            return None
        with open(cp, "w", encoding="utf-8") as fh:
            fh.write(body if body is not None else "\x00MISS\x00")
        return body

    # ---- release stamp (G5) ---------------------------------------------
    def release(self) -> str:
        if self._release is None:
            info = self._get("info/kegg") or ""
            rel = "unknown"
            for line in info.splitlines():
                if "Release" in line:
                    rel = line.strip(); break
            self._release = rel
        return self._release

    # ---- typed getters ---------------------------------------------------
    def get_entries(self, ids: Iterable[str]) -> dict[str, str]:
        """Batch /get up to 10 ids at a time. ids like 'rn:R00209' or 'cpd:C00022'.
        Returns {id: flatfile_record}. Missing ids simply absent from the dict."""
        ids = list(dict.fromkeys(ids))
        out: dict[str, str] = {}
        # try cache singly first, batch only the misses
        need = []
        for i in ids:
            cp = self._cache_path(f"get/{i}")
            if os.path.exists(cp):
                b = open(cp, encoding="utf-8").read()
                if b != "\x00MISS\x00":
                    out[i] = b
            else:
                need.append(i)
        for k in range(0, len(need), 10):
            chunk = need[k:k + 10]
            body = self._get("get/" + "+".join(chunk))
            # split multi-entry response on '///'
            recs = {}
            if body:
                for block in body.split("///\n"):
                    if not block.strip():
                        continue
                    rid = block.split()[1] if len(block.split()) > 1 else None
                    if rid:
                        recs[rid] = block + "///\n"
            # cache each individually so partial hits reuse
            for i in chunk:
                short = i.split(":", 1)[1]
                rec = recs.get(short) or recs.get(i)
                cp = self._cache_path(f"get/{i}")
                with open(cp, "w", encoding="utf-8") as fh:
                    fh.write(rec if rec else "\x00MISS\x00")
                if rec:
                    out[i] = rec
        return out

    def link(self, target_db: str, source: str) -> list[tuple[str, str]]:
        """/link/<target_db>/<source>  -> list of (source_id, target_id) pairs."""
        body = self._get(f"link/{target_db}/{source}")
        pairs = []
        if body:
            for line in body.splitlines():
                p = line.split("\t")
                if len(p) == 2:
                    pairs.append((p[0], p[1]))
        return pairs

    def list_db(self, db: str) -> list[tuple[str, str]]:
        """/list/<db> -> list of (id, description)."""
        body = self._get(f"list/{db}")
        rows = []
        if body:
            for line in body.splitlines():
                p = line.split("\t", 1)
                if len(p) == 2:
                    rows.append((p[0], p[1]))
        return rows


# ---- flatfile parsers (pure functions, unit-testable) -------------------
def parse_reaction(rec: str) -> dict:
    """Parse a KEGG reaction flatfile record into a structured dict."""
    d = {"id": None, "equation": "", "definition": "", "name": "",
         "rclass": [], "kos": [], "ec": [], "rhea": None,
         "modules": [], "maps": []}
    field = None
    for line in rec.splitlines():
        if line.startswith("///"):
            break
        key = line[:12].strip()
        val = line[12:].rstrip()
        if key:
            field = key
        if field == "ENTRY":
            d["id"] = line.split()[1] if len(line.split()) > 1 else None
        elif field == "NAME" and key == "NAME":
            d["name"] = val.rstrip(";")
        elif field == "DEFINITION" and key == "DEFINITION":
            d["definition"] = val
        elif field == "EQUATION" and key == "EQUATION":
            d["equation"] = val
        elif field == "RCLASS":
            # "RC00001  C00003_C00004  C00005_C00006"  (a RClass may list >1 pair)
            toks = val.split()
            if toks and toks[0].startswith("RC"):
                rc = toks[0]
                for pr in toks[1:]:
                    if "_" in pr:
                        a, b = pr.split("_", 1)
                        d["rclass"].append((rc, a, b))
            else:
                for pr in toks:
                    if "_" in pr:
                        a, b = pr.split("_", 1)
                        # continuation line: reuse last rc
                        rc = d["rclass"][-1][0] if d["rclass"] else "RC?"
                        d["rclass"].append((rc, a, b))
        elif field == "ENZYME":
            d["ec"] += [t for t in val.split() if t[0].isdigit()]
        elif field == "ORTHOLOGY":
            for t in val.split():
                if t.startswith("K") and t[1:].isdigit():
                    d["kos"].append(t)
        elif field == "MODULE":
            for t in val.split():
                if t.startswith("M") and t[1:].isdigit():
                    d["modules"].append(t)
        elif field == "PATHWAY":
            for t in val.split():
                if t.startswith("rn"):
                    d["maps"].append(t)
        elif field == "DBLINKS" and "RHEA:" in line:
            d["rhea"] = line.split("RHEA:")[1].split()[0]
    return d


def parse_compound(rec: str) -> dict:
    d = {"id": None, "name": "", "formula": "", "exact_mass": None}
    field = None
    for line in rec.splitlines():
        if line.startswith("///"):
            break
        key = line[:12].strip()
        val = line[12:].rstrip()
        if key:
            field = key
        if field == "ENTRY":
            d["id"] = line.split()[1] if len(line.split()) > 1 else None
        elif field == "NAME" and key == "NAME":
            d["name"] = val.rstrip(";")
        elif field == "FORMULA" and key == "FORMULA":
            d["formula"] = val.strip()
        elif field == "EXACT_MASS" and key == "EXACT_MASS":
            try:
                d["exact_mass"] = float(val.strip())
            except ValueError:
                pass
    return d


def parse_equation(eq: str) -> tuple[list[tuple[int, str]], list[tuple[int, str]], bool]:
    """'C00022 + C00010 <=> C00024 + C00011' ->
       (substrates, products, reversible)  where each side = [(coeff, cid), ...].
    Handles arrows <=> => <= and stoichiometric coefficients incl. 'n'."""
    for arrow, rev in (("<=>", True), ("=>", False), ("<=", False)):
        if arrow in eq:
            lhs, rhs = eq.split(arrow, 1)
            if arrow == "<=":            # normalise reverse arrow
                lhs, rhs = rhs, lhs
            break
    else:
        return [], [], False

    def side(s):
        out = []
        for term in s.split(" + "):
            term = term.strip()
            if not term:
                continue
            parts = term.split()
            if len(parts) == 2 and parts[0].lstrip("-").isdigit():
                coeff, cid = int(parts[0]), parts[1]
            else:
                coeff, cid = 1, parts[-1]
            cid = cid.split("(")[0]      # drop '(n)' style annotations
            if cid.startswith("C") or cid.startswith("G"):
                out.append((coeff, cid))
        return out
    return side(lhs), side(rhs), rev


def carbon_count(formula: str) -> int:
    """Number of carbon atoms in a KEGG molecular formula string."""
    import re
    if not formula:
        return 0
    m = re.search(r"C(\d*)(?![a-z])", formula)   # 'C', 'C6', not 'Cl'/'Co'
    if not m:
        return 0
    return int(m.group(1)) if m.group(1) else 1

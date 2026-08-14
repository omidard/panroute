#!/usr/bin/env python3
"""PanRoute local web server (stdlib only — no external deps, runs internally).

Serves the frontend and streams TRUE live pipeline progress over Server-Sent Events.

    python -m server.app            # then open http://localhost:8000

Endpoints:
    GET /                       -> web/index.html
    GET /assets/... /web/...    -> static files (KEGG map png, layout.json, js, css)
    GET /api/resolve?q=acetate  -> [{cid,name}]  (KEGG compound name search)
    GET /api/run?start=&end=&feedstock=  -> text/event-stream of live events
"""
import os, sys, json, re, threading, subprocess, urllib.parse, urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from panroute.engine import run_query
from panroute.keggfetch import KeggClient

PORT = int(os.environ.get("PANROUTE_PORT", "8000"))
DOCS = os.path.join(ROOT, "docs")           # the current client (same tree deployed to github.io)
ENZ_DIR = os.path.join(DOCS, "data", "enzymes")
ENZ_PY = os.path.join(ROOT, "bin", "enzyme_characterize.py")
_enz_locks = {}                              # rid -> Lock (serialise concurrent runs of the same reaction)
_enz_locks_guard = threading.Lock()
MIME = {".html": "text/html", ".js": "text/javascript", ".css": "text/css",
        ".json": "application/json", ".png": "image/png", ".svg": "image/svg+xml",
        ".gz": "application/gzip", ".woff2": "font/woff2", ".ico": "image/x-icon"}


def _enz_lock(rid):
    with _enz_locks_guard:
        return _enz_locks.setdefault(rid, threading.Lock())


def kegg_find_compound(q):
    """Resolve a metabolite name (or a C##### id) to KEGG compound candidates."""
    q = q.strip()
    if q.upper().startswith("C") and q[1:].isdigit():
        cl = KeggClient(os.path.join(ROOT, "cache"))
        rec = cl.get_entries([f"cpd:{q.upper()}"]).get(f"cpd:{q.upper()}", "")
        nm = q.upper()
        for line in rec.splitlines():
            if line.startswith("NAME"):
                nm = line[12:].strip().rstrip(";"); break
        return [{"cid": q.upper(), "name": nm}]
    try:
        url = f"https://rest.kegg.jp/find/compound/{urllib.parse.quote(q)}"
        body = urllib.request.urlopen(url, timeout=20).read().decode()
    except Exception:
        return []
    out = []
    for line in body.splitlines():
        parts = line.split("\t")
        if len(parts) == 2:
            cid = parts[0].replace("cpd:", "")
            name = parts[1].split(";")[0]
            out.append({"cid": cid, "name": name})
        if len(out) >= 15:
            break
    return out


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def _send(self, code, ctype, body, headers=None):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        if headers:
            for k, v in headers.items():
                self.send_header(k, v)
        self.end_headers()
        if body is not None:
            self.wfile.write(body if isinstance(body, bytes) else body.encode())

    def do_GET(self):
        u = urllib.parse.urlparse(self.path)
        path, qs = u.path, urllib.parse.parse_qs(u.query)

        if path == "/api/resolve":
            q = (qs.get("q") or [""])[0]
            return self._send(200, "application/json", json.dumps(kegg_find_compound(q)),
                              {"Access-Control-Allow-Origin": "*"})
        if path == "/api/run":
            return self._sse(qs)
        if path == "/api/enzyme":
            return self._sse_enzyme(qs)
        if path == "/api/ping":                          # lets the client detect a live backend
            return self._send(200, "application/json", '{"live":true}', {"Access-Control-Allow-Origin": "*"})

        # everything else: serve the docs/ client tree (same files deployed to github.io)
        rel = "index.html" if path in ("/", "") else path.lstrip("/")
        target = os.path.normpath(os.path.join(DOCS, rel))
        if not target.startswith(DOCS):                  # path-traversal guard
            return self._send(403, "text/plain", "forbidden")
        return self._file(target)

    def _file(self, fp):
        if not os.path.isfile(fp):
            return self._send(404, "text/plain", "not found")
        ext = os.path.splitext(fp)[1]
        with open(fp, "rb") as f:
            self._send(200, MIME.get(ext, "application/octet-stream"), f.read())

    def _sse(self, qs):
        start = (qs.get("start") or [""])[0].upper()
        end = (qs.get("end") or [""])[0].upper()
        feedstock = (qs.get("feedstock") or [""])[0].upper() or None
        maxlen = int((qs.get("max_len") or ["5"])[0])
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Connection", "keep-alive")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        def emit(ev, data):
            try:
                self.wfile.write(f"event: {ev}\ndata: {json.dumps(data)}\n\n".encode())
                self.wfile.flush()
            except (BrokenPipeError, ConnectionResetError):
                raise
        try:
            for ev, payload in run_query(start, end, feedstock, max_len=maxlen):
                emit(ev, payload)
            emit("close", {})
        except (BrokenPipeError, ConnectionResetError):
            return
        except Exception as e:
            import traceback; traceback.print_exc()
            try:
                emit("error", {"message": str(e)})
            except Exception:
                pass


    def _sse_enzyme(self, qs):
        """On-demand enzyme characterisation: gather KEGG orthologues, cluster, run MPEK + TemStaPro
        LIVE for this reaction, streaming progress. The result is cached to docs/data/enzymes/<rid>.json
        so a second request (and the static github.io client) get it instantly — compute happens only
        when a step is actually requested, never batch-precomputed."""
        rid = (qs.get("rid") or [""])[0]
        ko = (qs.get("ko") or [""])[0]
        sub = (qs.get("sub") or [""])[0]
        name = (qs.get("name") or [""])[0]
        temps = (qs.get("temps") or ["37"])[0]
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()

        def emit(ev, data):
            self.wfile.write(f"event: {ev}\ndata: {json.dumps(data)}\n\n".encode()); self.wfile.flush()

        # validate ids (list-arg subprocess already avoids shell injection; this rejects junk early)
        if not re.fullmatch(r"[A-Za-z0-9]+", rid or "") or not re.fullmatch(r"[A-Za-z0-9,]+", ko or ""):
            return emit("error", {"message": "bad rid/ko"})
        outp = os.path.join(ENZ_DIR, f"{rid}.json")
        try:
            if os.path.exists(outp):
                emit("progress", {"msg": "cached result — loading"})
                return emit("done", json.load(open(outp)))
            lock = _enz_lock(rid)
            if not lock.acquire(blocking=False):
                emit("progress", {"msg": "another request is computing this enzyme — waiting…"})
                lock.acquire()                            # block until the in-flight run finishes
            try:
                if os.path.exists(outp):                  # produced while we waited
                    return emit("done", json.load(open(outp)))
                emit("progress", {"msg": f"starting live analysis for {rid} (KO {ko}) — this runs the models once, then caches"})
                cmd = [sys.executable, ENZ_PY, "--rid", rid, "--ko", ko, "--sub-cid", sub,
                       "--sub-name", name, "--temps", temps, "--max-genes", "150", "--max-reps", "40", "--stage", "all"]
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
                for line in proc.stdout:
                    line = line.rstrip()
                    if line:
                        emit("progress", {"msg": line})
                proc.wait()
                if os.path.exists(outp):
                    emit("done", json.load(open(outp)))
                else:
                    emit("error", {"message": "analysis produced no result (no sequences, or a tool failed) — see server log"})
            finally:
                lock.release()
        except (BrokenPipeError, ConnectionResetError):
            return
        except Exception as e:
            import traceback; traceback.print_exc()
            try: emit("error", {"message": str(e)})
            except Exception: pass


def main():
    srv = ThreadingHTTPServer(("0.0.0.0", PORT), Handler)
    print(f"PanRoute live UI  ->  http://localhost:{PORT}")
    print("Ctrl-C to stop.")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        srv.shutdown()


if __name__ == "__main__":
    main()

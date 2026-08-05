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
import os, sys, json, urllib.parse, urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from panroute.engine import run_query
from panroute.keggfetch import KeggClient

PORT = int(os.environ.get("PANROUTE_PORT", "8000"))
MIME = {".html": "text/html", ".js": "text/javascript", ".css": "text/css",
        ".json": "application/json", ".png": "image/png", ".svg": "image/svg+xml"}


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

        if path == "/" or path == "/index.html":
            return self._file(os.path.join(ROOT, "web", "index.html"))
        if path.startswith("/web/"):
            return self._file(os.path.join(ROOT, path[1:]))
        if path.startswith("/assets/"):
            return self._file(os.path.join(ROOT, path[1:]))

        if path == "/api/resolve":
            q = (qs.get("q") or [""])[0]
            return self._send(200, "application/json", json.dumps(kegg_find_compound(q)),
                              {"Access-Control-Allow-Origin": "*"})

        if path == "/api/run":
            return self._sse(qs)

        return self._send(404, "text/plain", "not found")

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

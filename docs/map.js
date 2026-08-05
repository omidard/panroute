/* PanRoute live KEGG map (map01100). The whole global map is shown GRAY on a transparent
   background; the two searched metabolites blink; then the real reaction polylines that
   connect them light up in blue, one after another, like a Google-Maps route search finding
   all routes between two points. Coordinates are the real KGML polylines KEGG drew. */
const SVGNS = "http://www.w3.org/2000/svg";
function el(tag, a) { const e = document.createElementNS(SVGNS, tag); for (const k in a) e.setAttribute(k, a[k]); return e; }
const sleepM = ms => new Promise(r => setTimeout(r, ms));

class MapView {
  constructor(svg) { this.svg = svg; this.W = 4961; this.H = 3199; this.reset(true); }

  reset(hard) {
    this.svg.innerHTML = "";
    this.gEdges = el("g", {}); this.gNodes = el("g", {});
    if (hard) {
      const bg = el("image", { class: "bg", x: 0, y: 0, width: this.W, height: this.H });
      bg.setAttributeNS("http://www.w3.org/1999/xlink", "href", "assets/map01100/map01100_gray.png");
      bg.setAttribute("href", "assets/map01100/map01100_gray.png");
      this.svg.appendChild(bg);
    }
    this.svg.appendChild(this.gEdges); this.svg.appendChild(this.gNodes);
    this.svg.setAttribute("viewBox", `0 0 ${this.W} ${this.H}`);
    this._bbox = null; this._seen = new Set(); this._paths = [];
  }

  _grow(pts) { for (const [x, y] of pts) {
    if (!this._bbox) this._bbox = [x, y, x, y];
    this._bbox[0] = Math.min(this._bbox[0], x); this._bbox[1] = Math.min(this._bbox[1], y);
    this._bbox[2] = Math.max(this._bbox[2], x); this._bbox[3] = Math.max(this._bbox[3], y); } }

  focus(animate = true) {
    if (!this._bbox) return;
    let [x0, y0, x1, y1] = this._bbox;
    const pad = Math.max(300, (x1 - x0) * .3, (y1 - y0) * .3);
    x0 -= pad; y0 -= pad; x1 += pad; y1 += pad;
    const vw = this.svg.clientWidth || 900, vh = this.svg.clientHeight || 600, ar = vw / vh;
    let w = x1 - x0, h = y1 - y0; if (w / h > ar) h = w / ar; else w = h * ar;
    const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2, vb = `${cx - w / 2} ${cy - h / 2} ${w} ${h}`;
    animate ? this._animView(vb) : this.svg.setAttribute("viewBox", vb);
  }
  _animView(t) {
    const cur = (this.svg.getAttribute("viewBox")).split(/\s+/).map(Number), to = t.split(/\s+/).map(Number);
    const t0 = performance.now(), d = 800;
    const step = (n) => { const k = Math.min(1, (n - t0) / d), e = k < .5 ? 2 * k * k : 1 - Math.pow(-2 * k + 2, 2) / 2;
      this.svg.setAttribute("viewBox", cur.map((c, i) => c + (to[i] - c) * e).join(" "));
      if (k < 1) requestAnimationFrame(step); };
    requestAnimationFrame(step);
  }

  endpoint(xy, name, cls) {
    if (!xy) return; this._grow([xy]);
    const halo = el("circle", { cx: xy[0], cy: xy[1], r: 34, class: "ephalo " + cls });
    const c = el("circle", { cx: xy[0], cy: xy[1], r: 16, class: "epnode " + cls });
    this.gNodes.appendChild(halo); this.gNodes.appendChild(c);
    const t = el("text", { x: xy[0] + 30, y: xy[1] + 12, class: "nodelabel" }); t.textContent = name;
    this.gNodes.appendChild(t);
  }

  setEndpoints(start, end) {
    this.endpoint(end && end.xy, end ? end.name : "", "product");
    this.endpoint(start && start.xy, start ? start.name : "", "feedstock");
    this.focus(false);
  }

  _drawEdge(coords, cls, routeId) {
    const key = coords.map(p => p.join(",")).join(";");
    if (this._seen.has(key)) return Promise.resolve(); this._seen.add(key);
    this._grow(coords);
    const path = el("path", { d: "M" + coords.map(p => p.join(",")).join("L"), class: "edge" + cls });
    if (routeId != null) path.dataset.rid = routeId;
    this.gEdges.appendChild(path); this._paths.push({ path, key });
    const len = path.getTotalLength(); path.style.strokeDasharray = len; path.style.strokeDashoffset = len;
    return new Promise(res => { const t0 = performance.now(), d = 460;
      const a = (t) => { const k = Math.min(1, (t - t0) / d); path.style.strokeDashoffset = len * (1 - k);
        if (k < 1) requestAnimationFrame(a); else res(); };
      requestAnimationFrame(a); });
  }

  /* trace ALL routes' real map edges, spreading from the product toward the feedstock */
  async traceRoutes(routes, endName, endOnMap) {
    const edges = []; let gateway = null;
    for (const r of routes) for (const st of r.map.steps) {
      if (st.kind === "polyline" || st.kind === "connector") {
        edges.push({ coords: st.coords, cls: st.kind === "connector" ? " connector" : "" });
        if (!gateway) gateway = st.from_xy || st.to_xy || st.coords[0];
      }
    }
    if (!endOnMap && gateway) {                       // pin off-core-map product at the gateway
      const gx = gateway[0] - 150, gy = gateway[1] - 210;
      await this._drawEdge([[gateway[0], gateway[1]], [gx, gy]], " dashed");
      this.endpoint([gx, gy], endName, "product offmap");
    }
    for (const e of edges) { await this._drawEdge(e.coords, e.cls); await sleepM(60); this.focus(true); }
  }

  highlightRoute(route) {
    this._paths.forEach(p => p.path.classList.remove("hot"));
    for (const st of route.map.steps) if (st.coords) {
      const key = st.coords.map(p => p.join(",")).join(";");
      this._paths.forEach(p => { if (p.key === key) p.path.classList.add("hot"); });
    }
  }
}
window.MapView = MapView;

/* PanRoute live KEGG map (map01100). Renders the REAL global metabolic map, inverted to a
   dark "night map", and traces the route with an animated cyan highlight along the real
   reaction polylines KEGG drew — product → feedstock, like a route on Google Maps. */
const SVGNS = "http://www.w3.org/2000/svg";
function el(tag, attrs) { const e = document.createElementNS(SVGNS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]); return e; }

class MapView {
  constructor(svg) { this.svg = svg; this.W = 4961; this.H = 3199; this.reset(true); }

  reset(hard) {
    this.svg.innerHTML = "";
    this.gEdges = el("g", {}); this.gNodes = el("g", {});
    if (hard) {
      this.bg = el("image", { class: "bg", x: 0, y: 0, width: this.W, height: this.H,
        href: "/assets/map01100/map01100.png" });
      this.bg.setAttributeNS("http://www.w3.org/1999/xlink", "href", "/assets/map01100/map01100.png");
      this.svg.appendChild(this.bg);
    }
    this.svg.appendChild(this.gEdges); this.svg.appendChild(this.gNodes);
    this.svg.setAttribute("viewBox", `0 0 ${this.W} ${this.H}`);
    this._bbox = null;
  }

  _grow(pts) {
    for (const [x, y] of pts) {
      if (!this._bbox) this._bbox = [x, y, x, y];
      this._bbox[0] = Math.min(this._bbox[0], x); this._bbox[1] = Math.min(this._bbox[1], y);
      this._bbox[2] = Math.max(this._bbox[2], x); this._bbox[3] = Math.max(this._bbox[3], y);
    }
  }

  focus(animate = true) {
    if (!this._bbox) return;
    let [x0, y0, x1, y1] = this._bbox;
    const pad = Math.max(340, (x1 - x0) * 0.25, (y1 - y0) * 0.25);
    x0 -= pad; y0 -= pad; x1 += pad; y1 += pad;
    // keep aspect ratio of the svg viewport
    const vw = this.svg.clientWidth || 900, vh = this.svg.clientHeight || 600;
    let w = x1 - x0, h = y1 - y0;
    const ar = vw / vh;
    if (w / h > ar) h = w / ar; else w = h * ar;
    const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
    const vb = `${cx - w / 2} ${cy - h / 2} ${w} ${h}`;
    if (animate) this._animView(vb); else this.svg.setAttribute("viewBox", vb);
  }

  _animView(target) {
    const cur = (this.svg.getAttribute("viewBox") || `0 0 ${this.W} ${this.H}`).split(/\s+/).map(Number);
    const to = target.split(/\s+/).map(Number);
    const t0 = performance.now(), dur = 700;
    const step = (t) => {
      const k = Math.min(1, (t - t0) / dur); const e = k < .5 ? 2 * k * k : 1 - Math.pow(-2 * k + 2, 2) / 2;
      const vb = cur.map((c, i) => c + (to[i] - c) * e);
      this.svg.setAttribute("viewBox", vb.join(" "));
      if (k < 1) requestAnimationFrame(step);
    };
    requestAnimationFrame(step);
  }

  setEndpoints(start, end) {
    if (end && end.xy) { this._grow([end.xy]); this._node(end.xy, end.name, "product", true); }
    if (start && start.xy) { this._grow([start.xy]); this._node(start.xy, start.name, "feedstock", true); }
    this.focus(false);
  }

  _node(xy, label, cls, pulse) {
    const c = el("circle", { cx: xy[0], cy: xy[1], r: 16, class: "node " + cls });
    if (pulse) c.classList.add("pulse");
    this.gNodes.appendChild(c);
    if (label) {
      const t = el("text", { x: xy[0] + 22, y: xy[1] + 10, class: "nodelabel" });
      t.textContent = label; this.gNodes.appendChild(t);
    }
  }

  _drawPath(coords, cls) {
    this._grow(coords);
    const d = "M" + coords.map(p => p.join(",")).join("L");
    const path = el("path", { d, class: "edge" + cls });
    this.gEdges.appendChild(path);
    const len = path.getTotalLength();
    path.style.strokeDasharray = len; path.style.strokeDashoffset = len;
    return new Promise(res => {
      const t0 = performance.now(), dur = 520;
      const anim = t => { const k = Math.min(1, (t - t0) / dur);
        path.style.strokeDashoffset = len * (1 - k);
        if (k < 1) requestAnimationFrame(anim); else res(); };
      requestAnimationFrame(anim);
    });
  }

  /* animate one resolved step (retro order); on-map = real polyline, off-map = peripheral chip */
  async drawStep(step) {
    if (step.kind !== "offmap" && step.coords) {
      const anchor = step.from_xy || step.to_xy;   // substrate side (retro)
      if (step.from_xy) this._node(step.from_xy, step.from_name || "", "", false);
      if (step.to_xy) this._node(step.to_xy, step.to_name || "", "", false);
      this._anchor = step.from_xy || this._anchor;  // move toward feedstock
      this.focus(true);
      await this._drawPath(step.coords, step.kind === "connector" ? " connector" : "");
      await this._flushChips();
      return;
    }
    // off-map: stash until we have an on-map anchor, then lay chips outward
    const offName = step.from_xy ? step.to_name : step.from_name;   // whichever is off-map
    this._chips = this._chips || [];
    this._chips.push(offName);
    if (this._anchor) await this._flushChips();
  }

  async _flushChips() {
    if (!this._chips || !this._chips.length || !this._anchor) return;
    let [ax, ay] = this._anchor, dx = -220, dy = -160;
    while (this._chips.length) {
      const name = this._chips.shift();
      const nx = ax + dx, ny = ay + dy;
      this._grow([[nx, ny], [ax, ay]]);
      await this._drawPath([[ax, ay], [nx, ny]], " dashed");
      const chip = el("rect", { x: nx - 4, y: ny - 34, rx: 8, width: Math.max(120, name.length * 15 + 24),
        height: 46, fill: "#0c1524", stroke: "#5a6b86" });
      const t = el("text", { x: nx + 14, y: ny - 2, class: "nodelabel" }); t.textContent = name;
      this.gNodes.appendChild(chip); this.gNodes.appendChild(t);
      ax = nx + Math.max(120, name.length * 15 + 24); ay = ny; dx = -80;
      this.focus(true);
    }
  }

  finalizeRoute(mapRoute) {
    // draw remaining nodes of the traced route + off-map chips
    let lastOn = null;
    (mapRoute.nodes || []).forEach(n => {
      if (n.xy) { this._node(n.xy, "", "", false); lastOn = n.xy; }
    });
  }
}
window.MapView = MapView;

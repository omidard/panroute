/* PanRoute client-side engine — runs the full retrosynthetic search, genome gating,
   feedstock direction-gating and map resolution IN THE BROWSER (no backend), over the
   bundled static KEGG data. Mirrors the Python engine (panroute/*.py). Emits the same
   events the UI already consumes, so any start->end pair runs live.

   Node-testable: pass a custom `loader` (path -> parsed JSON) to the constructor. */
(function (global) {
  "use strict";

  const REV_MARGIN = 30;   // kJ/mol reversibility window (when ΔG available)

  class PanRoute {
    constructor(opts) {
      opts = opts || {};
      this.base = opts.base || "data/";
      this.mapBase = opts.mapBase || "assets/map01100/";
      this.loader = opts.loader || (url => fetch(url).then(r => r.json()));
      this.net = null; this.layout = null; this.tax = null; this.feed = null; this.thermo = {};
      this.koCache = {};
    }

    async load() {
      if (this.net) return;
      [this.net, this.tax, this.layout, this.feed] = await Promise.all([
        this.loader(this.base + "network.json"),
        this.loader(this.base + "taxonomy.json"),
        this.loader(this.mapBase + "layout.json"),
        this.loader(this.base + "feedstock_rules.json").catch(() => ({ feedstocks: {} })),
      ]);
      this.thermo = await this.loader(this.base + "thermo.json").catch(() => ({}));
      // adjacency
      this.out = {}; this.inn = {};
      for (const [s, d, rid] of this.net.edges) {
        (this.out[s] = this.out[s] || []).push([d, rid]);
        (this.inn[d] = this.inn[d] || []).push([s, rid]);
      }
    }

    cname(c) { const v = this.net.compounds[c]; return v ? (v.n || c).split(";")[0] : c; }

    // ---- retro search ----
    revDist(end) {
      const dist = { [end]: 0 }, q = [end];
      for (let i = 0; i < q.length; i++) {
        const v = q[i];
        for (const [s] of (this.inn[v] || [])) if (!(s in dist)) { dist[s] = dist[v] + 1; q.push(s); }
      }
      return dist;
    }

    enumerate(start, end, maxLen, maxRoutes) {
      const dist = this.revDist(end);
      if (!(start in dist) || dist[start] > maxLen) return [];
      const routes = []; let exp = 0;
      const stack = [[start, [start], new Set([start])]];
      while (stack.length && routes.length < maxRoutes && exp < 400000) {
        const [node, path, onp] = stack.pop();
        if (node === end && path.length > 1) { routes.push(path); continue; }
        const rem = maxLen - (path.length - 1); if (rem <= 0) continue;
        const nxt = new Set();
        for (const [v] of (this.out[node] || [])) {
          if (onp.has(v)) continue;
          const dv = dist[v]; if (dv === undefined || dv > rem - 1) continue;
          nxt.add(v);
        }
        exp++;
        [...nxt].sort((a, b) => (dist[b] || 1e9) - (dist[a] || 1e9))
          .forEach(v => stack.push([v, path.concat(v), new Set(onp).add(v)]));
      }
      routes.sort((a, b) => a.length - b.length);
      return routes.slice(0, maxRoutes).map(p => this.buildRoute(p));
    }

    reactionsFor(u, v) {   // reactions realising step u->v (unique)
      const seen = new Set(), out = [];
      for (const [d, rid] of (this.out[u] || [])) if (d === v && !seen.has(rid)) {
        seen.add(rid);
        const r = this.net.rxn[rid] || {};
        out.push({ rid, kos: r.k || [], ec: r.e ? [r.e] : [] });
      }
      return out;
    }

    buildRoute(path) {
      const steps = [];
      for (let i = 0; i < path.length - 1; i++) {
        const rxns = this.reactionsFor(path[i], path[i + 1]);
        const ec = new Set(); rxns.forEach(r => r.ec.forEach(e => ec.add(e)));
        steps.push({ from: path[i], to: path[i + 1], reactions: rxns, enzymes: [...ec].slice(0, 2).join("/") });
      }
      return { length: path.length - 1,
        path: path.map(c => ({ cid: c, name: this.cname(c) })), steps };
    }

    // ---- map resolution (mirror mapviz.py) ----
    resolveStep(from, to, rid) {
      const L = this.layout;
      if (L.reactions[rid]) return { kind: "polyline", coords: L.reactions[rid], reaction: rid };
      const key = [from, to].sort().join("|");
      for (const r of (L.pair_rxn[key] || [])) if (L.reactions[r]) return { kind: "polyline", coords: L.reactions[r], reaction: r };
      const a = L.compounds[from], b = L.compounds[to];
      if (a && b) return { kind: "connector", coords: [a, b], reaction: rid };
      return { kind: "offmap", coords: null, reaction: rid };
    }
    xy(c) { return this.layout.compounds[c] || null; }

    // ---- feedstock direction gating (mirror feedstock.py) ----
    compileExpr(expr) {
      let w = expr.trim().replace(/\s*,\s*/g, ",").replace(/(K\d{5})/g, "H.has('$1')")
        .replace(/\+/g, " ").replace(/\s+/g, " ").trim().replace(/ /g, " && ").replace(/,/g, " || ");
      return new Function("H", "return (" + w + ");");
    }
    feedStatus(cid, H) {
      const f = (this.feed.feedstocks || {})[cid]; if (!f) return "n/a";
      try { if (f.uptake && this.compileExpr(f.uptake.definition)(H)) return "uptake"; } catch (e) {}
      try { if (f.overflow_ambiguous && this.compileExpr(f.overflow_ambiguous.definition)(H)) return "overflow_capable"; } catch (e) {}
      return "none";
    }
    feedKos(cid) {
      const f = (this.feed.feedstocks || {})[cid]; const s = new Set(); if (!f) return s;
      for (const k of ["uptake", "overflow_ambiguous"]) if (f[k]) (f[k].definition.match(/K\d{5}/g) || []).forEach(x => s.add(x));
      return s;
    }

    async koOrgs(ko) {
      if (this.koCache[ko]) return this.koCache[ko];
      const d = await this.loader(this.base + "ko/" + ko + ".json").catch(() => ({ orgs: [] }));
      return (this.koCache[ko] = new Set(d.orgs || []));
    }

    // ---- thermo feasibility of a route (ΔG if available, else KEGG arrow) ----
    routeFeasible(route) {
      let ok = true, dgs = [];
      for (const st of route.steps) {
        const rid = st.reactions[0].rid, dg = this.thermo[rid];
        if (typeof dg === "number") { dgs.push(dg); if (dg > REV_MARGIN) ok = false; }
      }
      return { feasible: ok, dG_sum: dgs.length ? Math.round(dgs.reduce((a, b) => a + b, 0) * 10) / 10 : null };
    }
    species(name) { const t = name.replace(/'/g, "").split(/\s+/); return t.length >= 2 ? t.slice(0, 2).join(" ") : t[0]; }

    // ---- full run: emits events (mirror engine.run_query) ----
    async run(start, end, feedstock, emit, opts) {
      opts = opts || {}; const maxLen = opts.maxLen || 5, maxRoutes = opts.maxRoutes || 60;
      await this.load();
      emit("phase", { msg: "searching routes product → feedstock", pct: 12 });
      emit("endpoints", { start: { cid: start, name: this.cname(start), xy: this.xy(start) },
        end: { cid: end, name: this.cname(end), xy: this.xy(end) }, map: this.layout.image });
      const routes = this.enumerate(start, end, maxLen, maxRoutes);
      if (!routes.length) { emit("done", { error: `no route from ${this.cname(start)} to ${this.cname(end)}` }); return; }
      routes.forEach((r, i) => { r.id = i; r.map = this.resolveRoute(r); });
      const shortest = routes.reduce((a, b) => b.length < a.length ? b : a, routes[0]);
      // animate shortest, retro order
      const steps = shortest.map.steps;
      for (let k = 0; k < steps.length; k++) {
        const st = steps[k];
        emit("explore", { step: Object.assign({}, st, {
          from_name: this.cname(st.from), to_name: this.cname(st.to), from_xy: this.xy(st.from), to_xy: this.xy(st.to) }),
          index: k, total: steps.length, pct: 12 + Math.round(28 * (k + 1) / steps.length) });
      }
      emit("routes", { routes, shortest_len: shortest.length, n_routes: routes.length,
        capped: routes.length >= maxRoutes });
      // thermo feasibility
      const feas = {};
      routes.forEach(r => { feas[r.id] = this.routeFeasible(r); emit("thermo", { route_id: r.id, feasible: feas[r.id].feasible, length: r.length, dG_sum: feas[r.id].dG_sum }); });

      if (opts.skipGating) {              // genome data not uploaded yet: routes+map+pathways only
        emit("done", { n_routes: routes.length, shortest: shortest.length,
          genome_pending: true, kegg_release: this.net.kegg_release });
        return;
      }
      emit("phase", { msg: "gating KEGG genomes on the routes", pct: 55 });
      const routeKos = new Set(); routes.forEach(r => r.steps.forEach(s => s.reactions.forEach(x => x.kos.forEach(k => routeKos.add(k)))));
      const fkos = feedstock ? this.feedKos(feedstock) : new Set();
      const allKos = new Set([...routeKos, ...fkos]);
      // fetch per-KO genome lists (on demand), build org->KOset
      const orgKO = {};
      await Promise.all([...allKos].map(async ko => {
        const orgs = await this.koOrgs(ko);
        orgs.forEach(o => (orgKO[o] = orgKO[o] || new Set()).add(ko));
      }));
      // gate
      const stepKos = r => r.steps.map(s => new Set(s.reactions.flatMap(x => x.kos)));
      const routeStepKos = routes.map(stepKos);
      const encodes = (H, ri) => routeStepKos[ri].every(ks => { for (const k of ks) if (H.has(k)) return true; return false; });
      const termKos = new Set(); routes.forEach(r => r.steps[r.steps.length - 1].reactions.forEach(x => x.kos.forEach(k => termKos.add(k))));
      const bySpecies = {}; let t0 = new Set();
      for (const code in orgKO) {
        const H = orgKO[code], t = this.tax[code]; if (!t) continue;
        const idx = []; routes.forEach(r => { if (encodes(H, r.id)) idx.push(r.id); });
        for (const k of termKos) if (H.has(k)) { t0.add(t[0]); break; }
        if (!idx.length) continue;
        const sp = t[0];
        const feasible = idx.some(i => feas[i].feasible);
        const cur = bySpecies[sp];
        if (!cur || idx.length > cur.n_routes)
          bySpecies[sp] = { species: sp, gram: t[1], domain: t[2], n_routes: idx.length, route_idx: idx,
            code, thermo_feasible: feasible };
      }
      const rows = Object.values(bySpecies).sort((a, b) => b.n_routes - a.n_routes);
      const routeGenomes = {};      // route id -> # species that encode it
      for (const row of rows) for (const i of row.route_idx) routeGenomes[i] = (routeGenomes[i] || 0) + 1;
      emit("route_genomes", routeGenomes);
      emit("phase", { msg: `found ${rows.length} species with a native route`, pct: 82 });
      const gc = { Gpos: 0, Gneg: 0, Arch: 0, Other: 0 };
      let up = 0, ov = 0;
      for (const r of rows) {
        gc[r.gram] = (gc[r.gram] || 0) + 1;
        let feed = "n/a";
        if (feedstock) { feed = this.feedStatus(feedstock, orgKO[r.code] || new Set()); up += feed === "uptake"; ov += feed === "overflow_capable"; }
        r.feedstock = feed;
        emit("organism", r);
      }
      emit("done", { n_routes: routes.length, shortest: shortest.length,
        T0: t0.size, T2: rows.length, T3: feedstock ? up : null, overflow_excluded: feedstock ? ov : null,
        gram: gc, kegg_release: this.net.kegg_release });
    }

    resolveRoute(route) {
      const nodes = route.path.map(p => ({ cid: p.cid, name: p.name, xy: this.xy(p.cid), onmap: !!this.xy(p.cid) }));
      const steps = route.steps.map(st => Object.assign({ from: st.from, to: st.to, enzymes: st.enzymes },
        this.resolveStep(st.from, st.to, st.reactions[0].rid)));
      return { nodes: nodes.reverse(), steps: steps.reverse() };
    }
  }

  global.PanRoute = PanRoute;
  if (typeof module !== "undefined" && module.exports) module.exports = PanRoute;
})(typeof window !== "undefined" ? window : globalThis);

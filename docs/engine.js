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
      this.aliases = await this.loader(this.base + "aliases.json").catch(() => ({}));
      this.qual = await this.loader(this.base + "rxnqual.json").catch(() => ({}));   // reaction curation 0..3
      this.rxnko = await this.loader(this.base + "rxnko.json").catch(() => ({}));    // multi-EC complex KO groups
      this.sides = await this.loader(this.base + "rxnsides.json").catch(() => ({})); // L/R sides for signed ΔG
      this.dgSuspect = new Set(await this.loader(this.base + "thermo_suspect.json").catch(() => [])); // |ΔG|>200 = unbalanced artifact

      // ---- collapse synonym compound IDs to one canonical node, IN THE GRAPH ----
      // KEGG fragments a metabolite across several compound IDs (e.g. 2-acetolactate C00900 /
      // (S)-C06010), which splits a single pathway across un-connected nodes. Merge each alias
      // group to one representative so the graph is connected.
      this.canon = {}; this.canonName = {}; this.canonXY = {};
      const seen = new Set();
      for (const cid in this.aliases) {
        const g = this.aliases[cid], key = g.join(",");
        if (seen.has(key)) continue; seen.add(key);
        const rep = g.find(c => this.layout.compounds[c]) || g[0];   // prefer one with map coords
        for (const c of g) this.canon[c] = rep;
        const names = g.map(c => (this.net.compounds[c] || {}).n || c);
        this.canonName[rep] = names.reduce((a, b) => (a && a.length <= b.length ? a : b));   // shortest = generic
        const xc = g.find(c => this.layout.compounds[c]);
        if (xc) this.canonXY[rep] = this.layout.compounds[xc];
      }
      const CN = c => this.canon[c] || c;

      // adjacency (over canonical nodes)
      this.out = {}; this.inn = {};
      for (const [s, d, rid] of this.net.edges) {
        const cs = CN(s), cd = CN(d);
        if (cs === cd) continue;                       // self-loop introduced by the merge
        (this.out[cs] = this.out[cs] || []).push([cd, rid]);
        (this.inn[cd] = this.inn[cd] || []).push([cs, rid]);
      }
    }

    cn(c) { return this.canon[c] || c; }
    cname(c) { const r = this.cn(c); if (this.canonName[r]) return this.canonName[r].split(";")[0];
      const v = this.net.compounds[r] || this.net.compounds[c]; return v ? (v.n || c).split(";")[0] : c; }

    // ---- retro search ----
    revDist(endSet) {
      const dist = {}, q = [];
      for (const e of endSet) if (!(e in dist)) { dist[e] = 0; q.push(e); }
      for (let i = 0; i < q.length; i++) {
        const v = q[i];
        for (const [s] of (this.inn[v] || [])) if (!(s in dist)) { dist[s] = dist[v] + 1; q.push(s); }
      }
      return dist;
    }

    enumerate(startSet, endSet, maxLen, maxRoutes) {
      const dist = this.revDist(endSet);
      const starts = [...startSet].filter(s => (s in dist) && dist[s] <= maxLen);
      if (!starts.length) return [];
      // Enumerate a LARGE candidate pool (not just maxRoutes), then rank by QUALITY so
      // canonical routes (no reaction-reuse, well-characterised, short) survive instead of
      // being crowded out of a small length-only cap by exotic carbon-skeleton shortcuts.
      const ENUM = Math.max(maxRoutes * 8, 600);
      const paths = []; let exp = 0;
      const stack = starts.map(s => [s, [s], new Set([s])]);
      while (stack.length && paths.length < ENUM && exp < 800000) {
        const [node, path, onp] = stack.pop();
        if (endSet.has(node) && path.length > 1) { paths.push(path); continue; }
        const rem = maxLen - (path.length - 1); if (rem <= 0) continue;
        const nxt = new Set();
        for (const [v] of (this.out[node] || [])) {
          if (onp.has(v)) continue;
          const dv = dist[v]; if (dv === undefined || dv > rem - 1) continue;
          nxt.add(v);
        }
        exp++;
        [...nxt].sort((a, b) => ((b in dist) ? dist[b] : 1e9) - ((a in dist) ? dist[a] : 1e9))
          .forEach(v => stack.push([v, path.concat(v), new Set(onp).add(v)]));
      }
      const built = paths.map(p => this.buildRoute(p));
      // A reaction that realises >1 step of a single linear route is a carbon-skeleton graph
      // ARTIFACT — one RCLASS atom-mapping matched to two different metabolite transitions
      // (e.g. glucose→6-acetyl-glucose→acetate both "via R00327"). A real pathway never reuses
      // one reaction along a simple path, and no genome encodes these — drop them outright so
      // the results are trustworthy (if that empties the set, the honest "can't trace this /
      // FBA territory" state is shown instead of biochemical nonsense).
      const clean = built.filter(r => r.repeats === 0);
      const use = clean.length ? clean : [];
      use.sort((a, b) => (a.uncur - b.uncur) || (a.length - b.length));
      return use.slice(0, maxRoutes);
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
      // #obscure reactions (no subsystem/BiGG/Rhea mapping) and reaction-reuse (the SAME
      // reaction realising >1 step is a carbon-skeleton graph artifact — one enzyme is not a
      // whole pathway) — both push a route down the ranking below genuine pathways.
      const uncur = steps.filter(st => ((this.qual || {})[st.reactions[0].rid] || 0) === 0).length;
      const prim = steps.map(st => st.reactions[0].rid);
      const repeats = prim.length - new Set(prim).size;
      return { length: path.length - 1, uncur, repeats,
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
    xy(c) { return this.canonXY[c] || this.layout.compounds[c] || (this.canon[c] ? this.layout.compounds[this.canon[c]] : null) || null; }

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
      // ΔG is stored as-written (left => right). A step that crosses a reaction from its
      // product side to its substrate side is a REVERSE traversal whose effective ΔG is -dg
      // (audit A5). Determine the direction from rxnsides and sign ΔG accordingly.
      let ok = true, dgs = [], unknown = 0;
      const CN = c => this.canon[c] || c;
      for (const st of route.steps) {
        const rid = st.reactions[0].rid, dg = this.thermo[rid];
        // |ΔG|>200 kJ/mol almost always means an unbalanced/generic-R-group equation, not real
        // thermodynamics — treat as unknown rather than let it flip feasibility (audit A5).
        if (typeof dg !== "number" || this.dgSuspect.has(rid)) { unknown++; continue; }
        let eff = dg;
        const sd = this.sides[rid];
        if (sd) {
          const L = new Set(sd.L.map(CN)), R = new Set(sd.R.map(CN));
          const f = CN(st.from), t = CN(st.to);
          if (R.has(f) && L.has(t) && !(L.has(f) && R.has(t))) eff = -dg;   // reverse traversal
          else if (L.has(f) && R.has(t)) eff = dg;                          // forward (as-written)
          // else ambiguous (shared compound / not found): keep as-written, best effort
        }
        dgs.push(eff);
        if (eff > REV_MARGIN) ok = false;                                   // strongly uphill as traversed
      }
      return { feasible: ok, dG_sum: dgs.length ? Math.round(dgs.reduce((a, b) => a + b, 0) * 10) / 10 : null,
        dG_known: dgs.length, dG_unknown: unknown };
    }
    species(name) { const t = name.replace(/'/g, "").split(/\s+/); return t.length >= 2 ? t.slice(0, 2).join(" ") : t[0]; }

    // ---- full run: emits events (mirror engine.run_query) ----
    async run(start, end, feedstock, emit, opts) {
      opts = opts || {}; const maxLen = opts.maxLen || 5, maxRoutes = opts.maxRoutes || 60;
      await this.load();
      emit("phase", { msg: "searching routes product → feedstock", pct: 12 });
      emit("endpoints", { start: { cid: start, name: this.cname(start), xy: this.xy(start) },
        end: { cid: end, name: this.cname(end), xy: this.xy(end) }, map: this.layout.image });
      // synonym IDs are already merged in the graph (this.canon); search canonical endpoints
      const startC = this.cn(start), endC = this.cn(end);
      const routes = this.enumerate(new Set([startC]), new Set([endC]), maxLen, maxRoutes);
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
      // reaction satisfied: a multi-EC complex needs ALL catalytic components (AND across EC
      // groups), OR within isozymes; a single-enzyme reaction = OR over its KOs. A step is
      // satisfied if ANY realising reaction is satisfied; a route is encoded if every step is.
      const RK = this.rxnko || {};
      const rxnSat = (H, rx) => { const g = RK[rx.rid];
        return g ? g.every(grp => grp.some(k => H.has(k))) : (rx.kos || []).some(k => H.has(k)); };
      const encodes = (H, route) => route.steps.every(s => s.reactions.some(rx => rxnSat(H, rx)));
      const termKos = new Set(); routes.forEach(r => r.steps[r.steps.length - 1].reactions.forEach(x => x.kos.forEach(k => termKos.add(k))));
      const bySpecies = {}; let t0 = new Set();
      for (const code in orgKO) {
        const H = orgKO[code], t = this.tax[code]; if (!t) continue;
        const idx = []; routes.forEach(r => { if (encodes(H, r)) idx.push(r.id); });
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

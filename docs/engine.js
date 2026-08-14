/* PanRoute client-side engine — runs the full retrosynthetic search, genome gating,
   feedstock direction-gating and map resolution IN THE BROWSER (no backend), over the
   bundled static KEGG data. Mirrors the Python engine (panroute/*.py). Emits the same
   events the UI already consumes, so any start->end pair runs live.

   Node-testable: pass a custom `loader` (path -> parsed JSON) to the constructor. */
(function (global) {
  "use strict";

  const REV_MARGIN = 30;   // kJ/mol reversibility window (when ΔG available)

  // binary min-heap keyed by `.w` (best-first search over weighted partial paths)
  class MinHeap {
    constructor() { this.a = []; }
    size() { return this.a.length; }
    push(x) { const a = this.a; a.push(x); let i = a.length - 1;
      while (i > 0) { const p = (i - 1) >> 1; if (a[p].w <= a[i].w) break; [a[p], a[i]] = [a[i], a[p]]; i = p; } }
    pop() { const a = this.a, top = a[0], last = a.pop();
      if (a.length) { a[0] = last; let i = 0; const n = a.length;
        for (;;) { let l = 2 * i + 1, r = l + 1, m = i;
          if (l < n && a[l].w < a[m].w) m = l; if (r < n && a[r].w < a[m].w) m = r;
          if (m === i) break; [a[m], a[i]] = [a[i], a[m]]; i = m; } }
      return top; } }

  class PanRoute {
    constructor(opts) {
      opts = opts || {};
      this.base = opts.base || "data/";
      this.mapBase = opts.mapBase || "assets/map01100/";
      this.loader = opts.loader || (url => fetch(url).then(r => r.json()));
      this.net = null; this.layout = null; this.tax = null; this.feed = null; this.thermo = {};
      this.koCache = {}; this.traits = {};
      // host-tractability tie-break for chassis/donor ranking. Native coverage is ALWAYS primary
      // (the point is to clone the fewest genes into an organism that already runs most of the
      // pathway); this only orders organisms that are otherwise equal. Curated, not exhaustive.
      this.PREF = { "Escherichia coli": 5, "Bacillus subtilis": 5, "Pseudomonas putida": 4.5,
        "Corynebacterium glutamicum": 4.5, "Lactococcus lactis": 3.5, "Cupriavidus necator": 3.5,
        "Vibrio natriegens": 3.5, "Clostridium acetobutylicum": 3, "Bacillus megaterium": 2.5,
        "Streptomyces coelicolor": 2.5, "Clostridium ljungdahlii": 2.5, "Pseudomonas fluorescens": 2.5 };
      this.PREF_GENUS = { Escherichia: 3, Bacillus: 3, Pseudomonas: 2.5, Corynebacterium: 2.5,
        Lactococcus: 2, Lactobacillus: 1.5, Streptomyces: 1.5, Clostridium: 1.5, Cupriavidus: 2,
        Rhodococcus: 1.5, Synechocystis: 1.5, Synechococcus: 1.5, Vibrio: 1.5 };
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
      this.excluded = await this.loader(this.base + "excluded.json").catch(() => ({}));   // currency/C1-sink metabolites + why
      this.traits = await this.loader(this.base + "traits.json").catch(() => ({}));   // species -> {ox,temp,safety} for chassis/donor ranking

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

      this.prev = await this.loader(this.base + "rxnprev.json").catch(() => ({})); // reaction genome prevalence

      // adjacency (over canonical nodes)
      this.out = {}; this.inn = {};
      for (const [s, d, rid] of this.net.edges) {
        const cs = CN(s), cd = CN(d);
        if (cs === cd) continue;                       // self-loop introduced by the merge
        (this.out[cs] = this.out[cs] || []).push([cd, rid]);
        (this.inn[cd] = this.inn[cd] || []).push([cs, rid]);
      }

      // ---- biological-plausibility edge WEIGHTS ----
      // Atom-conservation (RCLASS) admits carbon-skeleton edges that no organism runs (58% of
      // reactions are encoded by ZERO genomes). Hop-count search then finds absurd shortcuts.
      // Weight each edge by -log(genome prevalence of its BEST realising reaction): a real
      // 10-step pathway (all high-prevalence reactions) is cheaper than a 3-step path through
      // a rare one, so lowest-weight search recovers true metabolism instead of nonsense.
      const NGEN = this.net.n_genomes || 10513;
      this.NGEN = NGEN;
      const bestPrev = {};                             // "cs>cd" -> max reaction prevalence
      for (const [s, d, rid] of this.net.edges) {
        const cs = CN(s), cd = CN(d); if (cs === cd) continue;
        const key = cs + ">" + cd, p = this.prev[rid] || 0;
        if (bestPrev[key] === undefined || p > bestPrev[key]) bestPrev[key] = p;
      }
      this.outW = {};                                  // cs -> [[cd, weight], ...] (deduped)
      const STEP = 0.15;                               // small per-hop cost: breaks plausibility ties toward shorter
      for (const cs in this.out) {
        const seen = new Set(), arr = [];
        for (const [cd] of this.out[cs]) {
          if (seen.has(cd)) continue; seen.add(cd);
          const p = bestPrev[cs + ">" + cd] || 0;
          arr.push([cd, -Math.log((p + 1) / (NGEN + 1)) + STEP]);
        }
        this.outW[cs] = arr;
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
      // Weighted k-shortest-SIMPLE-paths (Dijkstra-style best-first). A min-heap always expands
      // the lowest-weight partial path, so complete routes emerge in ~increasing biological
      // IMPLAUSIBILITY: a true pathway made of high-prevalence reactions is popped long before
      // any shortcut through a rare/zero-genome reaction. Hop distance still prunes dead ends.
      const heap = new MinHeap();
      for (const s of starts) heap.push({ w: 0, node: s, path: [s], onp: new Set([s]) });
      const paths = [], seen = new Set();
      const POOL = Math.max(maxRoutes, 60);
      let pops = 0; const POP_CAP = 200000, HEAP_CAP = 400000;
      while (heap.size() && paths.length < POOL && pops < POP_CAP) {
        const cur = heap.pop(); pops++;
        if (endSet.has(cur.node) && cur.path.length > 1) {
          const sig = cur.path.join(","); if (!seen.has(sig)) { seen.add(sig); paths.push(cur.path); }
          continue;
        }
        const rem = maxLen - (cur.path.length - 1); if (rem <= 0) continue;
        const nbrs = this.outW[cur.node]; if (!nbrs) continue;
        for (const [v, w] of nbrs) {
          if (cur.onp.has(v)) continue;
          const dv = dist[v]; if (dv === undefined || dv > rem - 1) continue;
          if (heap.size() >= HEAP_CAP) break;
          const onp = new Set(cur.onp); onp.add(v);
          heap.push({ w: cur.w + w, node: v, path: cur.path.concat(v), onp });
        }
      }
      const built = paths.map(p => this.buildRoute(p));
      // A reaction that realises >1 step of a single linear route is a carbon-skeleton graph
      // ARTIFACT — one RCLASS atom-mapping matched to two different metabolite transitions
      // (e.g. glucose→6-acetyl-glucose→acetate both "via R00327"). A real pathway never reuses
      // one reaction along a simple path, and no genome encodes these — drop them outright so
      // the results are trustworthy (if that empties the set, the honest "can't trace this /
      // FBA territory" state is shown instead of biochemical nonsense).
      // `built` is already in ascending-weight (descending-plausibility) order from the heap —
      // do NOT re-sort by length, which is what buried real pathways before. Just drop the
      // reaction-reuse artifacts and keep the most plausible routes.
      const clean = built.filter(r => r.repeats === 0);
      return clean.slice(0, maxRoutes);
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

    // host-tractability score (chassis/donor tie-break only; native coverage is primary)
    hostScore(species, gram) {
      let s = 0;
      if (this.PREF[species] != null) s = this.PREF[species];
      else { const g = (species || "").split(" ")[0]; if (this.PREF_GENUS[g] != null) s = this.PREF_GENUS[g]; }
      const tr = (this.traits || {})[species] || {};
      if (tr.safety === "GRAS/QPS") s += 1.5; else if (tr.safety === "pathogen") s -= 2; else if (tr.safety === "opportunist") s -= 0.5;
      if (tr.temp === "mesophile") s += 0.4;
      if (tr.ox === "facultative" || tr.ox === "aerobe") s += 0.3;
      return s;
    }

    /* ---- heterologous-expression design ----
       Called only when NO genome encodes a full native route (T2 == 0). For each plausible
       route it finds the organism that natively runs the MOST steps (the best chassis / partial-
       pathway host), then mines the ENTIRE enzyme space — every KEGG genome that carries the gene
       — for donor organisms supplying each MISSING step. Result: a concrete "clone gene(s) X from
       donor A into chassis C" plan, per route, with the native vs heterologous split.
       Everything is derived from the real KEGG KO→genome mapping; nothing is invented. */
    heteroScenarios(routes, orgKO, rxnSat, feas) {
      const HET_ROUTES = 15, HET_SHOW = 6, MAX_DONORS = 6;
      const cand = routes.slice(0, HET_ROUTES);
      const codes = Object.keys(orgKO);
      const stepSat = (H, step) => step.reactions.some(rx => rxnSat(H, rx));

      // Pass A — best chassis per candidate route: maximise natively-encoded steps, break ties by
      // host tractability. An organism must carry ≥1 step to be a partial-pathway host.
      const best = {};                       // route.id -> {code, cover, sat:[stepIdx], hostScore}
      for (const route of cand) {
        let top = null;
        for (const code of codes) {
          const H = orgKO[code], t = this.tax[code]; if (!t) continue;
          const sat = [];
          for (let si = 0; si < route.steps.length; si++) if (stepSat(H, route.steps[si])) sat.push(si);
          if (!sat.length) continue;
          const hs = this.hostScore(t[0], t[1]);
          if (!top || sat.length > top.cover || (sat.length === top.cover && hs > top.hostScore))
            top = { code, cover: sat.length, sat, hostScore: hs };
        }
        if (top) best[route.id] = top;
      }

      // Pass B — for each chosen chassis, split native vs heterologous and mine donors for the gaps.
      const scen = [];
      for (const route of cand) {
        const b = best[route.id]; if (!b) continue;      // nobody carries any step of this route
        const chassisCode = b.code, t = this.tax[chassisCode];
        const satSet = new Set(b.sat), hetero = [];
        for (let si = 0; si < route.steps.length; si++) {
          if (satSet.has(si)) continue;
          const step = route.steps[si], rx = step.reactions[0] || {};
          const stepKOs = new Set(); step.reactions.forEach(r => (r.kos || []).forEach(k => stepKOs.add(k)));
          // donor mining across the whole enzyme space: every genome that satisfies this step
          const donorMap = {};                            // species -> {gram, kos:Set}
          for (const code of codes) {
            if (code === chassisCode) continue;
            const H = orgKO[code]; if (!stepSat(H, step)) continue;
            const tt = this.tax[code]; if (!tt) continue;
            const cur = donorMap[tt[0]] || (donorMap[tt[0]] = { gram: tt[1], kos: new Set() });
            for (const k of stepKOs) if (H.has(k)) cur.kos.add(k);
          }
          const donors = Object.entries(donorMap)
            .map(([sp, v]) => { const tr = this.traits[sp] || {};
              return { species: sp, gram: v.gram, kos: [...v.kos], ox: tr.ox, temp: tr.temp, safety: tr.safety, hs: this.hostScore(sp, v.gram) }; })
            .sort((a, c) => c.hs - a.hs || a.species.localeCompare(c.species))
            .slice(0, MAX_DONORS);
          // uncreditable = we can't name a bacterial gene to clone. Two very different reasons:
          //  no_ko            → no enzyme (KO) is mapped to this transition at all (often an RCLASS
          //                     graph shortcut) — the step may not be a real single reaction.
          //  no_bacterial_donor → the enzyme IS known (KO exists) but no bacterial genome carries it,
          //                     so the gene must come from a plant/fungal/engineered source.
          const reason = stepKOs.size === 0 ? "no_ko" : (donors.length === 0 ? "no_bacterial_donor" : null);
          hetero.push({ idx: si, from: step.from, from_name: this.cname(step.from), to: step.to, to_name: this.cname(step.to),
            rid: rx.rid, ec: (rx.ec || []).join("/"), enzymes: step.enzymes, kos: [...stepKOs],
            donors, n_donor_species: Object.keys(donorMap).length, uncreditable: !!reason, uncreditable_reason: reason });
        }
        if (!hetero.length) continue;
        const tr = this.traits[t[0]] || {};
        scen.push({ route_id: route.id, length: route.length,
          feasible: !!(feas[route.id] && feas[route.id].feasible),
          chassis: { code: chassisCode, species: t[0], gram: t[1], domain: t[2], cover: b.cover,
            ox: tr.ox, temp: tr.temp, safety: tr.safety, hostScore: b.hostScore },
          native_idx: [...satSet], hetero, n_clone: hetero.length,
          n_uncreditable: hetero.filter(h => h.uncreditable).length });
      }
      // realizable designs first: fewest steps with NO nameable bacterial gene (a design where
      // every gap can be filled from a real donor is actionable; one with an unfindable enzyme is
      // a research problem), then fewest genes to clone, thermo-feasible, better host, shorter.
      scen.sort((a, c) => a.n_uncreditable - c.n_uncreditable || a.n_clone - c.n_clone ||
        (c.feasible - a.feasible) || (c.chassis.hostScore - a.chassis.hostScore) || a.length - c.length);
      return scen.slice(0, HET_SHOW);
    }

    // ---- full run: emits events (mirror engine.run_query) ----
    async run(start, end, feedstock, emit, opts) {
      opts = opts || {}; const maxLen = opts.maxLen || 5, maxRoutes = opts.maxRoutes || 60;
      await this.load();
      emit("phase", { msg: "searching routes product → feedstock", pct: 12 });
      emit("endpoints", { start: { cid: start, name: this.cname(start), xy: this.xy(start) },
        end: { cid: end, name: this.cname(end), xy: this.xy(end) }, map: this.layout.image });
      // synonym IDs are already merged in the graph (this.canon); search canonical endpoints
      const startC = this.cn(start), endC = this.cn(end);
      // if an endpoint isn't in the carbon-skeleton network, say WHY (usually an excluded
      // currency / C1-sink metabolite like CO2) instead of a bare "no route".
      const inGraph = c => !!(this.out[c] || this.inn[c]);
      for (const [role, c, raw] of [["start (feedstock)", startC, start], ["product", endC, end]]) {
        if (inGraph(c)) continue;
        const ex = this.excluded[c] || this.excluded[raw];
        emit("done", { excluded: true, error: ex
          ? `${this.cname(raw)} can't be used as the ${role}. ${ex.why}`
          : `${this.cname(raw)} is not in the carbon-skeleton network — no atom-conserved (RCLASS) reactions connect it, so it can't start or end a route.` });
        return;
      }
      const routes = this.enumerate(new Set([startC]), new Set([endC]), maxLen, maxRoutes);
      if (!routes.length) { emit("done", { error: `No genome-independent carbon-skeleton route from ${this.cname(start)} to ${this.cname(end)} within ${maxLen} steps. They may be too far apart, or the conversion needs C–C cleavage/condensation the method can't trace.` }); return; }
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
      // no native producer? design heterologous-expression scenarios (chassis + donor genes)
      let hetero = null;
      if (rows.length === 0) {
        emit("phase", { msg: "no native producer — designing heterologous expression", pct: 90 });
        hetero = this.heteroScenarios(routes, orgKO, rxnSat, feas);
        emit("hetero", { scenarios: hetero });
      }
      emit("done", { n_routes: routes.length, shortest: shortest.length,
        T0: t0.size, T2: rows.length, T3: feedstock ? up : null, overflow_excluded: feedstock ? ov : null,
        gram: gc, hetero: hetero ? hetero.length : 0, kegg_release: this.net.kegg_release });
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

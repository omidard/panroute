/* PanRoute static site (github.io): loads a precomputed run bundle and REPLAYS the exact
   live experience client-side — real KEGG-map trace, live organism discovery, thermo
   feasibility, interactive results. No backend. Data is KEGG-derived (research use). */
const $ = s => document.querySelector(s);
const GCOL = { Gpos: "#3a8bff", Gneg: "#ff8a3a", Arch: "#22d18c", Other: "#8aa0bf" };
const GNAME = { Gpos: "Gram-positive", Gneg: "Gram-negative", Arch: "Archaea", Other: "other" };
const map = new MapView($("#map"));
let ST = null, RUNNING = false;
const sleep = ms => new Promise(r => setTimeout(r, ms));

const engine = new PanRoute();
let GENOME_READY = false, CPD = [];
fetch("data/ready.json").then(r => { if (r.ok) GENOME_READY = true; }).catch(() => {});
fetch("data/compounds.json").then(r => r.json()).then(d => { CPD = Object.entries(d).map(([cid, name]) => ({ cid, name }));
  document.dispatchEvent(new Event("cpd-ready")); }).catch(() => {});

/* ---- autocomplete over bundled compound names ---- */
// strip leading stereo/config descriptors so "glucose" matches "D-Glucose", "beta-D-Glucose"
const STEREO = /^(\(?[0-9rsezn,+\-' ]+\)?-|d-|l-|dl-|alpha-|beta-|gamma-|delta-|meso-|cis-|trans-|sn-|o-|m-|p-|n-)+/i;
function searchCpd(q) {
  q = q.trim().toLowerCase(); if (q.length < 2) return [];
  if (/^c\d{5}$/.test(q)) { const m = CPD.find(x => x.cid.toLowerCase() === q); return m ? [m] : []; }
  const scored = [];
  for (const x of CPD) {
    const n = x.name.toLowerCase(), base = n.replace(STEREO, "");
    let rank = -1;
    if (n === q || base === q) rank = 0;                 // exact (ignoring D-/beta-/… prefix)
    else if (base.startsWith(q) || n.startsWith(q)) rank = 1;   // starts with the query
    else if (new RegExp(`[ (\\-,]${q}`).test(n)) rank = 2;      // query at a word boundary
    else if (n.includes(q)) rank = 3;                    // anywhere
    if (rank >= 0) scored.push({ x, rank, len: n.length });
  }
  // best rank first; within a rank the SHORTER name (the base metabolite, not a derivative)
  scored.sort((a, b) => a.rank - b.rank || a.len - b.len);
  return scored.slice(0, 12).map(s => s.x);
}
function wireAC(inputId, acId, cidId) {
  const inp = $(inputId), ac = $(acId), cid = $(cidId); let items = [], sel = -1; inp.dataset.cid = "";
  inp.addEventListener("input", () => { cid.textContent = ""; inp.dataset.cid = "";
    items = searchCpd(inp.value); sel = -1;
    ac.innerHTML = items.map((x, i) => `<li data-i="${i}">${x.name}<span class="c">${x.cid}</span></li>`).join("");
    ac.classList.toggle("show", items.length > 0); });
  ac.addEventListener("mousedown", e => { const li = e.target.closest("li"); if (!li) return;
    const x = items[+li.dataset.i]; inp.value = x.name; inp.dataset.cid = x.cid; cid.textContent = x.cid; ac.classList.remove("show"); });
  inp.addEventListener("keydown", e => { if (!ac.classList.contains("show")) return; const lis = [...ac.children];
    if (e.key === "ArrowDown") { sel = Math.min(sel + 1, lis.length - 1); e.preventDefault(); }
    else if (e.key === "ArrowUp") { sel = Math.max(sel - 1, 0); e.preventDefault(); }
    else if (e.key === "Enter") { if (sel >= 0) { lis[sel].dispatchEvent(new Event("mousedown", { bubbles: true })); e.preventDefault(); } return; }
    lis.forEach((l, i) => l.classList.toggle("sel", i === sel)); });
  inp.addEventListener("blur", () => setTimeout(() => ac.classList.remove("show"), 160));
}
wireAC("#endInput", "#endAc", "#endCid"); wireAC("#startInput", "#startAc", "#startCid"); wireAC("#feedInput", "#feedAc", "#feedCid");
function resolveInput(id) { const inp = $(id); if (inp.dataset.cid) return inp.dataset.cid;
  const r = searchCpd(inp.value); if (r[0]) { inp.dataset.cid = r[0].cid; inp.value = r[0].name; return r[0].cid; } return ""; }

/* ---- example quick-picks (pinned to exact KEGG compound ids) ---- */
const EXAMPLES = [
  { p: ["succinate", "C00042"], s: ["pyruvate", "C00022"], f: ["", ""] },
  { p: ["L-lactate", "C00186"], s: ["pyruvate", "C00022"], f: ["", ""] },
  { p: ["acetoin", "C00466"], s: ["pyruvate", "C00022"], f: ["", ""] },
  { p: ["2,3-butanediol", "C03044"], s: ["pyruvate", "C00022"], f: ["acetate", "C00033"] },
];
$("#examples").innerHTML = "<span class='exlab'>try:</span>" + EXAMPLES.map((e, i) =>
  `<button class="ex" data-i="${i}">${e.p[0]} ← ${e.s[0]}</button>`).join("");
[...$("#examples").querySelectorAll(".ex")].forEach(b => b.onclick = () => { const e = EXAMPLES[+b.dataset.i];
  [["#endInput", "#endCid", e.p], ["#startInput", "#startCid", e.s], ["#feedInput", "#feedCid", e.f]].forEach(([iid, cid, v]) => {
    $(iid).value = v[0]; $(iid).dataset.cid = v[1]; $(cid).textContent = v[1]; });
  $("#query").dispatchEvent(new Event("submit")); });

$("#query").addEventListener("submit", e => { e.preventDefault();
  const end = resolveInput("#endInput"), start = resolveInput("#startInput"), feed = resolveInput("#feedInput");
  if (!start || !end) { $("#phase").textContent = "pick a valid product and start metabolite"; return; }
  if (!RUNNING) runLive(start, end, feed); });

/* ---- run the live in-browser search ---- */
async function runLive(start, end, feed) {
  RUNNING = true;
  ST = { routes: [], byId: {}, feas: {}, orgs: [], ep: null, query: { start, end, feed } };
  map.reset(true);
  $("#intro").classList.add("hidden"); $("#scrollcue").classList.add("hidden");
  $("#orglist").innerHTML = ""; $("#orgCount").textContent = "0"; $("#routeCount").textContent = "0"; $("#feasCount").textContent = "0";
  $("#results").classList.add("hidden"); $("#drawer").classList.add("hidden");
  $("#runBtn").disabled = true; $("#mapstatus").classList.remove("hidden"); $("#mapstatus").textContent = "searching…";
  try {
    await engine.run(start, end, feed || null, (t, d) => handleEvent(t, d), { skipGating: !GENOME_READY, maxLen: 13, maxRoutes: 60 });
  } catch (err) { $("#phase").textContent = "✕ " + (err.message || err); console.error(err); }
  RUNNING = false; $("#runBtn").disabled = false;
}

function handleEvent(event, data) {
  if (event === "phase") { $("#phase").textContent = "▸ " + data.msg; $("#pbar").style.width = (data.pct || 0) + "%"; }
  else if (event === "endpoints") { ST.ep = data; map.setEndpoints(data.start, data.end);
    $("#mapstatus").textContent = `tracing  ${data.end.name}  →  ${data.start.name}`; }
  else if (event === "explore") { /* map trace runs on the 'routes' event */ }
  else if (event === "routes") { ST.routes = data.routes; data.routes.forEach(r => ST.byId[r.id] = r);
    ST.capped = data.capped;
    $("#routeCount").textContent = data.capped ? data.n_routes + "+" : data.n_routes;
    $("#mapstatus").textContent = `connecting  ${ST.ep.end.name} → ${ST.ep.start.name}  via ${data.n_routes} real pathways`;
    map.traceRoutes(data.routes, ST.ep.end.name, !!(ST.ep.end && ST.ep.end.xy)); }
  else if (event === "thermo") { ST.feas[data.route_id] = data; $("#feasCount").textContent = Object.values(ST.feas).filter(x => x.feasible).length; }
  else if (event === "route_genomes") { ST.routeGenomes = data; }
  else if (event === "organism") { addOrganism(data); }
  else if (event === "done") { ST.done = data;
    if (data.error) { $("#phase").textContent = "✕ " + data.error; $("#mapstatus").textContent = data.error; return; }
    if (data.genome_pending) { $("#phase").textContent = "✓ routes found — organism results pending (genome data uploading)";
      $("#pbar").style.width = "100%"; $("#mapstatus").textContent = `${data.n_routes} routes found`; renderPathwaysOnly(); return; }
    $("#phase").textContent = "✓ complete"; $("#pbar").style.width = "100%";
    $("#mapstatus").textContent = `${data.T2.toLocaleString()} species · ${data.n_routes} routes`; renderResults(); }
}

/* results with routes+pathways only (genome bundle not uploaded yet) */
function renderPathwaysOnly() {
  $("#rtitle").innerHTML = `<i>${ST.ep.end.name}</i> <span style="color:var(--dim)">from</span> ${ST.ep.start.name}`;
  $("#herorow").innerHTML = `<div class="hero glass"><div class="ribbon" style="background:var(--amber)">${ST.routes.length} ROUTES</div>
    <div class="hero-head">${ST.ep.end.name}<span>from ${ST.ep.start.name}</span></div>
    <div class="hero-num">${ST.routes.length}<span>native routes found · shortest ${ST.done.shortest} steps</span></div>
    <div class="hero-route" style="color:var(--amber)">Which organisms encode these routes is computed once the genome data finishes uploading — the routes and map are live now.</div></div>`;
  COLS = [...ST.routes].sort((a, b) => a.length - b.length); COLIDX = {}; COLS.forEach((r, i) => COLIDX[r.id] = i);
  renderPathways(ST.routes);
  document.querySelector(".catalog").style.display = "none"; $("#funnelbar").style.display = "none";
  $("#results").classList.remove("hidden"); $("#scrollcue").classList.remove("hidden");
}

/* ---- live organism panel ---- */
function addOrganism(o) {
  ST.orgs.push(o);
  $("#orgCount").textContent = ST.orgs.length.toLocaleString();
  if ($("#orglist").children.length >= 200) return;
  const li = document.createElement("li");
  const fb = o.thermo_feasible ? `<span class="badge feas" title="a route is thermodynamically feasible">ΔG ✓</span>`
                               : `<span class="badge infeas" title="no feasible route direction">ΔG ✕</span>`;
  const feed = o.feedstock === "overflow_capable" ? ` <span class="badge infeas" title="overflow only, not uptake">overflow</span>` : "";
  li.innerHTML = `<span class="gdot" style="background:${GCOL[o.gram] || GCOL.Other};box-shadow:0 0 7px ${GCOL[o.gram] || GCOL.Other}"></span>
    <span class="sp">${o.species}</span><span class="badge routes">${o.n_routes}×</span>${fb}${feed}`;
  li.onclick = () => openOrganism(o);
  $("#orglist").appendChild(li);
}

/* ---- results: hero card(s) + species×routes catalog + funnel bar (mirrors bioconversion_overview) ---- */
let COLS = [], COLIDX = {};
function renderResults() {
  const d = ST.done;
  document.querySelector(".catalog").style.display = ""; $("#funnelbar").style.display = "";
  $("#rtitle").innerHTML = `<i>${ST.ep.end.name}</i> <span style="color:var(--dim)">from</span> ${ST.ep.start.name}`;
  const g = d.gram, tot = (g.Gpos + g.Gneg + g.Arch + g.Other) || 1;
  const gramBar = ["Gpos", "Gneg", "Arch", "Other"].filter(k => g[k]).map(k =>
    `<span class="gseg" title="${GNAME[k]} ${g[k]}" style="flex:${g[k]};background:${GCOL[k]}"></span>`).join("");
  const gramLeg = ["Gpos", "Gneg", "Arch", "Other"].filter(k => g[k]).map(k =>
    `<span><span class="gdot" style="background:${GCOL[k]};display:inline-block;margin-right:5px"></span>${GNAME[k]} <b>${g[k]}</b></span>`).join("");
  const shortest = ST.routes.reduce((a, b) => b.length < a.length ? b : a, ST.routes[0]);
  const enzymes = shortest ? [...new Set(shortest.steps.map(s => s.enzymes).filter(Boolean))].join(" → ") : "";
  const nFeas = Object.values(ST.feas).filter(x => x.feasible).length;

  // HERO ROW: big product card + stat cards
  $("#herorow").innerHTML = `
    <div class="hero glass">
      <div class="ribbon" style="${d.T2 > 0 ? "" : "background:var(--amber)"}">${d.T2 > 0 ? "NATIVE ROUTE" : (ST.routes.length ? "NO GENOME-ENCODED ROUTE" : "NO ROUTE FOUND")}</div>
      <div class="hero-head">${ST.ep.end.name}<span>from ${ST.ep.start.name}</span></div>
      <div class="hero-num">${d.T2.toLocaleString()}<span>${d.T2 > 0 ? "prokaryote species encode a full native route" : (ST.routes.length ? "genomes encode a route within the search horizon — carbon-skeleton routes exist but the native path may be longer (see below)" : "carbon-skeleton routes within the search horizon")}</span></div>
      <div class="hero-gram"><div class="gbar">${gramBar}</div><div class="gleg">${gramLeg}</div></div>
      <div class="hero-route">committed route&nbsp; <b>${enzymes || "—"}</b></div>
    </div>
    ${statCard(ST.routes.length, "native routes", `shortest ${d.shortest} steps`, "#2f7bff")}
    ${statCard(nFeas, "thermo-feasible routes", "eQuilibrator ΔG", "#22d18c")}
    ${d.T3 != null ? statCard(d.T3, "can take up the feedstock", `${d.overflow_excluded} overflow-only excluded`, "#39c0ff") : ""}`;

  // CATALOG columns = routes sorted by length
  COLS = [...ST.routes].sort((a, b) => a.length - b.length); COLIDX = {};
  COLS.forEach((r, i) => COLIDX[r.id] = i);
  renderPathways(ST.routes);
  renderCatalog("");
  $("#spFilter").oninput = e => renderCatalog(e.target.value.toLowerCase());
  $("#scrollcue").classList.remove("hidden");

  // FUNNEL BAR (clickable tiers)
  const tiers = [["terminal enzyme", d.T0, "#7d8ca6", "has only the LAST enzyme — overcounts"],
                 ["encodes full route", d.T2, "#2f7bff", "complete native route (honest headline)"]];
  if (d.T3 != null) tiers.push(["+ feedstock uptake", d.T3, "#22d18c", "can also take up the feedstock"]);
  const mx = Math.max(...tiers.map(t => t[1]), 1);
  $("#funnelbar").innerHTML = `<div class="fbtitle">having the last enzyme ≠ having the pathway</div>` +
    tiers.map(([lab, v, c], i) => `<div class="fseg" data-tier="${i}" style="--c:${c}">
      <div class="fbar" style="width:${Math.max(30, 190 * v / mx)}px"></div>
      <div class="fnum">${v.toLocaleString()}</div><div class="flab">${lab}</div></div>`).join(
      `<div class="fsep">▸</div>`);
  [...document.querySelectorAll(".fseg")].forEach(el => el.onclick = () => openTier(+el.dataset.tier, tiers[+el.dataset.tier]));
  $("#results").classList.remove("hidden");
}
function statCard(n, label, sub, c) {
  return `<div class="statcard glass" style="--c:${c}"><div class="sc-num">${(+n).toLocaleString()}</div>
    <div class="sc-lab">${label}</div><div class="sc-sub">${sub}</div></div>`;
}
function renderCatalog(filter) {
  const rows = ST.orgs.filter(o => o.species.toLowerCase().includes(filter))
    .sort((a, b) => b.n_routes - a.n_routes).slice(0, 400);
  const header = `<div class="cat-row cat-hd">
    <span class="c-sp">species (${rows.length})</span>
    <span class="c-badges">Gram · ΔG · feed</span>
    <span class="c-strip">routes by length →</span></div>`;
  $("#catalog").innerHTML = header + rows.map(o => {
    const set = new Set(o.route_idx || []);
    const strip = COLS.map(r => {
      const has = set.has(r.id);
      const f = ST.feas[r.id];
      const col = !has ? "transparent" : (f ? (f.feasible ? "var(--green)" : "var(--red)") : "var(--blue)");
      return `<span class="cell${has ? " on" : ""}" title="route ${r.length} steps${has ? "" : " — not encoded"}"
        style="background:${col}"></span>`;
    }).join("");
    return `<div class="cat-row" data-i="${ST.orgs.indexOf(o)}">
      <span class="c-sp"><span class="gdot" style="background:${GCOL[o.gram] || GCOL.Other}"></span><i>${o.species}</i></span>
      <span class="c-badges"><span class="badge routes">${o.n_routes}×</span>
        ${o.thermo_feasible ? '<span class="badge feas">✓</span>' : '<span class="badge infeas">✕</span>'}
        ${o.feedstock === "overflow_capable" ? '<span class="badge infeas" title="overflow only">ov</span>' :
          o.feedstock === "uptake" ? '<span class="badge feas" title="feedstock uptake">up</span>' : ''}</span>
      <span class="c-strip">${strip}</span></div>`;
  }).join("");
  [...$("#catalog").querySelectorAll(".cat-row:not(.cat-hd)")].forEach(r =>
    r.onclick = () => openOrganism(ST.orgs[+r.dataset.i]));
}

/* ---- distinct pathways: ranked by length, click opens a full detail PAGE ---- */
function renderPathways(routes) {
  const G = ST.routeGenomes || null;
  // when genomes are known, show only routes a genome actually encodes (drops spurious
  // carbon-skeleton shortcuts), ranked by length then by #species; else fall back to length.
  const encoded = G ? routes.filter(r => (G[r.id] || 0) > 0) : routes;
  const noneEncoded = G && encoded.length === 0;
  const useList = encoded.length ? encoded : routes;
  // rank: MOST GENOME-ENCODED first (the tool's whole premise — the canonical pathway is the
  // one thousands of organisms actually carry), then fewest obscure reactions, then shorter.
  // Length is a weak final tiebreak, NOT primary — that is what buried real long pathways.
  const list = [...useList].sort((a, b) =>
    (a.repeats || 0) - (b.repeats || 0) ||
    (G ? (G[b.id] || 0) - (G[a.id] || 0) : 0) ||
    (a.uncur || 0) - (b.uncur || 0) || a.length - b.length || a.id - b.id);
  if (noneEncoded) {
    $("#pwTitle").innerHTML = `Distinct pathways · <b>0</b> encoded by any KEGG genome
      <span class="hint">— ${routes.length} carbon-skeleton route(s) exist but no sequenced prokaryote encodes one within the ${ST.done ? ST.done.shortest : ""}-step search horizon. The native pathway may be longer than the search depth, or split across organisms. The skeleton routes are shown below for inspection.</span>`;
  } else {
    const hidden = G ? routes.length - encoded.length : 0;
    const note = G
      ? (hidden > 0 ? `<span class="hint">· ${hidden} carbon-skeleton route(s) no sequenced genome encodes are hidden</span>` : "")
      : (ST.capped ? `<span class="hint">(search cap — more exist)</span>` : "");
    $("#pwTitle").innerHTML = `Distinct pathways · <b>${list.length}${(!G && ST.capped) ? "+" : ""}</b>${G ? " encoded by genomes" : " found"} · ranked ${note}
      <span class="hint">click a pathway to open its full report (structures · reactions · cross-DB directionality · species)</span>`;
  }
  $("#pwList").innerHTML = list.map((r, rank) => {
    const f = ST.feas[r.id];
    const feas = f ? (f.feasible ? `<span class="badge feas">ΔG feasible ${f.dG_sum ?? ""}</span>`
                                 : `<span class="badge infeas">ΔG infeasible</span>`) : "";
    const gb = G ? `<span class="badge routes">${(G[r.id] || 0).toLocaleString()} species</span>` : "";
    const chain = r.path.map((p, i) => (i === 0 ? "" :
      `<span class="ar">→</span><span class="enz">${r.steps[i - 1].enzymes || "?"}</span><span class="ar">→</span>`)
      + `<span class="met">${p.name}</span>`).join(" ");
    const rx = r.steps.map(s => s.reactions[0].rid).join(" · ");
    return `<div class="pw" data-id="${r.id}"><div class="pw-top">
        <span><b class="pwrank">#${rank + 1}</b> · ${r.length} steps${rank === 0 ? " · top-ranked" : ""}${r.repeats ? ` <span class="hint" title="a reaction is reused across steps">· reuses a reaction</span>` : ""}</span>
        ${gb}${feas}<span class="pwopen">open report ↗</span></div>
      <div class="chain">${chain}</div><div class="rx">${rx}</div></div>`;
  }).join("");
  window.PW_ORDER = list.map(r => r.id);          // ranked order for prev/next on the report page
  [...$("#pwList").children].forEach(el => el.onclick = () => openRoutePage(ST.byId[+el.dataset.id]));
}

/* store route context + open the dedicated report page IN THE SAME TAB, so history.back()
   (and bfcache) return the user to the exact results view. A compact list of all ranked
   pathways travels along so the report can offer prev/next navigation without a round-trip. */
function openRoutePage(r) {
  const ranked = (window.PW_ORDER || ST.routes.map(x => x.id));
  const routes = ranked.map(id => ST.byId[id]).filter(Boolean);   // full objects, ranked
  const feas = {}; ranked.forEach(id => { if (ST.feas[id]) feas[id] = ST.feas[id]; });
  const genomes = ST.routeGenomes || {};
  // compact org list (species that encode ≥1 route) so the report can filter per pathway locally
  const orgs = ST.orgs.map(o => ({ species: o.species, gram: o.gram, domain: o.domain,
    n_routes: o.n_routes, feedstock: o.feedstock, route_idx: o.route_idx || [] }));
  sessionStorage.setItem("panroute_route", JSON.stringify({
    routes, feas, genomes, orgs, current: r.id,
    query: { start: ST.ep.start, end: ST.ep.end },
    genome_pending: !!(ST.done && ST.done.genome_pending),
  }));
  saveSession();                       // so Back restores results even if bfcache misses
  sessionStorage.setItem("panroute_return", "1");
  location.href = "route.html";
}

/* prefill the two inputs from stored compound ids (name looked up in the bundle) */
function prefillFromCids(start, end, feed) {
  const nm = cid => (CPD.find(x => x.cid === cid) || {}).name || cid;
  [["#endInput", "#endCid", end], ["#startInput", "#startCid", start], ["#feedInput", "#feedCid", feed]]
    .forEach(([iid, cidEl, c]) => { const el = $(iid);
      el.value = c ? nm(c) : ""; el.dataset.cid = c || ""; $(cidEl).textContent = c || ""; });
}

/* restore the results view when the user comes BACK from the route report (or reloads).
   bfcache usually restores the live DOM for free; this is the fallback when it doesn't. */
window.addEventListener("pageshow", ev => {
  if (ev.persisted) return;                       // bfcache already restored everything
  const ret = sessionStorage.getItem("panroute_return");
  const last = sessionStorage.getItem("panroute_last");
  if (ret && last && !RUNNING) {
    sessionStorage.removeItem("panroute_return");
    const q = JSON.parse(last);
    prefillFromCids(q.start, q.end, q.feed);
    if (CPD.length) runLive(q.start, q.end, q.feed || "");
    else document.addEventListener("cpd-ready", () => runLive(q.start, q.end, q.feed || ""), { once: true });
  }
});

/* persist the active query so the results view can be rebuilt on return / reload */
function saveSession() {
  if (!ST.ep) return;
  sessionStorage.setItem("panroute_last", JSON.stringify({
    start: ST.ep.start.cid, end: ST.ep.end.cid, feed: (ST.query && ST.query.feed) || "" }));
}

/* ---- drawer ---- */
function openOrganism(o) {
  $("#drawer").classList.remove("hidden");
  $("#dtitle").innerHTML = `<i>${o.species}</i>`;
  const routes = (o.route_idx || []).map(i => ST.byId[i]).filter(Boolean);
  $("#dbody").innerHTML =
    `<div style="color:var(--dim);font-size:12px;margin-bottom:10px">${GNAME[o.gram] || "—"} · ${o.domain} ·
      encodes <b>${o.n_routes}</b> native route(s) · feedstock: <b>${o.feedstock}</b>
      <span class="hint">· click a route for its full report</span></div>` +
    (routes.map(routeBox).join("") || `<p class="hint">route detail unavailable</p>`);
  [...$("#dbody").querySelectorAll(".routebox")].forEach(el =>
    el.onclick = () => { const rr = ST.byId[+el.dataset.id]; if (rr) openRoutePage(rr); });
}
function routeBox(r) {
  const f = ST.feas[r.id];
  const feas = f ? (f.feasible ? `<span class="badge feas">ΔG feasible ${f.dG_sum ?? ""} kJ/mol</span>` : `<span class="badge infeas">ΔG infeasible</span>`) : "";
  const chain = r.path.map((p, i) => (i === 0 ? "" : `<span class="ar">→</span><span class="enz">${r.steps[i - 1].enzymes || "?"}</span><span class="ar">→</span>`) + `<span class="met">${p.name}</span>`).join(" ");
  const rxns = r.steps.map(s => s.reactions[0].rid).join(" · ");
  return `<div class="routebox" data-id="${r.id}" style="cursor:pointer"><div class="rtop"><span>route · ${r.length} steps</span>${feas}<span class="pwopen">open report ↗</span></div>
    <div class="chain">${chain}</div><div class="rx" style="margin-top:6px">${rxns}</div></div>`;
}
function openTier(i, t) {
  $("#drawer").classList.remove("hidden"); $("#dtitle").textContent = t[0];
  const txt = ["Species that carry ONLY the last enzyme of the route — this overcounts, because having the final step says nothing about whether the cell can reach the precursor. This is the metric a naïve KEGG survey reports.",
    "Species that encode a COMPLETE native route from the start metabolite to the product. This is the honest headline count.",
    "Species that encode a full route AND can take up the feedstock in the consuming direction (overflow-only carriers excluded)."][i];
  $("#dbody").innerHTML = `<p style="font-size:13.5px;color:var(--ink)">${txt}</p>
    <p style="font-size:26px;font-weight:700;color:var(--neon)">${t[1].toLocaleString()} species</p>`;
}
$("#closeResults").onclick = () => $("#results").classList.add("hidden");
$("#closeDrawer").onclick = () => $("#drawer").classList.add("hidden");

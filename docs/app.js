/* PanRoute static site (github.io): loads a precomputed run bundle and REPLAYS the exact
   live experience client-side — real KEGG-map trace, live organism discovery, thermo
   feasibility, interactive results. No backend. Data is KEGG-derived (research use). */
const $ = s => document.querySelector(s);
const GCOL = { Gpos: "#3a8bff", Gneg: "#ff8a3a", Arch: "#22d18c", Other: "#8aa0bf" };
const GNAME = { Gpos: "Gram-positive", Gneg: "Gram-negative", Arch: "Archaea", Other: "other" };
const map = new MapView($("#map"));
let ST = null, RUNNING = false;
const sleep = ms => new Promise(r => setTimeout(r, ms));

/* ---- load run index ---- */
async function loadIndex() {
  try {
    const idx = await fetch("runs/index.json").then(r => r.json());
    const sel = $("#runSelect");
    idx.runs.forEach(r => {
      const o = document.createElement("option");
      o.value = r.slug;
      o.textContent = `${r.end_name} ← ${r.start_name}  ·  ${(r.T2 || 0).toLocaleString()} species`;
      sel.appendChild(o);
    });
    if (idx.runs[0]) sel.value = idx.runs[0].slug;
  } catch (e) { $("#phase").textContent = "no runs published yet"; }
}
loadIndex();

$("#runBtn").addEventListener("click", () => {
  const slug = $("#runSelect").value;
  if (slug && !RUNNING) replay(slug);
});

/* ---- replay a bundle ---- */
async function replay(slug) {
  RUNNING = true;
  ST = { routes: [], byId: {}, feas: {}, orgs: [], ep: null };
  map.reset(true);
  $("#intro").classList.add("hidden");
  $("#orglist").innerHTML = ""; $("#orgCount").textContent = "0";
  $("#routeCount").textContent = "0"; $("#feasCount").textContent = "0";
  $("#results").classList.add("hidden"); $("#drawer").classList.add("hidden");
  $("#runBtn").disabled = true; $("#mapstatus").classList.remove("hidden");
  $("#mapstatus").textContent = "loading…";

  let bundle;
  try { bundle = await fetch(`runs/${slug}.json`).then(r => r.json()); }
  catch (e) { $("#phase").textContent = "failed to load run"; RUNNING = false; $("#runBtn").disabled = false; return; }

  let orgShown = 0;
  for (const { event, data } of bundle.events) {
    if (event === "phase") { $("#phase").textContent = "▸ " + data.msg; $("#pbar").style.width = (data.pct || 0) + "%"; await sleep(120); }
    else if (event === "endpoints") { ST.ep = data; map.setEndpoints(data.start, data.end);
      $("#mapstatus").textContent = `tracing  ${data.end.name}  →  ${data.start.name}`; await sleep(300); }
    else if (event === "explore") { /* the map trace runs on the 'routes' event via traceRoutes */ }
    else if (event === "routes") { ST.routes = data.routes; data.routes.forEach(r => ST.byId[r.id] = r);
      $("#routeCount").textContent = data.n_routes;
      $("#mapstatus").textContent = `connecting  ${ST.ep.end.name} → ${ST.ep.start.name}  via ${data.n_routes} real pathways`;
      await map.traceRoutes(data.routes, ST.ep.end.name, !!(ST.ep.end && ST.ep.end.xy)); }
    else if (event === "thermo") { ST.feas[data.route_id] = data;
      $("#feasCount").textContent = Object.values(ST.feas).filter(x => x.feasible).length; }
    else if (event === "organism") { addOrganism(data);
      orgShown++; if (orgShown < 80) await sleep(24); else if (orgShown % 8 === 0) await sleep(4); }
    else if (event === "done") { ST.done = data;
      $("#phase").textContent = "✓ complete"; $("#pbar").style.width = "100%";
      $("#mapstatus").textContent = `${data.T2.toLocaleString()} species · ${data.n_routes} routes`;
      renderResults(); }
  }
  RUNNING = false; $("#runBtn").disabled = false;
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
      <div class="ribbon">${d.T2 > 0 ? "NATIVE ROUTE" : "ENGINEERED ONLY"}</div>
      <div class="hero-head">${ST.ep.end.name}<span>from ${ST.ep.start.name}</span></div>
      <div class="hero-num">${d.T2.toLocaleString()}<span>prokaryote species encode a full native route</span></div>
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

/* ---- distinct pathways (clickable, highlight on map) ---- */
function renderPathways(routes) {
  const list = [...routes].sort((a, b) => a.length - b.length);
  $("#pwTitle").innerHTML = `Distinct pathways · <b>${routes.length}</b> found
    <span class="hint">click one to highlight it on the map above</span>`;
  $("#pwList").innerHTML = list.map(r => {
    const f = ST.feas[r.id];
    const feas = f ? (f.feasible ? `<span class="badge feas">ΔG feasible ${f.dG_sum ?? ""}</span>`
                                 : `<span class="badge infeas">ΔG infeasible</span>`) : "";
    const chain = r.path.map((p, i) => (i === 0 ? "" :
      `<span class="ar">→</span><span class="enz">${r.steps[i - 1].enzymes || "?"}</span><span class="ar">→</span>`)
      + `<span class="met">${p.name}</span>`).join(" ");
    const rx = r.steps.map(s => s.reactions[0].rid).join(" · ");
    return `<div class="pw" data-id="${r.id}"><div class="pw-top"><span>pathway ${r.id + 1} · ${r.length} steps</span>${feas}</div>
      <div class="chain">${chain}</div><div class="rx">${rx}</div></div>`;
  }).join("");
  [...$("#pwList").children].forEach(el => el.onclick = () => {
    map.highlightRoute(ST.byId[+el.dataset.id]);
    $("#mapwrap").scrollIntoView({ behavior: "smooth", block: "center" });
  });
}

/* ---- drawer ---- */
function openOrganism(o) {
  $("#drawer").classList.remove("hidden");
  $("#dtitle").innerHTML = `<i>${o.species}</i>`;
  const routes = (o.route_idx || []).map(i => ST.byId[i]).filter(Boolean);
  $("#dbody").innerHTML =
    `<div style="color:var(--dim);font-size:12px;margin-bottom:10px">${GNAME[o.gram] || "—"} · ${o.domain} ·
      encodes <b>${o.n_routes}</b> native route(s) · feedstock: <b>${o.feedstock}</b></div>` +
    (routes.map(routeBox).join("") || `<p class="hint">route detail unavailable</p>`);
}
function routeBox(r) {
  const f = ST.feas[r.id];
  const feas = f ? (f.feasible ? `<span class="badge feas">ΔG feasible ${f.dG_sum ?? ""} kJ/mol</span>` : `<span class="badge infeas">ΔG infeasible</span>`) : "";
  const chain = r.path.map((p, i) => (i === 0 ? "" : `<span class="ar">→</span><span class="enz">${r.steps[i - 1].enzymes || "?"}</span><span class="ar">→</span>`) + `<span class="met">${p.name}</span>`).join(" ");
  const rxns = r.steps.map(s => s.reactions[0].rid).join(" · ");
  return `<div class="routebox"><div class="rtop"><span>route · ${r.length} steps</span>${feas}</div>
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

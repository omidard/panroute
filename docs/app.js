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
    else if (event === "explore") { $("#mapstatus").textContent = `tracing route on KEGG map · step ${data.index + 1}/${data.total}`;
      await map.drawStep(data.step); await sleep(160); }
    else if (event === "routes") { ST.routes = data.routes; data.routes.forEach(r => ST.byId[r.id] = r);
      $("#routeCount").textContent = data.n_routes;
      const sh = data.routes.reduce((a, b) => b.length < a.length ? b : a, data.routes[0]); map.finalizeRoute(sh.map); }
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

/* ---- results ---- */
function renderResults() {
  const d = ST.done;
  $("#rtitle").innerHTML = `${ST.ep.end.name} <span style="color:var(--dim)">from</span> ${ST.ep.start.name}
    &nbsp;·&nbsp; native capacity across KEGG genomes`;
  const tiers = [["terminal enzyme (last step)", d.T0, "#7d8ca6"], ["encodes a full route", d.T2, "#2f7bff"]];
  if (d.T3 != null) tiers.push(["+ feedstock uptake", d.T3, "#22d18c"]);
  const max = Math.max(...tiers.map(t => t[1]), 1);
  $("#funnel").innerHTML = tiers.map(([lab, v, c], i) =>
    `<div class="tierbar" data-tier="${i}"><div class="lab">${lab}</div>
      <div class="track"><div class="fillb" style="width:${Math.max(6, 100 * v / max)}%;
        background:linear-gradient(90deg,${c},${c}cc)">${v.toLocaleString()}</div></div></div>`).join("");
  const g = d.gram, tot = (g.Gpos + g.Gneg + g.Arch + g.Other) || 1;
  let acc = 0; const segs = ["Gpos", "Gneg", "Arch", "Other"].filter(k => g[k]).map(k => { const f = g[k] / tot, a0 = acc; acc += f; return { k, f, a0 }; });
  const R = 62, C = 2 * Math.PI * R;
  $("#donut").innerHTML = `<div class="donutwrap"><svg width="150" height="150" viewBox="0 0 150 150">
    ${segs.map(s => `<circle cx="75" cy="75" r="${R}" fill="none" stroke="${GCOL[s.k]}" stroke-width="20"
      stroke-dasharray="${s.f * C} ${C}" stroke-dashoffset="${-s.a0 * C}" transform="rotate(-90 75 75)"/>`).join("")}
    <text x="75" y="80" text-anchor="middle" fill="#dce6f5" font-size="22" font-weight="700">${d.T2}</text></svg>
    <div>${segs.map(s => `<div style="margin:4px 0"><span class="gdot" style="background:${GCOL[s.k]};display:inline-block;margin-right:6px"></span>${GNAME[s.k]} · <b>${g[s.k]}</b></div>`).join("")}</div></div>`;
  renderSpecies(""); $("#spFilter").oninput = e => renderSpecies(e.target.value.toLowerCase());
  $("#results").classList.remove("hidden");
  document.querySelectorAll(".tierbar").forEach(t => t.onclick = () => openTier(+t.dataset.tier, tiers[+t.dataset.tier]));
}
function renderSpecies(filter) {
  const rows = ST.orgs.filter(o => o.species.toLowerCase().includes(filter)).sort((a, b) => b.n_routes - a.n_routes);
  $("#sptable").innerHTML = rows.slice(0, 300).map(o =>
    `<li data-i="${ST.orgs.indexOf(o)}"><span class="gdot" style="background:${GCOL[o.gram] || GCOL.Other}"></span>
      <span class="sp">${o.species}</span><span class="badge routes">${o.n_routes} routes</span>
      ${o.thermo_feasible ? '<span class="badge feas">ΔG ✓</span>' : '<span class="badge infeas">ΔG ✕</span>'}</li>`).join("");
  [...$("#sptable").children].forEach(li => li.onclick = () => openOrganism(ST.orgs[+li.dataset.i]));
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

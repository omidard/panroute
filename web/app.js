/* PanRoute frontend controller: resolves metabolites, streams the live run over SSE,
   drives the KEGG-map animation, the live discovery panel, and the interactive results. */
const $ = s => document.querySelector(s);
const GCOL = { Gpos: "#3a8bff", Gneg: "#ff8a3a", Arch: "#22d18c", Other: "#8aa0bf" };
const GNAME = { Gpos: "Gram-positive", Gneg: "Gram-negative", Arch: "Archaea", Other: "other" };
const map = new MapView($("#map"));
let ST = null, evtSrc = null;

/* ---------------- autocomplete ---------------- */
function wireAC(inputId, acId, cidId) {
  const inp = $(inputId), ac = $(acId), cid = $(cidId);
  let t, items = [], sel = -1;
  inp.dataset.cid = "";
  inp.addEventListener("input", () => {
    cid.textContent = ""; inp.dataset.cid = "";
    clearTimeout(t);
    const q = inp.value.trim(); if (q.length < 2) { ac.classList.remove("show"); return; }
    t = setTimeout(async () => {
      const r = await fetch("/api/resolve?q=" + encodeURIComponent(q)).then(r => r.json()).catch(() => []);
      items = r; sel = -1;
      ac.innerHTML = r.map((x, i) => `<li data-i="${i}">${x.name}<span class="c">${x.cid}</span></li>`).join("");
      ac.classList.toggle("show", r.length > 0);
    }, 220);
  });
  ac.addEventListener("click", e => {
    const li = e.target.closest("li"); if (!li) return;
    const x = items[+li.dataset.i];
    inp.value = x.name; inp.dataset.cid = x.cid; cid.textContent = x.cid; ac.classList.remove("show");
  });
  inp.addEventListener("keydown", e => {
    if (!ac.classList.contains("show")) return;
    const lis = [...ac.children];
    if (e.key === "ArrowDown") { sel = Math.min(sel + 1, lis.length - 1); e.preventDefault(); }
    else if (e.key === "ArrowUp") { sel = Math.max(sel - 1, 0); e.preventDefault(); }
    else if (e.key === "Enter") { if (sel >= 0) { lis[sel].click(); e.preventDefault(); } return; }
    lis.forEach((l, i) => l.classList.toggle("sel", i === sel));
  });
  inp.addEventListener("blur", () => setTimeout(() => ac.classList.remove("show"), 180));
}
wireAC("#endInput", "#endAc", "#endCid");
wireAC("#startInput", "#startAc", "#startCid");
wireAC("#feedInput", "#feedAc", "#feedCid");

async function ensureCid(inputId) {
  const inp = $(inputId);
  if (inp.dataset.cid) return inp.dataset.cid;
  const q = inp.value.trim(); if (!q) return "";
  const r = await fetch("/api/resolve?q=" + encodeURIComponent(q)).then(r => r.json()).catch(() => []);
  if (r[0]) { inp.dataset.cid = r[0].cid; inp.value = r[0].name; return r[0].cid; }
  return "";
}

/* ---------------- run ---------------- */
$("#query").addEventListener("submit", async e => {
  e.preventDefault();
  const end = await ensureCid("#endInput"), start = await ensureCid("#startInput");
  const feed = await ensureCid("#feedInput");
  if (!start || !end) { $("#phase").textContent = "pick a valid product and start metabolite"; return; }
  startRun(start, end, feed);
});

let animQ = [], animRunning = false;
async function pump() {
  if (animRunning) return; animRunning = true;
  while (animQ.length) { const s = animQ.shift(); await map.drawStep(s); }
  animRunning = false;
}

function startRun(start, end, feed) {
  if (evtSrc) evtSrc.close();
  ST = { routes: [], byId: {}, feas: {}, orgs: [], start, end, feed, done: null };
  map.reset(true);
  $("#orglist").innerHTML = ""; $("#orgCount").textContent = "0";
  $("#routeCount").textContent = "0"; $("#feasCount").textContent = "0";
  $("#results").classList.add("hidden"); $("#drawer").classList.add("hidden");
  $("#runBtn").disabled = true; $("#mapstatus").classList.remove("hidden");
  $("#mapstatus").textContent = "initialising…";
  animQ = [];

  const url = `/api/run?start=${start}&end=${end}` + (feed ? `&feedstock=${feed}` : "");
  evtSrc = new EventSource(url);

  evtSrc.addEventListener("phase", e => { const d = JSON.parse(e.data);
    $("#phase").textContent = "▸ " + d.msg; $("#pbar").style.width = (d.pct || 0) + "%"; });

  evtSrc.addEventListener("endpoints", e => { const d = JSON.parse(e.data);
    ST.ep = d; map.setEndpoints(d.start, d.end);
    $("#mapstatus").textContent = `tracing  ${d.end.name}  →  ${d.start.name}`; });

  evtSrc.addEventListener("explore", e => { const d = JSON.parse(e.data);
    animQ.push(d.step); pump();
    $("#mapstatus").textContent = `tracing route on KEGG map · step ${d.index + 1}/${d.total}`; });

  evtSrc.addEventListener("routes", e => { const d = JSON.parse(e.data);
    ST.routes = d.routes; d.routes.forEach(r => ST.byId[r.id] = r);
    $("#routeCount").textContent = d.n_routes;
    const shortest = d.routes.reduce((a, b) => b.length < a.length ? b : a, d.routes[0]);
    map.finalizeRoute(shortest.map); });

  evtSrc.addEventListener("thermo", e => { const d = JSON.parse(e.data);
    ST.feas[d.route_id] = d;
    const n = Object.values(ST.feas).filter(x => x.feasible).length;
    $("#feasCount").textContent = n; });

  evtSrc.addEventListener("organism", e => { const d = JSON.parse(e.data); addOrganism(d); });

  evtSrc.addEventListener("done", e => { ST.done = JSON.parse(e.data);
    $("#phase").textContent = "✓ complete"; $("#pbar").style.width = "100%";
    $("#mapstatus").textContent = `${ST.done.T2} species · ${ST.done.n_routes} routes`;
    $("#runBtn").disabled = false; renderResults(); });

  evtSrc.addEventListener("error", e => { try { const d = JSON.parse(e.data);
    $("#phase").textContent = "✕ " + (d.message || "error"); } catch (_) {}
    $("#runBtn").disabled = false; });
  evtSrc.addEventListener("close", () => { evtSrc.close(); $("#runBtn").disabled = false; });
}

/* ---------------- live organism panel ---------------- */
function addOrganism(o) {
  ST.orgs.push(o);
  $("#orgCount").textContent = ST.orgs.length;
  const li = document.createElement("li");
  const feasBadge = o.thermo_feasible
    ? `<span class="badge feas" title="a route is thermodynamically feasible">ΔG ✓</span>`
    : `<span class="badge infeas" title="no feasible route direction">ΔG ✕</span>`;
  const feed = o.feedstock === "overflow_capable"
    ? ` <span class="badge infeas" title="overflow only, not uptake">overflow</span>` : "";
  li.innerHTML = `<span class="gdot" style="background:${GCOL[o.gram] || GCOL.Other};
      box-shadow:0 0 7px ${GCOL[o.gram] || GCOL.Other}"></span>
    <span class="sp">${o.species}</span>
    <span class="badge routes">${o.n_routes}×</span>${feasBadge}${feed}`;
  li.onclick = () => openOrganism(o);
  const ul = $("#orglist");
  if (ul.children.length < 200) ul.appendChild(li);
}

/* ---------------- results ---------------- */
function renderResults() {
  const d = ST.done;
  $("#rtitle").innerHTML = `${ST.ep.end.name} <span style="color:var(--dim)">from</span> ${ST.ep.start.name}
    &nbsp;·&nbsp; native capacity across KEGG genomes`;
  // funnel
  const tiers = [["terminal enzyme (last step)", d.T0, "#7d8ca6"],
                 ["encodes a full route", d.T2, "#2f7bff"]];
  if (d.T3 != null) tiers.push(["+ feedstock uptake", d.T3, "#22d18c"]);
  const max = Math.max(...tiers.map(t => t[1]), 1);
  $("#funnel").innerHTML = tiers.map(([lab, v, c], i) =>
    `<div class="tierbar" data-tier="${i}"><div class="lab">${lab}</div>
      <div class="track"><div class="fillb" style="width:${Math.max(6, 100 * v / max)}%;
        background:linear-gradient(90deg,${c},${c}cc)">${v.toLocaleString()}</div></div></div>`).join("");
  // donut
  const g = d.gram; const tot = (g.Gpos + g.Gneg + g.Arch + g.Other) || 1;
  let acc = 0; const segs = ["Gpos", "Gneg", "Arch", "Other"].filter(k => g[k]).map(k => {
    const frac = g[k] / tot, a0 = acc; acc += frac; return { k, frac, a0 }; });
  const R = 62, C = 2 * Math.PI * R;
  $("#donut").innerHTML = `<div class="donutwrap"><svg width="150" height="150" viewBox="0 0 150 150">
    ${segs.map(s => `<circle cx="75" cy="75" r="${R}" fill="none" stroke="${GCOL[s.k]}" stroke-width="20"
       stroke-dasharray="${s.frac * C} ${C}" stroke-dashoffset="${-s.a0 * C}"
       transform="rotate(-90 75 75)"/>`).join("")}
    <text x="75" y="80" text-anchor="middle" fill="#dce6f5" font-size="22" font-weight="700">${d.T2}</text></svg>
    <div>${segs.map(s => `<div style="margin:4px 0"><span class="gdot" style="background:${GCOL[s.k]};
       display:inline-block;margin-right:6px"></span>${GNAME[s.k]} · <b>${g[s.k]}</b></div>`).join("")}</div></div>`;
  // species table
  renderSpecies("");
  $("#spFilter").oninput = e => renderSpecies(e.target.value.toLowerCase());
  $("#results").classList.remove("hidden");
  document.querySelectorAll(".tierbar").forEach(t => t.onclick = () =>
    openTier(+t.dataset.tier, tiers[+t.dataset.tier]));
}

function renderSpecies(filter) {
  const rows = ST.orgs.filter(o => o.species.toLowerCase().includes(filter))
    .sort((a, b) => b.n_routes - a.n_routes);
  $("#sptable").innerHTML = rows.slice(0, 250).map((o, i) =>
    `<li data-i="${ST.orgs.indexOf(o)}"><span class="gdot" style="background:${GCOL[o.gram] || GCOL.Other}"></span>
      <span class="sp">${o.species}</span>
      <span class="badge routes">${o.n_routes} routes</span>
      ${o.thermo_feasible ? '<span class="badge feas">ΔG ✓</span>' : '<span class="badge infeas">ΔG ✕</span>'}</li>`).join("");
  [...$("#sptable").children].forEach(li => li.onclick = () => openOrganism(ST.orgs[+li.dataset.i]));
}

/* ---------------- detail drawer ---------------- */
function openOrganism(o) {
  const dr = $("#drawer"); dr.classList.remove("hidden");
  $("#dtitle").innerHTML = `<i>${o.species}</i>`;
  const routes = (o.route_idx || []).map(i => ST.byId[i]).filter(Boolean);
  $("#dbody").innerHTML =
    `<div style="color:var(--dim);font-size:12px;margin-bottom:10px">
      ${GNAME[o.gram] || "—"} · ${o.domain} · encodes <b>${o.n_routes}</b> native route(s) ·
      feedstock: <b>${o.feedstock}</b></div>` +
    routes.map(r => routeBox(r)).join("") ||
    `<p class="hint">route detail unavailable</p>`;
}
function routeBox(r) {
  const f = ST.feas[r.id];
  const feas = f ? (f.feasible ? `<span class="badge feas">ΔG feasible ${f.dG_sum ?? ""} kJ/mol</span>`
                               : `<span class="badge infeas">ΔG infeasible</span>`) : "";
  const chain = r.path.map((p, i) => {
    const enz = i < r.steps.length ? `<span class="ar">→</span><span class="enz">${r.steps[i].enzymes || "?"}</span><span class="ar">→</span>` : "";
    return (i === 0 ? "" : enz) + `<span class="met">${p.name}</span>`;
  }).join(" ");
  const rxns = r.steps.map(s => s.reactions[0].rid).join(" · ");
  return `<div class="routebox"><div class="rtop"><span>route · ${r.length} steps</span>${feas}</div>
    <div class="chain">${chain}</div><div class="rx" style="margin-top:6px">${rxns}</div></div>`;
}
function openTier(i, t) {
  const dr = $("#drawer"); dr.classList.remove("hidden");
  $("#dtitle").textContent = t[0];
  const txt = ["Species that carry ONLY the last enzyme of the route — this overcounts, because having the final step says nothing about whether the cell can reach the precursor. This is the metric a naïve KEGG survey reports.",
    "Species that encode a COMPLETE native route from the start metabolite to the product. This is the honest headline count.",
    "Species that encode a full route AND can take up the feedstock in the consuming direction (overflow-only carriers are excluded)."][i];
  $("#dbody").innerHTML = `<p style="font-size:13.5px;color:var(--ink)">${txt}</p>
    <p style="font-size:26px;font-weight:700;color:var(--neon)">${t[1].toLocaleString()} species</p>`;
}
$("#closeResults").onclick = () => $("#results").classList.add("hidden");
$("#closeDrawer").onclick = () => $("#drawer").classList.add("hidden");

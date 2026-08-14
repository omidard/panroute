/* PanRoute route report: a modularised metabolic map. Reactions are grouped into subsystem
   modules (KEGG pathway); each module is a box showing its metabolites (drawn from SMILES)
   connected by enzyme-labelled reaction arrows. Below: per-reaction cross-database
   directionality (KEGG/Rhea/MetaCyc) + our eQuilibrator ΔrG′° + xref links, and the species
   that encode the route. */
const GCOL = { Gpos: "#3a8bff", Gneg: "#ff8a3a", Arch: "#22d18c", Other: "#8aa0bf" };
const SUBCOL = ["#39c0ff", "#22d18c", "#ff8a3a", "#b07cff", "#ffd24a", "#ff6b9d", "#4defff", "#7ee787"];
const KEGG = r => `https://www.kegg.jp/entry/${r}`;
let INFO = {}, SMI = {};

function dbLink(db, id) {
  const u = { bigg: "http://bigg.ucsd.edu/universal/reactions/" + id.replace(/^R_/, ""),
    rhea: "https://www.rhea-db.org/rhea/" + id, seed: "https://modelseed.org/biochem/reactions/" + id,
    metacyc: "https://metacyc.org/META/NEW-IMAGE?object=" + id, kegg: KEGG(id) }[db];
  return `<a class="xref" href="${u}" target="_blank" rel="noopener"><span class="db">${db}</span>${id}</a>`;
}
const dgCls = dg => typeof dg === "number" ? (dg < -5 ? "fav" : dg > 5 ? "unfav" : "rev") : "rev";
function primarySub(rid) { const s = (INFO[rid] || {}).subsystems || []; return s[0] || "Other / unassigned"; }
const gdotR = g => `<span class="gdot" style="background:${GCOL[g] || GCOL.Other}"></span>`;
const UNCREDR = {
  no_ko: ["no enzyme mapped", "no KEGG ortholog is annotated for this exact transition — it may be a carbon-skeleton graph shortcut rather than a single reaction, so no gene can be named"],
  no_bacterial_donor: ["no bacterial donor", "the enzyme is known but no sequenced bacterium carries it — the gene must come from a plant, fungal or engineered source"],
};

let BUNDLE = null;
async function main() {
  const raw = sessionStorage.getItem("panroute_route");
  if (!raw) { document.getElementById("rpt").innerHTML = "<p>No route selected. Open a pathway from the main page.</p>"; return; }
  BUNDLE = JSON.parse(raw);
  if (BUNDLE.route && !BUNDLE.routes) { BUNDLE.routes = [BUNDLE.route]; BUNDLE.current = BUNDLE.route.id; }  // legacy bundle
  [INFO, SMI] = await Promise.all([
    fetch("data/rxninfo.json").then(r => r.json()).catch(() => ({})),
    fetch("data/smiles.json").then(r => r.json()).catch(() => ({})),
  ]);
  render(BUNDLE.current != null ? BUNDLE.current : (BUNDLE.routes[0] || {}).id);
}

function speciesFor(routeId) {
  if (BUNDLE.orgs) return BUNDLE.orgs.filter(o => (o.route_idx || []).includes(routeId)).sort((a, b) => b.n_routes - a.n_routes);
  return BUNDLE.species || [];   // legacy single-route bundle
}

function render(routeId) {
  const routes = BUNDLE.routes, q = BUNDLE.query;
  let idx = routes.findIndex(x => x.id === routeId); if (idx < 0) idx = 0;
  const r = routes[idx]; BUNDLE.current = r.id;
  const feasObj = (BUNDLE.feas && (BUNDLE.feas[r.id] || (BUNDLE.feas.feasible !== undefined ? BUNDLE.feas : null))) || null;
  // heterologous-design scenario for this route (present only when no genome makes the product)
  const het = (BUNDLE.hetero && BUNDLE.hetero[r.id]) || null;
  const heteroIdx = new Set(het ? het.hetero.map(h => h.idx) : []);

  // subsystem colour map for this route
  const subs = [...new Set(r.steps.map(s => primarySub(s.reactions[0].rid)))];
  const subColor = {}; subs.forEach((s, i) => subColor[s] = SUBCOL[i % SUBCOL.length]);

  // group consecutive steps into subsystem modules
  const modules = []; let cur = null;
  r.steps.forEach((st, i) => {
    const sub = primarySub(st.reactions[0].rid);
    if (!cur || cur.sub !== sub) { cur = { sub, steps: [], idx: [] }; modules.push(cur); }
    cur.steps.push(st); cur.idx.push(i);
  });

  const canvases = [];   // {id, smiles}
  function metNode(cid, name, role) {
    const smi = SMI[cid];
    const id = "s" + Math.random().toString(36).slice(2, 9);
    if (smi) canvases.push({ id, smiles: smi });
    const struct = smi ? `<canvas class="struct" id="${id}" width="118" height="98"></canvas>`
                       : `<div class="struct nostruct">${cid}<br>(no structure)</div>`;
    return `<div class="metnode ${role}">${struct}<div class="mname">${name}</div></div>`;
  }
  function rxStep(st, gi) {
    const rx = st.reactions[0], x = INFO[rx.rid] || {};
    const dg = x.our_dg;
    const dgchip = typeof dg === "number"
      ? `<span class="dgchip dg ${dgCls(dg)}" style="background:rgba(120,120,120,.15)">ΔG ${dg}</span>` : "";
    // in engineering mode, tag every step native (chassis) vs heterologous (clone)
    const engCls = het ? (heteroIdx.has(gi) ? " s-het" : " s-nat") : "";
    const flag = het ? `<span class="stepflag ${heteroIdx.has(gi) ? "het" : "nat"}">${heteroIdx.has(gi) ? "✚ clone" : "native"}</span>` : "";
    return `<div class="rxstep${engCls}"><div class="enz">${st.enzymes || x.ec || "?"}</div>
      <div class="rxarrow"></div>
      <div class="rxmeta"><a href="${KEGG(rx.rid)}" target="_blank">${rx.rid}</a>${x.ec ? " · EC " + x.ec : ""}</div>
      ${flag}${dgchip}</div>`;
  }

  const modHtml = modules.map(m => {
    const c = subColor[m.sub];
    const path = [r.path[m.idx[0]]].concat(m.idx.map(i => r.path[i + 1]));
    let flow = metNode(path[0].cid, path[0].name, m.idx[0] === 0 ? "startpt" : "");
    m.steps.forEach((st, k) => {
      flow += rxStep(st, m.idx[k]);
      const pi = m.idx[k] + 1;
      flow += metNode(r.path[pi].cid, r.path[pi].name, pi === r.path.length - 1 ? "endpt" : "");
    });
    return `<div class="module"><div class="modhead"><span class="mdot" style="background:${c}"></span>
      ${m.sub}<span class="msub">· ${m.steps.length} reaction${m.steps.length > 1 ? "s" : ""}</span></div>
      <div class="modflow">${flow}</div></div>`;
  }).join("");

  // reaction details (cross-DB directionality)
  const rxDetail = r.steps.map((st, i) => st.reactions.map(rx => {
    const x = INFO[rx.rid] || {};
    const rows = [["KEGG (arrow)", x.kegg_dir || "—"]];
    if (x.rhea_dir) rows.push(["Rhea", `<b>${x.rhea_dir}</b>`]);
    if (x.metacyc_dir) rows.push(["MetaCyc", `<b>${x.metacyc_dir}</b>`]);
    if (typeof x.our_dg === "number") rows.push(["our ΔrG′° (eQuilibrator)",
      `<span class="dg ${dgCls(x.our_dg)}">${x.our_dg} kJ/mol</span> · ${x.our_dir || ""}`]);
    const xrefs = [dbLink("kegg", rx.rid)].concat((x.bigg || []).map(b => dbLink("bigg", b)))
      .concat(x.rhea ? [dbLink("rhea", x.rhea)] : []).concat((x.seed || []).map(s => dbLink("seed", s)))
      .concat(x.metacyc ? [dbLink("metacyc", x.metacyc)] : []);
    return `<div class="rxcard"><div class="rxhead"><span class="rxid">${rx.rid}</span>
        <span class="rxec">EC ${x.ec || "—"} · step ${i + 1} · ${st.enzymes || "?"}${x.subsystems && x.subsystems.length ? " · " + x.subsystems.join(", ") : ""}</span></div>
      ${x.eq ? `<div class="rxeq">${x.eq}</div>` : ""}
      <div class="dirtable">${rows.map(([k, v]) => `<span class="k">${k}</span><span class="v">${v}</span>`).join("")}</div>
      <div class="xrefs">${xrefs.join("")}</div>
      <div class="enzrow"><button class="enzview${heteroIdx.has(i) ? " het" : ""}" data-rid="${rx.rid}" data-ko="${(rx.kos || []).join(",")}" data-sub="${st.from}" data-subname="${((r.path.find(p => p.cid === st.from) || {}).name || "").replace(/"/g, "&quot;")}">⚗ enzyme candidates</button><span class="enzann" id="enzann-${rx.rid}"></span></div>
      </div>`;
  }).join("")).join("");

  // species (filtered locally for this route)
  const species = speciesFor(r.id);
  let spHtml;
  if (BUNDLE.genome_pending) spHtml = `<div class="pending">Species that encode this route appear once the genome data finishes deploying — reload the report.</div>`;
  else if (!species.length) spHtml = `<p style="color:var(--dim)">No KEGG genome encodes this exact route.</p>`;
  else spHtml = `<p style="color:var(--dim);font-size:13px">${species.length} species encode this route</p>
      <div class="splist">${species.slice(0, 300).map(o => `<div class="spitem">
        <span class="gdot" style="background:${GCOL[o.gram] || GCOL.Other}"></span><i>${o.species}</i>
        ${o.thermo_feasible ? '<span class="badge feas">ΔG✓</span>' : ''}</div>`).join("")}</div>`;

  const feas = feasObj ? (feasObj.feasible ? `<span class="badge feas">ΔG feasible ${feasObj.dG_sum ?? ""} kJ/mol</span>`
                                           : `<span class="badge infeas">ΔG infeasible</span>`) : "";
  // heterologous-design plan (engineering mode)
  const planHtml = het ? `
    <div class="engbanner">
      <div class="eb-title">✚ Heterologous expression design — no genome makes this natively</div>
      <div class="eb-body">Proposed host: chassis ${gdotR(het.chassis.gram)}<b>${het.chassis.species}</b>, which natively runs
        <b>${het.chassis.cover}/${het.length}</b> steps. Clone <b>${het.n_clone}</b> gene${het.n_clone > 1 ? "s" : ""} for the
        missing (purple) step${het.n_clone > 1 ? "s" : ""}${het.n_uncreditable ? `; <b>${het.n_uncreditable}</b> step${het.n_uncreditable > 1 ? "s have" : " has"} no nameable bacterial gene` : ""}.
        <span style="color:var(--dim)">Genome potential — a strain-design starting point, not a guarantee of function.</span></div>
    </div>
    <div class="planlist">${het.hetero.map(h => {
      const head = `<span class="ps-idx">step ${h.idx + 1}</span> <b>${h.from_name} → ${h.to_name}</b>
        <span class="ps-rx"><a href="${KEGG(h.rid)}" target="_blank" rel="noopener">${h.rid}</a>${h.ec ? " · EC " + h.ec : ""}${h.kos.length ? " · KO " + h.kos.join(", ") : ""}</span>`;
      if (h.uncreditable) { const u = UNCREDR[h.uncreditable_reason] || ["unavailable", ""];
        return `<div class="planstep unc"><div class="ps-head">${head}</div><div class="ps-warn">⚠ <b>${u[0]}</b> — ${u[1]}</div></div>`; }
      const donors = h.donors.map(d => `<span class="donor" title="gene ${d.kos.join(", ")}">${gdotR(d.gram)}<i>${d.species}</i><span class="dko">${d.kos.join("/")}</span>${d.safety ? `<span class="dsafe ${d.safety === "pathogen" ? "bad" : d.safety === "GRAS/QPS" ? "good" : ""}">${d.safety}</span>` : ""}</span>`).join("");
      const more = h.n_donor_species > h.donors.length ? `<span class="dmore">+${(h.n_donor_species - h.donors.length).toLocaleString()} more donor species</span>` : "";
      return `<div class="planstep"><div class="ps-head">${head}</div><div class="ps-donors"><span class="pl">clone from:</span>${donors}${more}</div></div>`;
    }).join("")}</div>` : "";
  // navigation: back to the results view + move between ranked pathways without a round-trip
  const total = routes.length;
  const opts = routes.map((x, i) => `<option value="${x.id}"${x.id === r.id ? " selected" : ""}>#${i + 1} · ${x.length} steps · ${(BUNDLE.genomes && BUNDLE.genomes[x.id]) || 0} species</option>`).join("");
  const nav = `<div class="rnav">
      <a class="rback" href="#" id="toResults">← back to results</a>
      <div class="rpager">
        <button class="pgbtn" id="prevPw"${idx === 0 ? " disabled" : ""}>‹ prev</button>
        <select id="pwSelect" class="pwselect" title="jump to a pathway">${opts}</select>
        <button class="pgbtn" id="nextPw"${idx === total - 1 ? " disabled" : ""}>next ›</button>
        <span class="rpos">pathway ${idx + 1} of ${total}</span>
      </div></div>`;
  document.getElementById("rpt").innerHTML = nav + `
    <div class="rtitle"><i>${q.end.name}</i> from ${q.start.name}</div>
    <div class="rmeta">${r.length} steps · ${r.steps.reduce((a, s) => a + s.reactions.length, 0)} reactions · spans ${subs.length} subsystem${subs.length > 1 ? "s" : ""} ${feas}</div>
    <div class="subchips">${subs.map(s => `<span class="subchip" style="border-color:${subColor[s]}44;color:${subColor[s]}">${s}</span>`).join("")}</div>
    ${het ? `<h2 class="sec">Heterologous expression plan</h2>${planHtml}` : ""}
    <h2 class="sec">Pathway map · modularised by subsystem${het ? " <span class='hint' style='font-weight:400'>· green = native to chassis · purple = clone</span>" : ""}</h2>
    ${modHtml}
    <h2 class="sec">Enzyme engineering · variant kinetics &amp; thermostability</h2>
    <div class="engzone">
      <button id="enzAll">▶ Analyze all ${new Set(r.steps.flatMap(s => s.reactions.map(x => x.rid))).size} enzymes in this pathway</button>
      <span class="hint">One pass over the whole pathway: every step's KEGG orthologues → cluster at 80% → kcat/Km (MPEK) + thermostability (TemStaPro). The models load once for all steps (far faster than per-enzyme), and each reaction's result caches. Then open any step's <b>enzyme candidates</b> for its 3D landscape + ranked variants.</span>
      <div id="enzAllProg" class="el-prog"></div>
    </div>
    <h2 class="sec">Reaction details · directionality across databases</h2>
    <p style="color:var(--dim);font-size:12.5px;margin:-6px 0 14px">Directionality as reported by each source (KEGG, Rhea, MetaCyc) plus our own component-contribution ΔrG′° where computed. Cross-refs via MetaNetX; not every reaction is mapped in every database.</p>
    ${rxDetail}
    <h2 class="sec">Species that encode this route</h2>
    ${spHtml}
    <p style="color:var(--dim);font-size:11.5px;margin-top:24px">Genome <i>potential</i>, not proof of production. KEGG / MetaNetX-derived · research use.</p>`;

  // one button analyses every enzyme in the pathway; per-step links then view each result
  const allBtn = document.getElementById("enzAll");
  if (allBtn) allBtn.onclick = () => runPathwayEnzymes();
  [...document.querySelectorAll(".enzview")].forEach(b => b.onclick = () =>
    openEnzymeLab(b.dataset.rid, { ko: b.dataset.ko, sub: b.dataset.sub, name: b.dataset.subname }));

  // wire navigation
  const go = id => { window.scrollTo(0, 0); render(id); };
  document.getElementById("toResults").onclick = e => { e.preventDefault();
    // history.back() returns to the exact results view (bfcache); app.js re-runs as a fallback
    if (history.length > 1) history.back(); else location.href = "index.html"; };
  const prev = document.getElementById("prevPw"), next = document.getElementById("nextPw");
  if (prev) prev.onclick = () => { if (idx > 0) go(routes[idx - 1].id); };
  if (next) next.onclick = () => { if (idx < total - 1) go(routes[idx + 1].id); };
  document.getElementById("pwSelect").onchange = e => go(+e.target.value);
  // keyboard: ← / → move between pathways
  document.onkeydown = e => { if (e.key === "ArrowLeft" && idx > 0) go(routes[idx - 1].id);
    else if (e.key === "ArrowRight" && idx < total - 1) go(routes[idx + 1].id); };

  // draw structures
  if (window.SmilesDrawer) {
    const drawer = new SmilesDrawer.Drawer({ width: 116, height: 96, padding: 6, bondThickness: 1.1, compactDrawing: true });
    canvases.forEach(c => { try {
      SmilesDrawer.parse(c.smiles, tree => drawer.draw(tree, c.id, "light", false), () => {});
    } catch (e) {} });
  }
}
/* ================= enzyme lab: variant kinetics + thermostability deep-dive =================
   Loads a precomputed per-reaction bundle (data/enzymes/<rid>.json) — real KEGG orthologue
   sequences clustered at 80%, each scored by MPEK (kcat/Km) and TemStaPro (thermostability).
   Renders a 3D scatter (Km x kcat x thermostability) + a sortable ranked table + a sequence /
   source-organism viewer. Predictions are offline-precomputed (the models are too large to run
   in-browser); an un-computed reaction says so honestly. */
const ENZCACHE = {};
let EL_STATE = null;

function elClose() { document.getElementById("enzlab").classList.add("hidden");
  try { Plotly.purge("el-scatter"); } catch (e) {} }
document.addEventListener("click", e => { if (e.target && e.target.id === "el-close") elClose(); });
document.addEventListener("keydown", e => { if (e.key === "Escape") elClose(); });

async function openEnzymeLab(rid, opts) {
  opts = opts || {};
  const lab = document.getElementById("enzlab"); lab.classList.remove("hidden");
  document.getElementById("el-title").innerHTML = `Enzyme candidates · <span style="color:var(--neon)">${rid}</span>`;
  const body = document.getElementById("el-body");
  body.innerHTML = `<div class="el-loading">loading variant predictions…</div>`;
  let data = ENZCACHE[rid];
  if (!data) {
    try { const r = await fetch(`data/enzymes/${rid}.json`); if (!r.ok) throw 0; data = await r.json(); ENZCACHE[rid] = data; }
    catch (e) { data = null; }
  }
  if (data && data.variants && data.variants.length) {
    EL_STATE = { rid, data, sortKey: "rank", sortDir: 1, selected: 0 };
    return renderEnzLab();
  }
  // not cached -> direct the user to the ONE whole-pathway button (analysing per-enzyme reloads the
  // gigabyte models every time; the pathway button computes them all in a single model-load).
  body.innerHTML = `<div class="el-empty"><div class="ee-ic">⚗</div>
    <p><b>This enzyme hasn't been analysed yet.</b></p>
    <p>Use <b>▶ Analyze all enzymes in this pathway</b> at the top of the report — one pass characterises
    every step at once (kcat / Km via MPEK, thermostability via TemStaPro), which is much faster than
    running enzymes one by one, and each reaction's result is cached. Then reopen this step for its 3D
    landscape and ranked variant table.</p>
    <button id="el-goall" class="enzview" style="font-size:13px;padding:8px 16px">▶ Analyze the whole pathway now</button>
    <p class="hint" style="margin-top:14px">This runs on the analysis backend (<code>python -m server.app</code> → localhost:8000);
    the static site shows results once they're computed.</p></div>`;
  const gb = document.getElementById("el-goall");
  if (gb) gb.onclick = () => { elClose(); const a = document.getElementById("enzAll"); if (a) { a.scrollIntoView({ block: "center" }); a.click(); } };
}

/* Analyse EVERY enzyme in the pathway in one backend pass (one MPEK + one TemStaPro model-load for all
   steps), streaming progress. On completion each reaction has a cached bundle to open. */
let ENZ_ALL_RUNNING = false;
function runPathwayEnzymes() {
  if (ENZ_ALL_RUNNING) return;
  const steps = [...document.querySelectorAll(".enzview")].map(b => ({
    rid: b.dataset.rid, ko: b.dataset.ko, sub_cid: b.dataset.sub, sub_name: b.dataset.subname }));
  const uniq = []; const seen = new Set();
  for (const s of steps) { if (s.rid && s.ko && !seen.has(s.rid)) { seen.add(s.rid); uniq.push(s); } }
  const prog = document.getElementById("enzAllProg");
  const btn = document.getElementById("enzAll");
  const log = (m, cls) => { if (prog) { const d = document.createElement("div"); d.className = "pl" + (cls ? " " + cls : ""); d.textContent = m; prog.appendChild(d); prog.scrollTop = prog.scrollHeight; } };
  if (!uniq.length) { log("no reactions with a mapped enzyme (KO) to analyse.", "err"); return; }
  ENZ_ALL_RUNNING = true; prog.innerHTML = "";
  if (btn) { btn.disabled = true; btn.textContent = `analysing ${uniq.length} enzymes…`; }
  log(`requesting analysis of ${uniq.length} enzymes across the pathway…`);
  let es, gotAny = false;
  const q = new URLSearchParams({ steps: JSON.stringify(uniq) });
  try { es = new EventSource(`api/pathway_enzymes?${q}`); }
  catch (e) { log("no analysis backend reachable.", "err"); ENZ_ALL_RUNNING = false; return; }
  const finish = () => { ENZ_ALL_RUNNING = false; if (btn) { btn.disabled = false; btn.textContent = `▶ Analyze all ${uniq.length} enzymes in this pathway`; } };
  es.addEventListener("progress", e => { gotAny = true; try { log(JSON.parse(e.data).msg); } catch (x) {} });
  es.addEventListener("done", e => { gotAny = true; es.close(); finish();
    let d = {}; try { d = JSON.parse(e.data); } catch (x) {}
    const rids = d.rids || [];
    log(`✓ analysed ${rids.length}/${d.total || uniq.length} enzymes — open any step's “enzyme candidates”.`, "ok");
    rids.forEach(rid => { const b = document.querySelector(`.enzview[data-rid="${rid}"]`); if (b) { b.classList.add("done"); b.textContent = "⚗ enzyme candidates ✓"; } });
  });
  es.addEventListener("error", e => { let m = ""; try { m = JSON.parse(e.data).message; } catch (x) {} if (m) { log("error: " + m, "err"); es.close(); finish(); } });
  es.onerror = () => { es.close();
    if (!gotAny) { log("could not reach the analysis backend (this static site has none).", "err");
      log("run it live with:  python -m server.app   → open localhost:8000", "err"); }
    finish(); };
}

function renderEnzLab() {
  const { data } = EL_STATE;
  const body = document.getElementById("el-body");
  const kmVals = data.variants.map(v => v.km).filter(x => x != null);
  const hasKin = kmVals.length > 0;
  body.innerHTML = `
    <div class="el-sum">
      <b>${data.n_variants}</b> enzyme variants (clustered from <b>${data.n_sequences || data.n_ko_genes}</b> KEGG orthologues at 80% identity)
      · substrate <b>${data.substrate.name || data.substrate.cid || "—"}</b>
      · KO ${(data.ko || []).join(", ")}
      <div class="el-method">kcat/Km: ${data.method.kinetics} · thermostability: ${data.method.thermostability} · <span style="color:var(--amber)">model predictions — a ranked shortlist for cloning, not measured values</span></div>
    </div>
    <div class="el-grid">
      <div class="el-left">
        <div class="el-plottitle">Km × kcat × thermostability${hasKin ? "" : " (kinetics unavailable — no substrate SMILES)"}</div>
        <div id="el-scatter"></div>
        <div class="el-axnote">x = Km (mM, log) · y = kcat (1/s, log) · z = P(Tm &gt; 55 °C) · colour = rank · lower-left-high is best</div>
      </div>
      <div class="el-right">
        <div id="el-seq"></div>
      </div>
    </div>
    <div class="el-tablewrap"><table class="el-table" id="el-table"></table></div>`;
  drawScatter(); drawTable(); selectVariant(EL_STATE.selected);
}

function drawScatter() {
  const V = EL_STATE.data.variants;
  const pts = V.filter(v => v.km != null && v.kcat != null);
  const z = v => { const c = v.thermo || {}; return c["55"] != null ? c["55"] : (c["50"] != null ? c["50"] : 0); };
  const trace = {
    type: "scatter3d", mode: "markers",
    x: pts.map(v => v.km), y: pts.map(v => v.kcat), z: pts.map(z),
    text: pts.map(v => `#${v.rank} ${v.organisms[0] || v.rep_gene}<br>kcat ${fmtNum(v.kcat)} /s · Km ${fmtNum(v.km)} mM<br>P(Tm>55°)=${(z(v)).toFixed(2)} · ${v.thermo_label || ""}`),
    hoverinfo: "text",
    marker: { size: pts.map(v => 5 + Math.min(9, Math.log2((v.cluster_size || 1) + 1) * 2)),
      color: pts.map(v => v.rank), colorscale: "Viridis", reversescale: true, opacity: .9,
      line: { width: 0 } },
  };
  const dark = { paper_bgcolor: "rgba(0,0,0,0)", font: { color: "#9fb0cc", size: 11 },
    margin: { l: 0, r: 0, t: 6, b: 0 },
    scene: { xaxis: { title: "Km (mM)", type: "log", gridcolor: "#22304a", color: "#7d8ca6" },
      yaxis: { title: "kcat (1/s)", type: "log", gridcolor: "#22304a", color: "#7d8ca6" },
      zaxis: { title: "P(Tm>55°C)", gridcolor: "#22304a", color: "#7d8ca6" },
      bgcolor: "rgba(0,0,0,0)" } };
  try { Plotly.newPlot("el-scatter", pts.length ? [trace] : [], dark, { responsive: true, displayModeBar: false }); }
  catch (e) { document.getElementById("el-scatter").innerHTML = `<div class="hint" style="padding:20px">3D view unavailable</div>`; }
}

const COLS = [
  ["rank", "#", v => v.rank],
  ["organism", "top organism", v => v.organisms[0] || v.rep_gene],
  ["kcat", "kcat (1/s)", v => v.kcat],
  ["km", "Km (mM)", v => v.km],
  ["kcat_km", "kcat/Km", v => v.kcat_km],
  ["thermo", "P(Tm>55°)", v => (v.thermo || {})["55"]],
  ["thermo_label", "Tm range", v => v.thermo_label],
  ["length", "len", v => v.length],
  ["cluster_size", "orthologs", v => v.cluster_size],
];
function drawTable() {
  const S = EL_STATE, V = S.data.variants.slice();
  const col = COLS.find(c => c[0] === S.sortKey) || COLS[0];
  V.sort((a, b) => { const x = col[2](a), y = col[2](b);
    if (x == null) return 1; if (y == null) return -1;
    return (typeof x === "string" ? x.localeCompare(y) : x - y) * S.sortDir; });
  const head = `<thead><tr>${COLS.map(c =>
    `<th data-k="${c[0]}" class="${S.sortKey === c[0] ? "on" : ""}">${c[1]}${S.sortKey === c[0] ? (S.sortDir > 0 ? " ▲" : " ▼") : ""}</th>`).join("")}</tr></thead>`;
  const rows = V.map(v => `<tr data-id="${v.id}" class="${v.id === S.data.variants[S.selected].id ? "sel" : ""}">
    <td>${v.rank}</td><td class="org"><i>${v.organisms[0] || v.rep_gene}</i>${v.organisms.length > 1 ? ` <span class="hint">+${v.organisms.length - 1}</span>` : ""}</td>
    <td>${fmtNum(v.kcat)}</td><td>${fmtNum(v.km)}</td><td>${fmtNum(v.kcat_km)}</td>
    <td>${(v.thermo || {})["55"] != null ? (v.thermo["55"]).toFixed(2) : "—"}</td>
    <td>${v.thermo_label || "—"}</td><td>${v.length}</td><td>${v.cluster_size}</td></tr>`).join("");
  const t = document.getElementById("el-table");
  t.innerHTML = head + `<tbody>${rows}</tbody>`;
  [...t.querySelectorAll("th")].forEach(th => th.onclick = () => {
    const k = th.dataset.k; if (S.sortKey === k) S.sortDir *= -1; else { S.sortKey = k; S.sortDir = (k === "km") ? 1 : -1; }
    drawTable(); });
  [...t.querySelectorAll("tbody tr")].forEach(tr => tr.onclick = () => {
    const idx = S.data.variants.findIndex(v => v.id === tr.dataset.id); selectVariant(idx); });
}
function selectVariant(idx) {
  EL_STATE.selected = idx; const v = EL_STATE.data.variants[idx];
  const seq = (v.sequence || "").replace(/(.{60})/g, "$1\n");
  const curve = v.thermo || {};
  const bars = [40, 45, 50, 55, 60, 65].map(t => { const p = curve[t]; return p == null ? "" :
    `<div class="tb"><div class="tb-bar" style="height:${Math.round(p * 54)}px;background:${p > .5 ? "var(--green)" : "var(--amber)"}"></div><div class="tb-l">${t}</div></div>`; }).join("");
  document.getElementById("el-seq").innerHTML = `
    <div class="es-h">#${v.rank} · <i>${v.organisms[0] || v.rep_gene}</i></div>
    <div class="es-kin"><span>kcat <b>${fmtNum(v.kcat)}</b> /s</span><span>Km <b>${fmtNum(v.km)}</b> mM</span><span>kcat/Km <b>${fmtNum(v.kcat_km)}</b></span></div>
    <div class="es-therm"><div class="es-lab">thermostability — P(Tm &gt; T)</div><div class="tbars">${bars || '<span class="hint">no thermostability call</span>'}</div><div class="es-tl">${v.thermo_label || ""}</div></div>
    <div class="es-org"><div class="es-lab">clone from (${v.organisms.length} organism${v.organisms.length > 1 ? "s" : ""} share this variant cluster)</div>
      <div class="es-orgs">${v.organisms.slice(0, 12).map(o => `<span class="oc"><i>${o}</i></span>`).join("")}${v.organisms.length > 12 ? `<span class="hint">+${v.organisms.length - 12} more</span>` : ""}</div></div>
    <div class="es-seq"><div class="es-lab">representative sequence · ${v.length} aa · <span class="hint">${v.rep_gene}</span></div><pre>${seq}</pre></div>`;
  drawTable();
}
function fmtNum(x) { if (x == null) return "—"; const a = Math.abs(x);
  if (a === 0) return "0"; if (a < 0.01 || a >= 10000) return x.toExponential(2);
  return (+x).toPrecision(3).replace(/\.?0+$/, ""); }

main();

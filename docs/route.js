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
      <div class="xrefs">${xrefs.join("")}</div></div>`;
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
    <h2 class="sec">Reaction details · directionality across databases</h2>
    <p style="color:var(--dim);font-size:12.5px;margin:-6px 0 14px">Directionality as reported by each source (KEGG, Rhea, MetaCyc) plus our own component-contribution ΔrG′° where computed. Cross-refs via MetaNetX; not every reaction is mapped in every database.</p>
    ${rxDetail}
    <h2 class="sec">Species that encode this route</h2>
    ${spHtml}
    <p style="color:var(--dim);font-size:11.5px;margin-top:24px">Genome <i>potential</i>, not proof of production. KEGG / MetaNetX-derived · research use.</p>`;

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
main();

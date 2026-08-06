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

async function main() {
  const raw = sessionStorage.getItem("panroute_route");
  if (!raw) { document.getElementById("rpt").innerHTML = "<p>No route selected. Open a pathway from the main page.</p>"; return; }
  const D = JSON.parse(raw);
  [INFO, SMI] = await Promise.all([
    fetch("data/rxninfo.json").then(r => r.json()).catch(() => ({})),
    fetch("data/smiles.json").then(r => r.json()).catch(() => ({})),
  ]);
  const r = D.route, q = D.query;

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
  function rxStep(st) {
    const rx = st.reactions[0], x = INFO[rx.rid] || {};
    const dg = x.our_dg;
    const dgchip = typeof dg === "number"
      ? `<span class="dgchip dg ${dgCls(dg)}" style="background:rgba(120,120,120,.15)">ΔG ${dg}</span>` : "";
    return `<div class="rxstep"><div class="enz">${st.enzymes || x.ec || "?"}</div>
      <div class="rxarrow"></div>
      <div class="rxmeta"><a href="${KEGG(rx.rid)}" target="_blank">${rx.rid}</a>${x.ec ? " · EC " + x.ec : ""}</div>
      ${dgchip}</div>`;
  }

  const modHtml = modules.map(m => {
    const c = subColor[m.sub];
    const path = [r.path[m.idx[0]]].concat(m.idx.map(i => r.path[i + 1]));
    let flow = metNode(path[0].cid, path[0].name, m.idx[0] === 0 ? "startpt" : "");
    m.steps.forEach((st, k) => {
      flow += rxStep(st);
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

  // species
  let spHtml;
  if (D.genome_pending) spHtml = `<div class="pending">Species that encode this route appear once the genome data finishes deploying — reload the report.</div>`;
  else if (!D.species || !D.species.length) spHtml = `<p style="color:var(--dim)">No KEGG genome encodes this exact route.</p>`;
  else spHtml = `<p style="color:var(--dim);font-size:13px">${D.species.length} species encode this route</p>
      <div class="splist">${D.species.slice(0, 300).map(o => `<div class="spitem">
        <span class="gdot" style="background:${GCOL[o.gram] || GCOL.Other}"></span><i>${o.species}</i>
        ${o.thermo_feasible ? '<span class="badge feas">ΔG✓</span>' : ''}</div>`).join("")}</div>`;

  const feas = D.feas ? (D.feas.feasible ? `<span class="badge feas">ΔG feasible ${D.feas.dG_sum ?? ""} kJ/mol</span>`
                                         : `<span class="badge infeas">ΔG infeasible</span>`) : "";
  document.getElementById("rpt").innerHTML = `
    <a class="rback" href="index.html">← back to search</a>
    <div class="rtitle"><i>${q.end.name}</i> from ${q.start.name}</div>
    <div class="rmeta">${r.length} steps · ${r.steps.reduce((a, s) => a + s.reactions.length, 0)} reactions · spans ${subs.length} subsystem${subs.length > 1 ? "s" : ""} ${feas}</div>
    <div class="subchips">${subs.map(s => `<span class="subchip" style="border-color:${subColor[s]}44;color:${subColor[s]}">${s}</span>`).join("")}</div>
    <h2 class="sec">Pathway map · modularised by subsystem</h2>
    ${modHtml}
    <h2 class="sec">Reaction details · directionality across databases</h2>
    <p style="color:var(--dim);font-size:12.5px;margin:-6px 0 14px">Directionality as reported by each source (KEGG, Rhea, MetaCyc) plus our own component-contribution ΔrG′° where computed. Cross-refs via MetaNetX; not every reaction is mapped in every database.</p>
    ${rxDetail}
    <h2 class="sec">Species that encode this route</h2>
    ${spHtml}
    <p style="color:var(--dim);font-size:11.5px;margin-top:24px">Genome <i>potential</i>, not proof of production. KEGG / MetaNetX-derived · research use.</p>`;

  // draw structures
  if (window.SmilesDrawer) {
    const drawer = new SmilesDrawer.Drawer({ width: 116, height: 96, padding: 6, bondThickness: 1.1, compactDrawing: true });
    canvases.forEach(c => { try {
      SmilesDrawer.parse(c.smiles, tree => drawer.draw(tree, c.id, "light", false), () => {});
    } catch (e) {} });
  }
}
main();

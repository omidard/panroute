/* PanRoute route report page. Reads the clicked route from sessionStorage, fetches the
   cross-database reaction info (rxninfo.json), and renders: the pathway, each reaction with
   its cross-DB xrefs (BiGG/Rhea/SEED/MetaCyc) + directionality from every source + our own
   eQuilibrator ΔG, and the species that encode this route. */
const GCOL = { Gpos: "#3a8bff", Gneg: "#ff8a3a", Arch: "#22d18c", Other: "#8aa0bf" };
const KEGG = r => `https://www.kegg.jp/entry/${r}`;

function dbLink(db, id) {
  const u = {
    bigg: "http://bigg.ucsd.edu/universal/reactions/" + id.replace(/^R_/, ""),
    rhea: "https://www.rhea-db.org/rhea/" + id,
    seed: "https://modelseed.org/biochem/reactions/" + id,
    metacyc: "https://metacyc.org/META/NEW-IMAGE?object=" + id,
    kegg: KEGG(id),
  }[db];
  return `<a class="xref" href="${u}" target="_blank" rel="noopener"><span class="db">${db}</span>${id}</a>`;
}

function dgClass(dg, dir) {
  if (typeof dg === "number") return dg < -5 ? "fav" : dg > 5 ? "unfav" : "rev";
  return "rev";
}

async function main() {
  const raw = sessionStorage.getItem("panroute_route");
  if (!raw) { document.getElementById("rpt").innerHTML = "<p>No route selected. Open a pathway from the main page.</p>"; return; }
  const D = JSON.parse(raw);
  const info = await fetch("rxninfo.json").then(r => r.json()).catch(() => ({}));
  const r = D.route, q = D.query;
  const rankSteps = r.length;
  const chain = r.path.map((p, i) => {
    const cls = i === 0 ? "start" : (i === r.path.length - 1 ? "end" : "");
    const enz = i > 0 ? `<span class="ar">→</span><span class="enz">${r.steps[i - 1].enzymes || "?"}</span><span class="ar">→</span>` : "";
    return enz + `<span class="met ${cls}">${p.name}</span>`;
  }).join(" ");

  const feas = D.feas ? (D.feas.feasible
    ? `<span class="badge feas">ΔG feasible ${D.feas.dG_sum ?? ""} kJ/mol</span>`
    : `<span class="badge infeas">ΔG infeasible</span>`) : "";

  // reactions
  const rxHtml = r.steps.map((st, i) => st.reactions.map(rx => {
    const x = info[rx.rid] || {};
    const rows = [];
    rows.push(["KEGG (arrow)", x.kegg_dir || "—"]);
    if (x.rhea_dir) rows.push(["Rhea", `<b>${x.rhea_dir}</b>`]);
    if (x.metacyc_dir) rows.push(["MetaCyc", `<b>${x.metacyc_dir}</b>`]);
    if (typeof x.our_dg === "number")
      rows.push(["our ΔrG′° (eQuilibrator)", `<span class="dg ${dgClass(x.our_dg)}">${x.our_dg} kJ/mol</span> · ${x.our_dir || ""}`]);
    else if (x.our_dir)
      rows.push(["our direction", x.our_dir]);
    const xrefs = [dbLink("kegg", rx.rid)]
      .concat((x.bigg || []).map(b => dbLink("bigg", b)))
      .concat(x.rhea ? [dbLink("rhea", x.rhea)] : [])
      .concat((x.seed || []).map(s => dbLink("seed", s)))
      .concat(x.metacyc ? [dbLink("metacyc", x.metacyc)] : []);
    return `<div class="rxcard">
      <div class="rxhead"><span class="rxid">${rx.rid}</span>
        <span class="rxec">EC ${x.ec || rx.ec.join(", ") || "—"} · step ${i + 1} · enzyme <b style="color:var(--cyan)">${st.enzymes || "?"}</b></span></div>
      ${x.eq ? `<div class="rxeq">${x.eq}</div>` : ""}
      <div class="dirtable">${rows.map(([k, v]) => `<span class="k">${k}</span><span class="v">${v}</span>`).join("")}</div>
      <div class="xrefs">${xrefs.join("")}</div>
    </div>`;
  }).join("")).join("");

  // species
  let spHtml;
  if (D.genome_pending) {
    spHtml = `<div class="pending">Species that encode this route are computed once the genome data finishes uploading — reload the report then.</div>`;
  } else if (!D.species || !D.species.length) {
    spHtml = `<p style="color:var(--dim)">No KEGG genome encodes this exact route.</p>`;
  } else {
    spHtml = `<p style="color:var(--dim);font-size:13px">${D.species.length} species encode this route</p>
      <div class="splist">${D.species.slice(0, 300).map(o =>
        `<div class="spitem"><span class="gdot" style="background:${GCOL[o.gram] || GCOL.Other}"></span>
          <i>${o.species}</i>${o.thermo_feasible ? '<span class="badge feas">ΔG✓</span>' : ''}</div>`).join("")}</div>`;
  }

  document.getElementById("rpt").innerHTML = `
    <a class="rback" href="index.html">← back to search</a>
    <div class="rtitle"><i>${q.end.name}</i> from ${q.start.name}</div>
    <div class="rmeta">pathway · ${rankSteps} steps · ${r.steps.reduce((a, s) => a + s.reactions.length, 0)} reactions ${feas}</div>
    <div class="rchain">${chain}</div>
    <h2 class="sec">Reactions · directionality across databases</h2>
    <p style="color:var(--dim);font-size:12.5px;margin:-6px 0 14px">Directionality is shown as reported by each source (KEGG, Rhea, MetaCyc) plus our own component-contribution ΔrG′° where computed. Cross-references via MetaNetX; not every reaction is mapped in every database.</p>
    ${rxHtml}
    <h2 class="sec">Species that encode this route</h2>
    ${spHtml}
    <p style="color:var(--dim);font-size:11.5px;margin-top:24px">Genome <i>potential</i>, not proof of production. KEGG/MetaNetX-derived · research use.</p>`;
}
main();

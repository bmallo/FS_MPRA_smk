/* No network calls or external libraries: this viewer works from a local file. */
"use strict";
const Pilot = (() => {
  const tracks = ["tf_small", "msp", "nucleosome"];
  const finite = x => typeof x === "number" && Number.isFinite(x);
  const escape = x => String(x ?? "").replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const format = (x, digits=3) => finite(x) ? x.toLocaleString(undefined, {maximumFractionDigits:digits}) : "unmeasured";
  const pvalue = x => !finite(x) ? "unmeasured" : x === 0 ? "0 (reported)" : x < .00001 ? x.toExponential(2) : format(x,5);
  function index(data) {
    if (data?.schema_version !== "fsmpra.pilot.v1" || !data.sample?.id || !Array.isArray(data.variants) || !data.region?.sequence) throw Error("Expected a fsmpra.pilot.v1 measurement file.");
    const byId = new Map(), byCell = new Map(), length = data.region.sequence.length;
    for (const variant of data.variants) {
      if (byId.has(variant.variant_id) || !Number.isInteger(variant.n_reads) || variant.n_reads < 0) throw Error("Invalid variant identity or molecule count.");
      for (const track of tracks) {
        const group = variant.tracks?.[track];
        for (const field of ["occupancy", "delta_vs_WT", "eligible_reads", "called_reads"]) {
          if (!Array.isArray(group?.[field]) || group[field].length !== length || group[field].some(x => x !== null && !finite(x))) throw Error(`Invalid ${track} ${field} array.`);
        }
      }
      byId.set(variant.variant_id, variant);
      if (variant.class === "SNV") {
        const key = `${variant.oligo_pos1}:${variant.alt}`;
        if (byCell.has(key)) throw Error("Duplicate SNV coordinate/allele.");
        byCell.set(key, variant);
      }
    }
    const wt = data.variants.filter(v => v.class === "WT");
    if (wt.length !== 1) throw Error("Expected one WT sequence.");
    return {byId, byCell, wt:wt[0], ordered:[...data.variants].sort((a,b) => (a.oligo_pos1 ?? -1)-(b.oligo_pos1 ?? -1) || a.alt.localeCompare(b.alt))};
  }
  function metric(variant, name) { return name === "n_reads" ? variant.n_reads : variant.summary?.[name] ?? null; }
  function path(values, x, y) {
    let active = false, result = "";
    values.forEach((value, i) => {
      if (!finite(value)) { active = false; return; }
      result += `${active ? "L" : "M"}${x(i).toFixed(2)},${y(value).toFixed(2)} `;
      active = true;
    });
    return result;
  }
  function normalized(histogram) {
    const sum = histogram.reduce((a,b) => a+b, 0);
    return histogram.map(n => sum ? n/sum : null);
  }
  function color(value, limit, count=false) {
    if (!finite(value)) return "#33434b";
    const t = Math.min(1, Math.abs(value) / Math.max(limit, 1e-9));
    const a = count ? [27,46,57] : [213,221,224];
    const b = count ? [109,216,178] : value < 0 ? [42,121,181] : [205,87,68];
    return `rgb(${a.map((v,i) => Math.round(v+(b[i]-v)*t)).join(",")})`;
  }
  return {tracks, finite, escape, format, pvalue, index, metric, path, normalized, color};
})();
if (typeof module !== "undefined") module.exports = Pilot;

if (typeof document !== "undefined") (() => {
  const $ = id => document.getElementById(id), esc = Pilot.escape, fmt = Pilot.format;
  const datasets = new Map();
  let data, lookup, selectedId, geometry, comparison=null, scaleLimit=1;
  const names = {tf_small:"Small TF-track footprints", msp:"Accessible segments · MSP", nucleosome:"Nucleosome-sized footprints"};
  const colors = {tf_small:"#eab970", msp:"#65c7a5", nucleosome:"#ba9ddd"};
  function status(message, error=false) { $("status").textContent=message; $("status").classList.toggle("error", error); }
  function minimum() { return Math.max(1, Number($("minimum").value)||1); }
  function reaction(sample) { const r=sample.reaction; return r ? `${r.duration_minutes} min · ${r.hia5_volume_uL} µL Hia5 / ${fmt(r.cells/1e6)}M cells` : "Reaction not documented"; }
  function current() { return lookup.byId.get(selectedId) || lookup.wt; }
  function install(value) {
    if(value.schema_version==="fsmpra.collection.v1") {
      if(!Array.isArray(value.datasets)||!value.datasets.length)throw Error("Empty dataset collection");
      const ids=new Set();
      for(const entry of value.datasets) { Pilot.index(entry); if(ids.has(entry.sample.id))throw Error("Duplicate sample ID in collection"); ids.add(entry.sample.id); }
      datasets.clear(); comparison=value.comparison;
      for(const entry of value.datasets)datasets.set(entry.sample.id,entry);
      $("sample").replaceChildren(...[...datasets.keys()].map(id=>new Option(id,id)));
      const id=value.datasets[value.datasets.length-1].sample.id;
      $("sample").value=id; activate(id); return;
    }
    comparison=null;
    Pilot.index(value); datasets.set(value.sample.id, value);
    $("sample").replaceChildren(...[...datasets.keys()].map(id => new Option(id,id)));
    $("sample").value=value.sample.id; activate(value.sample.id);
  }
  function activate(id) {
    data=datasets.get(id); lookup=Pilot.index(data);
    const query = new URLSearchParams(location.hash.slice(1));
    const requested=query.get("variant");
    selectedId=(lookup.ordered.find(v => v.source_variant===requested) || lookup.wt).variant_id;
    $("search").value="";
    $("app").hidden=false;
    const n=lookup.ordered.filter(v=>v.class==="SNV" && v.n_reads>0).length;
    $("stats").innerHTML=[
      [data.sample.biological_replicate, "Biological preparation"],
      [fmt(data.counts.measured,0), "WT / SNV molecules measured"],
      [`${n} / ${lookup.byCell.size}`, "SNVs with sampled evidence"],
      [fmt(lookup.wt.n_reads,0), "Wild-type molecules"]
    ].map(([value,label])=>`<div class="stat"><strong>${esc(value)}</strong><span>${esc(label)}</span></div>`).join("");
    status(`${reaction(data.sample)} · ${data.sampling.fraction===1?"complete genotype ledger":fmt(data.sampling.fraction*100,2)+"% ledger sample"} · ${data.region.sequence.length}-bp mutable element · descriptive measurements; no significance tests`);
    $("methods").innerHTML=`<p>${esc(data.region.status)}</p><ul>${data.limitations.map(x=>`<li>${esc(x)}</li>`).join("")}</ul><pre>${esc(JSON.stringify({sample:data.sample,method:data.method,sampling:data.sampling,counts:data.counts,provenance:data.provenance},null,2))}</pre>`;
    renderOptions(); drawMap(); renderDetail(); renderComparison();
  }
  function renderComparison() {
    $("comparison-panel").hidden=!comparison;
    if(!comparison)return;
    $("agreement-summary").innerHTML=`<table><thead><tr><th>Runs</th><th>Reaction comparison</th><th>Shared SNVs ≥${comparison.display_minimum_reads}</th><th>MSP Δ Pearson r</th><th>Small-footprint Δ Pearson r</th></tr></thead><tbody>${comparison.pairwise.map(p=>`<tr><td>${esc(p.sample_a)} / ${esc(p.sample_b)}</td><td>${p.comparison_type==="same_stated_reaction"?"Same stated reaction":p.comparison_type==="different_reaction_conditions"?"Different reactions":"Unknown"}</td><td>${p.SNVs_above_minimum_in_both}</td><td>${fmt(p.msp_common_site_delta_pearson)}</td><td>${fmt(p.tf_small_common_site_delta_pearson)}</td></tr>`).join("")}</tbody></table><p>Descriptive correlations of region-wide variant-minus-WT summaries on common sites. These do not test localized footprint effects. Unequal depth and uncalibrated observability can affect agreement; no significance tests or confidence intervals are provided.</p>`;
    $("comparison-summary").innerHTML=`<table><thead><tr><th>Run / reaction</th><th>Ledger used</th><th>Full plasmids¹</th><th>WT molecules</th><th>SNVs observed</th><th>SNVs ≥${comparison.display_minimum_reads}</th><th>Median SNV n</th><th>WT MSP mean</th></tr></thead><tbody>${comparison.sample_qc.map(q=>`<tr><td>${esc(q.sample_id)}<br>${esc(reaction(datasets.get(q.sample_id).sample))}</td><td>${fmt(q.sampling_fraction*100,2)}%</td><td>${Pilot.finite(q.full_plasmid_fraction_of_metrics_matched_sample)?fmt(q.full_plasmid_fraction_of_metrics_matched_sample*100,1)+"%":"unmeasured"}</td><td>${fmt(q.WT_molecules,0)}</td><td>${q.SNVs_observed}/${q.SNVs_total}</td><td>${q.SNVs_above_display_minimum}</td><td>${fmt(q.SNV_median_molecules,0)}</td><td>${fmt(q.WT_msp_fraction_mean)}</td></tr>`).join("")}</tbody></table><p>¹ Fraction of metrics-matched molecules in the sampled genotype ledger; not total sequencing yield. The table uses a fixed ${comparison.display_minimum_reads}-molecule threshold. Runs are not pooled.</p>`;
    const recoveryRows=[];
    for(const [sampleId,evidence] of Object.entries(comparison.recovery||{})) {
      if(!evidence)continue;
      const summary=evidence.summary,total=Object.values(summary.by_plasmid_class).reduce((n,g)=>n+g.molecules,0);
      for(const [pc,g] of Object.entries(summary.by_plasmid_class)) {
        recoveryRows.push(`<tr><td>${esc(sampleId)}</td><td>${esc(pc)}</td><td>${fmt(g.molecules,0)}</td><td>${total?fmt(g.molecules/total*100,2)+"%":"unmeasured"}</td><td>${fmt(g.mean_read_length,0)}</td><td>${fmt((g.genotype_classes.WT||0)+(g.genotype_classes.SNV||0),0)}</td><td>${fmt(g.passes_below_10,0)}</td></tr>`);
      }
    }
    $("recovery-summary").innerHTML=recoveryRows.length?`<table><thead><tr><th>Run</th><th>Upstream class</th><th>Molecules</th><th>Of imported ledger</th><th>Mean length (bp)</th><th>Assigned WT/SNV</th><th>Passes &lt;10</th></tr></thead><tbody>${recoveryRows.join("")}</tbody></table><p>Counts before measurement selection. Full-plasmid classification depends on spanning both ends of the rotated reference; it differs from the ledger's length-based structure label. WT/SNV counts here have not yet passed the sequence and quality checks. Runs without a recovery summary are omitted.</p>`:'<p>No detailed recovery summary is available for these imports.</p>';
    $("comparison-methods").innerHTML=`<ul>${comparison.limitations.map(s=>`<li>${esc(s)}</li>`).join("")}</ul><pre>${esc(JSON.stringify({pairwise:comparison.pairwise,processing:comparison.processing_evidence,provenance:comparison.provenance},null,2))}</pre>`;
    const W=900,H=145,L=45,R=15,T=30,B=25,n=data.region.sequence.length,palette=["#eab970","#65c7a5","#ba9ddd"];
    $("comparison-profiles").innerHTML=["tf_small","msp"].map(track=>{
      const groups=[...datasets.values()].map((d,i)=>({d,wt:d.variants.find(v=>v.class==="WT"),color:palette[i%palette.length]}));
      const values=groups.flatMap(g=>g.wt.tracks[track].occupancy).filter(Pilot.finite);
      const max=track==="msp"?1:Math.max(.05,...values)*1.1;
      const x=i=>L+i/(n-1)*(W-L-R),y=v=>T+(1-v/max)*(H-T-B);
      let body=`<text x="${L}" y="15">WT · ${esc(names[track])}</text>`;
      for(const tick of [0,max/2,max])body+=`<line x1="${L}" x2="${W-R}" y1="${y(tick)}" y2="${y(tick)}" stroke="#33434b"/><text x="${L-5}" y="${y(tick)+3}" text-anchor="end">${fmt(tick*100,1)}%</text>`;
      for(let i=0;i<n;i+=40)body+=`<text x="${x(i)}" y="${H-5}">${data.region.oligo_start1+i}</text>`;
      groups.forEach((g,i)=>{body+=`<path d="${Pilot.path(g.wt.tracks[track].occupancy,x,y)}" fill="none" stroke="${g.color}" stroke-width="1.5"/><text x="${380+i*165}" y="15" style="fill:${g.color}">${esc(g.d.sample.id)}</text>`;});
      return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" role="img" aria-label="WT landscapes across runs">${body}</svg>`;
    }).join("");
  }
  function renderOptions() {
    const term=$("search").value.toLowerCase().trim();
    const options=lookup.ordered.filter(v=>v.source_variant.toLowerCase().includes(term));
    $("variant").replaceChildren(...options.map(v=>new Option(`${v.source_variant}   ·   n=${v.n_reads}`,v.variant_id)));
    $("variant").value=selectedId;
  }
  function choose(id) {
    if (!lookup.byId.has(id)) return;
    selectedId=id;
    if (![...$("variant").options].some(option=>option.value===id)) {
      $("search").value=""; renderOptions();
    }
    $("variant").value=id;
    const query=new URLSearchParams({sample:data.sample.id,variant:current().source_variant});
    try { history.replaceState(null,"",`#${query}`); } catch (_) { /* Embedded viewers may disallow history writes. */ }
    drawMap(); renderDetail();
  }
  function drawMap() {
    if (!data) return;
    const canvas=$("alleles"), dpr=window.devicePixelRatio||1;
    const width=canvas.getBoundingClientRect().width||1000, height=160;
    canvas.width=width*dpr; canvas.height=height*dpr;
    const ctx=canvas.getContext("2d"); ctx.scale(dpr,dpr);
    const left=26, right=12, top=10, row=27, cw=(width-left-right)/data.region.sequence.length;
    geometry={left,top,row,cw};
    const metric=$("metric").value, count=metric==="n_reads";
    const values=lookup.ordered.filter(v=>v.class==="SNV" && v.n_reads>=minimum()).map(v=>Pilot.metric(v,metric)).filter(Pilot.finite).map(Math.abs).sort((a,b)=>a-b);
    scaleLimit=Math.max(count?1:.005,values[Math.floor(.98*Math.max(0,values.length-1))]||0);
    ctx.font="11px system-ui"; ctx.textBaseline="middle";
    for (let j=0;j<4;j++) {
      ctx.fillStyle="#c5d2d9";ctx.fillText("ACGT"[j],5,top+j*row+row/2);
      for(let i=0;i<data.region.sequence.length;i++) {
        const variant=lookup.byCell.get(`${data.region.oligo_start1+i}:${"ACGT"[j]}`);
        const value=variant && variant.n_reads>=minimum()?Pilot.metric(variant,metric):null;
        ctx.fillStyle=Pilot.color(value,scaleLimit,count);
        ctx.fillRect(left+i*cw,top+j*row,Math.max(1,cw-.5),row-3);
        if(variant?.variant_id===selectedId) {ctx.strokeStyle="#fff";ctx.lineWidth=1.5;ctx.strokeRect(left+i*cw,top+j*row,Math.max(1,cw-.5),row-3);}
      }
    }
    ctx.fillStyle="#aebfc9";ctx.textAlign="center";
    for(let i=0;i<data.region.sequence.length;i+=20)ctx.fillText(String(data.region.oligo_start1+i),left+(i+.5)*cw,top+4*row+15);
    $("scale").textContent=count?`0 → ${fmt(scaleLimit,0)} molecules · color saturated above this`:`Blue: loss · red: gain · ±${fmt(scaleLimit)} fraction · gray: reference / low coverage / missing`;
  }
  function mapVariant(event) {
    if(!geometry)return null;
    const r=$("alleles").getBoundingClientRect(), {left,top,row,cw}=geometry;
    const i=Math.floor((event.clientX-r.left-left)/cw), j=Math.floor((event.clientY-r.top-top)/row);
    if(i<0||i>=data.region.sequence.length||j<0||j>3)return null;
    return lookup.byCell.get(`${data.region.oligo_start1+i}:${"ACGT"[j]}`);
  }
  function profileSVG(variant,track) {
    const W=900,H=130,L=46,R=16,T=26,B=25, n=data.region.sequence.length;
    const delta=$("mode").value==="delta", group=variant.tracks[track], wt=lookup.wt.tracks[track];
    const values=delta?group.delta_vs_WT:group.occupancy, baseline=delta?Array(n).fill(0):wt.occupancy;
    const all=[...values,...baseline].filter(Pilot.finite);
    const max=delta?Math.max(.02,...all.map(Math.abs))*1.08:track==="tf_small"?Math.max(.04,...all)*1.12:1;
    const min=delta?-max:0, x=i=>L+i/(n-1)*(W-L-R), y=v=>T+(max-v)/(max-min)*(H-T-B);
    let body=`<text x="${L}" y="14" class="track-title">${esc(names[track])}</text>`;
    for(const tick of [min,(min+max)/2,max])body+=`<line x1="${L}" x2="${W-R}" y1="${y(tick)}" y2="${y(tick)}" stroke="#33434b"/><text x="${L-8}" y="${y(tick)+3}" text-anchor="end">${fmt(tick*100,1)}%</text>`;
    for(let i=0;i<n;i+=40)body+=`<text x="${x(i)}" y="${H-6}" text-anchor="middle">${data.region.oligo_start1+i}</text>`;
    if(variant.class==="SNV") {const p=variant.oligo_pos1-data.region.oligo_start1;body+=`<line x1="${x(p)}" x2="${x(p)}" y1="${T}" y2="${H-B}" stroke="#c5d2d9" stroke-dasharray="2 3"/>`;}
    body+=`<path d="${Pilot.path(baseline,x,y)}" fill="none" stroke="#c3cdd3" stroke-width="1.4" stroke-dasharray="5 4"/><path d="${Pilot.path(values,x,y)}" fill="none" stroke="${colors[track]}" stroke-width="1.8"/>`;
    if(!values.some(Pilot.finite))body+=`<text x="${W/2}" y="70" text-anchor="middle">No measured track for this allele</text>`;
    return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(names[track])}">${body}</svg>`;
  }
  function histogramSVG(variant,key) {
    const vals=Pilot.normalized(variant[key]), wt=Pilot.normalized(lookup.wt[key]);
    const W=440,H=175,L=38,R=12,T=12,B=30,max=Math.max(.05,...[...vals,...wt].filter(Pilot.finite))*1.1;
    const cw=(W-L-R)/vals.length, y=v=>H-B-v/max*(H-T-B);
    let body="";
    for(const tick of [0,max/2,max])body+=`<line x1="${L}" x2="${W-R}" y1="${y(tick)}" y2="${y(tick)}" stroke="#33434b"/><text x="${L-6}" y="${y(tick)+3}" text-anchor="end">${fmt(tick*100,0)}%</text>`;
    vals.forEach((v,i)=>{if(Pilot.finite(v))body+=`<rect x="${L+i*cw}" y="${y(v)}" width="${cw-1}" height="${H-B-y(v)}" fill="#65c7a5" opacity=".65"/>`;});
    body+=`<path d="${Pilot.path(wt,i=>L+(i+.5)*cw,y)}" fill="none" stroke="#d7e0e6" stroke-width="1.5" stroke-dasharray="4 3"/>`;
    for(const v of [0,.25,.5,.75,1])body+=`<text x="${L+v*(W-L-R)}" y="${H-8}" text-anchor="middle">${v}</text>`;
    return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" role="img" aria-label="Per-molecule distribution">${body}</svg>`;
  }
  function moleculesSVG(variant) {
    const rows=[...lookup.wt.examples.map(v=>({...v,group:"WT"})),...(variant.class==="WT"?[]:variant.examples.map(v=>({...v,group:variant.source_variant})))];
    if(!rows.length)return '<p class="notice">No molecule examples available.</p>';
    const W=900,L=125,R=12,H=rows.length*32+32,n=data.region.sequence.length,x=p=>L+p/n*(W-L-R);
    let body="";
    rows.forEach((m,i)=>{
      const y=10+i*32;
      body+=`<g><title>${esc(m.read_id)} · strand ${esc(m.strand)} · global m6A/AT ${fmt(m.global_methylation_fraction)}</title><text x="2" y="${y+13}">${esc(m.group)} · ${esc(m.strand)}</text><rect x="${L}" y="${y}" width="${W-L-R}" height="20" fill="#293740"/>`;
      for(const track of ["msp","nucleosome","tf_small"]){
        if(m.tracks[track]===null)continue;
        for(const [a,b] of m.tracks[track])body+=`<rect x="${x(a)}" y="${y+(track==="tf_small"?6:2)}" width="${Math.max(1,x(b)-x(a))}" height="${track==="tf_small"?8:16}" fill="${colors[track]}" opacity="${track==="msp"?.5:.9}"/>`;
      }
      body+="</g>";
    });
    if(variant.class==="SNV"){const p=variant.oligo_pos1-data.region.oligo_start1;body+=`<line x1="${x(p+.5)}" x2="${x(p+.5)}" y1="5" y2="${H-25}" stroke="white" stroke-dasharray="2 3"/>`;}
    body+=`<text x="${L}" y="${H-5}">${data.region.oligo_start1}</text><text x="${W-R}" y="${H-5}" text-anchor="end">${data.region.oligo_start1+n-1} · oligo position</text>`;
    return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" role="img" aria-label="Individual molecules">${body}</svg>`;
  }
  function renderDetail() {
    const v=current(), mpra=data.mpra.filter(r=>r.variant_id===v.variant_id);
    $("detail-title").textContent=`${v.source_variant} · chromatin profiles`;
    $("detail-subtitle").textContent=`${fmt(v.n_reads,0)} measured molecules · ${fmt(v.n_forward,0)} forward / ${fmt(v.n_reverse,0)} reverse · oligo coordinates`;
    $("variant-info").innerHTML=`<strong>${esc(v.source_variant)}</strong><br>${v.class==="SNV"?`Oligo position ${v.oligo_pos1}<br>Genomic position ${v.genomic_pos1.toLocaleString()} (1-based)<br>`:"Reference sequence<br>"}MSP mean: ${fmt(v.summary.msp_fraction_mean)}<br>Common-site Δ MSP: ${fmt(v.summary.msp_common_site_delta)}${v.n_reads<minimum()?'<p class="warning">Below the display coverage threshold. Inspect descriptively.</p>':""}`;
    $("mpra").innerHTML=mpra.length?`<table><thead><tr><th>Source</th><th>log₂ effect</th><th>Reported p</th></tr></thead><tbody>${mpra.map(r=>`<tr><td>${esc(r.element)}</td><td>${fmt(r.effect_log2)}</td><td>${Pilot.pvalue(r.p_value)}</td></tr>`).join("")}</tbody></table><p>External assay; measurements kept separate.</p>`:'<p>No allele-specific published measurement.</p>';
    $("profiles").innerHTML=Pilot.tracks.map(t=>profileSVG(v,t)).join("");
    $("msp-distribution").innerHTML=histogramSVG(v,"msp_fraction_histogram");
    $("meth-distribution").innerHTML=histogramSVG(v,"methylation_histogram");
    $("molecules").innerHTML=moleculesSVG(v);
    const compatible=comparison?[...datasets.values()].filter(d=>d.library_id===data.library_id && d.reference_sequence_sha256===data.reference_sequence_sha256 && d.region.start0===data.region.start0 && d.region.end0===data.region.end0):[];
    $("across-runs").innerHTML=compatible.length<2?"":`<h3>This allele across runs</h3><table><thead><tr><th>Run / reaction</th><th>Molecules</th><th>Δ MSP</th><th>Δ small footprints</th></tr></thead><tbody>${compatible.map(d=>{
      const other=d.variants.find(o=>o.variant_id===v.variant_id);
      return `<tr><td><button data-sample="${esc(d.sample.id)}">${esc(d.sample.id)}</button><br>${esc(reaction(d.sample))}</td><td>${fmt(other?.n_reads,0)}</td><td>${fmt(other?.summary.msp_common_site_delta)}</td><td>${fmt(other?.summary.tf_small_common_site_delta)}</td></tr>`;
    }).join("")}</tbody></table><p>Each difference uses its own run's WT and common-site mask. Estimates are descriptive; reaction conditions and coverage differ.</p>`;
    $("across-runs").querySelectorAll("button").forEach(button=>button.addEventListener("click",()=>{const label=v.source_variant;$("sample").value=button.dataset.sample;activate(button.dataset.sample);const target=lookup.ordered.find(o=>o.source_variant===label);if(target)choose(target.variant_id);}));
  }
  function download(text,filename,type) {
    const url=URL.createObjectURL(new Blob([text],{type})), anchor=document.createElement("a");
    anchor.href=url;anchor.download=filename;anchor.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  $("file").addEventListener("change",async event=>{
    for(const file of event.target.files){try{install(JSON.parse(await file.text()));}catch(error){status(`${file.name}: ${error.message}`,true);}}
    event.target.value="";
  });
  $("sample").addEventListener("change",()=>activate($("sample").value));
  $("search").addEventListener("input",()=>{if(data)renderOptions();});
  $("variant").addEventListener("change",()=>choose($("variant").value));
  $("wildtype").addEventListener("click",()=>choose(lookup.wt.variant_id));
  for(const [id,delta] of [["previous",-1],["next",1]])$(id).addEventListener("click",()=>{
    const ids=[...$("variant").options].map(o=>o.value), i=ids.indexOf(selectedId);
    if(ids.length)choose(ids[Math.max(0,Math.min(ids.length-1,i+delta))]);
  });
  $("metric").addEventListener("change",drawMap);
  $("minimum").addEventListener("input",()=>{if(data){drawMap();renderDetail();}});
  $("mode").addEventListener("change",()=>{if(data)renderDetail();});
  $("alleles").addEventListener("click",event=>{const v=mapVariant(event);if(v)choose(v.variant_id);});
  $("alleles").addEventListener("mousemove",event=>{const v=mapVariant(event);$("alleles").title=v?`${v.source_variant} · n=${v.n_reads}\n${$("metric").selectedOptions[0].text}: ${fmt(Pilot.metric(v,$("metric").value))}`:"Reference allele or outside the design";});
  window.addEventListener("resize",drawMap);
  $("export").addEventListener("click",()=>{const v=current();download(JSON.stringify({schema_version:data.schema_version,sample:data.sample,region:data.region,method:data.method,sampling:data.sampling,limitations:data.limitations,provenance:data.provenance,variant:v,mpra:data.mpra.filter(r=>r.variant_id===v.variant_id)},null,2),`${data.sample.id}_${v.source_variant}.json`,"application/json");});
  $("export-summary").addEventListener("click",()=>{
    const rows=data.variants.map(v=>v.summary), fields=Object.keys(rows[0]);
    const quote=x=>'"'+String(x??"").replace(/"/g,'""')+'"';
    download([fields.map(quote).join(","),...rows.map(r=>fields.map(k=>quote(r[k])).join(","))].join("\n"),`${data.sample.id}_summary.csv`,"text/csv");
  });
  $("export-profile").addEventListener("click",()=>{
    const pieces=[...$("profiles").querySelectorAll("svg")].map((svg,i)=>`<g transform="translate(0,${i*130})">${svg.innerHTML}</g>`).join("");
    const svg=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 900 420"><rect width="900" height="420" fill="#182229"/><style>text{font:10px sans-serif;fill:#b6c5cf}.track-title{font-size:12px;fill:#e2ebef}</style>${pieces}<text x="46" y="413">${esc(data.sample.id)} · ${esc(current().source_variant)} · descriptive pilot; no significance tests</text></svg>`;
    download(svg,`${data.sample.id}_${current().source_variant}_profiles.svg`,"image/svg+xml");
  });
  try {const embedded=JSON.parse($("pilot-data").textContent);if(embedded)install(embedded);}catch(error){status(error.message,true);}
})();

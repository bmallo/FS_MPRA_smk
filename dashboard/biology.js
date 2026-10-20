'use strict';
(()=>{
const B=JSON.parse(document.getElementById('biology-data').textContent),origin=B.coordinates.genomic_start1;
const concordance=new Map(B.motif_concordance.map(r=>[r.hit_id,r]));
const byAudit=new Map(B.barcode_audits.map(r=>[r.variant_id+'|'+r.region_id,r]));
const fields=[['delta_pp','All eligible barcodes'],['split_0_delta_pp','Barcode half A'],['split_1_delta_pp','Barcode half B'],['molecule_weighted_delta_pp','Molecule-weighted'],['depth_ge3_delta_pp','Barcodes with ≥3 molecules'],['depth_ge5_delta_pp','Barcodes with ≥5 molecules'],['shared_at_radius0_delta_pp','Unchanged A/T bases; mutation excluded'],['shared_at_radius2_delta_pp','Unchanged A/T bases; mutation ±2 bp excluded'],['distal_m6a_standardized_delta_pp','Standardized distal m6A quartiles'],['distal_nuc_standardized_delta_pp','Standardized distal protection quartiles']];
for(const r of B.regions)$('biologyRegion').appendChild(option(r.id,r.id+' · '+r.genomic_start1+'–'+r.genomic_end1));
function focus(){
 const chosen=$('biologyRegion').value,v=selected(),p=v.relative_pos0;
 if(chosen==='whole_tf')return {id:'whole_tf',start0:0,end0:211};
 if(chosen!=='follow')return B.regions.find(r=>r.id===chosen);
 if(p===null)return {id:'whole_tf',start0:0,end0:211};
 return B.regions.find(r=>r.start0<=p&&p<r.end0)||{id:'whole_tf',start0:Math.max(0,p-18),end0:Math.min(211,p+19)};
}
function motifDelta(h,v){if(v.relative_pos0===null||v.relative_pos0<h.start0||v.relative_pos0>=h.end0)return null;const w=h.weights_genomic[v.relative_pos0-h.start0];return w['ACGT'.indexOf(v.alt)]-w['ACGT'.indexOf(v.ref)];}
function update(){
 const v=selected(),r=focus(),a=Math.max(0,r.start0-4),b=Math.min(211,r.end0+4),span=b-a;
 const query=$('biologyMotifQuery').value.toLowerCase(),db=$('biologyDatabase').value;
 const hits=B.motifs.filter(h=>h.start0<b&&h.end0>a&&(db==='all'||h.database===db)&&(h.name+' '+h.matrix_id).toLowerCase().includes(query));
 const c=$('biologySequence'),ctx=c.getContext('2d'),L=65,R=20,W=c.width-L-R,x=p=>L+(p-a+.5)/span*W;
 ctx.clearRect(0,0,c.width,c.height);ctx.font='12px system-ui';ctx.textAlign='center';
 ctx.fillStyle='#526674';
 for(let p=a;p<b;p++){if((p-a)%Math.max(1,Math.ceil(span/8))===0)ctx.fillText(String(origin+p),x(p),17);const base=D.sequence[p];ctx.fillStyle='AT'.includes(base)?'#167b77':'#6b7c87';if(span<=60)ctx.fillText(base,x(p),42);else ctx.fillRect(x(p)-1.5,29,3,'AT'.includes(base)?15:5);}
 if(v.relative_pos0!==null&&v.relative_pos0>=a&&v.relative_pos0<b){ctx.fillStyle='#da6838';ctx.fillText(v.alt,x(v.relative_pos0),62);ctx.strokeStyle='#da6838';ctx.setLineDash([3,3]);ctx.beginPath();ctx.moveTo(x(v.relative_pos0),22);ctx.lineTo(x(v.relative_pos0),c.height-5);ctx.stroke();ctx.setLineDash([]);}
 ctx.textAlign='left';
 for(const [i,h] of hits.slice(0,8).entries()){
  const y=78+i*25,left=x(Math.max(a,h.start0))-.45*W/span,right=x(Math.min(b,h.end0)-1)+.45*W/span;
  ctx.fillStyle=h.database==='JASPAR2024'?'#72529a':'#2376b3';ctx.fillRect(left,y,Math.max(2,right-left),5);
  ctx.fillText((h.strand==='+'?'→ ':'← ')+h.name+' · '+h.matrix_id,L,y+18);
 }
 $('biologySequenceInfo').textContent='GRCh38 chr19:'+String(origin+a)+'–'+String(origin+b-1)+' · '+r.id+' · '+(v.relative_pos0===null?'WT reference':variantLabel(v)+' · A/T '+v.at_change)+' · '+hits.length+' overlapping candidate matches.';
 $('biologyMotifCount').textContent=hits.length+' candidate matches · '+Math.min(8,hits.length)+' drawn';
 $('biologyMotifTable').innerHTML=table(['Collection / ID','Candidate','Genomic interval / strand','WT score (bits)','Uniform tail','GC-matched tail','Variant Δ score (bits)','ρ score / local footprint','ρ score / activity 1','ρ score / activity 2'],hits.map(h=>[h.database+' / '+h.matrix_id,h.name,h.genomic_start1+'–'+h.genomic_end1+' '+h.strand,fmt(h.score_bits),h.sequence_tail_p.toExponential(2),h.gc_matched_tail_p.toExponential(2),v.relative_pos0!==null&&h.start0<=v.relative_pos0&&v.relative_pos0<h.end0?fmt(motifDelta(h,v)):'Outside selected substitution',fmt(concordance.get(h.hit_id).tf_local_delta_pp_spearman),fmt(concordance.get(h.hit_id).published_LDLR_log2_spearman),fmt(concordance.get(h.hit_id).published_LDLR_2_log2_spearman)]));
 $('biologyLiterature').innerHTML=B.literature.map(l=>'Literature anchor: <a href="'+esc(l.url)+'">'+esc(l.name)+'</a> '+esc(l.sequence)+' at chr19:'+l.genomic_start1+'–'+l.genomic_end1+'. Exact sequence placement; binding in this preparation is unverified.').join(' ');
 const endpoint=$('biologyAuditEndpoint').value==='local_tf'?r.id:$('biologyAuditEndpoint').value;
 const audit=byAudit.get(v.variant_id+'|'+endpoint);
 $('biologyAuditTitle').textContent=variantLabel(v)+' · '+endpoint+' · sensitivity checks';
 $('biologyAudit').innerHTML=audit?table(['Comparison','Variant − WT (percentage points)'],fields.map(([key,label])=>[label,fmt(audit[key],3)])):'Select a measured SNV to inspect its sensitivity checks.';
 $('biologyAuditSupport').textContent=audit?audit.nonzero_barcodes+' / '+audit.n_barcodes+' variant barcodes have nonzero coverage in the audited region; largest single-barcode deletion changes its mean by '+fmt(audit.max_leave_one_barcode_change_pp,3)+' pp. Barcode counts: halves '+audit.split_0_barcodes+'/'+audit.split_1_barcodes+', depth ≥3 '+audit.depth_ge3_barcodes+', depth ≥5 '+audit.depth_ge5_barcodes+'. Shared A/T positions: '+audit.shared_at_radius0_positions+' (mutation excluded), '+audit.shared_at_radius2_positions+' (±2 bp excluded). Minimum variant support per distal bin: m6A '+audit.distal_m6a_min_bin_barcodes+', protection '+audit.distal_nuc_min_bin_barcodes+'.'+(endpoint.startsWith('whole_')?' Audit endpoint is the entire 211-bp element, even when the sequence view is zoomed.':''):'';
 const local=B.region_variants.filter(row=>row.region_id===r.id);
 $('biologyLocalVariants').innerHTML=local.length?table(['Variant','A/T change','Local small footprint Δ pp','Promoter accessibility Δ pp','Nucleosome-sized Δ pp','Published LDLR log₂','Published LDLR_2 log₂'],local.map(row=>[row.genomic_variant+' ('+row.source_variant+')',row.at_change,fmt(row.delta_pp),fmt(row.msp_delta_pp),fmt(row.nuc_85_250_delta_pp),fmt(row.published_LDLR_log2),fmt(row.published_LDLR_2_log2)])):'Choose a WT neighborhood or the adjacent-base region to see every local substitution.';
 c.setAttribute('aria-label',$('biologySequenceInfo').textContent);
}
$('biologyRegionSummary').innerHTML=table(['Region','Genomic interval','WT small-footprint coverage','Local variants measured','A/T positions','Lower / higher local footprint','Median accessibility Δ pp','Median nucleosome-sized Δ pp'],B.region_summaries.map(r=>[r.region_id,r.genomic_start1+'–'+r.genomic_end1,percent(r.wt_tf_coverage),r.n_measured+'/'+r.n_designed_snvs,r.at_positions,r.delta_pp_negative+' / '+r.delta_pp_positive,fmt(r.msp_delta_pp_median),fmt(r.nuc_85_250_delta_pp_median)]));
for(const id of ['biologyRegion','biologyDatabase','biologyAuditEndpoint'])$(id).addEventListener('change',update);
$('biologyMotifQuery').addEventListener('input',update);
for(const id of ['previousVariant','nextVariant','functionalLollipop','variantTable','heatmap'])$(id).addEventListener('click',update);
for(const id of ['variant','region'])$(id).addEventListener('change',update);
$('comparisonPanel').addEventListener('keydown',update);
update();
})();

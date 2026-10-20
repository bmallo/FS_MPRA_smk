'use strict';
const D=JSON.parse(document.getElementById('raw-data').textContent),$=id=>document.getElementById(id);
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt=x=>Number.isFinite(x)?x.toFixed(3):'Not measured';
const table=(headers,rows)=>'<table><thead><tr>'+headers.map(x=>'<th>'+esc(x)+'</th>').join('')+'</tr></thead><tbody>'+rows.map(r=>'<tr>'+r.map(x=>'<td>'+esc(x)+'</td>').join('')+'</tr>').join('')+'</tbody></table>';
const option=(v,t)=>{const o=document.createElement('option');o.value=v;o.textContent=t;return o;};
const records=D.comparisons.map(r=>Object.fromEntries(D.comparison_fields.map((k,i)=>[k,r[i]])));
const byComparison=new Map(records.map(r=>[[r.source_variant,r.scenario,r.region,r.metric].join('|'),r]));
const variants=D.variants.filter(v=>v.source_variant!=='WT');
for(const v of variants)$('variant').appendChild(option(v.source_variant,'chr19:'+v.genomic_pos1+':'+v.ref+'>'+v.alt+' ('+v.source_variant+')'));
for(const r of D.regions)$('region').appendChild(option(r.id,r.id+' · '+(11089285+r.start0)+'–'+(11089284+r.end0)));
for(const s of D.scenarios)$('scenario').appendChild(option(s.id,s.label));
function record(s,metric){return byComparison.get([$('variant').value,s,$('region').value,metric].join('|'));}
function draw(id,field,base=null){
 const canvas=$(id),ctx=canvas.getContext('2d'),region=D.regions.find(r=>r.id===$('region').value),v=variants.find(r=>r.source_variant===$('variant').value),p=v.oligo_pos1-21;
 const a=region.start0,b=region.end0,L=55,R=18,T=18,B=32,x=i=>L+(i-a)/Math.max(1,b-a-1)*(canvas.width-L-R),y=z=>canvas.height-B-z*(canvas.height-B-T);
 ctx.clearRect(0,0,canvas.width,canvas.height);ctx.font='12px system-ui';
 const left=Math.max(a,p-2),right=Math.min(b-1,p+2);if(left<=right){ctx.fillStyle='#eff0ed';ctx.fillRect(x(left),T,Math.max(2,x(right)-x(left)),canvas.height-B-T);}
 for(let j=0;j<=4;j++){ctx.strokeStyle='#dce4e8';ctx.beginPath();ctx.moveTo(L,y(j/4));ctx.lineTo(canvas.width-R,y(j/4));ctx.stroke();ctx.fillStyle='#526674';ctx.fillText((25*j)+'%',2,y(j/4)+4);}
 for(let i=a;i<b;i+=Math.max(1,Math.ceil((b-a)/7))){ctx.fillStyle='#526674';ctx.textAlign=i===a?'left':'center';ctx.fillText(String(11089285+i),x(i),canvas.height-10);}ctx.textAlign='left';
 for(const [name,color] of [['WT','#253d4b'],[v.source_variant,'#da6838']]){
  const profile=D.profiles[$('scenario').value+'|'+name][field];ctx.strokeStyle=color;ctx.fillStyle=color;ctx.lineWidth=1.5;ctx.beginPath();let previous=false;
  for(let i=a;i<b;i++){const z=profile[i];if(!Number.isFinite(z)||(base&&(D.sequence[i]!==base||Math.abs(i-p)<=2)))continue;
   if(!previous){ctx.moveTo(x(i),y(z));previous=true;}else ctx.lineTo(x(i),y(z));
  }ctx.stroke();
  if(base)for(let i=a;i<b;i++){const z=profile[i];if(!Number.isFinite(z)||D.sequence[i]!==base||Math.abs(i-p)<=2)continue;ctx.beginPath();ctx.arc(x(i),y(z),3,0,2*Math.PI);ctx.fill();}
 }
 canvas.setAttribute('aria-label',field+' for '+v.source_variant+' versus WT; '+region.id);
}
function update(){
 const scenario=$('scenario').value,raw=record(scenario,'retained_m6a_shared');
 const retention=D.retention.find(r=>r.scenario===scenario&&r.source_variant===$('variant').value),wt=D.retention.find(r=>r.scenario===scenario&&r.source_variant==='WT');
 $('support').textContent=retention.molecules+' selected variant molecules / '+retention.barcodes+' barcode labels; WT '+wt.molecules+' molecules / '+wt.barcodes+' barcode labels. Weighting: '+retention.weighting+'. Barcode labels in molecule-weighted groups are counts, not eligibility requirements.';
 draw('rawA','retained_m6a_fraction','A');draw('rawT','retained_m6a_fraction','T');draw('footprint',$('track').value);
 const metrics=[['retained_m6a_shared','Shared-site retained m6A'],['retained_m6a_shared_A','Shared reference A'],['retained_m6a_shared_T','Shared reference T'],['low_call_fraction_le25','Molecules with ≤25% retained calls'],['tf_small','Small footprints'],['msp','MSP accessibility'],['nuc_85_250','Nucleosome-sized'],['nuc_gt250','Long protection']];
 $('selected').innerHTML=table(['Measurement','Positions','WT fraction','Variant fraction','Difference (pp)','WT / variant callable molecules','WT / variant callable barcodes'],metrics.map(([key,label])=>{const r=record(scenario,key);return[label,r.shared_opportunities,fmt(r.wt_value),fmt(r.variant_value),fmt(r.delta_pp),r.wt_molecules+' / '+r.variant_molecules,r.wt_barcodes+' / '+r.variant_barcodes];}));
 $('sensitivity').innerHTML=table(['Comparison','Raw m6A Δ pp','Small footprint Δ pp','Accessibility Δ pp','Nucleosome-sized Δ pp','Callable raw WT / variant molecules'],D.scenarios.map(s=>{const r=record(s.id,'retained_m6a_shared');return[s.label,fmt(r.delta_pp),fmt(record(s.id,'tf_small').delta_pp),fmt(record(s.id,'msp').delta_pp),fmt(record(s.id,'nuc_85_250').delta_pp),r.wt_molecules+' / '+r.variant_molecules];}));
 const res=D.resampling.filter(r=>r.source_variant===$('variant').value&&r.region===$('region').value);
 $('resampling').innerHTML=res.length?table(['Method','Measurement','Repeats','Min Δ pp','Median Δ pp','Max Δ pp','Negative / positive'],res.map(r=>[r.kind,r.metric,r.n,fmt(r.minimum),fmt(r.median),fmt(r.maximum),r.negative+' / '+r.positive])):'Repeated-selection summaries are available for position 125 in WT_P3 and adjacent-region variants across the whole promoter.';
}
for(const id of ['variant','scenario','region','track'])$(id).addEventListener('change',update);
for(const [id,step] of [['previous',-1],['next',1]])$(id).addEventListener('click',()=>{const i=variants.findIndex(v=>v.source_variant===$('variant').value);$('variant').value=variants[(i+step+variants.length)%variants.length].source_variant;update();});
$('variant').value='125A>T';$('region').value='WT_P3';$('scenario').value='primary_barcodes';$('track').value='tf_small';update();

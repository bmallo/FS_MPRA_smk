'use strict';
const D=JSON.parse(document.getElementById('atlas-data').textContent);
const $=id=>document.getElementById(id), WT=D.variants[D.wt_index];
const mean=a=>a.length && a.every(Number.isFinite)?a.reduce((s,x)=>s+x,0)/a.length:null;
const fmt=(v,n=2)=>Number.isFinite(v)?v.toFixed(n):'Not measured';
const percent=v=>Number.isFinite(v)?fmt(v*100)+'%':'Not measured';
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const C=D.coordinates;
const baseLabel=p=>String(C.genomic_start1+p);
const variantLabel=v=>v.relative_pos0===null?'WT':v.genomic_variant+' (oligo '+v.source_variant+')';
const intervalLabel=r=>C.chrom+':'+r.genomic_start1+'–'+r.genomic_end1;
function tick(ctx,p,x,y){ctx.textAlign=p===0?'left':p===210?'right':'center';ctx.fillText(baseLabel(p),x,y);ctx.textAlign='left';}
const stateNames=['None of these three','Small TF only','Nucleosome-sized only','Small TF + nucleosome-sized','Long only','Small TF + long','Nucleosome-sized + long','All three'];
const option=(value,text)=>{const o=document.createElement('option');o.value=value;o.textContent=text;return o;};
const table=(headers,rows)=>'<table><thead><tr>'+headers.map(x=>'<th>'+esc(x)+'</th>').join('')+'</tr></thead><tbody>'+rows.map(row=>'<tr>'+row.map(x=>'<td>'+esc(x)+'</td>').join('')+'</tr>').join('')+'</tbody></table>';
function drawLines(id,series,mutation=null,window=null){
 const c=$(id),ctx=c.getContext('2d'),w=c.width,h=c.height,L=55,R=15,T=15,B=34;
 ctx.clearRect(0,0,w,h);ctx.font='12px system-ui';
 const values=series.flatMap(s=>s.values).filter(Number.isFinite),max=Math.max(.01,...values)*1.12;
 const x=p=>L+p/210*(w-L-R),y=v=>h-B-v/max*(h-T-B);
 if(window){ctx.fillStyle='#eaf0f4';ctx.fillRect(x(window.start0),T,x(window.end0-1)-x(window.start0),h-B-T);}
 for(let j=0;j<=4;j++){const v=j*max/4,yy=y(v);ctx.strokeStyle='#dce4e8';ctx.beginPath();ctx.moveTo(L,yy);ctx.lineTo(w-R,yy);ctx.stroke();ctx.fillStyle='#526674';ctx.fillText((v*100).toFixed(1)+'%',1,yy+4);}
 for(let p=0;p<=210;p+=30){ctx.fillStyle='#526674';tick(ctx,p,x(p),h-12);}
 for(const s of series){ctx.strokeStyle=s.color;ctx.lineWidth=2;ctx.beginPath();let connected=false;for(let p=0;p<s.values.length;p++){const v=s.values[p];if(!Number.isFinite(v)){connected=false;continue;}if(!connected){ctx.moveTo(x(p),y(v));connected=true;}else ctx.lineTo(x(p),y(v));}ctx.stroke();}
 if(mutation!==null){ctx.strokeStyle='#8d5c28';ctx.setLineDash([4,3]);ctx.beginPath();ctx.moveTo(x(mutation),T);ctx.lineTo(x(mutation),h-B);ctx.stroke();ctx.setLineDash([]);}
}
function metric(v){
 const id=$('region').value;
 if(id.startsWith('r'))return v.region_effects[Number(id.slice(1))];
 const t=Number(id.slice(1)),a=mean(v.coverage[t]),b=mean(WT.coverage[t]);return a===null?null:100*(a-b);
}
function selectedWindow(){const id=$('region').value;return id.startsWith('r')?D.regions[Number(id.slice(1))]:null;}
function color(v,max){if(!Number.isFinite(v))return '#e3e7ea';const a=Math.min(1,Math.abs(v)/max),target=v>=0?[185,53,54]:[38,105,164];return 'rgb('+target.map(x=>Math.round(250+(x-250)*a)).join(',')+')';}
let heatScale=1;
const cellMap=new Map(D.variants.filter(v=>v.relative_pos0!==null).map(v=>[v.relative_pos0+':'+v.alt,v]));
function heatmap(){
 const c=$('heatmap'),ctx=c.getContext('2d'),L=45,R=10,T=10,H=88,cw=(c.width-L-R)/211;
 ctx.clearRect(0,0,c.width,c.height);heatScale=Math.max(.01,...D.variants.map(metric).filter(Number.isFinite).map(Math.abs));
 for(let a=0;a<4;a++){ctx.fillStyle='#183445';ctx.font='13px system-ui';ctx.fillText('ACGT'[a],15,T+a*22+16);for(let p=0;p<211;p++){const v=cellMap.get(p+':'+'ACGT'[a]);ctx.fillStyle=color(v?metric(v):null,heatScale);ctx.fillRect(L+p*cw,T+a*22,cw+.1,21);}}
 ctx.fillStyle='#526674';for(let p=0;p<=210;p+=30)tick(ctx,p,L+(p+.5)*cw,T+H+19);
 ctx.fillText('Scale: −'+fmt(heatScale)+' to +'+fmt(heatScale)+' percentage points',L,T+H+42);
}
function heatVariant(e){const c=$('heatmap'),b=c.getBoundingClientRect(),x=(e.clientX-b.left)*c.width/b.width,y=(e.clientY-b.top)*c.height/b.height;const p=Math.floor((x-45)/((c.width-55)/211)),a=Math.floor((y-10)/22);return p>=0&&p<211&&a>=0&&a<4?cellMap.get(p+':'+'ACGT'[a]):null;}
function updateTable(){
 const query=$('query').value.toLowerCase();const rows=D.variants.filter(v=>v.source_variant!=='WT'&&variantLabel(v).toLowerCase().includes(query)&&($('support').value!=='floor'||v.coverage_floor_met));
 rows.sort((a,b)=>{const x=metric(a),y=metric(b);return (y===null?-1:Math.abs(y))-(x===null?-1:Math.abs(x))||a.oligo_pos1-b.oligo_pos1;});
 $('tableCount').textContent=rows.length+' variants · sorted by absolute observed difference, not confidence. All rows remain available by scrolling.';
 $('variantTable').innerHTML=rows.map(v=>'<tr class="selectable" data-id="'+esc(v.variant_id)+'"><td>'+esc(variantLabel(v))+'</td><td>'+fmt(metric(v))+'</td><td>'+v.n_barcodes.toLocaleString()+'</td><td>'+v.n_molecules.toLocaleString()+'</td><td>'+esc(v.at_change)+'</td><td>'+(v.coverage_floor_met?'Met':'Low / missing')+'</td></tr>').join('');
}
function selected(){return D.variants[Number($('variant').value)];}
const navigation=D.variants.map((v,i)=>({v,i})).filter(({v})=>v.relative_pos0!==null)
 .sort((a,b)=>a.v.relative_pos0-b.v.relative_pos0||'ACGT'.indexOf(a.v.alt)-'ACGT'.indexOf(b.v.alt)).map(({i})=>i);
function stepVariant(direction){
 const current=navigation.indexOf(Number($('variant').value));
 const next=current<0?(direction>0?0:navigation.length-1):(current+direction+navigation.length)%navigation.length;
 $('variant').value=navigation[next];updateVariant();
}
let lollipopPoints=[];
function drawLollipop(){
 const c=$('functionalLollipop'),ctx=c.getContext('2d'),source=$('functionalSource').value,v=selected();
 const L=55,R=15,T=20,B=34,w=c.width,h=c.height;
 const records=D.variants.filter(a=>a.relative_pos0!==null&&Number.isFinite(a.functional[source]));
 const limit=Math.max(.5,...records.map(a=>Math.abs(a.functional[source])))*1.12;
 const x=p=>L+p/210*(w-L-R),y=value=>T+(limit-value)/(2*limit)*(h-T-B);
 ctx.clearRect(0,0,w,h);ctx.font='12px system-ui';
 for(let j=-2;j<=2;j++){const value=limit*j/2,yy=y(value);ctx.strokeStyle=j===0?'#9caeb8':'#e3e9ed';ctx.beginPath();ctx.moveTo(L,yy);ctx.lineTo(w-R,yy);ctx.stroke();ctx.fillStyle='#526674';ctx.fillText(value.toFixed(1),3,yy+4);}
 for(let p=0;p<=210;p+=30){ctx.fillStyle='#526674';tick(ctx,p,x(p),h-12);}
 ctx.fillStyle='#526674';ctx.fillText('log₂ effect',3,12);
 const window=selectedWindow();
 if(window){ctx.fillStyle='rgba(80,110,130,0.06)';ctx.fillRect(x(window.start0),T,x(window.end0-1)-x(window.start0),h-T-B);}
 lollipopPoints=records.map(a=>({variant:a,x:x(a.relative_pos0),y:y(a.functional[source])}));
 function point(p,highlight){ctx.strokeStyle=highlight?'#da6838':'rgba(109,129,141,0.42)';ctx.lineWidth=highlight?2.5:1;ctx.beginPath();ctx.moveTo(p.x,y(0));ctx.lineTo(p.x,p.y);ctx.stroke();ctx.fillStyle=highlight?'#da6838':'#96a8b2';ctx.beginPath();ctx.arc(p.x,p.y,highlight?5.5:2,0,2*Math.PI);ctx.fill();if(highlight){ctx.strokeStyle='#fff';ctx.lineWidth=1.5;ctx.stroke();}}
 for(const p of lollipopPoints)point(p,false);
 const highlighted=lollipopPoints.find(p=>p.variant.variant_id===v.variant_id);
 if(v.relative_pos0!==null){ctx.strokeStyle='#8d5c28';ctx.setLineDash([4,3]);ctx.beginPath();ctx.moveTo(x(v.relative_pos0),T);ctx.lineTo(x(v.relative_pos0),h-B);ctx.stroke();ctx.setLineDash([]);}
 if(highlighted)point(highlighted,true);
 const value=v.functional[source];
 $('functionalSelected').textContent=v.relative_pos0===null?'WT reference selected; lollipops represent SNV effects relative to the published reference.':
  variantLabel(v)+' · '+source+' · '+(Number.isFinite(value)?(value>0?'+':'')+fmt(value,3)+' log₂ activity effect':'functional effect not measured; no lollipop drawn');
 $('functionalHover').textContent='Hover over a lollipop for its variant and effect; click it to select that variant.';
 c.setAttribute('aria-label',source+' published log2 effects. '+$('functionalSelected').textContent);
}
function lollipopAt(event){
 const c=$('functionalLollipop'),b=c.getBoundingClientRect(),x=(event.clientX-b.left)*c.width/b.width,y=(event.clientY-b.top)*c.height/b.height;
 let nearest=null,distance=100;
 for(const p of lollipopPoints){const d=(x-p.x)**2+(y-p.y)**2;if(d<distance){distance=d;nearest=p;}}
 return nearest;
}
function updateVariant(){
 const v=selected(),t=Number($('track').value),a=mean(v.coverage[t]),b=mean(WT.coverage[t]);
 const order=navigation.indexOf(Number($('variant').value));
 $('variantCounter').textContent=order<0?'WT reference · Next starts at the first SNV.':'Variant '+(order+1)+' of '+navigation.length+' · '+variantLabel(v);
 drawLollipop();
 $('selectedInfo').textContent=variantLabel(v)+' · '+v.n_barcodes.toLocaleString()+' barcodes · '+v.n_molecules.toLocaleString()+' molecules'+(v.oligo_pos1!==null?' · '+C.assembly+' forward strand · A/T '+v.at_change:'')+(v.coverage_floor_met?'':' · LOW OR MISSING COVERAGE');
 drawLines('comparison',[{values:WT.coverage[t],color:'#253d4b'},{values:v.coverage[t],color:'#da6838'}],v.relative_pos0,selectedWindow());
 $('selectedStats').textContent='Whole-element covered fraction: WT '+percent(b)+'; variant '+percent(a)+'; difference '+fmt(a===null?null:100*(a-b))+' pp. Any overlapping call: WT '+percent(WT.presence[t])+'; variant '+percent(v.presence[t])+'. Shading marks the selected region; dashed line marks the mutation.';
 const features=[['whole_m6a_fraction','Whole-plasmid observed m6A'],['distal_m6a_fraction','Distal observed m6A'],['distal_nucleosome_fraction','Distal nucleosome-class coverage']];
 $('context').innerHTML=table(['Measurement','WT','Variant'],features.map(([key,label])=>[label,percent(WT.context[key]),percent(v.context[key])]));
 $('functional').innerHTML=table(['Source','Published log₂ effect'],D.mpra_sources.map(s=>[s.id,fmt(v.functional[s.id],3)]));
}
function selectVariant(v){if(!v)return;$('variant').value=D.variants.indexOf(v);updateVariant();}
$('cohort').textContent=D.sample_id+' · '+C.assembly+' '+C.chrom+':'+C.genomic_start1+'–'+C.genomic_end1+' (+) · 211 bp · '+WT.n_barcodes.toLocaleString()+' WT barcodes / '+WT.n_molecules.toLocaleString()+' WT molecules · '+(D.variants.length-1)+' designed SNVs';
for(const [id,legend,tracks] of [['wtSmall','wtLegend1',[0,1,2]],['wtProtection','wtLegend2',[3,4,5,6]],['wtAccess','wtLegend3',[7]]]){
 $(legend).innerHTML=tracks.map(t=>'<span><i class="key" style="background:'+D.colors[t]+'"></i>'+esc(D.track_labels[t])+'</span>').join('');drawLines(id,tracks.map(t=>({values:WT.coverage[t],color:D.colors[t]})));
}
$('wtCategories').innerHTML=table(['Category','Mean coverage','Any call'],D.tracks.map((t,i)=>[D.track_labels[i],percent(mean(WT.coverage[i])),percent(WT.presence[i])]));
$('wtStates').innerHTML=table(['Combination','Molecule fraction'],stateNames.map((s,i)=>[s,percent(WT.states[i])]));
$('regionTable').innerHTML=table(['ID','Genomic peak (1-based)','Genomic interval (inclusive)','WT mean coverage'],D.regions.filter(r=>r.kind==='WT_peak_neighborhood').map(r=>[r.id,r.genomic_peak1,intervalLabel(r),percent(mean(WT.coverage[0].slice(r.start0,r.end0)))]));
D.regions.forEach((r,i)=>$('region').appendChild(option('r'+i,r.id+' · '+intervalLabel(r))));
for(const t of [0,7,3,4,5,6])$('region').appendChild(option('t'+t,'Whole element · '+D.track_labels[t]));
D.variants.forEach((v,i)=>$('variant').appendChild(option(i,variantLabel(v)+' · '+v.n_barcodes+' barcodes')));
D.tracks.forEach((t,i)=>$('track').appendChild(option(i,D.track_labels[i])));
D.mpra_sources.forEach(s=>$('functionalSource').appendChild(option(s.id,s.id)));
$('variant').value=navigation[0];$('track').value='0';
$('previousVariant').addEventListener('click',()=>stepVariant(-1));
$('nextVariant').addEventListener('click',()=>stepVariant(1));
$('comparisonPanel').addEventListener('keydown',event=>{
 if(event.altKey||event.ctrlKey||event.metaKey||event.shiftKey||['INPUT','SELECT','TEXTAREA','BUTTON'].includes(event.target.tagName)||event.target.isContentEditable)return;
 if(event.key==='ArrowLeft'||event.key==='ArrowRight'){event.preventDefault();stepVariant(event.key==='ArrowRight'?1:-1);}
});
$('functionalSource').addEventListener('change',drawLollipop);
$('functionalLollipop').addEventListener('mousemove',event=>{const p=lollipopAt(event);$('functionalHover').textContent=p?variantLabel(p.variant)+' · '+fmt(p.variant.functional[$('functionalSource').value],3)+' log₂ · '+p.variant.n_barcodes+' Fiber-seq barcodes':'Hover over a lollipop for its variant and effect; click it to select that variant.';});
$('functionalLollipop').addEventListener('click',event=>{const p=lollipopAt(event);if(p)selectVariant(p.variant);});
$('region').addEventListener('change',()=>{heatmap();updateTable();updateVariant();});
for(const id of ['query','support'])$(id).addEventListener(id==='query'?'input':'change',updateTable);
for(const id of ['variant','track'])$(id).addEventListener('change',updateVariant);
$('variantTable').addEventListener('click',e=>{const row=e.target.closest('tr[data-id]');if(row)selectVariant(D.variants.find(v=>v.variant_id===row.dataset.id));});
$('heatmap').addEventListener('mousemove',e=>{const v=heatVariant(e);$('hover').textContent=v?variantLabel(v)+' · '+fmt(metric(v))+' pp · '+v.n_barcodes+' barcodes / '+v.n_molecules+' molecules · A/T '+v.at_change:'WT base, missing cell, or outside heatmap';});
$('heatmap').addEventListener('click',e=>selectVariant(heatVariant(e)));
$('download').addEventListener('click',()=>{const v=selected(),url=URL.createObjectURL(new Blob([JSON.stringify({descriptive_only:true,coordinates:C,sequence:D.sequence,profile_positions1:Array.from({length:C.length},(_,i)=>C.genomic_start1+i),variant:v,wt:WT,tracks:D.tracks,regions:D.regions},null,2)],{type:'application/json'})),a=document.createElement('a');a.href=url;a.download='LDLR_'+v.source_variant.replace(/[^A-Za-z0-9]/g,'_')+'.json';a.click();URL.revokeObjectURL(url);});
heatmap();updateTable();updateVariant();

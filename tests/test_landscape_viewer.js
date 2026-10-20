'use strict';
// Exercise the real inline controller against a minimal DOM/canvas, including
// navigation into a variant with no measurements. No network or browser server.
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const file=process.argv[2];if(!file)throw Error('Pass atlas.json');
const data=JSON.parse(fs.readFileSync(file,'utf8'));
const nodes=new Map();
const drawnText=new Set();
const context=new Proxy({fillText(text){drawnText.add(String(text));}}, {get:(o,k)=>o[k]||(()=>{}),set:(o,k,v)=>{o[k]=v;return true;}});
function node(id){if(!nodes.has(id))nodes.set(id,{id,value:'',innerHTML:'',textContent:'',width:1100,height:id==='heatmap'?150:id==='functionalLollipop'?230:260,
 children:[],events:{},appendChild(o){this.children.push(o);if(this.children.length===1)this.value=String(o.value);},
 addEventListener(e,fn){const previous=this.events[e];this.events[e]=previous?(...args)=>{previous(...args);fn(...args);}:fn;},setAttribute(k,v){this[k]=v;},getContext(){return context;},getBoundingClientRect(){return{left:0,top:0,width:this.width,height:this.height};},click(){}});return nodes.get(id);}
node('atlas-data').textContent=JSON.stringify(data);node('support').value='all';
const sandbox={document:{getElementById:node,createElement:()=>node('created'+nodes.size)},console,Blob,URL};
vm.createContext(sandbox);vm.runInContext(fs.readFileSync('dashboard/landscape.js','utf8'),sandbox);
if(process.argv[3]){
 node('biology-data').textContent=fs.readFileSync(process.argv[3],'utf8');
 for(const [id,value] of [['biologyRegion','follow'],['biologyDatabase','all'],['biologyAuditEndpoint','local_tf']]){
  node(id).value=value;node(id).children=[{value}];
 }
 vm.runInContext(fs.readFileSync('dashboard/biology.js','utf8'),sandbox);
 assert.match(node('biologyLiterature').innerHTML,/11089391–11089400/);
 assert.match(node('biologyRegionSummary').innerHTML,/WT_P3/);
}
assert.match(node('cohort').textContent,/633 designed SNVs/);
assert.match(node('cohort').textContent,/GRCh38 chr19:11089285–11089495/);
assert.ok(drawnText.has('11089285'));assert.ok(drawnText.has('11089495'));
assert.ok(node('regionTable').innerHTML.includes('Genomic peak'));
assert.ok(node('wtCategories').innerHTML.includes('Nucleosome-sized'));
assert.equal(node('variant').children.length,634);
assert.match(node('variantCounter').textContent,/Variant 1 of 633/);
assert.equal(node('functionalSource').children.length,2);
const first=String(node('variant').value);
node('nextVariant').events.click();assert.match(node('variantCounter').textContent,/Variant 2 of 633/);
node('previousVariant').events.click();assert.equal(String(node('variant').value),first);
node('previousVariant').events.click();assert.match(node('variantCounter').textContent,/Variant 633 of 633/);
node('nextVariant').events.click();assert.equal(String(node('variant').value),first);
const visited=new Set();
for(let i=0;i<633;i++){visited.add(String(node('variant').value));node('nextVariant').events.click();}
assert.equal(visited.size,633);assert.equal(String(node('variant').value),first);
const target=data.variants.findIndex(v=>v.source_variant==='125A>T');
node('variant').value=String(target);node('variant').events.change();
assert.match(node('functionalSelected').textContent,/125A>T.*-0\.070/);
assert.match(node('functionalSelected').textContent,/chr19:11089389:A>T/);
if(process.argv[3]){
 assert.match(node('biologyAuditTitle').textContent,/11089389:A>T.*WT_P3/);
 assert.match(node('biologyAudit').innerHTML,/-3\.560/);
 node('biologyMotifQuery').value='SREBF';node('biologyMotifQuery').events.input();
 assert.match(node('biologyMotifTable').innerHTML,/SREBF1/);
 assert.ok(!node('biologyMotifTable').innerHTML.includes('SP1.H12CORE'));
 node('biologyDatabase').value='HOCOMOCOv12';node('biologyDatabase').events.change();
 assert.match(node('biologyMotifCount').textContent,/^0 candidate matches/);
 node('biologyMotifQuery').value='';node('biologyMotifQuery').events.input();
 node('biologyDatabase').value='all';node('biologyDatabase').events.change();
 node('biologyAuditEndpoint').value='whole_msp';node('biologyAuditEndpoint').events.change();
 assert.match(node('biologyAuditTitle').textContent,/whole_msp/);
 assert.match(node('biologyAudit').innerHTML,/3\.036/);
 node('biologyAuditEndpoint').value='local_tf';node('biologyAuditEndpoint').events.change();
 node('biologyRegion').value='adjacent_143_150';node('biologyRegion').events.change();
 assert.match(node('biologyLocalVariants').innerHTML,/11089411:T&gt;G/);
 assert.equal((node('biologyLocalVariants').innerHTML.match(/<tr>/g)||[]).length,25);
 node('biologyRegion').value='follow';node('biologyRegion').events.change();
}
assert.equal(vm.runInContext('lollipopPoints.find(p=>p.variant.variant_id===selected().variant_id).x',sandbox),55+(125-21)/210*(1100-55-15));
node('functionalSource').value='published_LDLR_2';node('functionalSource').events.change();
assert.match(node('functionalSelected').textContent,/125A>T.*0\.000/);
let prevented=false;
node('comparisonPanel').events.keydown({key:'ArrowRight',target:{tagName:'SECTION'},preventDefault(){prevented=true;}});
assert.ok(prevented);assert.notEqual(String(node('variant').value),String(target));
const afterKey=String(node('variant').value);
node('comparisonPanel').events.keydown({key:'ArrowRight',target:{tagName:'SELECT'},preventDefault(){throw Error('Native select key intercepted');}});
assert.equal(String(node('variant').value),afterKey);
const point=vm.runInContext('lollipopPoints[0]',sandbox);
node('functionalLollipop').events.click({clientX:point.x,clientY:point.y});
assert.equal(data.variants[Number(node('variant').value)].variant_id,point.variant.variant_id);
// Missing functional measurements must not be plotted as zero effects.
vm.runInContext('delete selected().functional[$("functionalSource").value];drawLollipop()',sandbox);
assert.match(node('functionalSelected').textContent,/functional effect not measured/);
assert.equal(vm.runInContext('lollipopPoints.some(p=>p.variant.variant_id===selected().variant_id)',sandbox),false);
const missing=data.variants.findIndex(v=>v.n_barcodes===0);assert.ok(missing>=0);
node('variant').value=String(missing);node('variant').events.change();
assert.match(node('selectedStats').textContent,/variant Not measured/);
assert.match(node('selectedInfo').textContent,/LOW OR MISSING COVERAGE/);
if(process.argv[3])assert.match(node('biologyAudit').innerHTML,/Not measured/);
node('region').value='t7';node('region').events.change();
node('query').value='160A>G';node('query').events.input();
assert.match(node('variantTable').innerHTML,/Not measured/);
node('support').value='floor';node('support').events.change();
assert.equal(node('variantTable').innerHTML,'');
node('query').value='';node('query').events.input();
assert.match(node('tableCount').textContent,/630 variants/);
node('query').value='chr19:11089389:A>T';node('query').events.input();
assert.match(node('tableCount').textContent,/1 variants/);
assert.match(node('variantTable').innerHTML,/125A&gt;T/);
assert.equal(vm.runInContext('color(null,1)===color(0,1)',sandbox),false);
console.log('Viewer checks passed: genomic axes/alleles/search, all 633 navigation steps, wrapping, keyboard focus, synchronized lollipop/source selection, point clicks, filters, and missingness.');

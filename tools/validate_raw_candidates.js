// DOM/canvas contract test with the real embedded dataset; not a browser render.
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const html=fs.readFileSync(process.argv[2],'utf8');
const payload=html.match(/<script id="raw-data" type="application\/json">(.*?)<\/script>/s)[1];
const app=html.match(/<script>(.*?)<\/script>/s)[1];
const context=new Proxy({},{get:(o,k)=>o[k]||(()=>{}),set:(o,k,v)=>(o[k]=v,true)});
const elements=new Map();
function element(id){if(!elements.has(id))elements.set(id,{value:'',textContent:'',innerHTML:'',width:1100,height:240,events:{},appendChild(){},addEventListener(k,f){this.events[k]=f;},getContext:()=>context,setAttribute(){}});return elements.get(id);}
element('raw-data').textContent=payload;
const document={getElementById:element,createElement:()=>({})};
const sandbox=vm.createContext({document,console});vm.runInContext(app,sandbox);
const data=JSON.parse(payload);let checks=0;
for(const v of data.variants.filter(v=>v.source_variant!=='WT'))for(const s of data.scenarios)for(const r of data.regions){
 element('variant').value=v.source_variant;element('scenario').value=s.id;element('region').value=r.id;
 element('variant').events.change();assert(element('selected').innerHTML.includes('<table>'));assert(!element('selected').innerHTML.includes('NaN'));checks++;
}
element('next').events.click();element('previous').events.click();
for(const track of ['tf_small','msp','nuc_85_250','nuc_gt250']){element('track').value=track;element('track').events.change();}
console.log(JSON.stringify({complete:true,variant_scenario_region_combinations:checks,canvas:'mocked; no browser visual QA'}));

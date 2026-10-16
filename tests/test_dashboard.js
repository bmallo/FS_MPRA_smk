"use strict";
const assert = require("node:assert/strict");
const Pilot = require("../dashboard/app.js");

function variant(id, cls, count, values) {
  const tracks = Object.fromEntries(Pilot.tracks.map(name => [name, {
    occupancy:values, delta_vs_WT:values, eligible_reads:[count,count], called_reads:[0,0]
  }]));
  return {variant_id:id, class:cls, n_reads:count, alt:cls==="WT"?"":"T", oligo_pos1:cls==="WT"?null:21,
    tracks, summary:{msp_common_site_delta:count?0:null}};
}
const data = {schema_version:"fsmpra.pilot.v1", sample:{id:"test"}, region:{sequence:"AC"},
  variants:[variant("wt","WT",5,[0,0]), variant("snv","SNV",0,[null,null])]};
const lookup = Pilot.index(data);
assert.equal(lookup.byCell.get("21:T").variant_id,"snv");
assert.equal(Pilot.metric(lookup.wt,"msp_common_site_delta"),0);
assert.equal(Pilot.metric(lookup.byId.get("snv"),"msp_common_site_delta"),null);
assert.notEqual(Pilot.color(null,1),Pilot.color(0,1));
assert.deepEqual(Pilot.normalized([0,0]),[null,null]);
assert.deepEqual(Pilot.normalized([1,3]),[.25,.75]);
assert.equal(Pilot.path([0,null,1],i=>i,x=>x),"M0.00,0.00 M2.00,1.00 ");
assert.equal(Pilot.escape('<script>"&'),"&lt;script&gt;&quot;&amp;");
assert.equal(Pilot.pvalue(1e-12),"1.00e-12");
assert.equal(Pilot.pvalue(0),"0 (reported)");
assert.throws(()=>Pilot.index({...data,variants:[...data.variants,data.variants[0]]}),/identity/);
assert.throws(()=>Pilot.index({...data,variants:[data.variants[1]]}),/one WT/);
assert.throws(()=>Pilot.index({...data,schema_version:"unknown"}),/Expected/);
const invalid=JSON.parse(JSON.stringify(data));
invalid.variants[1].tracks.msp.occupancy=[null];
assert.throws(()=>Pilot.index(invalid),/array/);
// Validate actual output bundles too when paths are supplied on the command line.
for (const file of process.argv.slice(2)) {
  const actual = JSON.parse(require("node:fs").readFileSync(file,"utf8"));
  for (const data of actual.datasets || [actual]) {
    const indexed=Pilot.index(data);
    assert.equal(indexed.byId.size,data.variants.length);
    console.log(`Validated ${data.variants.length} alleles: ${data.sample.id}`);
  }
}
console.log("Dashboard checks passed: missingness, zero, gaps, joins, escaping, and schema rejection.");

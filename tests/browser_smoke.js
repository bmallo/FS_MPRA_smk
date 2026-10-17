/* Append to an exported pilot HTML for a real-browser interaction smoke test.
 * Uses only the local DOM; no automation server or network is required.
 */
(() => {
  const $ = id => document.getElementById(id);
  const check = (condition, message) => { if (!condition) throw Error(message); };
  const change = (id,value,event="change") => { $(id).value=value; $(id).dispatchEvent(new Event(event)); };
  const banner = document.createElement("div");
  banner.style="position:fixed;top:0;left:0;right:0;z-index:1000;padding:12px;background:#135e35;color:white;font:16px sans-serif";
  document.body.prepend(banner);
  try {
    const bundle=JSON.parse($("pilot-data").textContent);
    const data=bundle.datasets?bundle.datasets[bundle.datasets.length-1]:bundle;
    const snv=data.variants.find(v=>v.class==="SNV" && v.n_reads>0);
    const missing=data.variants.find(v=>v.class==="SNV" && !v.n_reads);
    check(!$("app").hidden,"App did not initialize");
    change("search",snv.source_variant,"input");
    change("variant",snv.variant_id);
    check($("detail-title").textContent.startsWith(snv.source_variant),"Selection did not update detail");
    change("mode","delta");
    check($("profiles").querySelectorAll("svg").length===3,"Missing profile panels");
    check(!/NaN|undefined/.test($("profiles").innerHTML),"Invalid profile coordinates");
    change("minimum","1000000","input");
    check($("variant-info").textContent.includes("Below the display"),"Coverage warning missing");
    $("wildtype").click();
    check($("search").value==="" && $("variant").selectedOptions[0].text.startsWith("WT"),"WT navigation left a stale search");
    $("next").click();
    check(!$("detail-title").textContent.startsWith("WT"),"Next allele failed");
    $("previous").click();
    check($("detail-title").textContent.startsWith("WT"),"Previous allele failed");
    if(missing) {
      change("variant",missing.variant_id);
      check($("profiles").textContent.includes("No measured track"),"Missing data shown as zero");
    }
    change("metric","n_reads");
    check($("scale").textContent.includes("molecules"),"Overview metric did not update");
    change("minimum","10","input"); change("mode","occupancy"); change("variant",snv.variant_id);
    if(bundle.datasets) {
      check(!$("comparison-panel").hidden,"Comparison panel hidden");
      check($("comparison-profiles").querySelectorAll("svg").length===2,"WT comparison profiles missing");
      check($("agreement-summary").textContent.includes("Pearson"),"Agreement summary missing");
      if(bundle.comparison.recovery && Object.values(bundle.comparison.recovery).some(Boolean)) {
        check($("recovery-summary").textContent.includes("plasmid_partial"),"Recovery classes missing");
        check($("comparison-summary").textContent.includes("100%"),"Full-ledger fraction missing");
      }
      check($("across-runs").querySelectorAll("button").length===bundle.datasets.length,"Allele comparison incomplete");
      $("across-runs").querySelector("button").click();
      check($("sample").value===bundle.datasets[0].sample.id,"Run switch failed");
      check($("detail-title").textContent.startsWith(snv.source_variant),"Run switch lost selected allele");
      check($("status").textContent.includes("10 min"),"Reaction condition missing");
      change("sample",data.sample.id);
    }
    banner.textContent="PASS · initialization, allele search/selection, WT navigation, next/previous, delta profiles, coverage threshold, missing data, overview metric";
  } catch(error) {
    banner.style.background="#a32424";
    banner.textContent="FAIL · "+error.message;
  }
})();

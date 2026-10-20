"""Reconcile the landscape with independently generated prior summaries."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from fsmpra.io import atomic_output, open_text, source_info, write_json
from fsmpra.processing import check_sources


def run(root,output):
    root=Path(root)
    old_profiles=Path('results/rep3_statistics_20261017/inference/profiles.tsv.gz')
    old_heterogeneity=Path('results/rep3_heterogeneity_20261017/summary/variant_heterogeneity.tsv.gz')
    sources=[source_info(p,True) for p in [root/'landscape/landscape.npz',root/'atlas/atlas.json',old_profiles,old_heterogeneity,Path(__file__)]]
    for stage in ['landscape','atlas','patterns']:
        a=json.loads((root/stage/'audit.json').read_text());assert a['complete'];check_sources(a['sources'])
    with np.load(root/'landscape/landscape.npz',allow_pickle=False) as a:data={k:a[k] for k in a.files}
    idx={v:i for i,v in enumerate(data['variant_ids'])}
    error=0.;rows=missing=0
    with open_text(old_profiles) as h:
        for r in csv.DictReader(h,delimiter='\t'):
            i=idx[r['variant_id']];t={'tf_small':0,'msp':7}[r['track']];p=int(r['relative_pos'])
            actual=data['coverage'][i,t,p]
            if not r['mean_equal_barcode']:
                assert np.isnan(actual);missing+=1
            else:error=max(error,abs(float(r['mean_equal_barcode'])-actual))
            rows+=1
    assert rows==634*2*211 and missing==422 and error<1e-10,(rows,missing,error)
    nuc_error=0.;nuc_rows=0
    with open_text(old_heterogeneity) as h:
        for r in csv.DictReader(h,delimiter='\t'):
            if r['cohort']=='footprint_analysis_eligible_with_barcode' and r['feature']=='element_nucleosome_fraction' and r['barcode_mean']:
                i=idx[r['variant_id']]
                nuc_error=max(nuc_error,abs(float(r['barcode_mean'])-data['coverage'][i,3:7].sum(axis=0).mean()))
                nuc_rows+=1
    assert nuc_rows==633 and nuc_error<1e-10,(nuc_rows,nuc_error)
    bundle=json.loads((root/'atlas/atlas.json').read_text())
    assert len(bundle['variants'])==634
    missing_variants=[]
    for v in bundle['variants']:
        i=idx[v['variant_id']]
        if v['n_barcodes']==0:
            missing_variants.append(v['source_variant'])
            assert all(x is None for track in v['coverage'] for x in track)
            assert all(x is None for x in v['region_effects'])
        else:
            assert abs(sum(v['states'])-1)<1e-6
            assert np.max(np.abs(np.array(v['coverage'])-data['coverage'][i]))<5.01e-8
        assert len(v['coverage'])==8 and all(len(track)==211 for track in v['coverage'])
    assert missing_variants==['160A>G']
    html=(root/'atlas/index.html').read_text()
    assert 'DATA_PLACEHOLDER' not in html and 'APP_PLACEHOLDER' not in html
    result=dict(status='passed',prior_profile_rows=rows,missing_profile_rows=missing,
                maximum_profile_error=error,nucleosome_union_rows=nuc_rows,
                maximum_nucleosome_union_error=nuc_error,missing_variants=missing_variants,
                source_audits_verified=['landscape','atlas','patterns'],sources=sources)
    check_sources(sources)
    with atomic_output(output) as out:write_json(out/'audit.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='sources'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();run(a.root,a.output)

"""Check raw audit conservation, independent aggregate agreement and packaging."""
import argparse
import re
from collections import defaultdict
from pathlib import Path
import numpy as np
from fsmpra.io import read_json,write_json
from fsmpra.processing import check_sources
from tools.render_raw_candidates import read_table


def run(viewer,output):
    viewer=Path(viewer);d=read_json(viewer/'viewer_data.json')
    check_sources(read_json(viewer/'audit.json')['sources'])
    check_sources(read_json(viewer/'summary_audit.json')['sources'])
    comparisons=read_table(viewer/'comparisons.tsv');profiles=read_table(viewer/'profiles.tsv.gz');retention=read_table(viewer/'retention.tsv')
    assert len(comparisons)==27*18*4*8 and len(profiles)==28*18*211 and len(retention)==28*18
    seen=set()
    for r in comparisons:
        key=tuple(r[k] for k in ['source_variant','scenario','region','metric']);assert key not in seen;seen.add(key)
        if r['delta_pp'] is not None:assert np.isclose(r['delta_pp'],100*(r['variant_value']-r['wt_value']))
        else:assert r['variant_value'] is None or r['wt_value'] is None
        for group in ['wt','variant']:
            assert r[group+'_nonzero_units']<=r[group+'_molecules']
            if not r[group+'_molecules']:assert r[group+'_value'] is None
    old=read_json('results/rep3_landscape_20261018/atlas_genomic_20261019/atlas.json')
    original={v['variant_id']:v for v in old['variants']};error=0.
    for r in profiles:
        p=r['genomic_pos1']-11089285;assert 0<=p<211
        if r['base'] not in 'AT':assert r['retained_m6a_fraction'] is None and r['callable_units']==0
        if r['scenario']=='primary_barcodes' and r['source_variant']!='WT':
            for field,track in [('tf_small',0),('msp',7)]:
                error=max(error,abs(r[field]-original[r['variant_id']]['coverage'][track][p]))
    assert error<1.01e-7,(error,'variant profiles differ from original full-data aggregates')
    for r in retention:
        assert r['footprint_molecules']<=r['molecules']
        if r['scenario']=='primary_barcodes' and r['source_variant']!='WT':
            v=original[r['variant_id']];assert r['barcodes']==v['n_barcodes'] and r['molecules']==v['n_molecules']
    groups=defaultdict(list)
    for r in read_table(viewer/'joint_molecule_states.tsv'):groups[r['source_variant'],r['region'],r['group']].append(r['fraction'])
    assert len(groups)==27*4*2
    for values in groups.values():
        assert len(values)==16
        assert all(v is None for v in values) or np.isclose(sum(values),1)
    html=(viewer/'index.html').read_text()
    for target in re.findall(r'href="([^"]+)"',html.split('</main>')[0]):assert (viewer/target).is_file(),target
    import json
    assert json.loads(re.search(r'<script id="raw-data" type="application/json">(.*?)</script>',html,re.S).group(1))==d
    write_json(output,dict(complete=True,comparisons=len(comparisons),profiles=len(profiles),retention_rows=len(retention),joint_groups=len(groups),max_variant_profile_rounding_error=error))
    print('Validated counts, missingness, deltas, joint mass, original variant aggregates and packaged exports.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--viewer',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();run(a.viewer,a.output)

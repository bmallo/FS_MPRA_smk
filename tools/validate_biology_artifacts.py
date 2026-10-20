"""Validate biological annotations against the preserved atlas and exports."""
import argparse
import csv
import re
from pathlib import Path
import numpy as np
from fsmpra.io import atomic_output,open_text,read_json,source_info,write_json
from fsmpra.processing import check_sources


def run(atlas,previous,output):
    atlas,previous=map(Path,(atlas,previous))
    audit=read_json(atlas/'audit.json');check_sources(audit['sources'])
    check_sources(read_json(atlas/'biology_audit.json')['sources'])
    for p in previous.iterdir():
        if p.name=='index.html':continue
        target=atlas/('base_atlas_audit.json' if p.name=='audit.json' else p.name)
        if source_info(p,True)['sha256']!=source_info(target,True)['sha256']:
            raise ValueError('Original artifact changed: '+p.name)
    d=read_json(atlas/'atlas.json');b=read_json(atlas/'biology.json')
    variants={v['variant_id']:v for v in d['variants']};wt=d['variants'][d['wt_index']]
    hits={h['hit_id']:h for h in b['motifs']}
    if len(hits)!=len(b['motifs']):raise ValueError('Duplicate motif hits')
    for h in hits.values():
        if h['sequence']!=d['sequence'][h['start0']:h['end0']] or h['genomic_start1']!=d['genomic_start1']+h['start0'] or h['genomic_end1']!=d['genomic_start1']+h['end0']-1:
            raise ValueError('Motif sequence or genomic position differs')
        if h['sequence_tail_p']>b['scan']['threshold']:raise ValueError('Unexpected motif cutoff')
    effects=0;seen=set()
    with open_text(atlas/'motif_variant_effects.tsv.gz') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            h=hits[row['hit_id']];v=variants[row['variant_id']];key=(h['hit_id'],v['variant_id'])
            if key in seen:raise ValueError('Duplicate motif/variant effect')
            seen.add(key)
            weights=np.array(h['weights_genomic']);sequence=list(h['sequence'])
            sequence[v['relative_pos0']-h['start0']]=v['alt']
            score=sum(weights[i,'ACGT'.index(base)] for i,base in enumerate(sequence))
            if not np.isclose(score,float(row['variant_score_bits']),atol=1e-10):raise ValueError('Variant motif score differs')
            for source in ('published_LDLR','published_LDLR_2'):
                if float(row[source+'_log2'])!=v['functional'][source]:raise ValueError('Functional effect changed')
            effects+=1
    if effects!=sum(3*(h['end0']-h['start0']) for h in hits.values()):raise ValueError('Incomplete motif substitution grid')
    seen=set();max_error=0.
    for r in b['barcode_audits']:
        key=r['variant_id'],r['region_id']
        if key in seen:raise ValueError('Duplicate barcode audit')
        seen.add(key);v=variants[r['variant_id']]
        if r['n_barcodes']!=v['n_barcodes'] or r['split_0_barcodes']+r['split_1_barcodes']!=r['n_barcodes']:
            raise ValueError('Barcode split conservation failed')
        if not 0<=r['nonzero_barcodes']<=r['n_barcodes'] or not 0<=r['depth_ge5_barcodes']<=r['depth_ge3_barcodes']<=r['n_barcodes']:
            raise ValueError('Invalid support counts')
        if not r['n_barcodes']:
            if r['delta_pp'] is not None:raise ValueError('Missing measurements replaced with zero')
            continue
        track=7 if r['region_id']=='whole_msp' else 0;a,z=r['start0'],r['end0']
        expected=100*np.mean(np.array(v['coverage'][track][a:z])-np.array(wt['coverage'][track][a:z]))
        max_error=max(max_error,abs(expected-r['delta_pp']))
    if len(seen)!=633*7 or max_error>1.01e-5:raise ValueError('Barcode effects differ from original profiles')
    html=(atlas/'index.html').read_text()
    for identifier,data in [('atlas-data',d),('biology-data',b)]:
        import json
        match=re.search(r'<script id="'+identifier+r'" type="application/json">(.*?)</script>',html,re.S)
        if not match or json.loads(match.group(1))!=data:raise ValueError('Embedded data differs')
    for target in re.findall(r'href="([^"]+)"',html.split('</main>')[0]):
        if not target.startswith(('#','http')) and not (atlas/target).is_file():raise ValueError('Broken export: '+target)
    with atomic_output(output) as out:
        write_json(out/'audit.json',dict(complete=True,original_artifacts_unchanged=True,motif_hits=len(hits),
                   motif_variant_pairs=effects,barcode_audit_rows=len(seen),max_rounding_difference_pp=max_error,
                   sources=[source_info(p,True) for p in [atlas/'audit.json',atlas/'atlas.json',atlas/'biology.json',atlas/'index.html',Path(__file__)]]))
    print(f'Validated {len(hits)} motif hits, {effects} motif/variant pairs, {len(seen)} barcode audits; original artifacts unchanged.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--atlas',required=True)
    p.add_argument('--previous',default='results/rep3_landscape_20261018/atlas_genomic_20261019')
    p.add_argument('--output',required=True);a=p.parse_args();run(a.atlas,a.previous,a.output)

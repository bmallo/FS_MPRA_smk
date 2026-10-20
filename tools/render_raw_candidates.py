"""Build a self-contained, offline raw-call audit and scientific figures."""
import argparse
import csv
import json
import shutil
from collections import defaultdict
from pathlib import Path
import numpy as np
from fsmpra.io import atomic_output, open_text, read_json, source_info, write_json, write_table
from fsmpra.processing import check_sources


def read_table(path):
    with open_text(path) as h:
        rows=list(csv.DictReader(h,delimiter='\t'))
    for row in rows:
        for k,v in row.items():
            if v=='':row[k]=None
            else:
                try:row[k]=float(v) if any(c in v for c in '.eE') else int(v)
                except ValueError:pass
    return rows


def figures(out,data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,1,figsize=(12,9),sharex=True)
    x=np.arange(93,113)+11089285
    for name,color in [('WT','#253d4b'),('125A>C','#467aaf'),('125A>G','#a6568e'),('125A>T','#cc7b32')]:
        p=data['profiles']['primary_barcodes|'+name]
        for ax,base in zip(axes[:2],'AT'):
            positions=[i for i in range(93,113) if data['sequence'][i]==base and abs(i-104)>2 and p['retained_m6a_fraction'][i] is not None]
            ax.plot(np.array(positions)+11089285,[100*p['retained_m6a_fraction'][i] for i in positions],'o-',color=color,label=name)
            ax.set(ylabel=f'Ref {base}: called fraction (%)',ylim=(0,100))
        axes[2].plot(x,np.array(p['tf_small'][93:113])*100,color=color,label=name)
    for ax in axes:
        ax.axvspan(11089387,11089391,color='grey',alpha=.12);ax.grid(alpha=.2)
        ax.ticklabel_format(axis='x',style='plain',useOffset=False)
    axes[0].legend(ncol=4);axes[2].set(ylabel='Small-footprint coverage (%)',xlabel='GRCh38 chr19 (1-based)',ylim=(0,None))
    fig.suptitle('Position 125: retained calls versus existing small footprints\nGray = mutation ±2 bp; equal barcode weight, sampled WT; lines guide the eye')
    fig.tight_layout()
    for ext in ('png','svg'):fig.savefig(out/f'position125_raw.{ext}',dpi=150)
    plt.close(fig)
    rows=data['comparisons'];keys=data['comparison_fields'];rows=[dict(zip(keys,r)) for r in rows]
    variants=[v for v in data['variants'] if v['source_variant']!='WT']
    scenarios=data['scenarios'];lookup={(r['source_variant'],r['scenario'],r['region'],r['metric']):r['delta_pp'] for r in rows}
    fig,axes=plt.subplots(1,3,figsize=(18,10),sharey=True)
    for ax,metric in zip(axes,['retained_m6a_shared','tf_small','msp']):
        mat=np.array([[lookup[v['source_variant'],s['id'],'WT_P3' if v['oligo_pos1']==125 else 'promoter',metric] for s in scenarios] for v in variants],dtype=float)
        lim=max(1,np.nanmax(np.abs(mat)));im=ax.imshow(np.ma.masked_invalid(mat),cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto')
        ax.set(xticks=range(len(scenarios)),xticklabels=[s['id'] for s in scenarios],yticks=range(len(variants)),yticklabels=[v['source_variant'] for v in variants],title=metric)
        ax.tick_params(axis='x',rotation=90,labelsize=8);fig.colorbar(im,ax=ax,label='Variant − WT (pp)',shrink=.6)
    fig.suptitle('Sensitivity comparisons: position 125 local WT_P3; positions 143–150 whole promoter\nEach panel has its own symmetric color scale. Empty cells are unavailable, not zero.')
    fig.tight_layout()
    for ext in ('png','svg'):fig.savefig(out/f'sensitivity.{ext}',dpi=150)
    plt.close(fig)


def run(extracted,summary,report,output):
    extracted,summary,report=map(Path,(extracted,summary,report))
    audit=read_json(summary/'audit.json');check_sources(audit['sources'])
    assets=Path('dashboard');paths=list(summary.iterdir())+[extracted/'variants.json',report,Path(__file__),assets/'raw_candidates.html',assets/'raw_candidates.js']
    sources=[source_info(p,True) for p in paths]
    rows=read_table(summary/'comparisons.tsv');profiles={}
    for r in read_table(summary/'profiles.tsv.gz'):
        key=r['scenario']+'|'+r['source_variant']
        if key not in profiles:profiles[key]={k:[] for k in ['retained_m6a_fraction','tf_small','msp','nuc_85_250','nuc_gt250']}
        for k in profiles[key]:profiles[key][k].append(r[k])
    grouped=defaultdict(list)
    for r in read_table(summary/'resampling.tsv'):
        if r['delta_pp'] is not None:grouped[tuple(r[k] for k in ['source_variant','region','metric','kind'])].append(r['delta_pp'])
    repeats=[]
    for key,values in grouped.items():
        a=np.array(values);r=dict(zip(['source_variant','region','metric','kind'],key))
        r.update(n=len(a),minimum=float(a.min()),median=float(np.median(a)),maximum=float(a.max()),negative=int((a<0).sum()),positive=int((a>0).sum()));repeats.append(r)
    variants=read_json(extracted/'variants.json');fields=list(rows[0])
    data=dict(variants=variants,sequence=variants[0]['sequence'],regions=audit['regions'],scenarios=audit['scenarios'],comparison_fields=fields,
              comparisons=[[r[k] for k in fields] for r in rows],profiles=profiles,retention=read_table(summary/'retention.tsv'),resampling=repeats)
    choices={'WT':'Reference architecture and reaction control', '125A>C':'A/T opportunity loss at local footprint candidate',
             '125A>G':'A/T opportunity loss at local footprint candidate','125A>T':'Preserves A/T opportunity; switches adenine strand',
             '143C>G':'Adjacent-region change preserving GC status','146C>G':'Adjacent-region change preserving GC status',
             '150C>G':'Adjacent-region change preserving GC status','144T>A':'Adjacent-region adenine-strand switch',
             '147T>A':'Adjacent-region adenine-strand switch','147T>G':'Functional discordance / A/T opportunity loss'}
    panel=[dict(source_variant=v['source_variant'],variant_id=v['variant_id'],genomic_pos1=v['genomic_pos1'],ref=v['ref'],alt=v['alt'],mutable_sequence=v['sequence'],rationale=choices[v['source_variant']],status='Proposed; not a new construct order') for v in variants if v['source_variant'] in choices]
    if len(panel)!=len(choices):raise ValueError('Validation construct absent from catalog')
    with atomic_output(output) as out:
        for p in summary.iterdir():shutil.copy2(p,out/('summary_audit.json' if p.name=='audit.json' else p.name))
        shutil.copy2(report,out/'report.md');write_table(out/'validation_panel.tsv',panel,list(panel[0]))
        write_json(out/'viewer_data.json',data)
        payload=json.dumps(data,separators=(',',':'),allow_nan=False).replace('&','\\u0026').replace('<','\\u003c').replace('>','\\u003e')
        html=(assets/'raw_candidates.html').read_text().replace('DATA_PLACEHOLDER',payload).replace('/* RAW_APP_PLACEHOLDER */',(assets/'raw_candidates.js').read_text())
        (out/'index.html').write_text(html);figures(out,data)
        check_sources(sources);write_json(out/'audit.json',dict(schema_version='fsmpra.raw_candidate_viewer.v1',complete=True,sources=sources,variants=27,scenarios=18))
    print(output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--extracted',default='results/rep3_raw_candidates_20261020/extracted');p.add_argument('--summary',default='results/rep3_raw_candidates_20261020/summary')
    p.add_argument('--report',default='docs/rep3_raw_candidates_20261020.md');p.add_argument('--output',required=True)
    a=p.parse_args();run(a.extracted,a.summary,a.report,a.output)

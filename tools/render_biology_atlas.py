"""Package biological annotations with the preserved genomic atlas, offline."""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np
from fsmpra.io import atomic_output,read_json,source_info,write_json
from fsmpra.processing import check_sources


def figures(out,atlas,biology):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d=read_json(atlas/'atlas.json');wt=d['variants'][d['wt_index']];origin=d['genomic_start1']
    fig,axes=plt.subplots(3,1,figsize=(14,9),sharex=True,gridspec_kw={'height_ratios':[3,1,2]})
    x=np.arange(211)+origin
    axes[0].plot(x,100*np.array(wt['coverage'][0]),color='#167b77',label='WT small footprints (10–60 bp)')
    for r in biology['regions']:
        if r['id'].startswith('WT_'):
            axes[0].axvspan(r['genomic_start1'],r['genomic_end1'],alpha=.08,color='#72529a')
            axes[0].text((r['genomic_start1']+r['genomic_end1'])/2,9.7,r['id'],ha='center')
    axes[0].set(ylabel='Covered fraction (%)',ylim=(0,11));axes[0].legend(loc='upper right')
    at=[origin+i for i,b in enumerate(d['sequence']) if b in 'AT']
    axes[1].vlines(at,0,1,color='#167b77');axes[1].set(yticks=[],ylabel='A/T sites')
    # Explicitly chosen context models, not the highest ranked identities.
    examples=[('JASPAR2024','MA0595.1',11089391),('JASPAR2024','MA0596.1',11089391),
              ('HOCOMOCOv12','SP1.H12CORE.2.P.B',11089397),('HOCOMOCOv12','WT1.H12CORE.0.P.B',11089407)]
    labels=[]
    for db,mid,start in examples:
        h=next(h for h in biology['motifs'] if h['database']==db and h['matrix_id']==mid and h['genomic_start1']==start)
        j=len(labels);labels.append(h['name']+' · '+mid)
        axes[2].plot([h['genomic_start1'],h['genomic_end1']],[j,j],lw=6,color='#72529a' if db=='JASPAR2024' else '#2376b3')
        axes[2].plot(h['genomic_end1'] if h['strand']=='+' else h['genomic_start1'],j,marker='>' if h['strand']=='+' else '<',color='#183445')
    axes[2].set(yticks=range(len(labels)),yticklabels=labels,ylim=(-.7,len(labels)-.3),
                xlabel='GRCh38 chr19 genomic position (1-based; forward axis)',xlim=(origin,origin+210))
    axes[2].ticklabel_format(axis='x',style='plain',useOffset=False)
    fig.suptitle('WT architecture, A/T opportunities and selected sequence candidates\nIllustrative SREBF / SP1 / WT1 matches; 114 total matches are available in the atlas')
    fig.tight_layout()
    for ext in ('png','svg'):fig.savefig(out/f'annotated_wt.{ext}',dpi=150)
    plt.close(fig)
    fields=['delta_pp','split_0_delta_pp','split_1_delta_pp','depth_ge5_delta_pp','distal_m6a_standardized_delta_pp','shared_at_radius2_delta_pp']
    names=['All barcodes','Half A','Half B','Depth ≥5','Distal m6A\nstandardized','Shared A/T\nexclude ±2 bp']
    fig,ax=plt.subplots(figsize=(11,4.8))
    for name,color in [('125A>C','#467aaf'),('125A>G','#a6568e'),('125A>T','#cc7b32')]:
        row=next(r for r in biology['barcode_audits'] if r['source_variant']==name and r['region_id']=='WT_P3')
        ax.plot(range(len(fields)),[row[k] for k in fields],'o-',label=row['genomic_variant'],color=color)
    ax.axhline(0,color='grey',lw=.8);ax.set(xticks=range(len(fields)),xticklabels=names,ylabel='Variant − WT coverage (percentage points)',
             title='chr19:11089389 substitutions · WT_P3 local footprint loss\nDescriptive sensitivity comparisons; no confidence intervals or new significance tests')
    ax.legend();ax.grid(axis='y',alpha=.2);fig.tight_layout()
    for ext in ('png','svg'):fig.savefig(out/f'position_125_sensitivity.{ext}',dpi=150)
    plt.close(fig)


def run(atlas,biology,report,output):
    atlas,biology,report=map(Path,(atlas,biology,report));assets=Path('dashboard')
    check_sources(read_json(atlas/'audit.json')['sources'])
    check_sources(read_json(biology/'audit.json')['sources'])
    paths=list(atlas.iterdir())+list(biology.iterdir())+[report,Path(__file__),assets/'biology.js',assets/'biology_sequence.html',assets/'biology_panel.html']
    sources=[source_info(p,True) for p in paths]
    data=read_json(biology/'biology.json');original=read_json(atlas/'atlas.json')
    if data['coordinates']!=original['coordinates']:raise ValueError('Annotation coordinates differ')
    with atomic_output(output) as out:
        for p in atlas.iterdir():shutil.copy2(p,out/('base_atlas_audit.json' if p.name=='audit.json' else p.name))
        for p in biology.iterdir():shutil.copy2(p,out/('biology_audit.json' if p.name=='audit.json' else p.name))
        shutil.copy2(report,out/'candidate_report.md')
        html=(atlas/'index.html').read_text()
        marker='<h3>Functional effect across the element</h3>'
        if html.count(marker)!=1:raise ValueError('Cannot locate sequence insertion point')
        html=html.replace(marker,(assets/'biology_sequence.html').read_text()+marker)
        marker='<section class="card"><h2>Interpretation and exports</h2>'
        if html.count(marker)!=1:raise ValueError('Cannot locate annotation panel insertion point')
        html=html.replace(marker,(assets/'biology_panel.html').read_text()+marker)
        payload=json.dumps(data,separators=(',',':'),allow_nan=False).replace('&','\\u0026').replace('<','\\u003c').replace('>','\\u003e')
        html=html.replace('</body>',f'<script id="biology-data" type="application/json">{payload}</script><script>'+(assets/'biology.js').read_text()+'</script></body>')
        (out/'index.html').write_text(html)
        figures(out,atlas,data)
        if (out/'atlas.json').read_bytes()!=(atlas/'atlas.json').read_bytes():raise ValueError('Original measurements changed')
        write_json(out/'audit.json',dict(schema_version='fsmpra.biology_atlas.v1',complete=True,original_atlas_unchanged=True,sources=sources))
        check_sources(sources)
    print(output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--atlas',default='results/rep3_landscape_20261018/atlas_genomic_20261019')
    p.add_argument('--biology',default='results/rep3_landscape_20261018/biology_20261020_complete')
    p.add_argument('--report',default='docs/rep3_candidate_biology_20261020.md')
    p.add_argument('--output',required=True)
    a=p.parse_args();run(a.atlas,a.biology,a.report,a.output)

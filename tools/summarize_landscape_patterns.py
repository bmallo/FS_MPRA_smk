"""Export explicitly post hoc candidate patterns from the descriptive atlas."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from fsmpra.io import atomic_output, source_info, write_json, write_table
from fsmpra.processing import check_sources


def run(atlas,output):
    atlas=Path(atlas)
    sources=[source_info(p,True) for p in [atlas/'atlas.json',atlas/'variant_summary.tsv',atlas/'audit.json',Path(__file__)]]
    bundle=json.loads((atlas/'atlas.json').read_text())
    check_sources(json.loads((atlas/'audit.json').read_text())['sources'])
    with open(atlas/'variant_summary.tsv') as h:rows=list(csv.DictReader(h,delimiter='\t'))
    by={r['source_variant']:r for r in rows}
    selected=[r for r in rows if r['oligo_pos1'] and 143<=int(r['oligo_pos1'])<=150]
    if len(selected)!=24:raise ValueError('Expected 24 variants in highlighted region')
    metrics=['msp_delta_pp','nuc_85_250_delta_pp','tf_small_delta_pp','published_LDLR_log2','published_LDLR_2_log2']
    summaries={}
    for key in metrics:
        values=np.array([float(r[key]) for r in selected])
        summaries[key]=dict(n=24,negative=int((values<0).sum()),positive=int((values>0).sum()),
                            median=float(np.median(values)),minimum=float(values.min()),maximum=float(values.max()))
    coordinates=bundle['coordinates']
    origin=coordinates['genomic_start1']
    xlabel=f"{coordinates['assembly']} {coordinates['chrom']} genomic position (1-based; forward strand)"
    wt=bundle['variants'][bundle['wt_index']]
    variants={v['source_variant']:v for v in bundle['variants']}
    region=next(r for r in bundle['regions'] if r['id']=='WT_P3')
    a,b=region['start0'],region['end0']
    local=[]
    for name in ['125A>C','125A>G','125A>T']:
        v=variants[name]
        local.append(dict(source_variant=name,genomic_variant=v['genomic_variant'],genomic_pos1=v['genomic_pos1'],assembly=coordinates['assembly'],chrom=coordinates['chrom'],strand=coordinates['strand'],region_genomic_start1=region['genomic_start1'],region_genomic_end1=region['genomic_end1'],n_barcodes=v['n_barcodes'],n_molecules=v['n_molecules'],at_change=v['at_change'],
                          wt_coverage=float(np.mean(wt['coverage'][0][a:b])),variant_coverage=float(np.mean(v['coverage'][0][a:b])),
                          published_LDLR_log2=v['functional']['published_LDLR'],published_LDLR_2_log2=v['functional']['published_LDLR_2']))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    check_sources(sources)
    with atomic_output(output) as out:
        write_table(out/'region_143_150.tsv',selected,list(selected[0]))
        write_table(out/'position_125.tsv',local,list(local[0]))
        write_json(out/'patterns.json',dict(selection='Post hoc descriptive examples selected after viewing atlas; not independent validation.',
                    coordinates=coordinates,selected_region=dict(oligo_start1=143,oligo_end1=150,genomic_start1=origin+122,genomic_end1=origin+129),region_143_150=summaries,position_125=local))
        fig,ax=plt.subplots(figsize=(11,4.5))
        x=np.arange(211)+origin
        ax.plot(x,100*np.array(wt['coverage'][0]),label='WT',color='#273f4b',lw=2)
        for name,color in [('125A>C','#467aaf'),('125A>G','#a6568e'),('125A>T','#cc7b32')]:
            ax.plot(x,100*np.array(variants[name]['coverage'][0]),label=variants[name]['genomic_variant'],color=color,lw=1.5)
        ax.axvspan(a+origin,b+origin-1,color='grey',alpha=.08);ax.axvline(origin+104,color='grey',ls='--',lw=1)
        ax.set(xlabel=xlabel,ylabel='Small-footprint covered fraction (%)',xlim=(origin+79,origin+124),
               title=f'All three substitutions at {origin+104} (oligo 125) reduce the observed local signal\nA→T retains A/T count; strand/context changes remain unresolved')
        ax.ticklabel_format(axis='x',style='plain',useOffset=False)
        ax.legend();ax.grid(alpha=.15);fig.tight_layout()
        for ext in ['png','svg']:fig.savefig(out/f'position_125.{ext}',dpi=160)
        plt.close(fig)
        fields=[('msp_delta_pp','Accessibility change (pp)'),('nuc_85_250_delta_pp','Nucleosome-sized coverage change (pp)'),
                ('published_LDLR_log2','Published LDLR activity (log₂)'),('published_LDLR_2_log2','Published LDLR_2 activity (log₂)')]
        fig,axes=plt.subplots(4,1,figsize=(9,10),sharex=True)
        for ax,(field,label) in zip(axes,fields):
            mat=np.full((4,8),np.nan)
            for r in selected:mat['ACGT'.index(r['alt']),int(r['oligo_pos1'])-143]=float(r[field])
            limit=float(np.nanmax(np.abs(mat)));im=ax.imshow(mat,aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
            for ai in range(4):
                for pi in range(8):
                    if np.isfinite(mat[ai,pi]):ax.text(pi,ai,f'{mat[ai,pi]:+.1f}',ha='center',va='center',fontsize=9,color='white' if abs(mat[ai,pi])>.65*limit else 'black')
            ax.set(yticks=range(4),yticklabels=list('ACGT'),ylabel='Alternate base',title=label)
            fig.colorbar(im,ax=ax)
        axes[-1].set(xticks=range(8),xticklabels=range(origin+122,origin+130),xlabel=xlabel+' (WT sequence CTCCTCCC)')
        axes[-1].tick_params(axis='x',labelrotation=35)
        fig.suptitle(f'Candidate region {origin+122}–{origin+129} (oligo 143–150) · all 24 substitutions shown\nSelected after inspection; published sources are not independent Fiber-seq validation')
        fig.tight_layout()
        for ext in ['png','svg']:fig.savefig(out/f'region_143_150.{ext}',dpi=160)
        plt.close(fig)
        write_json(out/'audit.json',dict(schema_version='fsmpra.landscape_patterns.v1',complete=True,post_hoc=True,sources=sources))
        check_sources(sources)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--atlas',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();run(a.atlas,a.output)

"""Package descriptive landscape tables, scientific figures and a local viewer."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks
from scipy.stats import pearsonr, spearmanr

from .io import atomic_output, open_text, read_json, source_info, write_json, write_table
from .processing import check_sources
from .landscape import TRACKS, START, LENGTH, number
from .genomic_coordinates import interval_fields, mapping_for

DEFAULT_COORDINATE_AUDIT = 'results/rep3_landscape_20261018/genomic_coordinates_20261019/audit.json'

LABELS = ['Small footprints (10–60 bp)', 'Other TF-class calls', 'Nuc-class <85 bp',
          'Nucleosome-sized (85–250 bp)', 'Protection 251–500 bp',
          'Protection 501–1,000 bp', 'Protection >1,000 bp', 'MSP accessibility']
COLORS = ['#167b77','#65a5a1','#b7a0c9','#72529a','#ba77a1','#c38644','#8f4a30','#2376b3']


def regions(wt):
    smooth = np.convolve(np.pad(wt, (2,2), mode='edge'), np.ones(5)/5, mode='valid')
    peaks, _ = find_peaks(smooth, prominence=.005, distance=10)
    result = []
    for i, pos in enumerate(peaks):
        a = min(max(0, int(pos)-10), LENGTH-20)
        result.append(dict(id=f'WT_P{i+1}', kind='WT_peak_neighborhood', start0=a, end0=a+20,
                           peak0=int(pos), oligo_peak1=int(pos)+21,
                           wt_peak_coverage=float(wt[pos])))
    for a in list(range(0,192,10))+[191]:
        result.append(dict(id=f'fixed_{a:03d}_{a+20:03d}',kind='fixed_window',start0=a,end0=a+20,
                           peak0=None,oligo_peak1=None,wt_peak_coverage=None))
    return result


def at_change(ref, alt):
    before, after = ref in 'AT', alt in 'AT'
    return 'loss' if before and not after else 'gain' if after and not before else 'retained_count'


def lists(array):
    array = np.asarray(array)
    if array.ndim == 0:
        return round(float(array),7) if np.isfinite(array) else None
    return [lists(x) for x in array]


def run(config, output=None, mpra_catalog='results/pilot_20261014/catalog/catalog.json',
        coordinate_audit=DEFAULT_COORDINATE_AUDIT):
    config=Path(config)
    cfg_source=source_info(config,True)
    cfg=read_json(config)
    root=Path(cfg['output_dir'])
    land=root/'landscape'
    output=Path(output or root/'atlas')
    assets=Path(__file__).resolve().parent.parent/'dashboard'
    paths=[config,land/'audit.json',land/'landscape.npz',land/'variants.tsv',
           land/'wt_call_lengths.tsv',land/'wt_tf_boundaries.tsv',Path(cfg['catalog']),
           Path(cfg['heterogeneity']),Path(mpra_catalog),Path(__file__),
           Path(__file__).with_name('landscape.py'),Path(__file__).with_name('io.py'),
           Path(__file__).with_name('processing.py'),Path(__file__).with_name('genomic_coordinates.py'),
           Path(coordinate_audit),assets/'landscape.html',assets/'landscape.js']
    sources=[source_info(p,True) for p in paths]
    if sources[0]!=cfg_source:raise ValueError('Config changed while reading')
    audit=read_json(land/'audit.json')
    if audit.get('schema_version')!='fsmpra.landscape.v1' or not audit.get('complete'):
        raise ValueError('Complete landscape required')
    check_sources(audit['sources'])
    with np.load(land/'landscape.npz',allow_pickle=False) as f:
        data={k:f[k] for k in f.files}
    with open(land/'variants.tsv') as h: variants=list(csv.DictReader(h,delimiter='\t'))
    if [v['variant_id'] for v in variants]!=list(data['variant_ids']):raise ValueError('Variant order differs')
    if list(data['tracks'])!=list(TRACKS):raise ValueError('Track order differs')
    wi=next(i for i,v in enumerate(variants) if v['genotype_class']=='WT')
    wt=data['coverage'][wi]
    library=read_json(cfg['catalog'])['library']
    coordinate_check=read_json(coordinate_audit)
    mapping=mapping_for(library)
    if (coordinate_check.get('schema_version')!='fsmpra.genomic_coordinates.v1'
            or not coordinate_check.get('complete') or coordinate_check.get('mapping')!=mapping
            or coordinate_check.get('mutable_sequence')!=library['mutable_sequence']
            or coordinate_check.get('library_id')!=library['id']):
        raise ValueError('Complete matching genomic coordinate audit required')
    if mapping['reference_start0']!=START or mapping['length']!=LENGTH:
        raise ValueError('Landscape and genomic element differ')
    sources.extend(coordinate_check['sources'])
    check_sources(sources)
    region_rows=regions(wt[0])
    for r in region_rows:
        r.update(assembly=mapping['assembly'],chrom=mapping['chrom'],strand=mapping['strand'],
                 **interval_fields(mapping,r['start0'],r['end0']),
                 genomic_peak1=mapping['genomic_start1']+r['peak0'] if r['peak0'] is not None else None)
    context={v['variant_id']:{} for v in variants}
    with open_text(cfg['heterogeneity']) as h:
        for row in csv.DictReader(h,delimiter='\t'):
            if row['cohort']=='footprint_analysis_eligible_with_barcode' and row['variant_id'] in context:
                context[row['variant_id']][row['feature']]=float(row['barcode_mean']) if row['barcode_mean'] else None
    mpra=read_json(mpra_catalog)
    sources.extend(x['source'] for x in mpra['mpra_audit'])
    check_sources(sources)
    previous=next(x for x in mpra['libraries'] if x['id']==library['id'])
    if previous['mutable_sequence']!=library['mutable_sequence'] or previous['genomic_mutable_start1']!=library['genomic_mutable_start1']:
        raise ValueError('Published MPRA mapping library differs')
    oldseq={v['variant_id']:v['sequence'] for v in previous['variants']}
    for v in library['variants']:
        if v['variant_id'] in context and oldseq.get(v['variant_id'])!=v['sequence']:
            raise ValueError('Published MPRA mapping sequence differs')
    functional={v['variant_id']:{} for v in variants}
    for row in mpra['mpra']:
        if row.get('join_status')=='matched' and row['variant_id'] in functional:
            target=functional[row['variant_id']]
            if row['mpra_source_id'] in target:raise ValueError('Duplicate functional mapping')
            target[row['mpra_source_id']]=float(row['effect_log2'])
    summary_rows, effect_rows, bundle_variants=[],[],[]
    for i,v in enumerate(variants):
        supported=bool(data['n_barcodes'][i]>=50 and data['n_molecules'][i]>=100)
        pos=int(v['ref_pos0'])-START if v['ref_pos0'] else None
        base=dict(variant_id=v['variant_id'],source_variant=v['source_variant'],relative_pos0=pos,
                  oligo_pos1=pos+21 if pos is not None else None,
                  genomic_pos1=library['genomic_mutable_start1']+pos if pos is not None else None,
                  ref=v['ref'],alt=v['alt'],at_change=at_change(v['ref'],v['alt']) if pos is not None else 'WT',
                  n_barcodes=int(data['n_barcodes'][i]),n_molecules=int(data['n_molecules'][i]),coverage_floor_met=supported)
        base.update(assembly=mapping['assembly'],chrom=mapping['chrom'],strand=mapping['strand'],
                    reference_pos0=START+pos if pos is not None else None,
                    genomic_variant=f"{mapping['chrom']}:{base['genomic_pos1']}:{v['ref']}>{v['alt']}" if pos is not None else 'WT')
        row=dict(base)
        for t,track in enumerate(TRACKS):
            row[track+'_coverage']=number(data['coverage'][i,t].mean())
            row[track+'_delta_pp']=number(100*(data['coverage'][i,t].mean()-wt[t].mean()))
            row[track+'_presence']=number(data['presence'][i,t])
        for feature in ['distal_m6a_fraction','distal_nucleosome_fraction','whole_m6a_fraction']:
            value=context[v['variant_id']].get(feature)
            reference=context[variants[wi]['variant_id']].get(feature)
            row[feature]=value
            row[feature+'_delta_pp']=100*(value-reference) if value is not None and reference is not None else None
        for source in mpra['mpra_sources']:row[source['id']+'_log2']=functional[v['variant_id']].get(source['id'])
        effects=[]
        for region in region_rows:
            a,b=region['start0'],region['end0']
            value=number(data['coverage'][i,0,a:b].mean())
            delta=100*(value-float(wt[0,a:b].mean())) if value is not None else None
            effects.append(delta)
            effect_rows.append(dict(base,region_id=region['id'],region_kind=region['kind'],
                                    start0=a,end0=b,**interval_fields(mapping,a,b),wt_coverage=float(wt[0,a:b].mean()),
                                    variant_coverage=value,delta_pp=delta))
        summary_rows.append(row)
        bundle_variants.append(dict(base,coverage=lists(data['coverage'][i]),presence=lists(data['presence'][i]),
                                    states=lists(data['states'][i]),region_effects=effects,
                                    context=context[v['variant_id']],functional=functional[v['variant_id']]))
    correlations=[]
    for source in mpra['mpra_sources']:
        for metric in ['tf_small_delta_pp','msp_delta_pp','nuc_85_250_delta_pp']:
            rows=[r for r in summary_rows if r['source_variant']!='WT' and r['coverage_floor_met'] and
                  r[metric] is not None and r[source['id']+'_log2'] is not None]
            x=[r[metric] for r in rows];y=[r[source['id']+'_log2'] for r in rows]
            correlations.append(dict(source=source['id'],metric=metric,n=len(rows),
                                     pearson=float(pearsonr(x,y)[0]),spearman=float(spearmanr(x,y)[0])))
    bundle=dict(schema_version='fsmpra.landscape_atlas.v2',sample_id=cfg['sample_id'],
                coordinates=dict(mapping,verified=True,display='genomic 1-based inclusive',interval_exports='0-based half-open and 1-based inclusive'),
                sequence=library['mutable_sequence'],element_reference_start0=START,
                genomic_start1=library['genomic_mutable_start1'],tracks=list(TRACKS),track_labels=LABELS,
                colors=COLORS,wt_index=wi,regions=region_rows,variants=bundle_variants,
                mpra_sources=mpra['mpra_sources'],descriptive_only=True)
    check_sources(sources)
    with atomic_output(output) as out:
        write_table(out/'variant_summary.tsv',summary_rows,list(summary_rows[0]))
        write_table(out/'region_effects.tsv.gz',effect_rows,list(effect_rows[0]))
        write_table(out/'wt_regions.tsv',region_rows,list(region_rows[0]))
        write_table(out/'functional_correlations.tsv',correlations,list(correlations[0]))
        write_coordinate_exports(out,mapping,library['mutable_sequence'],data,variants,land,region_rows)
        write_json(out/'coordinate_audit.json',coordinate_check)
        write_json(out/'atlas.json',bundle)
        text=json.dumps(bundle,separators=(',',':'),allow_nan=False).replace('&','\\u0026').replace('<','\\u003c').replace('>','\\u003e')
        html=(assets/'landscape.html').read_text().replace('/* APP_PLACEHOLDER */',(assets/'landscape.js').read_text())
        html=html.replace('DATA_PLACEHOLDER',text)
        (out/'index.html').write_text(html)
        plots(out,data,wi,region_rows,summary_rows,land,mapping)
        write_json(out/'audit.json',dict(schema_version='fsmpra.landscape_atlas.v2',complete=True,
                   coordinates=mapping,coordinate_audit=source_info(coordinate_audit,True),
                   descriptive_only=True,variants=len(variants),wt_peak_neighborhoods=sum(r['kind']=='WT_peak_neighborhood' for r in region_rows),
                   fixed_windows=21,mpra_relationship='external published assays, retained separately; not matched preparations',
                   sources=sources))
        check_sources(sources)
    return dict(output=str(output),regions=len(region_rows),variants=len(variants))


def write_coordinate_exports(out,mapping,sequence,data,variants,land,region_rows):
    bases=[dict(assembly=mapping['assembly'],chrom=mapping['chrom'],strand=mapping['strand'],
                genomic_pos1=mapping['genomic_start1']+p,relative_pos0=p,
                reference_pos0=START+p,oligo_pos1=mapping['oligo_start1']+p,ref=ref)
           for p,ref in enumerate(sequence)]
    write_table(out/'base_coordinates.tsv',bases,list(bases[0]))
    def profiles():
        for i,v in enumerate(variants):
            for t,track in enumerate(TRACKS):
                for p,base in enumerate(bases):
                    yield dict(base,variant_id=v['variant_id'],source_variant=v['source_variant'],track=track,
                               mean_equal_barcode=number(data['coverage'][i,t,p]),
                               mean_molecule_weighted=number(data['molecule_coverage'][i,t,p]),
                               n_barcodes=int(data['n_barcodes'][i]),n_molecules=int(data['n_molecules'][i]))
    fields=list(bases[0])+['variant_id','source_variant','track','mean_equal_barcode',
                          'mean_molecule_weighted','n_barcodes','n_molecules']
    write_table(out/'genomic_profiles.tsv.gz',profiles(),fields)
    with open(land/'wt_tf_boundaries.tsv') as h:
        boundaries=[dict(r,assembly=mapping['assembly'],chrom=mapping['chrom'],strand=mapping['strand'],
                         **interval_fields(mapping,int(r['relative_start0']),int(r['relative_end0'])))
                    for r in csv.DictReader(h,delimiter='\t')]
    write_table(out/'wt_tf_boundaries_genomic.tsv',boundaries,list(boundaries[0]))
    with open(out/'wt_regions.bed','w') as h:
        for r in region_rows:
            h.write(f"{mapping['chrom']}\t{r['genomic_start0']}\t{r['genomic_end0']}\t{r['id']}\t0\t{mapping['strand']}\n")


def plots(out,data,wi,regions,rows,land,mapping):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    origin=mapping['genomic_start1']
    x=np.arange(LENGTH)+origin
    xlabel=f"{mapping['assembly']} {mapping['chrom']} genomic position (1-based; forward strand)"
    wt=data['coverage'][wi]
    fig,axes=plt.subplots(3,1,figsize=(13,10),sharex=True)
    for t in [0,1,2]:axes[0].plot(x,100*wt[t],label=LABELS[t],color=COLORS[t])
    for t in [3,4,5,6]:axes[1].plot(x,100*wt[t],label=LABELS[t],color=COLORS[t])
    axes[2].plot(x,100*wt[7],label=LABELS[7],color=COLORS[7])
    for r in regions:
        if r['kind']=='WT_peak_neighborhood':
            axes[0].axvspan(r['genomic_start1'],r['genomic_end1'],color='grey',alpha=.08)
            axes[0].annotate(r['id'],(r['genomic_peak1'],100*wt[0,r['peak0']]),xytext=(0,10),textcoords='offset points',ha='center',fontsize=8)
    for ax in axes:
        ax.set(ylabel='Covered fraction (%)');ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.15)
    axes[2].set(xlabel=xlabel,xlim=(origin,origin+LENGTH-1))
    axes[2].ticklabel_format(axis='x',style='plain',useOffset=False)
    fig.suptitle('WT observed footprint landscape · equal barcode weighting\nOriginal call lengths; tracks can overlap; peaks are descriptive neighborhoods')
    fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(out/f'wt_landscape.{ext}',dpi=160)
    plt.close(fig)
    with open(land/'wt_call_lengths.tsv') as h:length_rows=list(csv.DictReader(h,delimiter='\t'))
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    for kind,color in [('tf',COLORS[0]),('nuc',COLORS[3])]:
        selected=[r for r in length_rows if r['caller_track']==kind]
        axes[0].plot([int(r['original_query_length']) for r in selected],
                     [float(r['weighted_calls_per_barcode']) for r in selected],'.',ms=3,label=kind,color=color)
    axes[0].set(xscale='log',yscale='log',xlabel='Original call length (bp; log)',ylabel='Weighted calls per barcode (log)')
    axes[0].legend();axes[0].set_title('Calls overlapping the element')
    labels=['None of these','TF','Nuc-sized','TF + nuc-sized','Long','TF + long','Nuc-sized + long','All three']
    axes[1].barh(labels,100*data['states'][wi],color='#72529a');axes[1].invert_yaxis()
    axes[1].set(xlabel='Barcode-weighted molecule fraction (%)',title='Co-occurring call categories anywhere in element')
    fig.suptitle('Nuc-sized: 85–250 bp · Long: >250 bp · “None” does not mean protein-free')
    fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(out/f'wt_lengths_states.{ext}',dpi=160)
    plt.close(fig)
    metrics=['tf_small_delta_pp','msp_delta_pp','nuc_85_250_delta_pp']
    fig,axes=plt.subplots(3,1,figsize=(15,6),sharex=True)
    for ax,metric in zip(axes,metrics):
        mat=np.full((4,LENGTH),np.nan)
        for r in rows:
            if r['relative_pos0'] is not None and r[metric] is not None:
                mat['ACGT'.index(r['alt']),r['relative_pos0']]=r[metric]
        limit=max(.01,float(np.nanmax(np.abs(mat))))
        im=ax.imshow(mat,aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit,extent=[origin-.5,origin+LENGTH-.5,3.5,-.5])
        ax.set(yticks=range(4),yticklabels=list('ACGT'),ylabel='Alternate base',title=metric.replace('_delta_pp','')+' · whole-element coverage difference')
        fig.colorbar(im,ax=ax,label='Percentage points')
    axes[-1].set_xlabel(xlabel+'; blank = WT base or no eligible data')
    axes[-1].ticklabel_format(axis='x',style='plain',useOffset=False)
    fig.suptitle('All SNVs, including low support · descriptive differences, no significance labels')
    fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(out/f'variant_effects.{ext}',dpi=160)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--output')
    p.add_argument('--mpra-catalog',default='results/pilot_20261014/catalog/catalog.json')
    p.add_argument('--coordinate-audit',default=DEFAULT_COORDINATE_AUDIT)
    a=p.parse_args();print(json.dumps(run(a.config,a.output,a.mpra_catalog,a.coordinate_audit)))


if __name__=='__main__':main()

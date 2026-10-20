"""Build sequence, motif and barcode-sensitivity annotations for a frozen atlas."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from . import biology_audit, biology_motifs
from .io import atomic_output,read_json,source_info,write_json,write_table
from .processing import check_sources


def run(atlas,output):
    atlas=Path(atlas);bundle=read_json(atlas/'atlas.json')
    if bundle.get('schema_version')!='fsmpra.landscape_atlas.v2' or not bundle['coordinates']['verified']:
        raise ValueError('Verified genomic atlas required')
    check_sources(read_json(atlas/'audit.json')['sources'])
    motif_files=[('extras/data/JASPAR2024_CORE_vertebrates_nr.jaspar','JASPAR2024'),
                 ('extras/data/HOCOMOCOv12_H12CORE_pcm.txt','HOCOMOCOv12')]
    sources=[source_info(p,True) for p in [atlas/'atlas.json',atlas/'audit.json',Path(__file__),
             Path(biology_motifs.__file__),Path(biology_audit.__file__),*[p for p,_ in motif_files]]]
    matrices=[];db_counts={}
    for path,db in motif_files:
        records=biology_motifs.read_matrices(path,db);db_counts[db]=len(records);matrices+=records
    if db_counts!={'JASPAR2024':879,'HOCOMOCOv12':1443}:raise ValueError('Unexpected versioned collection sizes')
    hits,settings=biology_motifs.scan(bundle['sequence'],matrices,bundle['genomic_start1'])
    print(f"Scanned {len(matrices)} matrices: {len(hits)} WT candidate matches",flush=True)
    audits,checks=biology_audit.run(bundle,'results/rep3_statistics_20261017/aggregates',
                                 'results/rep3_heterogeneity_20261017/summary/barcode_heterogeneity.tsv.gz')
    sources.extend(checks.pop('sources'))
    by_audit={(r['variant_id'],r['region_id']):r for r in audits}
    wt=bundle['variants'][bundle['wt_index']]
    regions=[dict(r) for r in bundle['regions'] if r['kind']=='WT_peak_neighborhood']
    regions.append(dict(id='adjacent_143_150',kind='post_hoc_adjacent_pattern',start0=122,end0=130,
                        genomic_start1=bundle['genomic_start1']+122,genomic_end1=bundle['genomic_start1']+129))
    summaries=[];region_variants=[];motif_effects=[];concordance=[]
    for r in regions:
        a,b=r['start0'],r['end0']
        candidates=[v for v in bundle['variants'] if v['relative_pos0'] is not None and a<=v['relative_pos0']<b]
        local_hits=[h for h in hits if h['start0']<b and h['end0']>a]
        r['sequence']=bundle['sequence'][a:b];r['at_positions']=sum(base in 'AT' for base in r['sequence'])
        r['wt_tf_coverage']=float(np.mean(wt['coverage'][0][a:b]));r['motif_hit_ids']=[h['hit_id'] for h in local_hits]
        summary=dict(region_id=r['id'],genomic_start1=r['genomic_start1'],genomic_end1=r['genomic_end1'],sequence=r['sequence'],
                     n_designed_snvs=len(candidates),n_measured=sum(v['n_barcodes']>0 for v in candidates),
                     wt_tf_coverage=r['wt_tf_coverage'],at_positions=r['at_positions'],motif_hits=len(local_hits),
                     motif_hits_both_backgrounds=sum(h['passes_both_backgrounds'] for h in local_hits))
        for v in candidates:
            row=dict(by_audit[v['variant_id'],r['id']])
            for track in (7,3):
                row[['msp_delta_pp','nuc_85_250_delta_pp'][track==3]]=100*(float(np.mean(v['coverage'][track]))-float(np.mean(wt['coverage'][track]))) if v['n_barcodes'] else None
            row.update({source['id']+'_log2':v['functional'].get(source['id']) for source in bundle['mpra_sources']})
            region_variants.append(row)
        for field in ('delta_pp','msp_delta_pp','nuc_85_250_delta_pp','published_LDLR_log2','published_LDLR_2_log2'):
            vals=[row[field] for row in region_variants if row['region_id']==r['id'] and row[field] is not None]
            summary[field+'_median']=float(np.median(vals)) if vals else None
            summary[field+'_negative']=sum(x<0 for x in vals);summary[field+'_positive']=sum(x>0 for x in vals)
        summaries.append(summary)
    for hit in hits:
        local_effects=[]
        for v in bundle['variants']:
            delta=biology_motifs.variant_delta(hit,v['relative_pos0'],v['ref'],v['alt'])
            if delta is not None:
                row=dict(hit_id=hit['hit_id'],variant_id=v['variant_id'],genomic_variant=v['genomic_variant'],
                                         source_variant=v['source_variant'],delta_score_bits=delta,
                                         wt_score_bits=hit['score_bits'],variant_score_bits=hit['score_bits']+delta,
                                         n_barcodes=v['n_barcodes'],at_change=v['at_change'],
                                         tf_local_delta_pp=100*float(np.mean(np.array(v['coverage'][0][hit['start0']:hit['end0']])-np.array(wt['coverage'][0][hit['start0']:hit['end0']]))) if v['n_barcodes'] else None,
                                         msp_delta_pp=by_audit[v['variant_id'],'whole_msp']['delta_pp'],
                                         published_LDLR_log2=v['functional'].get('published_LDLR'),
                                         published_LDLR_2_log2=v['functional'].get('published_LDLR_2'))
                motif_effects.append(row);local_effects.append(row)
        # These descriptive correlations retain every overlapping measured SNV;
        # they are not tests and are not used to select the motif matches.
        c=dict(hit_id=hit['hit_id'],database=hit['database'],matrix_id=hit['matrix_id'],name=hit['name'],
               genomic_start1=hit['genomic_start1'],genomic_end1=hit['genomic_end1'],n_designed=len(local_effects))
        for field in ('tf_local_delta_pp','msp_delta_pp','published_LDLR_log2','published_LDLR_2_log2'):
            measured=[r for r in local_effects if r[field] is not None]
            x=[r['delta_score_bits'] for r in measured];y=[r[field] for r in measured]
            c[field+'_n']=len(measured)
            c[field+'_spearman']=float(spearmanr(x,y)[0]) if len(x)>=3 and len(set(x))>1 and len(set(y))>1 else None
        concordance.append(c)
    # Independently published sequence anchor; no motif scan threshold is used
    # to decide its placement or to assign measured binding to a protein.
    sre='ATCACCCCAC';pos=bundle['sequence'].find(sre)
    if pos<0 or bundle['sequence'].find(sre,pos+1)>=0:raise ValueError('Expected unique literature SRE-1 sequence')
    literature=[dict(name='Published LDLR SRE-1 sequence',sequence=sre,start0=pos,end0=pos+len(sre),
                     genomic_start1=bundle['genomic_start1']+pos,genomic_end1=bundle['genomic_start1']+pos+len(sre)-1,
                     url='https://pmc.ncbi.nlm.nih.gov/articles/PMC20191/',
                     interpretation='Exact sequence anchor from published LDLR promoter work; not evidence that SREBP occupies these plasmids.')]
    result=dict(schema_version='fsmpra.biology.v1',descriptive_only=True,coordinates=bundle['coordinates'],
                databases=db_counts,scan=settings,motifs=hits,literature=literature,regions=regions,
                region_summaries=summaries,region_variants=region_variants,barcode_audits=audits,motif_concordance=concordance,
                checks=checks,limitations=[
                    'Motif redundancy and shared DNA preferences prevent unique TF assignment; database agreement is not independent replication.',
                    'No HepG2 expression or locus-specific ChIP evidence has been added.',
                    'Shared A/T checks reuse original calls and cannot remove sequence-dependent calling or methylation bias.',
                    'Distal-feature standardization may condition on biology or a mediator; it is a sensitivity comparison, not a causal correction.',
                    'One preparation and post hoc candidate selection; no new variant significance, confidence intervals or activator/repressor labels.'])
    check_sources(sources)
    with atomic_output(output) as out:
        write_json(out/'biology.json',result)
        for name,rows in [('candidate_region_summary.tsv',summaries),('candidate_region_variants.tsv',region_variants),
                          ('barcode_sensitivity.tsv',audits),('motif_variant_effects.tsv.gz',motif_effects),('motif_concordance.tsv',concordance)]:
            write_table(out/name,rows,list(rows[0]))
        flat=[{k:v for k,v in h.items() if k!='weights_genomic'} for h in hits]
        write_table(out/'motif_hits.tsv',flat,list(flat[0]))
        write_json(out/'audit.json',dict(schema_version='fsmpra.biology.v1',complete=True,sources=sources,
                   counts=dict(motifs=len(hits),motif_variant_effects=len(motif_effects),barcode_sensitivity_rows=len(audits)),checks=checks))
        check_sources(sources)
    return dict(output=str(output),hits=len(hits),regions=len(regions))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--atlas',default='results/rep3_landscape_20261018/atlas_genomic_20261019')
    p.add_argument('--output',required=True)
    a=p.parse_args();print(json.dumps(run(a.atlas,a.output)))

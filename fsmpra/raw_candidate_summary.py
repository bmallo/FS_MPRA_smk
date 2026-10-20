"""Predefined raw-call and read-filter sensitivity comparisons, without tests.

WT barcode and read samples remain distinct. No filter is chosen for agreement
with the published activity data, and original FiberHMM calls are not recalled.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .io import atomic_output,open_text,read_json,source_info,write_json,write_table
from .processing import check_sources

REGIONS=[('WT_P3',93,113),('adjacent_143_150',122,130),('adjacent_flanks',112,148),('promoter',0,211)]
SCENARIOS=[
 ('primary_barcodes','barcode','Supported barcodes, equal barcode weight',128,20),
 ('primary_molecules','molecule','Same supported-barcode cohort, equal molecule weight',128,20),
 ('genotype_only','molecule','Direct genotype; circular/eligible; no barcode requirement',128,20),
 ('genotype_all_full_element','molecule','Direct genotype; full element; circularity requirement removed',128,20),
 ('genotype_noncircular','molecule','Direct genotype; noncircular subset; full element',128,20),
 ('rq_9999','barcode','Supported cohort; read quality >=0.9999',128,20),
 ('passes20','barcode','Supported cohort; >=20 passes',128,20),
 ('passes40','barcode','Supported cohort; >=40 passes',128,20),
 ('mapq60','barcode','Supported cohort; mapping quality >=60',128,20),
 ('length_1pct','barcode','Supported cohort; read length within 1% of 4,647 bp',128,20),
 ('alignment_plus','barcode','Supported cohort; forward alignment',128,20),
 ('alignment_minus','barcode','Supported cohort; reverse alignment',128,20),
 ('barcode_depth3','barcode','Supported cohort; >=3 eligible molecules per barcode',128,20),
 ('barcode_depth5','barcode','Supported cohort; >=5 eligible molecules per barcode',128,20),
 ('baseq30','barcode','Supported cohort; callable positions require base quality >=30',128,30),
 ('retained_ml250','barcode','Supported cohort; retained ML threshold >=250 (not original probability calibration)',250,20),
 ('whole_m6a_below90','barcode','Supported cohort; whole-plasmid m6A fraction <0.90 (phenotype-dependent)',128,20),
 ('distal_m6a_10_80','barcode','Supported cohort; distal m6A 0.10–0.80 (phenotype-dependent)',128,20)]


def nmean(x,axis=0):
    valid=np.isfinite(x);n=valid.sum(axis=axis);s=np.where(valid,x,0).sum(axis=axis)
    return np.divide(s,n,out=np.full(np.shape(s),np.nan,dtype=float),where=n>0)


def barcode_means(values,barcodes):
    """Mean of observed molecules within barcode; missing outcomes stay missing."""
    order=np.argsort(barcodes,kind='stable');keys=barcodes[order]
    starts=np.r_[0,np.flatnonzero(keys[1:]!=keys[:-1])+1]
    if not len(keys):return np.empty((0,)+values.shape[1:]),np.array([],dtype=str)
    arr=values[order];valid=np.isfinite(arr)
    numer=np.add.reduceat(np.where(valid,arr,0),starts,axis=0);denom=np.add.reduceat(valid.astype(np.int32),starts,axis=0)
    return np.divide(numer,denom,out=np.full(numer.shape,np.nan),where=denom>0),keys[starts]


def scalar(values,barcodes,weighting):
    if not len(values):return dict(value=None,molecules=0,barcodes=0,nonzero_units=0)
    valid=np.isfinite(values)
    units=barcode_means(values,barcodes)[0] if weighting=='barcode' else values
    usable=units[np.isfinite(units)]
    return dict(value=float(usable.mean()) if len(usable) else None,molecules=int(valid.sum()),
                barcodes=len(set(barcodes[valid])- {''}),nonzero_units=int((usable>0).sum()))


def shared(sequence,start,end,position,base=None,radius=2):
    return np.array([p for p in range(start,end) if sequence[p] in (base or 'AT') and (position is None or abs(p-position)>radius)],dtype=int)


def methylation_outcome(scores,known,quality,positions,threshold=128,baseq=20):
    if not len(positions):return np.full(len(scores),np.nan)
    valid=(known[:,positions]&(quality[:,positions]>=baseq)).all(1)
    result=(scores[:,positions]>=threshold).mean(1)
    result[~valid]=np.nan
    return result


def resample_one_per_barcode(values,barcodes,rng,repeats=50):
    valid=np.isfinite(values)&(barcodes!='');values=values[valid];barcodes=barcodes[valid]
    if not len(values):return np.full(repeats,np.nan)
    order=np.argsort(barcodes,kind='stable');values=values[order];barcodes=barcodes[order]
    starts=np.r_[0,np.flatnonzero(barcodes[1:]!=barcodes[:-1])+1];counts=np.diff(np.r_[starts,len(values)])
    return np.array([values[starts+(rng.random(len(starts))*counts).astype(int)].mean() for _ in range(repeats)])


def joint_distribution(raw,coverage,footprints,barcodes):
    valid=np.isfinite(raw)&footprints&(barcodes!='')
    if not valid.any():return np.full(16,np.nan),0,0
    _,inverse,counts=np.unique(barcodes[valid],return_inverse=True,return_counts=True)
    weights=1/counts[inverse]/len(counts)
    protection=((coverage[valid]&4)>0).any(1).astype(int)+2*((coverage[valid]&8)>0).any(1).astype(int)
    methyl=np.searchsorted([.25,.5,.75],raw[valid],side='left')
    fractions=np.bincount(4*methyl+protection,weights=weights,minlength=16)
    if not np.isclose(fractions.sum(),1):raise ValueError('Joint molecule state mass lost')
    return fractions,int(valid.sum()),len(counts)


def run(extracted,output):
    extracted=Path(extracted);audit=read_json(extracted/'audit.json');check_sources(audit['sources'])
    paths=[extracted/'audit.json',extracted/'molecules.npz',extracted/'molecules.tsv.gz',extracted/'variants.json',Path(__file__)]
    sources=[source_info(p,True) for p in paths]
    with open_text(extracted/'molecules.tsv.gz') as h:meta=list(csv.DictReader(h,delimiter='\t'))
    with np.load(extracted/'molecules.npz') as f:arrays={k:f[k] for k in f.files}
    variants=read_json(extracted/'variants.json');sequence=variants[0]['sequence'];n=len(meta)
    if n!=audit['molecules'] or arrays['scores'].shape!=(n,211):raise ValueError('Extraction dimensions differ')
    vi=np.array([int(r['variant_index']) for r in meta]);bc=np.array([r['corrected_barcode'] for r in meta]);readids=np.array([r['read_id'] for r in meta])
    primary=np.array([r['footprint_analysis_eligible_with_barcode']=='1' for r in meta])
    barcode_sample=(vi!=0)|np.array([r['wt_barcode_sample']=='1' for r in meta]);read_sample=(vi!=0)|np.array([r['wt_genotype_sample']=='1' for r in meta])
    primary &= barcode_sample
    eligible=np.array([r['analysis_eligible']=='1' for r in meta])
    circular=np.array([r['circular']=='1' for r in meta])
    def field(name):return np.array([float(r[name]) if r[name] else np.nan for r in meta])
    depths={}
    for b in bc[primary]:depths[b]=depths.get(b,0)+1
    counts=np.array([depths.get(b,0) for b in bc])
    masks=dict(primary_barcodes=primary,primary_molecules=primary,genotype_only=read_sample&eligible,
               genotype_all_full_element=read_sample,genotype_noncircular=read_sample&~circular,
               rq_9999=primary&(field('read_quality')>=.9999),passes20=primary&(field('passes')>=20),
               passes40=primary&(field('passes')>=40),mapq60=primary&(field('mapq')>=60),
               length_1pct=primary&(np.abs(field('read_length')-4647)<=46.47),
               alignment_plus=primary&np.array([r['strand']=='+' for r in meta]),
               alignment_minus=primary&np.array([r['strand']=='-' for r in meta]),
               barcode_depth3=primary&(counts>=3),barcode_depth5=primary&(counts>=5),
               baseq30=primary,retained_ml250=primary,whole_m6a_below90=primary&(field('whole_m6a_fraction')<.9),
               distal_m6a_10_80=primary&(field('distal_m6a_fraction')>=.1)&(field('distal_m6a_fraction')<=.8))
    metrics=['retained_m6a_shared','retained_m6a_shared_A','retained_m6a_shared_T','low_call_fraction_le25',
             'tf_small','msp','nuc_85_250','nuc_gt250']
    rows=[];profiles=[];retention=[];repetitions=[];joint=[]
    for scenario,weighting,label,threshold,baseq in SCENARIOS:
        print('Summarizing '+scenario,flush=True)
        control_cache={}
        selected={j:np.flatnonzero(masks[scenario]&(vi==j)) for j in range(len(variants))}
        for j,v in enumerate(variants):
            idx=selected[j];barcodes=bc[idx];S=arrays['scores'][idx];K=arrays['known'][idx];Q=arrays['qualities'][idx];C=arrays['coverage'][idx];F=arrays['footprints_present'][idx].astype(bool)
            retention.append(dict(scenario=scenario,source_variant=v['source_variant'],variant_id=v['variant_id'],weighting=weighting,
                                  molecules=len(idx),barcodes=len(set(barcodes)-{''}),footprint_molecules=int(F.sum()),
                                  supported_molecules=int(np.array([meta[i]['footprint_analysis_eligible_with_barcode']=='1' for i in idx]).sum())))
            valid=K&(Q>=baseq)
            m=np.where(valid,S>=threshold,np.nan)
            fp=[np.where(F[:,None],(C&flag)>0,np.nan) for flag in [1,2,4,8]]
            if weighting=='barcode':
                pm,pb=barcode_means(m,barcodes);mp=nmean(pm);support=np.isfinite(pm).sum(0)
                fps=[nmean(barcode_means(x,barcodes)[0]) for x in fp]
            else:mp=nmean(m);support=valid.sum(0);fps=[nmean(x) for x in fp]
            edge=nmean(np.where(F[:,None],arrays['edges'][idx],np.nan))
            for p in range(211):
                profiles.append(dict(scenario=scenario,source_variant=v['source_variant'],variant_id=v['variant_id'],genomic_pos1=11089285+p,
                                     base=v['sequence'][p],adenine_strand='+' if v['sequence'][p]=='A' else '-' if v['sequence'][p]=='T' else '',
                                     retained_m6a_fraction=float(mp[p]) if np.isfinite(mp[p]) else None,callable_units=int(support[p]),callable_molecules=int(valid[:,p].sum()),
                                     tf_small=float(fps[0][p]) if np.isfinite(fps[0][p]) else None,msp=float(fps[1][p]) if np.isfinite(fps[1][p]) else None,
                                     nuc_85_250=float(fps[2][p]) if np.isfinite(fps[2][p]) else None,nuc_gt250=float(fps[3][p]) if np.isfinite(fps[3][p]) else None,
                                     tf_boundary_molecule_fraction=float(edge[p]) if np.isfinite(edge[p]) else None))
            if j==0:continue
            wtidx=selected[0];wbc=bc[wtidx];WS=arrays['scores'][wtidx];WK=arrays['known'][wtidx];WQ=arrays['qualities'][wtidx];WC=arrays['coverage'][wtidx];WF=arrays['footprints_present'][wtidx].astype(bool)
            pos=v['oligo_pos1']-21
            for region,a,b in REGIONS:
                values={};wtvalues={};opps={}
                for name,base in [('retained_m6a_shared',None),('retained_m6a_shared_A','A'),('retained_m6a_shared_T','T')]:
                    positions=shared(sequence,a,b,pos,base)
                    values[name]=methylation_outcome(S,K,Q,positions,threshold,baseq)
                    wtvalues[name]=methylation_outcome(WS,WK,WQ,positions,threshold,baseq);opps[name]=len(positions)
                name='low_call_fraction_le25'
                values[name]=np.where(np.isfinite(values['retained_m6a_shared']),values['retained_m6a_shared']<=.25,np.nan)
                wtvalues[name]=np.where(np.isfinite(wtvalues['retained_m6a_shared']),wtvalues['retained_m6a_shared']<=.25,np.nan)
                opps[name]=opps['retained_m6a_shared']
                for name,flag in [('tf_small',1),('msp',2),('nuc_85_250',4),('nuc_gt250',8)]:
                    values[name]=np.where(F,((C[:,a:b]&flag)>0).mean(1),np.nan)
                    wtvalues[name]=np.where(WF,((WC[:,a:b]&flag)>0).mean(1),np.nan);opps[name]=b-a
                if scenario=='primary_barcodes':
                    for group,raw,coverage,present,labels in [('variant',values['retained_m6a_shared'],C[:,a:b],F,barcodes),('WT',wtvalues['retained_m6a_shared'],WC[:,a:b],WF,wbc)]:
                        fractions,nm,nb=joint_distribution(raw,coverage,present,labels)
                        for state,fraction in enumerate(fractions):
                            joint.append(dict(source_variant=v['source_variant'],region=region,group=group,
                                              methylation_bin=['0_to_25','gt25_to_50','gt50_to_75','gt75_to_100'][state//4],
                                              protection_state=['neither_category','nuc_85_250_only','long_only','both_categories'][state%4],
                                              fraction=float(fraction) if np.isfinite(fraction) else None,molecules=nm,barcodes=nb))
                for name in metrics:
                    obs=scalar(values[name],barcodes,weighting)
                    # Many variants share the identical unchanged-site mask.
                    mask_key=tuple(shared(sequence,a,b,pos,{'retained_m6a_shared_A':'A','retained_m6a_shared_T':'T'}.get(name))) if name.startswith('retained_') or name=='low_call_fraction_le25' else ()
                    key=(region,name,mask_key)
                    if key not in control_cache:control_cache[key]=scalar(wtvalues[name],wbc,weighting)
                    control=control_cache[key]
                    delta=100*(obs['value']-control['value']) if obs['value'] is not None and control['value'] is not None else None
                    rows.append(dict(scenario=scenario,label=label,weighting=weighting,source_variant=v['source_variant'],variant_id=v['variant_id'],
                                     genomic_variant=f"chr19:{v['genomic_pos1']}:{v['ref']}>{v['alt']}",region=region,genomic_start1=11089285+a,genomic_end1=11089285+b-1,
                                     metric=name,shared_opportunities=opps[name],variant_value=obs['value'],wt_value=control['value'],delta_pp=delta,
                                     variant_molecules=obs['molecules'],wt_molecules=control['molecules'],variant_barcodes=obs['barcodes'],wt_barcodes=control['barcodes'],
                                     variant_nonzero_units=obs['nonzero_units'],wt_nonzero_units=control['nonzero_units']))
                    if scenario=='primary_barcodes' and ((region=='WT_P3' and v['oligo_pos1']==125) or (region=='promoter' and 143<=v['oligo_pos1']<=150)) and name in ('retained_m6a_shared','tf_small','msp','nuc_85_250'):
                        rng=np.random.default_rng(20261020+j)
                        vv=resample_one_per_barcode(values[name],barcodes,rng)
                        ww=resample_one_per_barcode(wtvalues[name],wbc,rng)
                        differences=100*(vv-ww)
                        for rep,difference in enumerate(differences):repetitions.append(dict(source_variant=v['source_variant'],region=region,metric=name,kind='one_molecule_per_barcode',repeat=rep,delta_pp=float(difference) if np.isfinite(difference) else None))
                        # Same-size WT subsets supply a descriptive sampling scale,
                        # not calibrated p-values or confidence intervals.
                        wb,_=barcode_means(wtvalues[name],wbc);wb=wb[np.isfinite(wb)]
                        vb,_=barcode_means(values[name],barcodes);size=int(np.isfinite(vb).sum())
                        if 0<size<len(wb):
                            for rep in range(100):
                                chosen=rng.choice(len(wb),size=size,replace=False);total=wb[chosen].sum()
                                difference=100*(total/size-(wb.sum()-total)/(len(wb)-size))
                                repetitions.append(dict(source_variant=v['source_variant'],region=region,metric=name,kind='WT_subset_sampling_scale',repeat=rep,delta_pp=float(difference)))
    check_sources(sources)
    with atomic_output(output) as out:
        write_table(out/'comparisons.tsv',rows,list(rows[0]));write_table(out/'profiles.tsv.gz',profiles,list(profiles[0]))
        write_table(out/'retention.tsv',retention,list(retention[0]));write_table(out/'resampling.tsv',repetitions,list(repetitions[0]))
        write_table(out/'joint_molecule_states.tsv',joint,list(joint[0]))
        write_json(out/'audit.json',dict(schema_version='fsmpra.raw_candidate_summary.v1',complete=True,sources=sources,
                   scenarios=[dict(id=s,weighting=w,label=l,ml_threshold=t,baseq=q) for s,w,l,t,q in SCENARIOS],
                   regions=[dict(id=n,start0=a,end0=b) for n,a,b in REGIONS],rows=len(rows),profiles=len(profiles),
                   raw_definition='Fraction of unchanged reference A/T positions, excluding mutation +/-2 bp, with a retained m6A call. All selected positions must be callable in a molecule; otherwise the endpoint is missing.',
                   strand_definition='Reference A = adenine on genomic forward strand; reference T = adenine on genomic reverse strand. Alignment strand is a separate sensitivity split.',
                   missingness='Per-barcode outcomes average available molecules, then barcodes equally. Missing modification channels and skipped question-mark sites are never zeros.',
                   sampling='WT primary: fixed 10,000-barcode sample/all its eligible molecules. Genotype-only WT: independent deterministic 1/32 read sample. All candidate molecules retained.',
                   inference='Post hoc descriptive sensitivity analyses, no new significance or confidence intervals. WT subset distributions describe sampling scale only.'))
        check_sources(sources)
    print(f'Complete: {len(rows)} comparisons, {len(profiles)} profile rows',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--extracted',default='results/rep3_raw_candidates_20261020/extracted');p.add_argument('--output',required=True)
    a=p.parse_args();run(a.extracted,a.output)

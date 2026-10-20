"""Barcode-level descriptive sensitivity checks for frozen footprint profiles."""
import csv
import hashlib
from pathlib import Path
import numpy as np
from .io import open_text,read_json,source_info
from .processing import check_sources


def finite(value):
    return float(value) if np.isfinite(value) else None


def shared_positions(sequence,start,end,mutation,radius=0):
    # Same reference bases in both groups; exclude the altered base even for A>T.
    return [p for p in range(start,end) if sequence[p] in 'AT' and (mutation is None or abs(p-mutation)>radius)]


def standardized(values,indices,feature,reference):
    """Equal-barcode WT-quartile standardization, requiring support in all bins."""
    edges,weights,means=reference
    vb=np.searchsorted(edges,feature[indices],side='right')
    nbins=len(edges)+1
    valid=np.isfinite(feature[indices])
    vn=np.array([np.sum((vb==b)&valid) for b in range(nbins)])
    if (vn<10).any():return None,int(vn.min()),edges.tolist()
    delta=sum(weights[b]*(values[indices[(vb==b)&valid]].mean()-means[b]) for b in range(nbins))
    return 100*float(delta),int(vn.min()),edges.tolist()


def run(bundle,aggregates,heterogeneity):
    agg=Path(aggregates);heterogeneity=Path(heterogeneity)
    paths=[agg/'audit.json',agg/'barcode_metadata.tsv.gz',agg/'barcode_counts.npy',agg/'barcode_molecules.npy',
           heterogeneity,heterogeneity.parent/'summary.json',Path(__file__)]
    sources=[source_info(p,p.suffix!='.npy') for p in paths]
    aa=read_json(agg/'audit.json')
    check_sources(aa['sources'])
    check_sources(read_json(heterogeneity.parent/'summary.json')['sources'])
    if aa['counts']['barcodes']!=568354:raise ValueError('Unexpected barcode cohort')
    variants=bundle['variants'];byid={v['variant_id']:i for i,v in enumerate(variants)}
    counts=np.load(agg/'barcode_counts.npy',mmap_mode='r');depth=np.load(agg/'barcode_molecules.npy',mmap_mode='r')
    n=len(depth);groups=np.empty(n,dtype=int);halves=np.empty(n,dtype=np.int8);barcodes={}
    with open_text(agg/'barcode_metadata.tsv.gz') as h:
        for i,r in enumerate(csv.DictReader(h,delimiter='\t')):
            if int(r['barcode_index'])!=i or int(r['molecules'])!=depth[i] or r['barcode'] in barcodes:raise ValueError('Barcode metadata differs')
            groups[i]=byid[r['variant_id']];barcodes[r['barcode']]=i
            halves[i]=hashlib.sha256(r['barcode'].encode()).digest()[0]%2
    if len(barcodes)!=n:raise ValueError('Incomplete barcode metadata')
    endpoints=[dict(r) for r in bundle['regions'] if r['kind']=='WT_peak_neighborhood']
    endpoints += [dict(id='adjacent_143_150',start0=122,end0=130),dict(id='whole_tf',start0=0,end0=211),dict(id='whole_msp',start0=0,end0=211)]
    values=np.zeros((n,len(endpoints)))
    for a in range(0,n,4096):
        b=min(a+4096,n);block=np.asarray(counts[a:b],dtype=float)/depth[a:b,None,None]
        for j,r in enumerate(endpoints):values[a:b,j]=block[:,int(r['id']=='whole_msp'),r['start0']:r['end0']].mean(1)
    # The heterogeneity table was built from precisely this supported-barcode cohort.
    features=np.full((n,2),np.nan);seen=np.zeros(n,dtype=bool)
    with open_text(heterogeneity) as h:
        for r in csv.DictReader(h,delimiter='\t'):
            i=barcodes[r['corrected_barcode']]
            if seen[i] or int(r['molecules'])!=depth[i] or r['variant_id']!=variants[groups[i]]['variant_id']:raise ValueError('Heterogeneity cohort differs')
            seen[i]=True
            for j,name in enumerate(['distal_m6a_fraction','distal_nucleosome_fraction']):
                if int(r[name+'_valid_molecules'])==depth[i]:features[i,j]=float(r[name+'_mean'])
    if not seen.all():raise ValueError('Missing barcode heterogeneity')
    wi=bundle['wt_index'];wtidx=np.flatnonzero(groups==wi);wt=variants[wi]
    profiles={};idxs={i:np.flatnonzero(groups==i) for i in range(len(variants))}
    references={};backgrounds={}
    for j in range(len(endpoints)):
        y=values[:,j];backgrounds[j]=dict(mean=float(y[wtidx].mean()),weighted=float(np.average(y[wtidx],weights=depth[wtidx])),
            halves=[float(y[wtidx[halves[wtidx]==h]].mean()) for h in (0,1)],
            depth={k:float(y[wtidx[depth[wtidx]>=k]].mean()) for k in (3,5)})
        r=endpoints[j];track=7 if r['id']=='whole_msp' else 0
        if abs(backgrounds[j]['mean']-np.mean(wt['coverage'][track][r['start0']:r['end0']]))>5.1e-8:
            raise ValueError('WT barcode endpoint differs from atlas')
        for fi in range(2):
            f=features[wtidx,fi];valid=np.isfinite(f);edges=np.unique(np.quantile(f[valid],[.25,.5,.75]))
            bins=np.searchsorted(edges,f,side='right')
            masks=[(bins==b)&valid for b in range(len(edges)+1)]
            weights=np.array([mask.sum() for mask in masks],dtype=float);weights/=weights.sum()
            references[j,fi]=(edges,weights,[float(y[wtidx[mask]].mean()) for mask in masks])
    rows=[];max_reconcile=0.
    for vi,v in enumerate(variants):
        idx=idxs[vi]
        if len(idx)!=v['n_barcodes'] or int(depth[idx].sum())!=v['n_molecules']:raise ValueError('Atlas barcode counts differ')
        if len(idx):
            # Only the small variant groups need full per-base arrays.
            if vi==wi:profile=np.array(wt['coverage'][0])
            else:profile=(np.asarray(counts[idx,0,:],dtype=float)/depth[idx,None]).mean(0)
            error=float(np.max(np.abs(profile-np.array(v['coverage'][0]))));max_reconcile=max(max_reconcile,error)
            if error>5.1e-8:raise ValueError('Barcode profile differs from atlas rounding')
            profiles[vi]=profile
        if vi==wi:continue
        for j,r in enumerate(endpoints):
            y=values[:,j];baseline=backgrounds[j]['mean'];a,b=r['start0'],r['end0']
            result=dict(variant_id=v['variant_id'],source_variant=v['source_variant'],genomic_variant=v['genomic_variant'],
                        region_id=r['id'],start0=a,end0=b,genomic_start1=bundle['genomic_start1']+a,genomic_end1=bundle['genomic_start1']+b-1,
                        mutation_in_region=a<=v['relative_pos0']<b,n_barcodes=len(idx),n_molecules=v['n_molecules'],
                        wt_coverage=baseline,variant_coverage=finite(y[idx].mean()) if len(idx) else None,
                        delta_pp=100*float(y[idx].mean()-baseline) if len(idx) else None,
                        nonzero_barcodes=int((y[idx]>0).sum()),at_change=v['at_change'])
            result['molecule_weighted_delta_pp']=100*float(np.average(y[idx],weights=depth[idx])-backgrounds[j]['weighted']) if len(idx) else None
            result['max_leave_one_barcode_change_pp']=100*float(np.max(np.abs(y[idx]-y[idx].mean()))/(len(idx)-1)) if len(idx)>1 else None
            for half in (0,1):
                ii=idx[halves[idx]==half]
                result[f'split_{half}_barcodes']=len(ii)
                result[f'split_{half}_delta_pp']=100*float(y[ii].mean()-backgrounds[j]['halves'][half]) if len(ii) else None
            for minimum in (3,5):
                ii=idx[depth[idx]>=minimum]
                result[f'depth_ge{minimum}_barcodes']=len(ii)
                result[f'depth_ge{minimum}_delta_pp']=100*float(y[ii].mean()-backgrounds[j]['depth'][minimum]) if len(ii) else None
            for fi,name in enumerate(['distal_m6a','distal_nuc']):
                delta,support,edges=standardized(y,idx,features[:,fi],references[j,fi]) if len(idx) else (None,0,[])
                result[name+'_standardized_delta_pp']=delta;result[name+'_min_bin_barcodes']=support
            for radius in (0,2):
                positions=shared_positions(bundle['sequence'],a,b,v['relative_pos0'],radius)
                result[f'shared_at_radius{radius}_positions']=len(positions)
                track=int(r['id']=='whole_msp')
                if len(idx) and positions:
                    variant_profile=profiles[vi] if track==0 else np.array(v['coverage'][7])
                    result[f'shared_at_radius{radius}_delta_pp']=100*float((variant_profile[positions]-np.array(wt['coverage'][track*7])[positions]).mean())
                else:result[f'shared_at_radius{radius}_delta_pp']=None
            rows.append(result)
    check_sources(sources)
    return rows,dict(barcodes=n,molecules=int(depth.sum()),max_profile_rounding_error=max_reconcile,
                    endpoints=endpoints,split='SHA256 barcode first byte parity; one deterministic split, not independent replication',
                    standardization='WT barcode-mean distal-feature quartiles, >=10 variant barcodes in every bin; sensitivity only, possible mediators',
                    shared_at='Unchanged WT A/T bases excluding mutation +/- radius; uses existing calls, not a caller-bias correction',sources=sources)

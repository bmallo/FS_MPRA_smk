"""Bounded raw MM/ML extraction for predefined replicate-3 candidates.

Never calls pysam.modified_bases: mixed modification tags crash that accessor
in this environment. Scores here are retained input calls, not raw kinetics.
"""
import argparse
import csv
import hashlib
import json
import time
from collections import Counter
from pathlib import Path
import numpy as np
from .footprint_join import region_blocks,region_calls
from .footprint_quality import select_barcodes
from .heterogeneity import MM_HEADER,decode_m6a
from .io import atomic_output,open_text,read_json,source_info,write_json,write_table
from .processing import check_sources

FIELDS=['read_id','observed_variant_id','corrected_barcode','assigned_variant_id','barcode_assignment_status',
        'passes','read_quality','mapq','read_length','strand','plasmid_class','circular','analysis_eligible',
        'footprint_analysis_eligible','footprint_analysis_eligible_with_barcode','m6a_status','footprint_status',
        'region_complete','whole_m6a_fraction','distal_m6a_fraction','distal_nucleosome_fraction',
        'element_tf_small_fraction','element_msp_fraction','element_nucleosome_fraction']


def sampled_read(read_id):
    return int.from_bytes(hashlib.sha256(('raw-audit-20261020|'+read_id).encode()).digest()[:8],'big')%32==0


def decode_scores(read):
    """Return scores and known opportunities in stored-query orientation.

    Missing channels and skipped '?' sites remain unavailable. Unlisted '.'
    positions mean no retained call, not a measured probability of zero.
    """
    forward=read.get_forward_sequence();n=len(forward)
    scores=np.zeros(n,dtype=np.uint8);known=np.zeros(n,dtype=bool)
    if not read.has_tag('MM') or not read.has_tag('ML'):return scores,known
    if read.has_tag('MN') and read.get_tag('MN')!=n:raise ValueError('MN length mismatch')
    ml=read.get_tag('ML')
    if getattr(ml,'typecode',None)!='B':raise ValueError('ML must be bytes')
    cursor=0
    for section in read.get_tag('MM').split(';'):
        if not section:continue
        header,*tokens=section.split(',');match=MM_HEADER.fullmatch(header)
        if not match:raise ValueError('Unsupported MM header')
        base,strand,codes,omitted=match.groups();codes=[codes] if codes.isdigit() else list(codes)
        skips=np.array([int(x) for x in tokens],dtype=int)
        if (skips<0).any():raise ValueError('Negative MM skip')
        positions=np.arange(n) if base=='N' else np.flatnonzero(np.frombuffer(forward.encode(),dtype=np.uint8)==ord(base))
        idx=np.cumsum(skips+1)-1
        if len(idx) and idx[-1]>=len(positions):raise ValueError('MM positions outside sequence')
        end=cursor+len(tokens)*len(codes)
        if end>len(ml):raise ValueError('Truncated ML')
        values=np.asarray(ml[cursor:end],dtype=np.uint8).reshape(len(tokens),len(codes))
        for j,code in enumerate(codes):
            if (base,strand,code) not in [('A','+','a'),('T','-','a')]:continue
            available=positions[idx] if omitted=='?' else positions
            called=positions[idx]
            if read.is_reverse:available=n-1-available;called=n-1-called
            known[available]=True
            np.maximum.at(scores,called,values[:,j])
        cursor=end
    if cursor!=len(ml):raise ValueError('Unused ML values')
    return scores,known


def project_read(read,expected,start=1234,end=1445):
    blocks=region_blocks(read,start,end);qmap=np.full(end-start,-1,dtype=int)
    for qa,qb,ra in blocks:
        lo,hi=ra-start,ra-start+qb-qa
        if (qmap[lo:hi]>=0).any():raise ValueError('Repeated reference bases')
        qmap[lo:hi]=np.arange(qa,qb)
    if (qmap<0).any():raise ValueError('Selected molecule does not span element')
    sequence=''.join(read.query_sequence[q] for q in qmap)
    if sequence!=expected:raise ValueError('Observed genotype differs from BAM sequence')
    scores,known=decode_scores(read)
    status,_,m128,_=decode_m6a(read)
    if not np.array_equal(np.flatnonzero(scores>=128),m128):raise ValueError('Independent m6A decoders disagree')
    coverage=np.zeros(end-start,dtype=np.uint8);edges=np.zeros(end-start,dtype=np.uint8)
    if read.has_tag('MA'):
        calls,_=region_calls(read,blocks)
        for track,items in calls.items():
            for call in items:
                length=call['query_length']
                flag=1 if track=='tf' and 10<=length<=60 else 2 if track=='msp' else 4 if track=='nuc' and 85<=length<=250 else 8 if track=='nuc' and length>250 else 0
                for lo,hi in call['ref_intervals0']:
                    coverage[lo-start:hi-start]|=flag
                    if flag==1:edges[lo-start]=1;edges[hi-start-1]=1
    qualities=np.array(read.query_qualities,dtype=np.uint8)[qmap]
    return scores[qmap],known[qmap],qualities,coverage,edges,int(read.has_tag('MA'))


def run(output):
    import pysam
    started=time.monotonic();cfg=read_json('config/ldlr_footprint_quality.json')
    feature_path=Path(cfg['features']);fa=read_json(cfg['feature_audit'])
    check_sources(fa['sources'])
    lib=read_json(cfg['catalog'])['library'];wt=next(v for v in lib['variants'] if v['class']=='WT')
    candidates=[v for v in lib['variants'] if v['class']=='SNV' and (v['oligo_pos1']==125 or 143<=v['oligo_pos1']<=150)]
    if len(candidates)!=27:raise ValueError('Expected 27 predefined substitutions')
    variants=[wt]+candidates;byid={v['variant_id']:i for i,v in enumerate(variants)}
    bcs,bc_counts=select_barcodes(cfg['barcode_metadata'],wt['variant_id'],10000,20261018)
    bam_path=Path(read_json(cfg['footprint_contract'])['bam'])
    sources=[source_info(p,p.suffix not in ('.bam','.gz')) for p in [Path('config/ldlr_footprint_quality.json'),Path(cfg['catalog']),Path(cfg['feature_audit']),feature_path,Path(cfg['barcode_metadata']),Path(cfg['barcode_audit']),Path(cfg['footprint_contract']),bam_path,Path(str(bam_path)+'.bai'),Path(__file__)]]
    sources += [source_info(Path(__file__).with_name(name),True) for name in ['heterogeneity.py','footprint_quality.py','footprint_join.py','io.py','processing.py']]
    frozen=next(s for s in fa['sources'] if s['path']==str(bam_path.resolve()))
    if source_info(bam_path)!=frozen:raise ValueError('Annotation BAM changed')
    selected={};meta=[];retention=Counter();bc_seen=Counter()
    with open_text(feature_path) as h:
        for r in csv.DictReader(h,delimiter='\t'):
            vi=byid.get(r['observed_variant_id'])
            if vi is None:continue
            retention[(variants[vi]['source_variant'],'sequence_qualified')]+=1
            if r['region_complete']!='1' or r['m6a_status']!='ok':continue
            primary=r['footprint_analysis_eligible_with_barcode']=='1'
            in_barcode_sample=vi==0 and r['corrected_barcode'] in bcs and primary
            in_genotype_sample=vi==0 and sampled_read(r['read_id'])
            if vi==0 and not (in_barcode_sample or in_genotype_sample):continue
            if r['read_id'] in selected:raise ValueError('Duplicate selected read')
            item={k:r[k] for k in FIELDS}
            item.update(variant_index=vi,source_variant=variants[vi]['source_variant'],wt_barcode_sample=int(in_barcode_sample),wt_genotype_sample=int(in_genotype_sample))
            selected[r['read_id']]=len(meta);meta.append(item)
            if in_barcode_sample:bc_seen[r['corrected_barcode']]+=1
            retention[(variants[vi]['source_variant'],'extracted')]+=1
            if len(meta)>300000:raise ValueError('Bounded extraction exceeds 300,000 molecules')
    if dict(bc_seen)!={k:v['eligible_depth'] for k,v in bcs.items()}:raise ValueError('WT barcode sample depth mismatch')
    for vi,v in enumerate(variants):
        if vi==0:continue
        expected=next(x for x in read_json('results/rep3_landscape_20261018/atlas_genomic_20261019/atlas.json')['variants'] if x['variant_id']==v['variant_id'])
        count=sum(r['variant_index']==vi and r['footprint_analysis_eligible_with_barcode']=='1' for r in meta)
        if count!=expected['n_molecules']:raise ValueError('Candidate primary cohort incomplete')
    n=len(meta);print(f'Selected {n:,} molecules; WT barcode sample {sum(bc_seen.values()):,}; scanning BAM',flush=True)
    arrays=dict(scores=np.zeros((n,211),dtype=np.uint8),known=np.zeros((n,211),dtype=bool),
                qualities=np.zeros((n,211),dtype=np.uint8),coverage=np.zeros((n,211),dtype=np.uint8),
                edges=np.zeros((n,211),dtype=np.uint8),footprints_present=np.zeros(n,dtype=np.uint8))
    seen=np.zeros(n,dtype=bool);total=0;hist=np.zeros(256,dtype=np.int64)
    with pysam.AlignmentFile(str(bam_path),'rb',threads=1) as bam:
        for read in bam.fetch(lib['contig']):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:continue
            total+=1;i=selected.get(read.query_name)
            if i is not None:
                if seen[i]:raise ValueError('Duplicate primary BAM read')
                seen[i]=True;r=meta[i];vi=r['variant_index']
                if read.query_length!=int(r['read_length']) or ('-' if read.is_reverse else '+')!=r['strand']:raise ValueError('BAM/features identity differs')
                values=project_read(read,variants[vi]['sequence'])
                for key,value in zip(arrays,values):arrays[key][i]=value
                hist+=np.bincount(values[0][values[0]>0],minlength=256)
                if r['footprint_status']=='evaluated':
                    for flag,field in [(1,'element_tf_small_fraction'),(2,'element_msp_fraction')]:
                        if abs(np.mean((arrays['coverage'][i]&flag)>0)-float(r[field]))>1e-10:raise ValueError('Footprint re-projection mismatch')
            if total%1000000==0:print(f'BAM {total:,}; matched {int(seen.sum()):,}/{n:,}',flush=True)
    if not seen.all() or total!=fa['counts']['primary_rows']:raise ValueError('BAM scan incomplete')
    check_sources(sources)
    with atomic_output(output) as out:
        np.savez_compressed(out/'molecules.npz',**arrays)
        write_table(out/'molecules.tsv.gz',meta,list(meta[0]))
        write_json(out/'variants.json',variants)
        write_table(out/'retention.tsv',[dict(source_variant=v,stage=s,molecules=c) for (v,s),c in retention.items()],['source_variant','stage','molecules'])
        write_json(out/'audit.json',dict(schema_version='fsmpra.raw_candidate_extract.v1',complete=True,sources=sources,
                   molecules=n,bam_primary_rows=total,wt_barcode_selection=bc_counts,retained_positive_ml_histogram=hist.tolist(),seconds=time.monotonic()-started,
                   selection='All sequence-qualified exact candidate genotypes spanning the element with known m6A status. WT: union of 10,000 barcode sample and independent SHA256 read sample (1/32); analysis keeps the two selections separate.',
                   limitations='Available retained MM/ML calls, not original kinetics/probabilities. No barcode or genotype re-calling. Sequence qualification (passes>=10,rq>=.99,mapq>=20) remains required.'))
        check_sources(sources)
    print(f'Complete: {n:,} molecules in {time.monotonic()-started:.1f}s',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();run(a.output)

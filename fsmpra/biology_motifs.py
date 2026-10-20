"""Versioned, descriptive PWM annotations; no inference of protein binding."""
import re
import numpy as np

BASES = 'ACGT'


def read_matrices(path, database):
    records=[]
    with open(path) as h:
        lines=[s.strip() for s in h if s.strip()]
    i=0
    while i<len(lines):
        if not lines[i].startswith('>'):
            raise ValueError('Expected motif header')
        header=lines[i][1:].split(None,1);mid=header[0];i+=1;rows=[]
        while i<len(lines) and not lines[i].startswith('>'):
            rows.append(lines[i]);i+=1
        if database=='JASPAR2024':
            if len(rows)!=4 or [r[0] for r in rows]!=list(BASES):
                raise ValueError('Invalid JASPAR base order')
            matrix=np.array([[float(x) for x in re.sub(r'[\[\]]','',r[1:]).split()] for r in rows]).T
            name=header[1] if len(header)>1 else mid
        else:
            matrix=np.array([[float(x) for x in r.split()] for r in rows]);name=mid.split('.')[0]
        if matrix.ndim!=2 or matrix.shape[1]!=4 or not np.isfinite(matrix).all() or (matrix<0).any() or (matrix.sum(1)<=0).any():
            raise ValueError('Invalid motif counts')
        records.append(dict(database=database,matrix_id=mid,name=name,counts=matrix))
    if len({r['matrix_id'] for r in records})!=len(records):raise ValueError('Duplicate motif IDs')
    return records


def pwm(counts,background):
    # One total pseudocount distributed according to the stated background.
    probabilities=(counts+background)/(counts.sum(1,keepdims=True)+1.)
    return np.log2(probabilities/background)


def score_distribution(weights,background,scale=100):
    """Exact tail probabilities for scores rounded to 0.01-bit contributions.

    This is a transparent local scanner, not FIMO. Returned probabilities are
    per-window sequence-match tails under an independent-base background.
    """
    quantized=np.rint(weights*scale).astype(int)
    shifts=quantized.min(1);q=quantized-shifts[:,None]
    distribution=np.array([1.])
    for row in q:
        nxt=np.zeros(len(distribution)+int(row.max()))
        for base,shift in enumerate(row):
            nxt[shift:shift+len(distribution)]+=background[base]*distribution
        distribution=nxt
    if not np.isclose(distribution.sum(),1.,atol=1e-10):raise ValueError('PWM probability mass lost')
    return quantized, int(shifts.sum()), np.cumsum(distribution[::-1])[::-1]


def scan(sequence,matrices,origin,threshold=1e-4):
    seq=np.array([BASES.index(b) for b in sequence]);uniform=np.full(4,.25)
    gc=(sequence.count('G')+sequence.count('C'))/len(sequence)
    matched=np.array([(1-gc)/2,gc/2,gc/2,(1-gc)/2])
    hits=[];tested=0
    for m in matrices:
        length=len(m['counts'])
        if length>len(seq):continue
        w=pwm(m['counts'],uniform)
        q,offset,tail=score_distribution(w,uniform)
        windows=np.lib.stride_tricks.sliding_window_view(seq,length)
        local=[]
        for strand,window in [('+',windows),('-',3-windows[:,::-1])]:
            totals=q[np.arange(length),window].sum(1)
            pvalues=tail[totals-offset];tested+=len(window)
            for start in np.flatnonzero(pvalues<=threshold):
                fragment=window[start];score=float(w[np.arange(length),fragment].sum())
                local.append((int(start),strand,fragment,score,float(pvalues[start])))
        if not local:continue
        wg=pwm(m['counts'],matched);qg,og,tg=score_distribution(wg,matched)
        for start,strand,fragment,score,pvalue in local:
            # Store genomic-forward matrices so variant scores use genomic ALT
            # directly even for reverse-strand motif hits.
            oriented=w if strand=='+' else w[::-1,::-1]
            pg=float(tg[int(qg[np.arange(length),fragment].sum())-og])
            hit=dict(hit_id=f"{m['database']}:{m['matrix_id']}:{start}:{strand}",
                     database=m['database'],matrix_id=m['matrix_id'],name=m['name'],
                     start0=start,end0=start+length,strand=strand,
                     genomic_start1=origin+start,genomic_end1=origin+start+length-1,
                     sequence=sequence[start:start+length],motif_oriented_sequence=''.join(BASES[b] for b in fragment),
                     score_bits=score,sequence_tail_p=pvalue,gc_matched_tail_p=pg,
                     passes_both_backgrounds=pg<=threshold,weights_genomic=oriented.tolist())
            hits.append(hit)
    hits.sort(key=lambda h:(h['sequence_tail_p'],h['database'],h['matrix_id'],h['start0'],h['strand']))
    return hits,dict(matrices=len(matrices),tested_windows_both_strands=tested,hit_count=len(hits),
                     threshold=threshold,background_uniform=uniform.tolist(),background_element_gc=matched.tolist(),
                     pseudocount_total=1.,score_grid_bits=.01,
                     interpretation='Unadjusted per-window sequence-match tails, not binding or variant-effect p-values. WT sites only; no de novo variant-created sites scanned.')


def variant_delta(hit,position,ref,alt):
    if position is None or not hit['start0']<=position<hit['end0']:return None
    weights=hit['weights_genomic'][position-hit['start0']]
    return float(weights[BASES.index(alt)]-weights[BASES.index(ref)])

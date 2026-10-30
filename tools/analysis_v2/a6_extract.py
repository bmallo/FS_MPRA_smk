"""Step A6 (part 1): per-molecule m6A calls at every reference A/T for WT molecules.

Same read filters and subsample as A5 (a5_m6a_extract.py), plus >= 10 passes and WT only (the --wt-region matches
the reference exactly). Without --wt-region every read counts (the naked pEF-GFP control). Per molecule and per
reference A/T: -1 = not observed (read != reference or not aligned), 0 = no call, 1 = call with ML >= 128.
Distal m6A fraction as in A5 (outside element +/-100 bp and barcode +/-25 bp; the whole plasmid without them).

Usage: python tools/analysis_v2/a6_extract.py --bam <bam> --contig <name> --reference <fa> \
           [--element 1234:1445 --barcode 3294:3309 --wt-region 1234:1445] --fraction 0.12 --output <new dir>
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
import pysam

_spec = importlib.util.spec_from_file_location('a5x', Path(__file__).with_name('a5_m6a_extract.py'))
a5x = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a5x)
MIN_PASSES = 10
ML_MIN = a5x.THRESHOLDS[0]


def molecule(rec, ref_arr, at_idx, L, dmask, wt_region):
    """(status, record): record holds the call vector over at_idx (-1/0/1) and per-molecule values."""
    if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
        return 'not_primary', None
    if rec.mapping_quality < a5x.MIN_MAPQ:
        return 'mapq', None
    if not a5x.full_length(rec, L):
        return 'not_full_length', None
    if rec.get_tag('rq') < a5x.MIN_RQ:
        return 'read_quality', None
    if int(rec.get_tag('np')) < MIN_PASSES:
        return 'low_passes', None
    ml = a5x.ml_array(rec)
    if ml is None:
        return 'no_m6a_tags', None
    pairs = rec.get_aligned_pairs(matches_only=True)
    q = np.fromiter((p[0] for p in pairs), np.int64, count=len(pairs))
    r = np.fromiter((p[1] for p in pairs), np.int64, count=len(pairs))
    qseq = np.frombuffer(rec.query_sequence.upper().encode(), np.uint8)
    if wt_region is not None:
        quals = np.asarray(rec.query_qualities) if rec.query_qualities is not None else None
        if not a5x.wt_exact(q, r, qseq, quals, ref_arr, wt_region):
            return 'not_wt', None
    ok = qseq[q] == ref_arr[r]
    full = np.full(L, -1, np.int8)
    full[r[ok]] = (ml[q[ok]] >= ML_MIN).astype(np.int8)
    vec = full[at_idx]
    d = vec[dmask[at_idx]]
    distal_at = int((d >= 0).sum())
    if distal_at < a5x.MIN_DISTAL_AT:
        return 'few_distal_at', None
    return 'ok', dict(vec=vec, distal_frac=float((d == 1).sum()) / distal_at, passes=int(rec.get_tag('np')),
                      strand=1 if rec.is_reverse else 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--contig', required=True)
    ap.add_argument('--reference', required=True)
    ap.add_argument('--element', help='start0:end0 (half-open)')
    ap.add_argument('--barcode', help='start0:end0 (half-open)')
    ap.add_argument('--wt-region', help='start0:end0 that must match the reference; omit to keep every read')
    ap.add_argument('--fraction', type=float, required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    span = lambda s: tuple(int(x) for x in s.split(':')) if s else None
    element, barcode, wt_region = span(a.element), span(a.barcode), span(a.wt_region)
    ref = a5x.ref_seq(a.reference)
    L = len(ref)
    ref_arr = np.frombuffer(ref.encode(), np.uint8)
    at_idx = np.flatnonzero((ref_arr == ord('A')) | (ref_arr == ord('T')))
    dmask = a5x.distal_mask(L, element, barcode) if element and barcode else np.ones(L, bool)
    status, vecs, rows, seen, t0 = {}, [], [], 0, time.time()
    with pysam.AlignmentFile(a.bam) as f:
        for rec in f.fetch(a.contig):
            if rec.is_secondary or rec.is_supplementary or not a5x.keep_read(rec.query_name, a.fraction):
                continue
            seen += 1
            st, x = molecule(rec, ref_arr, at_idx, L, dmask, wt_region)
            status[st] = status.get(st, 0) + 1
            if x is None:
                continue
            vecs.append(x['vec'])
            rows.append((x['distal_frac'], x['passes'], x['strand']))
            if seen % 100000 == 0:
                print(f'{seen} sampled reads, {len(rows)} kept, {time.time() - t0:.0f} s', flush=True)
    X = np.array(vecs, np.int8).reshape(len(vecs), at_idx.size)
    rows = np.array(rows, float).reshape(len(rows), 3)
    np.savez_compressed(out / 'molecules.npz', pos=at_idx, observed=np.packbits(X >= 0, axis=1),
                        methylated=np.packbits(X == 1, axis=1), distal_frac=rows[:, 0],
                        passes=rows[:, 1].astype(np.int32), strand=rows[:, 2].astype(np.int8))
    json.dump(dict(bam=a.bam, contig=a.contig, reference=a.reference, reference_length=L, element=element,
                   barcode=barcode, wt_region=wt_region, fraction=a.fraction, sampled_primary=seen, status=status,
                   kept=len(rows), n_at=int(at_idx.size), min_passes=MIN_PASSES, ml_min=ML_MIN,
                   seconds=round(time.time() - t0)), open(out / 'extract.json', 'w'), indent=1)
    print(json.dumps(status), len(rows))


def load(path):
    """(pos, X int8 molecules x positions with -1/0/1, distal_frac, passes, strand) from molecules.npz."""
    z = np.load(Path(path) / 'molecules.npz')
    P = z['pos'].size
    obs = np.unpackbits(z['observed'], axis=1, count=P).astype(bool)
    meth = np.unpackbits(z['methylated'], axis=1, count=P).astype(np.int8)
    X = np.where(obs, meth, np.int8(-1)).astype(np.int8)
    return z['pos'], X, z['distal_frac'], z['passes'], z['strand']


if __name__ == '__main__':
    main()

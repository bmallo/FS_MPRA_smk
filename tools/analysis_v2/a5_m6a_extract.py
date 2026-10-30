"""Step A5 (part 1): per-run m6A call tallies for methylation-call QC.

Streams one aligned BAM once and keeps a deterministic subsample of primary plasmid reads (CRC32 of the read
name). Read filters follow rules.json (MAPQ >= 20, full length, rq >= 0.99) except passes, which is kept as a
dimension because A5 tests m6A reliability by pass count. Per read:
  observed m6A  = a listed A+a / T-a call with ML >= threshold at an aligned reference A/T where read == reference
                  (the plasmid-LLR convention; MM '.' mode, so unlisted adenines are unmethylated)
  distal m6A    = that fraction outside the element +/-100 bp and the barcode +/-25 bp (the rules.json classes:
                  U <= 0.01, M >= 0.60, C in between; ML >= 128)
  WT            = the --wt-region matches the reference exactly (no mismatch or indel, every base >= Q20)
Counts per reference position are kept by class x alignment strand x pass bin, at ML >= 128 and >= 204.
Positions inside --wt-region are counted only from WT reads; all others from every read. The ML histogram covers
every listed call at counted positions, with the number of counted positions (unlisted = no call).

Usage: python tools/analysis_v2/a5_m6a_extract.py --bam <bam> --contig <name> --reference <fa> \
           --element 1234:1445 --barcode 3294:3309 --wt-region 1234:1445 --fraction 0.12 --output <new dir>
"""
import argparse
import gzip
import json
import time
import zlib
from pathlib import Path

import numpy as np
import pysam

CLASSES = ['U', 'C', 'M']
PASS_BINS = [0, 10, 15, 20, 30]          # bin k = [PASS_BINS[k], PASS_BINS[k+1]); last is 30+
PASS_LABELS = ['<10', '10-14', '15-19', '20-29', '30+']
THRESHOLDS = [128, 204]
U_MAX, M_MIN, MIN_DISTAL_AT = 0.01, 0.60, 500
MIN_MAPQ, MIN_RQ, MIN_BQ, END_SLACK, MAX_CLIP = 20, 0.99, 20, 50, 100


def ref_seq(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>')).upper()


def distal_mask(L, element, barcode):
    m = np.ones(L, bool)
    m[max(0, element[0] - 100):min(L, element[1] + 100)] = False
    m[max(0, barcode[0] - 25):min(L, barcode[1] + 25)] = False
    return m


def pass_bin(p):
    return int(np.searchsorted(PASS_BINS, p, side='right') - 1)


def klass(frac):
    return 0 if frac <= U_MAX else 2 if frac >= M_MIN else 1


def keep_read(name, fraction):
    return zlib.crc32(name.encode()) % 1_000_000 < fraction * 1_000_000


def full_length(rec, L):
    c = rec.cigartuples
    clip5 = c[0][1] if c[0][0] == 4 else 0
    clip3 = c[-1][1] if c[-1][0] == 4 else 0
    return (rec.reference_start <= END_SLACK and rec.reference_end >= L - END_SLACK and max(clip5, clip3) <= MAX_CLIP
            and rec.query_length <= 1.1 * L)


def wt_exact(q, r, qseq, quals, ref_arr, region):
    """Region [s, e) aligned base-for-base with no indel, every base equal to the reference and >= Q20."""
    s, e = region
    i, j = np.searchsorted(r, [s, e])
    if j - i != e - s or r[i] != s or q[j - 1] - q[i] != e - s - 1:
        return False
    qq = q[i:j]
    return bool(np.all(qseq[qq] == ref_arr[s:e]) and (quals is None or np.all(quals[qq] >= MIN_BQ)))


def ml_array(rec):
    """ML per stored-SEQ position for m6A calls (-1 = unlisted); None if the read has no m6A channel.
    A channel listed in MM with no calls (an unlabeled molecule) is kept: pysam omits it from modified_bases."""
    mods = rec.modified_bases
    keys = [k for k in (mods or {}) if k[2] == 'a' and k[0] in 'AT']
    out = np.full(rec.query_length, -1, np.int16)
    if not keys:
        mm = rec.get_tag('MM') if rec.has_tag('MM') else rec.get_tag('Mm') if rec.has_tag('Mm') else ''
        return out if any(p.startswith(('A+a', 'T-a')) for p in mm.split(';')) else None
    for k in keys:
        for pos, qual in mods[k]:
            if 0 <= pos < out.size:
                out[pos] = max(out[pos], qual)
    return out


def process(rec, ref_arr, at_mask, L, dmask, wt_region):
    """(status, record) for one alignment; record holds the per-read values and counted positions."""
    if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
        return 'not_primary', None
    if rec.mapping_quality < MIN_MAPQ:
        return 'mapq', None
    if not full_length(rec, L):
        return 'not_full_length', None
    if rec.get_tag('rq') < MIN_RQ:
        return 'read_quality', None
    ml = ml_array(rec)
    if ml is None:
        return 'no_m6a_tags', None
    pairs = rec.get_aligned_pairs(matches_only=True)
    q = np.fromiter((p[0] for p in pairs), np.int64, count=len(pairs))
    r = np.fromiter((p[1] for p in pairs), np.int64, count=len(pairs))
    qseq = np.frombuffer(rec.query_sequence.upper().encode(), np.uint8)
    quals = np.asarray(rec.query_qualities) if rec.query_qualities is not None else None
    wt = wt_exact(q, r, qseq, quals, ref_arr, wt_region)
    keep = at_mask[r] & (qseq[q] == ref_arr[r])
    if not wt:
        keep &= ~((r >= wt_region[0]) & (r < wt_region[1]))
    pos, mlv = r[keep], ml[q[keep]]
    d = dmask[pos]
    distal_at = int(d.sum())
    if distal_at < MIN_DISTAL_AT:
        return 'few_distal_at', None
    distal_frac = float((mlv[d] >= THRESHOLDS[0]).sum()) / distal_at
    return 'ok', dict(pos=pos, ml=mlv, wt=wt, distal_at=distal_at, distal_frac=distal_frac,
                      passes=int(rec.get_tag('np')), rq=float(rec.get_tag('rq')), strand=1 if rec.is_reverse else 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--contig', required=True)
    ap.add_argument('--reference', required=True)
    ap.add_argument('--element', required=True, help='start0:end0 (half-open)')
    ap.add_argument('--barcode', required=True, help='start0:end0 (half-open)')
    ap.add_argument('--wt-region', required=True, help='start0:end0 that must match the reference for WT')
    ap.add_argument('--fraction', type=float, required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    span = lambda s: tuple(int(x) for x in s.split(':'))
    element, barcode, wt_region = span(a.element), span(a.barcode), span(a.wt_region)
    ref = ref_seq(a.reference)
    L = len(ref)
    ref_arr = np.frombuffer(ref.encode(), np.uint8)
    at_mask = (ref_arr == ord('A')) | (ref_arr == ord('T'))
    dmask = distal_mask(L, element, barcode)
    counts = np.zeros((3, 2, len(PASS_LABELS), 1 + len(THRESHOLDS), L), np.int64)   # class, strand, pass bin, obs/m128/m204, pos
    ml_hist = np.zeros((3, len(PASS_LABELS), 256), np.int64)
    n_at = np.zeros((3, len(PASS_LABELS)), np.int64)
    status, rows, seen, t0 = {}, [], 0, time.time()
    el = (np.arange(L) >= element[0]) & (np.arange(L) < element[1])
    with pysam.AlignmentFile(a.bam) as f:
        for rec in f.fetch(a.contig):
            if rec.is_secondary or rec.is_supplementary or not keep_read(rec.query_name, a.fraction):
                continue
            seen += 1
            st, x = process(rec, ref_arr, at_mask, L, dmask, wt_region)
            status[st] = status.get(st, 0) + 1
            if x is None:
                continue
            c, b, s = klass(x['distal_frac']), pass_bin(x['passes']), x['strand']
            counts[c, s, b, 0] += np.bincount(x['pos'], minlength=L)
            for k, t in enumerate(THRESHOLDS):
                counts[c, s, b, 1 + k] += np.bincount(x['pos'][x['ml'] >= t], minlength=L)
            listed = x['ml'][x['ml'] >= 0]
            ml_hist[c, b] += np.bincount(listed, minlength=256)
            n_at[c, b] += x['ml'].size
            ein = el[x['pos']]
            rows.append((x['passes'], round(x['rq'], 5), '-' if s else '+', int(x['wt']), x['distal_at'],
                         round(x['distal_frac'], 5), CLASSES[c], int(ein.sum()),
                         int((x['ml'][ein] >= THRESHOLDS[0]).sum()), int((x['ml'][ein] >= THRESHOLDS[1]).sum())))
            if seen % 100000 == 0:
                print(f'{seen} sampled reads, {len(rows)} kept, {time.time() - t0:.0f} s', flush=True)
    np.savez_compressed(out / 'counts.npz', counts=counts, ml_hist=ml_hist, n_at=n_at,
                        axes=np.array(['class', 'strand', 'pass_bin', 'obs_m128_m204', 'ref_pos0']),
                        classes=np.array(CLASSES), pass_labels=np.array(PASS_LABELS), thresholds=np.array(THRESHOLDS))
    with gzip.open(out / 'reads.tsv.gz', 'wt') as fh:
        fh.write('passes\trq\tstrand\twt\tdistal_at\tdistal_m6a\tclass\telement_at\telement_m6a_128\telement_m6a_204\n')
        for row in rows:
            fh.write('\t'.join(map(str, row)) + '\n')
    json.dump(dict(bam=a.bam, contig=a.contig, reference=a.reference, reference_length=L, element=element,
                   barcode=barcode, wt_region=wt_region, fraction=a.fraction, sampled_primary=seen, status=status,
                   kept=len(rows), filters=dict(min_mapq=MIN_MAPQ, min_rq=MIN_RQ, min_bq_wt=MIN_BQ, end_slack=END_SLACK,
                   max_clip=MAX_CLIP, min_distal_at=MIN_DISTAL_AT, u_max=U_MAX, m_min=M_MIN, class_ml=THRESHOLDS[0]),
                   seconds=round(time.time() - t0)), open(out / 'extract.json', 'w'), indent=1)
    print(json.dumps(status), len(rows))


if __name__ == '__main__':
    main()

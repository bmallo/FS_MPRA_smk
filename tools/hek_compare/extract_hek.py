"""HEK293T (older LDLR library) WT extraction for the plasmid-LLR comparison with HepG2 rep3.

Applies the replicate-3 analysis filters to the existing plasmid-aligned HEK293T BAM (the raw
HiFi BAM is not on disk; genomic reads were never aligned):
  primary, MAPQ >= 20, np >= 10, rq >= 0.99, full length (ends within 50 bp of the plasmid
  ends, soft clips <= 100 bp, query <= 1.1x plasmid), m6A callable (A+a / T-a present).
WT = exact reference match, no indels, every base >= Q20 over the whole 318-bp mutable
promoter of the older library (ref0 3184-3502), which contains the HepG2 211-bp element on
the opposite strand (ref0 3236-3447).
Distal m6A = retained calls with ML >= 128 / aligned A/T, outside the shared element +/-100
bp and the barcode N run +/-25 bp, as in the replicate-3 site-methylation classes. Then
U (<= 0.01), M (>= 0.60) and C read sets are drawn exactly as tools/llr_compare/select_reads.py.

Usage: python tools/hek_compare/extract_hek.py --workdir results/hek293t_llr_20261028 --bam <aligned.bam> --procs N
"""
import argparse
import gzip
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pysam

L = 4718
MUT = (3184, 3502)             # older library mutable promoter, ref0 half-open
ELEM = (3236, 3447)            # HepG2 211-bp element (reverse strand here)
DISTAL = [(0, 1318), (1382, 3136), (3547, L)]   # minus element +/-100 and barcode (1343-1357) +/-25
U_MAX, M_MIN, ML_MIN = 0.01, 0.60, 128
SEED = 20261026
COLS = ['read_id', 'strand', 'passes', 'read_quality', 'mapq', 'status', 'wt', 'whole_at', 'whole_m6a',
        'distal_at', 'distal_m6a', 'element_at', 'element_m6a']


def ref_seq(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>')).upper()


def classify(rec, ref):
    if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
        return None
    row = dict(read_id=rec.query_name, strand='-' if rec.is_reverse else '+', passes=rec.get_tag('np'),
               read_quality=round(rec.get_tag('rq'), 5), mapq=rec.mapping_quality, wt='')
    c = rec.cigartuples
    clip5 = c[0][1] if c[0][0] == 4 else 0
    clip3 = c[-1][1] if c[-1][0] == 4 else 0
    if rec.mapping_quality < 20:
        return dict(row, status='mapq')
    if row['passes'] < 10:
        return dict(row, status='passes')
    if row['read_quality'] < 0.99:
        return dict(row, status='read_quality')
    if (rec.reference_start > 50 or rec.reference_end < L - 50 or max(clip5, clip3) > 100
            or rec.query_length > 1.1 * L):
        return dict(row, status='not_full_length')
    mods = rec.modified_bases
    if not mods or not any(k[0] in 'AT' and k[2] == 'a' for k in mods):
        return dict(row, status='no_m6a_tags')
    pairs = np.array(rec.get_aligned_pairs(matches_only=True), dtype=np.int64)
    q, r = pairs[:, 0], pairs[:, 1]
    seq = rec.query_sequence
    qual = np.asarray(rec.query_qualities)
    # WT over the mutable promoter: every ref base aligned once, no insertion, exact, >= Q20.
    inm = (r >= MUT[0]) & (r < MUT[1])
    qm = q[inm]
    n = MUT[1] - MUT[0]
    wt = (len(qm) == n and qm[-1] - qm[0] == n - 1
          and ''.join(seq[i] for i in qm) == ref[MUT[0]:MUT[1]] and qual[qm].min() >= 20)
    row['wt'] = int(wt)
    qa = np.frombuffer(seq.encode(), dtype=np.uint8)
    at = (qa[q] == ord('A')) | (qa[q] == ord('T'))
    called = np.zeros(len(seq), bool)
    for (base, _, code), calls in mods.items():
        if base in 'AT' and code == 'a':
            for pos, ml in calls:
                if ml >= ML_MIN:
                    called[pos] = True
    m = called[q] & at

    def frac(mask):
        return int(at[mask].sum()), int(m[mask].sum())
    row['whole_at'], row['whole_m6a'] = frac(np.ones(len(q), bool))
    dmask = np.zeros(len(q), bool)
    for lo, hi in DISTAL:
        dmask |= (r >= lo) & (r < hi)
    row['distal_at'], row['distal_m6a'] = frac(dmask)
    row['element_at'], row['element_m6a'] = frac((r >= ELEM[0]) & (r < ELEM[1]))
    return dict(row, status='eligible')


def worker(args):
    bam, refpath, k, nproc, out = args
    ref = ref_seq(refpath)
    with pysam.AlignmentFile(bam, 'rb') as af, gzip.open(out, 'wt') as h:
        for i, rec in enumerate(af.fetch(until_eof=True)):
            if i % nproc != k:
                continue
            row = classify(rec, ref)
            if row is not None:
                h.write('\t'.join(str(row.get(c, '')) for c in COLS) + '\n')
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workdir', required=True)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--n-control', type=int, default=93030)
    ap.add_argument('--n-chrom', type=int, default=60000)
    a = ap.parse_args()
    wd = Path(a.workdir)
    refpath = str(wd / 'ref/LDLR_HEK293T.fa')
    parts = [str(wd / f'data/molecules.part{k}.tsv.gz') for k in range(a.procs)]
    with mp.Pool(a.procs) as pool:
        pool.map(worker, [(a.bam, refpath, k, a.procs, p) for k, p in enumerate(parts)])
    import pandas as pd
    df = pd.concat([pd.read_csv(p, sep='\t', names=COLS, dtype={'wt': 'float'}) for p in parts], ignore_index=True)
    df = df.sort_values('read_id').reset_index(drop=True)
    for c in ('whole', 'distal', 'element'):
        df[f'{c}_m6a_fraction'] = df[f'{c}_m6a'] / df[f'{c}_at'].where(df[f'{c}_at'] > 0)
    df.to_csv(wd / 'data/molecules.tsv.gz', sep='\t', index=False)
    for p in parts:
        Path(p).unlink()

    wt = df[(df.status == 'eligible') & (df.wt == 1)]
    rng = np.random.default_rng(SEED)
    sets = {'wt_unmeth': wt[wt.distal_m6a_fraction <= U_MAX], 'wt_meth': wt[wt.distal_m6a_fraction >= M_MIN],
            'wt_chrom': wt[(wt.distal_m6a_fraction > U_MAX) & (wt.distal_m6a_fraction < M_MIN)]}
    n = {'wt_unmeth': a.n_control, 'wt_meth': a.n_control, 'wt_chrom': a.n_chrom}
    summary = {}
    for name, d in sets.items():
        take = d.read_id.values if len(d) <= n[name] else rng.choice(d.read_id.values, n[name], replace=False)
        with open(wd / 'config' / f'{name}.read_ids.txt', 'w') as h:
            h.write('\n'.join(sorted(take)) + '\n')
        summary[name] = dict(available=int(len(d)), selected=int(len(take)))
    with open(wd / 'config' / 'samples.tsv', 'w') as h:
        h.write('sample\trole\nwt_meth\tcontrol\nwt_unmeth\tcontrol\nwt_chrom\tchromatin\n')
    elig = df[df.status == 'eligible']
    funnel = dict(primary_records=int(len(df)), status=df.status.value_counts().to_dict(),
                  eligible=int(len(elig)), eligible_wt=int(len(wt)),
                  median_whole_m6a_fraction=float(elig.whole_m6a_fraction.median()),
                  wt_class_fractions={k: len(v) / max(len(wt), 1) for k, v in sets.items()})
    json.dump(dict(bam=a.bam, reference=refpath, mutable_ref0=MUT, element_ref0=ELEM, distal_ref0=DISTAL,
                   u_max=U_MAX, m_min=M_MIN, ml_min=ML_MIN, seed=SEED, funnel=funnel, sets=summary),
              open(wd / 'config' / 'selection.json', 'w'), indent=1)
    print(json.dumps(dict(funnel=funnel, sets=summary), indent=1))


if __name__ == '__main__':
    main()

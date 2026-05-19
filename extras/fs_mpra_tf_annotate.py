#!/usr/bin/env python
"""P4.1 — nominate candidate TFs for Stage-3 motif calls.

Post-hoc (reads the Stage-3 HDF5; does NOT touch the pipeline / its
deps). For each cross-variant motif call it scores the motif's
reference DNA (both strands, every offset) against JASPAR + HOCOMOCO
PWMs with a pure-numpy log-odds model, AND cross-checks the match
against the motif's per-position/per-base sensitivity profile — the
MPRA-specific corroboration: for a TRUE TF the positions whose
mutations most reduce occupancy should coincide with the PWM's
high-information-content positions. Emits <sample>_tf_candidates.tsv.

Usage:
  python extras/fs_mpra_tf_annotate.py --h5 results/.../sample.h5 \
      --jaspar extras/data/JASPAR2024_CORE_vertebrates_nr.jaspar \
      --hocomoco extras/data/HOCOMOCOv12_H12CORE_pcm.txt \
      [--top 5] [--gc 0.41] [--pseudocount 0.8]
"""
import argparse
import csv
import os
import sys

import numpy as np
import h5py

BASES = 'ACGT'
BIDX = {b: i for i, b in enumerate(BASES)}
_COMP = str.maketrans('ACGTN', 'TGCAN')


# ----------------------------- PWM I/O --------------------------------

def _pfm_records_jaspar(path):
    """Yield (matrix_id, tf, pfm[L,4]) from a concatenated JASPAR file:
    '>ID NAME' then 4 lines 'A  [ c c c ]' (rows=bases A,C,G,T)."""
    with open(path) as fh:
        lines = [ln.rstrip('\n') for ln in fh]
    i = 0
    while i < len(lines):
        if not lines[i].startswith('>'):
            i += 1
            continue
        hdr = lines[i][1:].split(None, 1)
        mid = hdr[0]
        tf = hdr[1].strip() if len(hdr) > 1 else mid
        rows = []
        for j in range(1, 5):
            seg = lines[i + j]
            seg = seg[seg.index('[') + 1: seg.index(']')] \
                if '[' in seg else seg.split(None, 1)[1]
            rows.append([float(x) for x in seg.split()])
        i += 5
        yield mid, tf, np.asarray(rows, float).T          # [L,4] ACGT


def _pfm_records_hocomoco(path):
    """Yield (matrix_id, tf, pfm[L,4]) from the concatenated HOCOMOCO
    file: '>NAME' then L lines of 4 tab counts (A C G T), repeated."""
    with open(path) as fh:
        lines = [x.rstrip('\n') for x in fh]
    mid, tf, rows = None, None, []
    for ln in lines:
        if ln.startswith('>'):
            if mid is not None and rows:
                yield mid, tf, np.asarray(rows, float)     # [L,4] ACGT
            mid = ln[1:].strip()
            tf = mid.split('.')[0]
            rows = []
        elif ln.strip():
            rows.append([float(v) for v in ln.split()])
    if mid is not None and rows:
        yield mid, tf, np.asarray(rows, float)


def load_pwms(jaspar, hocomoco, pseudocount, bg):
    """Build log2-odds matrices [L,4] + per-position info content from
    count matrices: M = log2((c + pc*bg)/(colN + pc)/bg)."""
    out = []
    srcs = []
    if jaspar and os.path.exists(jaspar):
        srcs.append(('JASPAR', _pfm_records_jaspar(jaspar)))
    if hocomoco and os.path.exists(hocomoco):
        srcs.append(('HOCOMOCO', _pfm_records_hocomoco(hocomoco)))
    bg = np.asarray(bg, float)
    for db, gen in srcs:
        for mid, tf, pfm in gen:
            if pfm.ndim != 2 or pfm.shape[1] != 4 or pfm.shape[0] < 4:
                continue
            n = pfm.sum(1, keepdims=True)
            n[n == 0] = 1.0
            p = (pfm + pseudocount * bg) / (n + pseudocount)
            lod = np.log2(p / bg)
            ic = float(np.sum(p * np.log2(p / bg)) / pfm.shape[0])
            ic_pos = np.sum(p * np.log2(p / bg), axis=1)        # [L]
            out.append({'db': db, 'id': mid, 'tf': tf,
                        'lod': lod.astype(np.float64),
                        'ic_pos': ic_pos.astype(np.float64),
                        'len': pfm.shape[0], 'ic_mean': ic})
    return out


# --------------------------- scoring ----------------------------------

def _seq_idx(seq):
    return np.array([BIDX.get(b, -1) for b in seq.upper()], int)


def _best_placement(lod, sidx):
    """Best (score, offset) of an L-wide PWM over a sequence (one
    strand). Windows containing an N (idx<0) are skipped."""
    L = lod.shape[0]
    n = sidx.size
    best_s, best_o = -1e18, -1
    for o in range(0, n - L + 1):
        w = sidx[o:o + L]
        if (w < 0).any():
            continue
        s = lod[np.arange(L), w].sum()
        if s > best_s:
            best_s, best_o = float(s), o
    return best_s, best_o


def annotate(h5_path, pwms, top, out_tsv):
    with h5py.File(h5_path, 'r') as f:
        if 'motifs' not in f or 'calls' not in f['motifs']:
            print('no motif calls in HDF5 — nothing to annotate')
            open(out_tsv, 'w').close()
            return 0
        calls = f['motifs']['calls']
        motifs = []
        for k in sorted(calls, key=lambda x: int(x[1:])):
            g = calls[k]
            seq = g.attrs['ref_sequence']
            seq = seq.decode() if isinstance(seq, bytes) else seq
            sc = g['sensitivity_count'][:]                  # [L,4]
            sd = g['sensitivity_mean_signed_delta'][:]       # [L,4]
            b = g.attrs['bin']
            motifs.append({
                'idx': int(k[1:]),
                'bin': b.decode() if isinstance(b, bytes) else b,
                'abs_start': int(g.attrs['abs_start']),
                'abs_end': int(g.attrs['abs_end']),
                'seq': seq,
                # per-position observed disruption magnitude:
                # Σ_base count·|signed Δ| (how much mutating that
                # column actually moved occupancy)
                'sens': np.sum(sc * np.abs(sd), axis=1),
            })

    rows = []
    for m in motifs:
        seq = m['seq']
        if not seq or len(seq) < 5:
            continue
        fwd = _seq_idx(seq)
        rev = _seq_idx(seq[::-1].translate(_COMP))
        sens = m['sens']
        sens_z = ((sens - sens.mean()) / sens.std()
                  if sens.std() > 0 else np.zeros_like(sens))
        cands = []
        for w in pwms:
            if w['len'] > len(seq):
                continue
            sf, of = _best_placement(w['lod'], fwd)
            sr, orr = _best_placement(w['lod'], rev)
            if sf < -1e17 and sr < -1e17:
                continue
            if sf >= sr:
                strand, score, off = '+', sf, of
            else:
                strand, score, off = '-', sr, orr
            # max attainable for this PWM (per-col best) -> normalized
            mx = w['lod'].max(1).sum()
            mn = w['lod'].min(1).sum()
            norm = (score - mn) / (mx - mn) if mx > mn else 0.0
            # MPRA corroboration: does PWM info-content track the
            # motif's per-position mutation sensitivity at the match?
            ic = w['ic_pos']
            if strand == '-':
                ic = ic[::-1]
            seg = sens_z[off:off + w['len']]
            if seg.size == w['len'] and seg.std() > 0 \
                    and ic.std() > 0:
                conc = float(np.corrcoef(
                    ic, (seg - seg.mean()) / seg.std())[0, 1])
            else:
                conc = float('nan')
            cands.append((norm, score, conc, strand, off, w))
        # rank: sequence fit first, concordance as corroboration
        cands.sort(key=lambda c: (c[0], (c[2] if c[2] == c[2]
                                         else -1)), reverse=True)
        for norm, score, conc, strand, off, w in cands[:top]:
            rows.append({
                'motif_idx': m['idx'], 'bin': m['bin'],
                'abs_start': m['abs_start'], 'abs_end': m['abs_end'],
                'motif_seq': seq, 'db': w['db'], 'tf': w['tf'],
                'matrix_id': w['id'], 'pwm_len': w['len'],
                'strand': strand, 'offset': off,
                'logodds': round(score, 3),
                'score_norm': round(norm, 4),
                'pwm_ic_mean': round(w['ic_mean'], 3),
                'ic_sensitivity_concordance':
                    ('' if conc != conc else round(conc, 4)),
            })

    hdr = ['motif_idx', 'bin', 'abs_start', 'abs_end', 'motif_seq',
           'db', 'tf', 'matrix_id', 'pwm_len', 'strand', 'offset',
           'logodds', 'score_norm', 'pwm_ic_mean',
           'ic_sensitivity_concordance']
    with open(out_tsv, 'w', newline='') as fo:
        w = csv.DictWriter(fo, fieldnames=hdr, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f'{len(rows)} candidate rows for {len(motifs)} motif(s) '
          f'-> {out_tsv}')
    return len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--h5', required=True)
    ap.add_argument('--jaspar', default=None)
    ap.add_argument('--hocomoco', default=None,
                    help='concatenated HOCOMOCO PCM file')
    ap.add_argument('--out', default=None)
    ap.add_argument('--top', type=int, default=5,
                    help='top candidates per motif (default 5)')
    ap.add_argument('--gc', type=float, default=0.41,
                    help='background GC fraction (LDLR plasmid ~0.41); '
                         'sets bg=[A,C,G,T]')
    ap.add_argument('--pseudocount', type=float, default=0.8)
    a = ap.parse_args()
    gc = a.gc
    bg = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    pwms = load_pwms(a.jaspar, a.hocomoco, a.pseudocount, bg)
    if not pwms:
        print('no PWMs loaded — give --jaspar and/or --hocomoco')
        sys.exit(2)
    print(f'loaded {len(pwms)} PWMs '
          f'(JASPAR+HOCOMOCO), bg GC={gc}')
    out = a.out or a.h5.rsplit('.h5', 1)[0] + '_tf_candidates.tsv'
    annotate(a.h5, pwms, a.top, out)


if __name__ == '__main__':
    main()

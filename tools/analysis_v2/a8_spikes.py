"""Step A8 (part 2): spike-in footprints of known size and position into WT molecules (rules.json call_confidence).

A spike makes every observed A/T in [start, start + L) unmarked (perfect protection; the leaky variant keeps each
marked adenine marked with probability 0.25). Each molecule is re-called with plasmid-LLR's Kadane scan (min_llr 5,
drawdown 25, min_opps 3) and compared with itself unspiked. Two emission sets in the same scan: the cell's
position-specific fold-A table ('position') and FiberHMM's hia5_pacbio hexamer table ('kmer'; A/T whose hexamer
contains an N get no FiberHMM context and are treated as not observed, which affects backbone positions only).

Detection: a call covering the spike centre and overlapping >= 50% of the spike. Writes spikes.tsv (emission x
mode x size x centre x group: detection spiked/unspiked, net, size-correct, length error, collateral), and for
HepG2 with --positions, floor.tsv (smallest size with net detection >= 0.5 in C vs A4's prediction).

Usage: python tools/analysis_v2/a8_spikes.py --extract <a6 dir> --element 1234:1445 [--reverse] \
           --emissions <fold-A emissions.tsv.gz> [--positions <a4 positions.tsv>] --procs 4 --output <new dir>
"""
import argparse
import importlib.util
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from numba import njit

_spec = importlib.util.spec_from_file_location('a7', Path(__file__).with_name('a7_snv_edits.py'))
a7 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a7)
LLR = '/mmfs1/gscratch/stergachislab/bmallo/large_home/code/mine/plasmid-llr'
_scan = a7._scan

SIZES = [10, 15, 20, 30, 40, 60, 80, 147]
LEAKY_SIZES, LEAKY_STEP, LEAKY_KEEP = [20, 40, 147], 5, 0.25
N_C = 20000
LAB_BINS = [(0.01, 0.2), (0.2, 0.4), (0.4, 0.6)]
OVERLAP = 0.5
COLLATERAL_GAP = 10
WIDTH, OLIGO0 = 211, 21
SEED = 20261029


def size_bin(length):
    for k, (_n, lo, hi) in enumerate(a7.BINS):
        if lo <= length <= hi:
            return k
    return -1


@njit
def spike_scan(X, hit, miss, pos, usable, c_lo, c_hi, keep, u, s0, e0, centre, L, el_s, el_e, reverse, b_lo, b_hi,
               spiked):
    """Per molecule: detected (call covers centre and overlaps >= 50% of [s0, e0)), size-correct, its length; and
    coverage counts per size bin over the element (oligo order). Columns c_lo..c_hi-1 are spiked when spiked=True:
    unmarked, except a marked adenine stays marked when u[r, k] < keep."""
    n, P = X.shape
    nb = b_lo.size
    det = np.zeros(n, np.bool_)
    sok = np.zeros(n, np.bool_)
    clen = np.zeros(n, np.int64)
    cov = np.zeros((nb, el_e - el_s), np.int64)
    step = np.zeros(P)
    ih = np.zeros(P, np.bool_)
    ii = np.zeros(P, np.bool_)
    sb = -1
    for q in range(nb):
        if b_lo[q] <= L <= b_hi[q]:
            sb = q
    for r in range(n):
        lo = -1
        hi = 0
        for j in range(P):
            x = X[r, j] if usable[j] else -1
            if spiked and c_lo <= j < c_hi and x == 1 and not (u[r, j - c_lo] < keep):
                x = 0
            ii[j] = x >= 0
            ih[j] = x == 1
            step[j] = hit[j] if x == 1 else (miss[j] if x == 0 else 0.0)
            if x >= 0:
                if lo < 0:
                    lo = j
                hi = j + 1
        if lo < 0:
            continue
        s, e, _l, _o, _la, _ra = _scan(step, ih, ii, lo, hi, a7.MIN_LLR, a7.MIN_OPPS, a7.DRAWDOWN)
        for c in range(s.size):
            a = pos[s[c]]
            b = pos[e[c] - 1] + 1
            length = b - a
            ov = min(b, e0) - max(a, s0)
            if a <= centre < b and ov >= 0.5 * L:
                det[r] = True
                clen[r] = length
                if sb >= 0 and b_lo[sb] <= length <= b_hi[sb]:
                    sok[r] = True
            for q in range(nb):
                if b_lo[q] <= length <= b_hi[q]:
                    for rr in range(max(a, el_s), min(b, el_e)):
                        cov[q, (el_e - 1 - rr) if reverse else (rr - el_s)] += 1
    return det, sok, clen, cov


@njit
def base_segments(X, hit, miss, pos, usable):
    """Unspiked calls of every molecule as flat [start, end) reference arrays with row offsets."""
    n, P = X.shape
    cap = 64 * n
    ss = np.empty(cap, np.int64)
    ee = np.empty(cap, np.int64)
    off = np.zeros(n + 1, np.int64)
    k = 0
    step = np.zeros(P)
    ih = np.zeros(P, np.bool_)
    ii = np.zeros(P, np.bool_)
    for r in range(n):
        lo = -1
        hi = 0
        for j in range(P):
            x = X[r, j] if usable[j] else -1
            ii[j] = x >= 0
            ih[j] = x == 1
            step[j] = hit[j] if x == 1 else (miss[j] if x == 0 else 0.0)
            if x >= 0:
                if lo < 0:
                    lo = j
                hi = j + 1
        if lo >= 0:
            s, e, _l, _o, _la, _ra = _scan(step, ih, ii, lo, hi, a7.MIN_LLR, a7.MIN_OPPS, a7.DRAWDOWN)
            for c in range(s.size):
                if k >= cap:
                    raise ValueError('segment buffer full')
                ss[k] = pos[s[c]]
                ee[k] = pos[e[c] - 1] + 1
                k += 1
        off[r + 1] = k
    return ss[:k], ee[:k], off


@njit
def detect(ss, ee, off, s0, e0, centre, L, sb_lo, sb_hi):
    """Detection and size-correct flags from stored calls (same rule as spike_scan)."""
    n = off.size - 1
    det = np.zeros(n, np.bool_)
    sok = np.zeros(n, np.bool_)
    for r in range(n):
        for c in range(off[r], off[r + 1]):
            a = ss[c]
            b = ee[c]
            if a <= centre < b and min(b, e0) - max(a, s0) >= 0.5 * L:
                det[r] = True
                if sb_lo <= b - a <= sb_hi:
                    sok[r] = True
    return det, sok


def kmer_steps(reference, pos):
    """FiberHMM hia5_pacbio hexamer LLR steps at each A/T column (call.py kmer_arrays without its A/T assertion,
    which fails where a hexamer contains an N); usable = FiberHMM assigns a context."""
    import fiberhmm
    from fiberhmm.core.bam_reader import encode_from_query_sequence
    from fiberhmm.core.model_io import load_model
    from fiberhmm.inference.tf_recaller import N_CTX, build_llr_tables
    seq = ''.join(l.strip() for l in open(reference) if not l.startswith('>')).upper()
    tab_hit, tab_miss = build_llr_tables(load_model(str(Path(fiberhmm.__file__).parent / 'models/hia5_pacbio.json')))
    codes = encode_from_query_sequence(seq, set(range(len(seq))), edge_trim=0, mode='pacbio-fiber', context_size=3)
    c = codes[pos]
    usable = (c >= 0) & (c < N_CTX)
    hit = np.where(usable, tab_hit[np.where(usable, c, 0)], 0.0)
    miss = np.where(usable, tab_miss[np.where(usable, c, 0)], 0.0)
    return hit, miss, usable


_G = {}


def _work(task):
    emis, mode, L, centre = task
    g = _G
    hit, miss, usable = g['emis'][emis]
    el_s, el_e, rev = g['el_s'], g['el_e'], g['reverse']
    cref = (el_e - 1 - centre) if rev else (el_s + centre)
    s0 = cref - (L - 1) // 2 if not rev else cref - L // 2
    e0 = s0 + L
    pos = g['pos']
    c_lo, c_hi = int(np.searchsorted(pos, s0)), int(np.searchsorted(pos, e0))
    keep = LEAKY_KEEP if mode == 'leaky' else 0.0
    rows = []
    for grp, X in g['X'].items():
        rng = np.random.default_rng([SEED, L, centre, len(grp)])
        width = max(c_hi - c_lo, 1)
        u = rng.random((X.shape[0], width)) if mode == 'leaky' else np.ones((X.shape[0], width))
        d1, k1, l1, cov1 = spike_scan(X, hit, miss, pos, usable, c_lo, c_hi, keep, u, s0, e0, cref, L, el_s, el_e,
                                      rev, g['b_lo'], g['b_hi'], True)
        ss, ee, off, cov0 = g['base'][(emis, grp)]
        sb = size_bin(L)
        d0, k0 = detect(ss, ee, off, s0, e0, cref, L, g['b_lo'][sb], g['b_hi'][sb])
        n = X.shape[0]
        dc = (cov1 - cov0) / n
        idx = np.arange(WIDTH)
        lo_i, hi_i = centre - (L - 1) // 2, centre - (L - 1) // 2 + L      # spike in element index (approximate)
        far = (idx < lo_i - COLLATERAL_GAP) | (idx >= hi_i + COLLATERAL_GAP)
        new = d1 & ~d0
        for sub, sel in g['groups'][grp]:
            m = sel
            nm = int(m.sum())
            if nm == 0:
                continue
            nw = new[m]
            rows.append(dict(emissions=emis, mode=mode, size=L, centre_oligo=centre + OLIGO0, group=sub, n=nm,
                             n_at_spiked=c_hi - c_lo, det_spiked=d1[m].mean(), det_unspiked=d0[m].mean(),
                             net=d1[m].mean() - d0[m].mean(), size_ok_spiked=k1[m].mean(),
                             size_ok_unspiked=k0[m].mean(), net_size_ok=k1[m].mean() - k0[m].mean(),
                             length_error_median=float(np.median(l1[m][nw] - L)) if nw.any() else np.nan,
                             collateral=float(np.abs(dc[:, far]).max()) if sub in ('C', 'M') and far.any()
                             else np.nan))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--extract', required=True)
    ap.add_argument('--element', required=True)
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--emissions', required=True)
    ap.add_argument('--positions', help='A4 positions.tsv (HepG2) for the floor check')
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--n-c', type=int, default=N_C)
    ap.add_argument('--centres', default='', help='comma-separated element indices (tests); default all')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    t0 = time.time()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    el_s, el_e = map(int, a.element.split(':'))
    pos, X, distal, _, _ = a7.a6x.load(a.extract)
    ref = json.load(open(Path(a.extract) / 'extract.json'))['reference']
    hit, miss = a7.emission_steps(a.emissions, pos)
    kh, km, kuse = kmer_steps(ref, pos)
    cls = a7.classes(distal)
    rng = np.random.default_rng(SEED)
    c_rows = np.sort(rng.choice(np.flatnonzero(cls == 'C'), size=min(a.n_c, int((cls == 'C').sum())), replace=False))
    m_rows = np.flatnonzero(cls == 'M')
    dc = distal[c_rows]
    groups = {'C': [('C', np.ones(c_rows.size, bool))] +
              [(f'C {lo:.2f}-{hi:.2f}', (dc > lo) & (dc <= hi)) for lo, hi in LAB_BINS],
              'M': [('M', np.ones(m_rows.size, bool))]}
    _G.update(X={'C': np.ascontiguousarray(X[c_rows]), 'M': np.ascontiguousarray(X[m_rows])}, pos=pos.astype(np.int64),
              emis={'position': (hit, miss, np.ones(pos.size, np.bool_)), 'kmer': (kh, km, kuse)},
              el_s=el_s, el_e=el_e, reverse=a.reverse, groups=groups,
              b_lo=np.array([b[1] for b in a7.BINS], np.int64), b_hi=np.array([b[2] for b in a7.BINS], np.int64))
    del X
    base = {}
    for e, (h_, m_, us_) in _G['emis'].items():
        for grp, Xg in _G['X'].items():
            ss, ee, off = base_segments(Xg, h_, m_, _G['pos'], us_)
            one = np.ones((Xg.shape[0], 1))
            _, _, _, cov0 = spike_scan(Xg, h_, m_, _G['pos'], us_, 0, 0, 0.0, one, 0, 1, 0, 1, el_s, el_e,
                                       a.reverse, _G['b_lo'], _G['b_hi'], False)
            base[(e, grp)] = (ss, ee, off, cov0)
    _G['base'] = base
    centres = [int(c) for c in a.centres.split(',')] if a.centres else list(range(WIDTH))
    tasks = [(e, 'perfect', L, c) for e in ('position', 'kmer') for L in SIZES for c in centres]
    tasks += [(e, 'leaky', L, c) for e in ('position', 'kmer') for L in LEAKY_SIZES for c in centres
              if c % LEAKY_STEP == 0]
    _work(tasks[0])                                       # compile before forking
    print(f'{len(tasks)} spikes; C {c_rows.size}, M {m_rows.size}; kmer usable {int(kuse.sum())}/{kuse.size}',
          flush=True)
    rows = []
    with mp.get_context('fork').Pool(a.procs) as pool:
        for i, r in enumerate(pool.imap_unordered(_work, tasks, chunksize=4)):
            rows += r
            if i % 500 == 0:
                print(f'  {i}/{len(tasks)}, {time.time() - t0:.0f} s', flush=True)
    t = pd.DataFrame(rows).sort_values(['emissions', 'mode', 'size', 'centre_oligo', 'group'])
    t.to_csv(out / 'spikes.tsv', sep='\t', index=False)
    summary = dict(molecules=dict(C=int(c_rows.size), M=int(m_rows.size)), spikes=len(tasks),
                   kmer_unusable_columns=[int(p) for p in pos[~kuse]],
                   settings=dict(sizes=SIZES, leaky_sizes=LEAKY_SIZES, leaky_keep=LEAKY_KEEP, overlap=OVERLAP,
                                 lab_bins=LAB_BINS, seed=SEED, emissions=a.emissions, extract=a.extract,
                                 element=[el_s, el_e], reverse=a.reverse))
    if a.positions:
        p = pd.read_csv(a.positions, sep='\t')
        c = t[(t.emissions == 'position') & (t['mode'] == 'perfect') & (t.group == 'C')]
        fl = []
        for oligo, g in c.groupby('centre_oligo'):
            ok = g[g.net >= 0.5]
            fl.append(dict(centre_oligo=oligo, smallest_detected=int(ok['size'].min()) if len(ok) else np.nan))
        fl = pd.DataFrame(fl).merge(p[['oligo_pos1', 'base', 'min_fp_llr_bp']], left_on='centre_oligo',
                                    right_on='oligo_pos1', how='left').drop(columns='oligo_pos1')
        fl.to_csv(out / 'floor.tsv', sep='\t', index=False)
    summary['seconds'] = round(time.time() - t0)
    json.dump(summary, open(out / 'summary.json', 'w'), indent=1)
    print(json.dumps(summary)[:2000])


if __name__ == '__main__':
    main()

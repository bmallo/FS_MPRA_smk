"""Footprints defined from the data, and every variant's effect on each (HepG2 run 3; after Phase A, exploratory).

1. Calls. Every molecule with labeling >= 0.2 is called with plasmid-LLR (min LLR 5, drawdown 25, run-3 fold-A
   emissions), keeping each call's midpoint, length and chance probability (docs/footprint_confidence_design.md:
   log10 P = dose * sum log(1 - p_a) over its unmarked adenines / ln 10; confident if P < 1e-3).
2. WT footprint map (variant-blind). Calls of 10-90 bp from WT molecules are binned by midpoint (oligo) x length
   per molecule, for moderately labeled molecules (0.2-0.6, "C") and heavily labeled ones (>= 0.6, "M"). Peaks of
   the smoothed C density become footprints (a midpoint window and a length range) by the v1 rule and settings
   (fsmpra.footprint_sizes.find_features, config/ldlr_footprint_sizes.json 'map').
3. Variant tests. Occupancy of a footprint = the molecule has a call inside its box (all calls, or confident calls
   only). Each SNV against WT with the A9 random-effects test (labeling-matched, clone-aware), with no mask, the A7
   mask and the +/-5 bp mask (masked adenines unobserved before calling, in the SNV's molecules and WT alike).
   BH at 5% across footprints within each SNV, mask and call set (descriptive, not the A10 hierarchy).

Writes calls.npz, footprints.tsv, wt_map.npz, effects.parquet, summary.json.

Usage: python tools/analysis_v2/fp_features.py --root results/analysis_v2_20261028 --catalog <catalog.json> \
           --procs 4 --output <new dir>
"""
import argparse
import importlib.util
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from numba import njit

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from fsmpra.footprint_sizes import find_features  # noqa: E402

_spec = importlib.util.spec_from_file_location('a10', HERE / 'a10_effects.py')
a10 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a10)
a9f, a7, a6x = a10.a9f, a10.a7, a10.a6x
_scan = a7._scan

EL = (1234, 1445)
WIDTH, OLIGO0 = 211, 21
LEN_RANGE = (10, 90)
KMAX = 32
REF_LABELING = 0.70
CONF_LOG10 = -3.0
C_RANGE, M_MIN = (0.2, 0.6), 0.6
Q = 0.05
EMISSIONS = 'results/llr_vs_fiberhmm_20261026/results/emissions.tsv.gz'


@njit
def call_table(X, hit, miss, lp, pos, el_s, el_e, kmax, min_llr, min_opps, drawdown):
    """Per molecule, calls with midpoint in the element: midpoint (oligo), length, start/end (oligo, unclipped) and
    S = sum of log(1 - p_a) over the call's unmarked adenines."""
    n, P = X.shape
    mid = np.full((n, kmax), -1, np.int64)
    ln = np.zeros((n, kmax), np.int64)
    st = np.zeros((n, kmax), np.int64)
    en = np.zeros((n, kmax), np.int64)
    S = np.zeros((n, kmax))
    step = np.zeros(P)
    ih = np.zeros(P, np.bool_)
    ii = np.zeros(P, np.bool_)
    for r in range(n):
        lo = -1
        hi = 0
        for j in range(P):
            x = X[r, j]
            ii[j] = x >= 0
            ih[j] = x == 1
            step[j] = hit[j] if x == 1 else (miss[j] if x == 0 else 0.0)
            if x >= 0:
                if lo < 0:
                    lo = j
                hi = j + 1
        if lo < 0:
            continue
        s, e, _l, _o, _la, _ra = _scan(step, ih, ii, lo, hi, min_llr, min_opps, drawdown)
        k = 0
        for c in range(s.size):
            a = pos[s[c]]
            b = pos[e[c] - 1] + 1
            m = (a + b - 1) // 2
            if m < el_s or m >= el_e or k >= kmax:
                continue
            tot = 0.0
            for j in range(s[c], e[c]):
                if X[r, j] == 0:
                    tot += lp[j]
            mid[r, k] = m - el_s + 21
            ln[r, k] = b - a
            st[r, k] = a - el_s + 21
            en[r, k] = b - 1 - el_s + 21
            S[r, k] = tot
            k += 1
    return mid, ln, st, en, S


def dose(labeling):
    L = np.clip(labeling, 0.01, 0.99)
    return np.log(1 - L) / np.log(1 - REF_LABELING)


def occupancy(tab, distal, feats):
    """Molecules x footprints: 1 if any call falls in the box, 2 if any confident call does (bit flags)."""
    mid, ln, _, _, S = tab
    conf = (dose(distal)[:, None] * S / np.log(10)) < CONF_LOG10
    out = np.zeros((mid.shape[0], len(feats)), np.uint8)
    for f, b in enumerate(feats):
        inbox = (mid >= b['mid_lo'] + OLIGO0) & (mid <= b['mid_hi'] + OLIGO0) & (ln >= b['len_lo']) & (ln <= b['len_hi'])
        out[:, f] = inbox.any(1) * 1 + (inbox & conf).any(1) * 2
    return out


_G = {}


def _tab(args):
    rows, cols = args
    g = _G
    X = g['X'][rows].copy()
    if cols is not None and cols.size:
        X[:, cols] = -1
    return call_table(X, g['hit'], g['miss'], g['lp'], g['pos'], EL[0], EL[1], KMAX, a7.MIN_LLR, a7.MIN_OPPS, a7.DRAWDOWN)


def tables(rows, cols, pool, procs):
    parts = pool.map(_tab, [(c, cols) for c in np.array_split(rows, max(1, procs * 8))])
    return tuple(np.concatenate([p[i] for p in parts]) for i in range(5))


def wt_map(tab, distal, rows, mp_cfg):
    """Calls per molecule by midpoint x length for C and M WT molecules; footprints from the C density."""
    mid, ln, *_ = tab
    lo, hi = LEN_RANGE
    L = hi - lo + 1

    def dens(sel):
        m, l = mid[sel], ln[sel]
        ok = (m >= OLIGO0) & (m < OLIGO0 + WIDTH) & (l >= lo) & (l <= hi)
        idx = (m[ok] - OLIGO0) * L + (l[ok] - lo)
        return np.bincount(idx, minlength=WIDTH * L).reshape(WIDTH, L) / max(sel.sum(), 1)
    d = distal[rows]
    C = dens((d >= C_RANGE[0]) & (d < C_RANGE[1]))
    M = dens(d >= M_MIN)
    feats, sm, bg = find_features(C, M, mp_cfg, LEN_RANGE)
    return feats, C, M, sm, bg


def _snv(task):
    """One SNV and mask: indicators for its molecules and (masked) WT, then the A9 contrast."""
    v, mask, cols = task
    g = _G
    rv = g['snv_rows'][v]
    if cols is None:
        yw, yv = g['occ_wt'], g['occ_all'][rv]
        Sy, So = g['ref_all'].Sy, g['ref_all'].So
        Syc, Soc = g['ref_conf'].Sy, g['ref_conf'].So
    else:
        tw = call_table(_masked(g['X'][g['wt2']], cols), g['hit'], g['miss'], g['lp'], g['pos'], EL[0], EL[1], KMAX,
                        a7.MIN_LLR, a7.MIN_OPPS, a7.DRAWDOWN)
        yw = occupancy(tw, g['distal'][g['wt2']], g['feats'])
        Sy, So = a10.stratum_sums((yw & 1).astype(float), g['st'][g['wt2']])
        Syc, Soc = a10.stratum_sums(((yw & 2) > 0).astype(float), g['st'][g['wt2']])
        tv = call_table(_masked(g['X'][rv], cols), g['hit'], g['miss'], g['lp'], g['pos'], EL[0], EL[1], KMAX,
                        a7.MIN_LLR, a7.MIN_OPPS, a7.DRAWDOWN)
        yv = occupancy(tv, g['distal'][rv], g['feats'])
    out = []
    for kind, bit, sy, so, ref in (('all', 1, Sy, So, g['ref_all']), ('confident', 2, Syc, Soc, g['ref_conf'])):
        Y = (yv & bit) > 0
        p = np.divide(sy, so, out=np.full_like(sy, np.nan), where=so > 0)
        r = a9f.contrast(Y, np.ones_like(Y), g['st'][rv], g['clone'][rv], p, so,
                         SimpleNamespace(V=ref.V, W2=ref.W2, tau2=ref.tau2))
        out.append(dict(snv=v, mask=mask, calls=kind, occ=Y.mean(0), **{k: r[k] for k in ('n', 'C', 'D', 'se', 'p', 'expected')}))
    return out


def _masked(X, cols):
    X = X.copy()
    X[:, cols] = -1
    return X


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    t0 = time.time()
    root = Path(a.root)
    ext = root / 'a9_statistics/extract/HepG2_run3'
    pos, X, distal, _, _ = a6x.load(ext)
    meta = np.load(ext / 'meta.npz')
    label, cls, clone = meta['label'], meta['cls'], meta['clone']
    hit, miss = a7.emission_steps(EMISSIONS, pos)
    t = pd.read_csv(EMISSIONS, sep='\t')
    t = t[t.fold == 'A']
    pa = dict(zip(t.ref_pos - 1, t.p_a))
    lp = np.log(1 - np.clip([pa.get(int(p), 0.0) for p in pos], 0, 0.9999))
    k2 = np.flatnonzero(distal >= a9f.E2_MIN)
    wt2 = k2[cls[k2] == 'WT']
    st = a9f.strata(distal)
    _G.update(X=X, hit=hit, miss=miss, lp=lp, pos=pos)
    with mp.get_context('fork').Pool(a.procs) as pool:
        tab = tables(k2, None, pool, a.procs)
    print(f'calls for {k2.size} molecules, {time.time() - t0:.0f} s', flush=True)
    row = np.full(distal.size, -1)
    row[k2] = np.arange(k2.size)
    np.savez_compressed(out / 'calls.npz', rows=k2, mid=tab[0], length=tab[1], start=tab[2], end=tab[3], S=tab[4])

    mp_cfg = json.load(open(ROOT / 'config/ldlr_footprint_sizes.json'))['map']
    wtab = tuple(x[row[wt2]] for x in tab)
    feats, C, M, sm, bg = wt_map(wtab, distal, wt2, mp_cfg)
    occ_all = occupancy(tab, distal[k2], feats)
    occ_full = np.zeros((distal.size, len(feats)), np.uint8)
    occ_full[k2] = occ_all
    ow = occ_full[wt2]
    for f, b in zip(range(len(feats)), feats):
        b['oligo_mid'] = b['mid_peak'] + OLIGO0
        b['oligo_lo'], b['oligo_hi'] = b['mid_lo'] + OLIGO0, b['mid_hi'] + OLIGO0
        b['span_lo'] = int(round(b['oligo_mid'] - b['len_peak'] / 2))
        b['span_hi'] = int(round(b['oligo_mid'] + b['len_peak'] / 2))
        b['wt_occupancy'] = float((ow[:, f] & 1).mean())
        b['wt_occupancy_confident'] = float(((ow[:, f] & 2) > 0).mean())
        b['confident_share'] = b['wt_occupancy_confident'] / b['wt_occupancy'] if b['wt_occupancy'] > 0 else np.nan
        for lo_, hi_ in ((0.2, 0.4), (0.4, 0.6), (0.6, 1.01)):
            sel = (distal[wt2] >= lo_) & (distal[wt2] < hi_)
            b[f'wt_occ_{lo_:.1f}'] = float((ow[sel, f] & 1).mean())
    fe = pd.DataFrame(feats)
    fe.to_csv(out / 'footprints.tsv', sep='\t', index=False)
    np.savez_compressed(out / 'wt_map.npz', C=C, M=M, smoothed=sm, background=bg)
    print(fe[['feature', 'oligo_mid', 'len_peak', 'oligo_lo', 'oligo_hi', 'len_lo', 'len_hi', 'wt_occupancy',
              'confident_share', 'c_over_m', 'flag']].round(4).to_string(), flush=True)

    # variant tests
    snv = np.flatnonzero(cls == 'SNV')
    snv2 = snv[distal[snv] >= a9f.E2_MIN]
    ref_all = a9f.Reference((ow & 1) > 0, np.ones(ow.shape, bool), st[wt2], clone[wt2])
    ref_conf = a9f.Reference((ow & 2) > 0, np.ones(ow.shape, bool), st[wt2], clone[wt2])
    masks = {m: a10.snv_masks(a.catalog, 'results/analysis_v2_20261028/a4_informative/visibility/variants_visibility.tsv',
                              pos, EL[0], EL[1], False, pm5=(m == 'pm5')) for m in ('A7', 'pm5')}
    snvs = sorted(set(label[snv2]) & set(masks['A7']))
    _G.update(distal=distal, st=st, clone=clone, wt2=wt2, feats=feats, occ_wt=ow, occ_all=occ_full,
              ref_all=ref_all, ref_conf=ref_conf, snv_rows={v: snv2[label[snv2] == v] for v in snvs})
    tasks = [(v, 'none', None) for v in snvs]
    tasks += [(v, m, masks[m][v]['mask']) for m in ('A7', 'pm5') for v in snvs if masks[m][v]['mask'].size]
    tasks += [(v, m, None) for m in ('A7', 'pm5') for v in snvs if not masks[m][v]['mask'].size]
    rows_out = []
    with mp.get_context('fork').Pool(a.procs) as pool:
        for i, res in enumerate(pool.imap_unordered(_snv, tasks, chunksize=2)):
            rows_out += res
            if i % 200 == 0:
                print(f'{i + 1}/{len(tasks)} variant x mask tasks, {time.time() - t0:.0f} s', flush=True)
    recs = []
    for r in rows_out:
        for f in range(len(feats)):
            recs.append(dict(snv=r['snv'], mask=r['mask'], calls=r['calls'], footprint=feats[f]['feature'],
                             occ=r['occ'][f], n=r['n'][f], C=r['C'][f], D=r['D'][f], se=r['se'][f], p=r['p'][f],
                             expected=r['expected'][f]))
    eff = pd.DataFrame(recs)
    eff['sig'] = False
    for _, g in eff.groupby(['snv', 'mask', 'calls']):
        eff.loc[g.index, 'sig'] = a10.bh(g.p.to_numpy(), Q)
    eff.to_parquet(out / 'effects.parquet', index=False)
    summ = dict(footprints=len(feats), molecules=int(k2.size), wt_molecules=int(wt2.size), snvs=len(snvs),
                tasks=len(tasks), sig_by_mask=eff.groupby(['mask', 'calls']).sig.sum().astype(int).to_dict().__repr__(),
                settings=dict(len_range=LEN_RANGE, c_range=C_RANGE, m_min=M_MIN, conf_log10=CONF_LOG10,
                              ref_labeling=REF_LABELING, map=mp_cfg), seconds=round(time.time() - t0))
    json.dump(summ, open(out / 'summary.json', 'w'), indent=1, default=float)
    print(json.dumps(summ, default=float)[:1500])


if __name__ == '__main__':
    main()

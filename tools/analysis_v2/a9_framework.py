"""Step A9 (part 2): calibrate the variant-vs-WT test on WT-only pseudo-variants; covariate balance; MDE per SNV.

Settings: rules.json statistics. Uses WT chromatin data and SNV metadata only (clone sizes, labeling, passes,
strand, read quality); no SNV chromatin endpoint is read.

Endpoints. E1: element adenine unmarked, molecules with labeling > 0.01. E2: element position covered by a
plasmid-LLR call of size bin b, molecules with labeling >= 0.2.
Test. D = mean over the variant's molecules of (y - p_WT(stratum)). Primary variance: random effects, binomial plus
the WT-estimated between-clone variance tau2 times sum_c (n_c/n)^2 (both sides), normal p. Reported comparators:
clone sandwich (CR1, t with C - 1 df; the prespecified test, which fails calibration) and naive binomial.
Pseudo-variants. Per SNV, --draws draws of WT clones matched to the SNV's clone sizes (smallest available WT clone
of size >= k, subsampled to k); reference = the other WT molecules.

Writes calibration.parquet (per draw x feature), calibration.tsv (summary), balance.tsv, mde.parquet,
mde_summary.tsv, summary.json.

Usage: python tools/analysis_v2/a9_framework.py --extract <a9 extract dir> --emissions <fold-A emissions.tsv.gz> \
           --element 1234:1445 [--reverse] --spikes <a8 spikes.tsv> --procs 4 --output <new dir>
"""
import argparse
import importlib.util
import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats

_spec = importlib.util.spec_from_file_location('a8r', Path(__file__).with_name('a8_repro.py'))
a8r = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a8r)
a7, a6x = a8r.a7, a8r.a7.a6x

N_STRATA = a8r.N_STRATA
E1_MIN, E2_MIN = 0.01, 0.2
MIN_N, MIN_C = 10, 3
MIN_EXPECTED = 5                 # expected minority count under WT (amendment 2)
DRAWS = 4
SEED = 20261029
ALPHA, POWER_Z = 0.05, stats.norm.ppf(0.8)
TIERS = [(1, 2), (3, 9), (10, 29), (30, 10 ** 9)]
SPIKE_SIZE = {'sub_TF': 15, 'TF': 30, 'PIC': 60, 'NUC': 147}
CALL_GROUPS = [('C 0.20-0.40', 0.2, 0.4), ('C 0.40-0.60', 0.4, 0.6), ('M', 0.6, 1.01)]
OLIGO0 = 21


def strata(distal):
    return np.minimum((distal / a7.STRATUM_WIDTH).astype(np.int64), N_STRATA - 1)


def tier(c):
    for lo, hi in TIERS:
        if lo <= c <= hi:
            return f'{lo}-{hi}' if hi < 10 ** 9 else f'{lo}+'
    return 'none'


def clone_matrix(clone):
    u, inv = np.unique(clone, return_inverse=True)
    return sparse.csr_matrix((np.ones(clone.size), (inv, np.arange(clone.size))), shape=(u.size, clone.size))


class Reference:
    """WT sums per labeling stratum and feature; the clone-robust (sandwich) variance of each stratum rate; the
    between-clone variance tau2 per feature (one-way random-effects moment estimator on stratum-centred residuals)
    and, per stratum, sum_c (n_c / N_s)^2 for the random-effects variance of the WT rate."""

    def __init__(self, Y, O, st, clone):
        F = Y.shape[1]
        self.Sy = np.zeros((N_STRATA, F))
        self.So = np.zeros((N_STRATA, F))
        for s in np.unique(st):
            m = st == s
            self.Sy[s] = (Y[m] & O[m]).sum(0)
            self.So[s] = O[m].sum(0)
        p = np.divide(self.Sy, self.So, out=np.zeros_like(self.Sy), where=self.So > 0)
        self.V = np.zeros_like(p)
        self.W2 = np.zeros_like(p)
        for s in np.unique(st):
            m = np.flatnonzero(st == s)
            A = clone_matrix(clone[m])
            cs = A @ ((Y[m] - p[s]) * O[m])
            nc = A @ O[m].astype(float)
            with np.errstate(divide='ignore', invalid='ignore'):
                self.V[s] = np.where(self.So[s] > 0, (cs ** 2).sum(0) / self.So[s] ** 2, 0.0)
                self.W2[s] = np.where(self.So[s] > 0, (nc ** 2).sum(0) / self.So[s] ** 2, 0.0)
        self.tau2 = tau2(Y, O, p[st], clone)

    def rates(self, Y=None, O=None, st=None):
        """WT rates per stratum, with the given molecules' contributions removed (the pseudo-variant's own)."""
        Sy, So = self.Sy, self.So
        if Y is not None:
            Sy, So = Sy.copy(), So.copy()
            np.subtract.at(Sy, st, Y * O)
            np.subtract.at(So, st, O.astype(float))
        p = np.divide(Sy, So, out=np.full_like(Sy, np.nan), where=So > 0)
        return p, So


def tau2(Y, O, pm, clone):
    """Between-clone variance of clone-mean residuals per feature (unbalanced one-way random-effects moments):
    Q = sum_c n_c d_c^2 - N dbar^2; E[Q] = tau2 (N - sum n_c^2 / N) + sum_c S_c / n_c - S / N, with S_c the summed
    binomial variances p(1-p) of clone c's molecules. Truncated at 0."""
    A = clone_matrix(clone)
    r = (Y - pm) * O
    sig = pm * (1 - pm) * O
    cs, nc, sc = A @ r, A @ O.astype(float), A @ sig
    N, S = nc.sum(0), sc.sum(0)
    with np.errstate(divide='ignore', invalid='ignore'):
        Q = np.where(nc > 0, cs ** 2 / nc, 0).sum(0) - cs.sum(0) ** 2 / N
        e0 = np.where(nc > 0, sc / nc, 0).sum(0) - S / N
        den = N - (nc ** 2).sum(0) / N
        t = (Q - e0) / den
    return np.where(np.isfinite(t), np.maximum(t, 0.0), 0.0)


def contrast(Y, O, st, clone, p, So, ref):
    """Variant-vs-WT test for each feature (columns of Y). Y, O: molecules x features (0/1; O = observed).
    p, So: WT rate and observation count per stratum x feature; ref: the Reference (V, W2, tau2).
    Three variances for the same D: random-effects (primary: binomial + tau2 * sum_c w_c^2, both sides), clone
    sandwich (CR1, t with C - 1 df) and naive binomial."""
    pm = p[st]
    use = O & np.isfinite(pm)
    pm = np.where(use, pm, 0.0)
    r = np.where(use, Y - pm, 0.0)
    n = use.sum(0).astype(float)
    D = np.divide(r.sum(0), n, out=np.full(n.size, np.nan), where=n > 0)
    A = clone_matrix(clone)
    cs = A @ (r - np.where(use, D, 0.0))
    nc = A @ use.astype(float)
    C = (nc > 0).sum(0).astype(float)
    with np.errstate(divide='ignore', invalid='ignore'):
        ns = np.zeros((N_STRATA, Y.shape[1]))
        np.add.at(ns, st, use.astype(float))
        w = ns / n
        binom_v = (pm * (1 - pm) * use).sum(0) / n ** 2
        binom_wt = (w ** 2 * np.nan_to_num(p * (1 - p) / So)).sum(0)
        v_naive = binom_v + binom_wt
        v_re = v_naive + ref.tau2 * ((nc ** 2).sum(0) / n ** 2 + (w ** 2 * ref.W2).sum(0))
        v_sand = C / (C - 1) * (cs ** 2).sum(0) / n ** 2 + (w ** 2 * np.nan_to_num(ref.V)).sum(0)
        se, se_sand, se_naive = np.sqrt(v_re), np.sqrt(v_sand), np.sqrt(v_naive)
        p_re = 2 * stats.norm.sf(np.abs(D / se))
        p_sand = 2 * stats.t.sf(np.abs(D / se_sand), np.maximum(C - 1, 1))
        p_naive = 2 * stats.norm.sf(np.abs(D / se_naive))
    expected = np.minimum(pm * use, (1 - pm) * use).sum(0)
    ok = (n >= MIN_N) & (C >= MIN_C) & (expected >= MIN_EXPECTED)
    return dict(n=n, C=C, D=D, expected=expected, se=se, se_sandwich=se_sand, se_naive=se_naive, p=np.where(ok, p_re, np.nan),
                p_sandwich=np.where(ok, p_sand, np.nan), p_naive=np.where(ok, p_naive, np.nan))


def draw_clones(sizes, pools, rng):
    """Molecule rows of WT clones matched to the target clone sizes (largest first): for size k, a random unused WT
    clone of the smallest available size >= k, subsampled to k; the largest unused clone if none is that big.
    pools: dict size -> list of row arrays. Each draw starts every size pool at a random offset and walks it, so
    no clone is used twice within a draw. Returns rows, clone labels and the number of short clones."""
    keys = np.array(sorted(pools))
    start, used = {}, {}
    rows, labels, short = [], [], 0

    def take_from(k):
        if k not in start:
            start[k], used[k] = int(rng.integers(len(pools[k]))), 0
        src = pools[k][(start[k] + used[k]) % len(pools[k])]
        used[k] += 1
        return src

    def free(k):
        return used.get(k, 0) < len(pools[k])

    for c, k in enumerate(sorted(sizes, reverse=True)):
        src = None
        for x in keys[np.searchsorted(keys, k):]:
            if free(x):
                src = take_from(x)
                break
        if src is None:
            x = next(x for x in keys[::-1] if free(x))
            src = take_from(x)
            short += 1
        take = src if src.size <= k else rng.choice(src, k, replace=False)
        rows.append(take)
        labels.append(np.full(take.size, c))
    return np.concatenate(rows), np.concatenate(labels), short


def smd(a, b):
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return np.nan
    s = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return (a.mean() - b.mean()) / s if s > 0 else 0.0


def bh_any(p, q=0.05):
    p = np.sort(p[np.isfinite(p)])
    m = p.size
    return bool(m and (p <= q * np.arange(1, m + 1) / m).any())


def _calls(args):
    X, hit, miss, pos, el_s, el_e, rev = args
    b_lo = np.array([b[1] for b in a7.BINS], np.int64)
    b_hi = np.array([b[2] for b in a7.BINS], np.int64)
    return a8r.call_bits(X, hit, miss, pos, el_s, el_e, rev, b_lo, b_hi)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--extract', required=True)
    ap.add_argument('--emissions', required=True)
    ap.add_argument('--element', required=True)
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--spikes', required=True)
    ap.add_argument('--draws', type=int, default=DRAWS)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    t0 = time.time()
    el_s, el_e = (int(x) for x in a.element.split(':'))
    pos, X, distal, passes, strand = a6x.load(a.extract)
    meta = np.load(Path(a.extract) / 'meta.npz')
    label, cls, clone, rq = meta['label'], meta['cls'], meta['clone'], meta['read_quality']
    wt = np.flatnonzero(cls == 'WT')
    Xw = X[wt]
    del X                                                     # SNV chromatin data are not used in A9
    dw, cw, stw = distal[wt], clone[wt], strata(distal[wt])
    print(f'{wt.size} WT molecules, {np.unique(cw).size} WT clones, {time.time() - t0:.0f} s', flush=True)

    # E1: element adenines, unmarked
    ecol = np.flatnonzero((pos >= el_s) & (pos < el_e))
    e1_ref = pos[ecol]
    e1_oligo = (el_e - 1 - e1_ref + OLIGO0) if a.reverse else (e1_ref - el_s + OLIGO0)
    k1 = dw > E1_MIN
    Y1 = (Xw[:, ecol] == 0)
    O1 = (Xw[:, ecol] >= 0) & k1[:, None]
    # E2: calls on WT molecules with labeling >= 0.2
    k2 = np.flatnonzero(dw >= E2_MIN)
    hit, miss = a7.emission_steps(a.emissions, pos)
    chunks = np.array_split(k2, max(1, a.procs * 8))
    with mp.get_context('fork').Pool(a.procs) as pool:
        bits = pool.map(_calls, [(Xw[c], hit, miss, pos, el_s, el_e, a.reverse) for c in chunks])
    nb, W = len(a7.BINS), el_e - el_s
    Y2 = np.zeros((wt.size, nb * W), bool)
    Y2[k2] = np.concatenate(bits).reshape(k2.size, nb * W)
    O2 = np.zeros_like(Y2)
    O2[k2] = True
    e2_bin = np.repeat([b[0] for b in a7.BINS], W)
    e2_oligo = np.tile(np.arange(W) + OLIGO0, nb)
    print(f'E1 {ecol.size} adenines; E2 {k2.size} molecules called, {time.time() - t0:.0f} s', flush=True)
    ref1, ref2 = Reference(Y1, O1, stw, cw), Reference(Y2, O2, stw, cw)

    # WT clone pools by size, SNV clone sizes
    wclone_u, wclone_inv = np.unique(cw, return_inverse=True)
    order = np.argsort(wclone_inv, kind='stable')
    bounds = np.cumsum(np.bincount(wclone_inv))
    groups = np.split(order, bounds[:-1])
    rng = np.random.default_rng(SEED)
    pools = {}
    for g in groups:
        pools.setdefault(g.size, []).append(g)
    for k in pools:
        rng.shuffle(pools[k])
    snv = np.flatnonzero(cls == 'SNV')
    snvs = sorted(set(label[snv]))
    print(f'{len(snvs)} SNVs, {snv.size} SNV molecules', flush=True)

    cal_rows, bal_rows, se_rows = [], [], []
    cov = dict(distal=distal, passes=passes.astype(float), strand=strand.astype(float), read_quality=rq)
    cov_wt = {k: x[wt] for k, x in cov.items()}
    for vi, v in enumerate(snvs):
        rows_v = snv[label[snv] == v]
        sizes = np.bincount(np.unique(clone[rows_v], return_inverse=True)[1])
        C_v = sizes.size
        b = dict(snv=v, molecules=rows_v.size, clones=C_v, tier=tier(C_v), draw=-1)
        for name, x in cov.items():
            b[f'smd_{name}'] = smd(x[rows_v], cov_wt[name])
        bal_rows.append(b)
        for d in range(a.draws):
            rr, cl, short = draw_clones(sizes, pools, rng)
            b = dict(snv=v, molecules=rr.size, clones=C_v, tier=tier(C_v), draw=d, short_clones=short)
            for name, x in cov_wt.items():
                b[f'smd_{name}'] = smd(x[rr], x)
            bal_rows.append(b)
            for ep, Y, O, ref in (('E1', Y1, O1, ref1), ('E2', Y2, O2, ref2)):
                p, So = ref.rates(Y[rr], O[rr], stw[rr])
                r = contrast(Y[rr], O[rr], stw[rr], cl, p, So, ref)
                f = pd.DataFrame(r)
                f['feature'] = np.arange(f.shape[0])
                f['endpoint'], f['snv'], f['draw'], f['clones_snv'] = ep, v, d, C_v
                cal_rows.append(f)
        if vi % 50 == 0:
            print(f'{vi + 1}/{len(snvs)} SNVs, {time.time() - t0:.0f} s', flush=True)

    cal = pd.concat(cal_rows, ignore_index=True)
    cal['tier'] = cal.clones_snv.map(tier)
    cal.to_parquet(out / 'calibration.parquet', index=False)
    bal = pd.DataFrame(bal_rows)
    bal.to_csv(out / 'balance.tsv', sep='\t', index=False)

    # calibration summary
    summ = []
    for ep, g0 in cal.groupby('endpoint'):
        for tr, g in list(g0.groupby('tier')) + [('all', g0)]:
            for test, col in (('random_effects', 'p'), ('sandwich', 'p_sandwich'), ('naive', 'p_naive')):
                pv = g[col].dropna().to_numpy()
                if pv.size == 0:
                    continue
                chi = stats.chi2.isf(np.clip(pv, 1e-300, 1), 1)
                anyd = g.groupby(['snv', 'draw'])[col].apply(lambda s: bh_any(s.to_numpy()))
                summ.append(dict(endpoint=ep, tier=tr, test=test, tests=pv.size, share_p05=float((pv < 0.05).mean()),
                                 lam=float(np.median(chi) / stats.chi2.ppf(0.5, 1)),
                                 bh_any=float(anyd.mean()), draws=int(anyd.size)))
    summ = pd.DataFrame(summ)
    summ.to_csv(out / 'calibration.tsv', sep='\t', index=False)

    # balance flags: SNV |SMD| above the 99th percentile of pseudo-variant |SMD| in its tier
    flags = {}
    real, pseudo = bal[bal.draw == -1], bal[bal.draw >= 0]
    for name in cov:
        col = f'smd_{name}'
        thr = pseudo.groupby('tier')[col].apply(lambda s: np.nanpercentile(np.abs(s), 99))
        f = np.abs(real[col]) > real.tier.map(thr)
        flags[name] = dict(threshold_by_tier=thr.round(3).to_dict(), flagged=int(f.sum()), of=int(real[col].notna().sum()))

    # MDE per SNV x feature: median SE over its draws
    se = cal.groupby(['snv', 'endpoint', 'feature']).agg(se=('se', 'median'), C=('C', 'median'), n=('n', 'median'))
    se = se.reset_index()
    m_tests = cal.dropna(subset=['p']).groupby(['snv', 'draw', 'endpoint']).size().groupby(['snv', 'endpoint']).median()
    se['m'] = [m_tests.get((s, e), np.nan) for s, e in zip(se.snv, se.endpoint)]
    se['mde'] = (stats.norm.ppf(1 - ALPHA / 2) + POWER_Z) * se.se
    se['mde_bonf'] = (stats.norm.ppf(1 - ALPHA / 2 / se.m.fillna(1)) + POWER_Z) * se.se
    se.loc[(se.C < MIN_C) | (se.n < MIN_N), ['mde', 'mde_bonf']] = np.nan
    e1 = se.endpoint == 'E1'
    se['oligo'] = np.where(e1, e1_oligo[np.minimum(se.feature, e1_oligo.size - 1)], e2_oligo[se.feature])
    se['bin'] = np.where(e1, '', e2_bin[se.feature])
    # occupancy MDE for E2: divide by net spike detection at that position and size (A8)
    sp = pd.read_csv(a.spikes, sep='\t')
    sp = sp[(sp.emissions == 'position') & (sp['mode'] == 'perfect')]
    share = np.array([((dw >= lo) & (dw < hi)).sum() for _, lo, hi in CALL_GROUPS], float)
    share /= share.sum()
    det = {}
    for bn, size in SPIKE_SIZE.items():
        s = sp[sp['size'] == size].pivot_table(index='centre_oligo', columns='group', values='net')
        det[bn] = sum(w * s[g] for (g, _, _), w in zip(CALL_GROUPS, share))
    se['detection'] = [det[b].get(o, np.nan) if b else np.nan for b, o in zip(se['bin'], se.oligo)]
    se['mde_occupancy'] = se.mde / se.detection
    se.to_parquet(out / 'mde.parquet', index=False)
    ms = []
    for (ep, bn), g in se.groupby(['endpoint', 'bin']):
        per = g.groupby('snv')[['mde', 'mde_bonf', 'mde_occupancy']].median()
        ms.append(dict(endpoint=ep, bin=bn or 'adenine', snvs_tested=int(per.mde.notna().sum()),
                       mde_median=per.mde.median(), mde_q10=per.mde.quantile(0.1), mde_q90=per.mde.quantile(0.9),
                       mde_bonf_median=per.mde_bonf.median(), mde_occ_median=per.mde_occupancy.median()))
    ms = pd.DataFrame(ms)
    ms.to_csv(out / 'mde_summary.tsv', sep='\t', index=False)
    info = dict(settings=dict(extract=a.extract, emissions=a.emissions, element=[el_s, el_e], reverse=a.reverse,
                              spikes=a.spikes, draws=a.draws, seed=SEED, min_n=MIN_N, min_clones=MIN_C,
                              min_expected=MIN_EXPECTED,
                              call_group_shares=dict(zip([g[0] for g in CALL_GROUPS], share.round(4).tolist()))),
                wt_molecules=int(wt.size), wt_clones=int(wclone_u.size), e2_molecules=int(k2.size),
                snvs=len(snvs), snv_molecules=int(snv.size), short_clones=int(bal.short_clones.fillna(0).sum()),
                balance=flags, seconds=round(time.time() - t0))
    json.dump(info, open(out / 'summary.json', 'w'), indent=1, default=float)
    print(summ.to_string())
    print(ms.to_string())
    print(json.dumps(info, default=float))


if __name__ == '__main__':
    main()

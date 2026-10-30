"""Step A8 (part 3): reproducibility of WT profiles between clone halves and between runs (rules.json call_confidence).

Per run: emissions rebuilt from the run's own A6 WT molecules (fold A = CRC32(name) % 2 == 0; p_a from M, p_p from
U, Jeffreys (k + 0.5) / (n + 1)), then every C molecule is called (plasmid-LLR scan, min_llr 5, drawdown 25,
min_opps 3). Two endpoints per element position, standardised over 0.05 distal strata with fixed weights:
  unmarked   caller-free: fraction of C molecules with the adenine unmarked (A/T positions only)
  cov_<bin>  fraction of C molecules with an LLR call of that size bin covering the position
Clone = exact barcode (unreadable barcode: the molecule alone). Halves = CRC32(barcode) % 2. Uncertainty: Poisson
bootstrap over clones (200 draws). z = difference / sqrt(se1^2 + se2^2).

Usage: python tools/analysis_v2/a8_repro.py --run HepG2_run3=<a6 dir>:<ids dir>:1234:1445 \
           --run HepG2_run1=...:1234:1445 --run HEK293T_run1=...:3236:3447:rev \
           --weights HepG2_run1=HepG2_run3 --compare HepG2_run1:HepG2_run3 --output <new dir>
"""
import argparse
import importlib.util
import json
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from numba import njit
from scipy import sparse

_spec = importlib.util.spec_from_file_location('a7', Path(__file__).with_name('a7_snv_edits.py'))
a7 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a7)
_scan = a7._scan

N_BOOT = 200
SEED = 20261029
WIDTH, OLIGO0 = 211, 21
N_STRATA = int(round(1 / a7.STRATUM_WIDTH))


def emissions(X, distal, fold):
    """Per column (p_a, p_p) from fold-A M and U molecules, Jeffreys pseudocount."""
    cls = a7.classes(distal)
    out = []
    for c in ('M', 'U'):
        sel = (cls == c) & (fold == 0)
        Y = X[sel]
        k = (Y == 1).sum(axis=0)
        n = (Y >= 0).sum(axis=0)
        out.append((k + 0.5) / (n + 1.0))
    return out[0], out[1], int(((cls == 'M') & (fold == 0)).sum()), int(((cls == 'U') & (fold == 0)).sum())


@njit
def call_bits(X, hit, miss, pos, el_s, el_e, reverse, b_lo, b_hi):
    """Per molecule x size bin x element position (oligo order): covered by a call of that bin."""
    n, P = X.shape
    nb = b_lo.size
    out = np.zeros((n, nb, el_e - el_s), np.bool_)
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
        s, e, _l, _o, _la, _ra = _scan(step, ih, ii, lo, hi, a7.MIN_LLR, a7.MIN_OPPS, a7.DRAWDOWN)
        for c in range(s.size):
            a = pos[s[c]]
            b = pos[e[c] - 1] + 1
            for q in range(nb):
                if b_lo[q] <= b - a <= b_hi[q]:
                    for rr in range(max(a, el_s), min(b, el_e)):
                        out[r, q, (el_e - 1 - rr) if reverse else (rr - el_s)] = True
    return out


def clone_ids(barcodes, name_crc):
    """Integer clone id per molecule (exact barcode; unreadable = its own clone) and its half."""
    bc = np.asarray(barcodes)
    readable = bc != ''
    key = np.where(readable, bc, np.char.add('#', name_crc.astype(str)))
    uniq, inv = np.unique(key, return_inverse=True)
    half_u = np.array([zlib.crc32(k.encode()) % 2 if not k.startswith('#') else int(k[1:]) % 2 for k in uniq])
    return inv, half_u[inv]


def standardised(num, den, strata, clones, weights, n_boot, rng):
    """Point estimate and bootstrap draws of sum_s w_s * (sum num / sum den) within stratum s.
    num, den: molecules x features (den may be 1-D for a shared denominator)."""
    uc, ci = np.unique(clones, return_inverse=True)
    W = rng.poisson(1.0, size=(n_boot, uc.size)).astype(np.float32)
    F = num.shape[1]
    est = np.zeros(F)
    boot = np.zeros((n_boot, F))
    wsum = 0.0
    for s in range(N_STRATA):
        sel = np.flatnonzero(strata == s)
        if sel.size == 0 or weights[s] == 0:
            continue
        A = sparse.csr_matrix((np.ones(sel.size, np.float32), (ci[sel], np.arange(sel.size))), shape=(uc.size, sel.size))
        cn = np.asarray(A @ num[sel].astype(np.float32))
        cd = np.asarray(A @ (den[sel].astype(np.float32) if den.ndim == 2 else den[sel, None].astype(np.float32)))
        with np.errstate(invalid='ignore', divide='ignore'):
            est += weights[s] * cn.sum(0) / cd.sum(0)
            boot += weights[s] * (W @ cn) / (W @ cd)
        wsum += weights[s]
    return est / wsum, boot / wsum


def load_run(spec):
    label, rest = spec.split('=', 1)
    parts = rest.split(':')
    return label, dict(a6=parts[0], ids=parts[1], element=(int(parts[2]), int(parts[3])),
                       reverse=len(parts) > 4 and parts[4] == 'rev')


def profile(label, cfg, out):
    """Per-run C-molecule matrices: unmarked (num, den) at element A/T and call bits; strata, clones, halves."""
    pos, X, distal, _, _ = a7.a6x.load(cfg['a6'])
    ids = np.load(Path(cfg['ids']) / 'ids.npz')
    if ids['name_crc'].size != distal.size:
        raise ValueError(f'{label}: ids do not line up with A6 molecules')
    fold = (ids['name_crc'] % 2).astype(int)
    p_a, p_p, n_m, n_u = emissions(X, distal, fold)
    hit, miss = np.log(p_p / p_a), np.log((1 - p_p) / (1 - p_a))
    el_s, el_e = cfg['element']
    cls = a7.classes(distal)
    C = np.flatnonzero(cls == 'C')
    XC = np.ascontiguousarray(X[C])
    bits = call_bits(XC, hit, miss, pos.astype(np.int64), el_s, el_e, cfg['reverse'],
                     np.array([b[1] for b in a7.BINS], np.int64), np.array([b[2] for b in a7.BINS], np.int64))
    el_cols = np.flatnonzero((pos >= el_s) & (pos < el_e))
    el_idx = (el_e - 1 - pos[el_cols]) if cfg['reverse'] else (pos[el_cols] - el_s)
    o = np.argsort(el_idx)
    el_cols, el_idx = el_cols[o], el_idx[o]
    V = XC[:, el_cols]
    clones, half = clone_ids(ids['barcode'][C], ids['name_crc'][C])
    strata = np.minimum((distal[C] / a7.STRATUM_WIDTH).astype(int), N_STRATA - 1)
    pd.DataFrame(dict(column=el_cols, ref_pos0=pos[el_cols], index0=el_idx, p_a=p_a[el_cols], p_p=p_p[el_cols])).to_csv(
        out / f'{label}_emissions_element.tsv', sep='\t', index=False)
    info = dict(molecules_C=int(C.size), emission_M_foldA=n_m, emission_U_foldA=n_u, clones=int(np.unique(clones).size),
                readable_barcode=float((ids['barcode'][C] != '').mean()), half_sizes=[int((half == h).sum()) for h in (0, 1)])
    return dict(unm_num=(V == 0), unm_den=(V >= 0), bits=bits.reshape(C.size, -1), strata=strata, clones=clones,
                half=half, el_idx=el_idx, info=info)


def estimates(P, rows, weights, rng):
    a, ab = standardised(P['unm_num'][rows], P['unm_den'][rows], P['strata'][rows], P['clones'][rows], weights,
                         N_BOOT, rng)
    b, bb = standardised(P['bits'][rows], np.ones(rows.size), P['strata'][rows], P['clones'][rows], weights, N_BOOT,
                         rng)
    return np.concatenate([a, b]), np.concatenate([ab, bb], axis=1)


def features(P):
    f = [('unmarked', int(i)) for i in P['el_idx']]
    f += [(f'cov_{b[0]}', j) for b in a7.BINS for j in range(WIDTH)]
    return pd.DataFrame(f, columns=['endpoint', 'index0'])


def compare(name, f, e1, b1, e2, b2):
    se1, se2 = np.nanstd(b1, axis=0, ddof=1), np.nanstd(b2, axis=0, ddof=1)
    d = e1 - e2
    z = d / np.sqrt(se1 ** 2 + se2 ** 2)
    t = f.copy()
    t['comparison'] = name
    t['oligo_pos1'] = t.index0 + OLIGO0
    t['est1'], t['est2'], t['se1'], t['se2'], t['diff'], t['z'] = e1, e2, se1, se2, d, z
    s = []
    for ep, g in t.groupby('endpoint', sort=False):
        ok = np.isfinite(g.z) & (g.se1 > 0) & (g.se2 > 0)
        g = g[ok]
        s.append(dict(comparison=name, endpoint=ep, positions=len(g), r=float(np.corrcoef(g.est1, g.est2)[0, 1]),
                      mean1=g.est1.mean(), mean2=g.est2.mean(), mean_diff=g['diff'].mean(),
                      mean_abs_diff=g['diff'].abs().mean(), sd_z=float(g.z.std()),
                      frac_abs_z_gt_1p96=float((g.z.abs() > 1.96).mean())))
    return t, pd.DataFrame(s)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', action='append', required=True)
    ap.add_argument('--weights', action='append', default=[], help='RUN=REF: use REF\'s C stratum distribution')
    ap.add_argument('--compare', action='append', default=[], help='RUN1:RUN2 (all molecules of each)')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    t0 = time.time()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    runs = dict(load_run(r) for r in a.run)
    P = {}
    for label, cfg in runs.items():
        P[label] = profile(label, cfg, out)
        print(label, json.dumps(P[label]['info']), f'{time.time() - t0:.0f} s', flush=True)
    wref = dict(w.split('=') for w in a.weights)
    W = {}
    for label in runs:
        ref = wref.get(label, label)
        W[label] = np.bincount(P[ref]['strata'], minlength=N_STRATA) / P[ref]['strata'].size
    rng = np.random.default_rng(SEED)
    tables, summaries, full = [], [], {}
    for label in runs:
        f = features(P[label])
        n = P[label]['strata'].size
        h = P[label]['half']
        e1, b1 = estimates(P[label], np.flatnonzero(h == 0), W[label], rng)
        e2, b2 = estimates(P[label], np.flatnonzero(h == 1), W[label], rng)
        t, s = compare(f'{label} halves', f, e1, b1, e2, b2)
        tables.append(t)
        summaries.append(s)
        full[label] = estimates(P[label], np.arange(n), W[label], rng)
        print(f'{label} halves done, {time.time() - t0:.0f} s', flush=True)
    for c in a.compare:
        r1, r2 = c.split(':')
        if not np.array_equal(W[r1], W[r2]):
            raise ValueError(f'{c}: runs must share stratum weights (--weights)')
        t, s = compare(f'{r1} vs {r2}', features(P[r1]), *full[r1], *full[r2])
        tables.append(t)
        summaries.append(s)
    pd.concat(tables).to_csv(out / 'positions.tsv.gz', sep='\t', index=False)
    S = pd.concat(summaries)
    S.to_csv(out / 'summary.tsv', sep='\t', index=False)
    json.dump(dict(runs={k: P[k]['info'] for k in runs}, weights={k: W[k].round(4).tolist() for k in runs},
                   n_boot=N_BOOT, seed=SEED, seconds=round(time.time() - t0)), open(out / 'summary.json', 'w'), indent=1)
    pd.set_option('display.width', 220)
    print(S.round(4).to_string(index=False))


if __name__ == '__main__':
    main()

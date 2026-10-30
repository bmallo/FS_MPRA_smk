"""Step A6 (part 2): what the in-sample M and U classes can and cannot tell us (settings: rules.json m_u_proxies).

Reads a6_extract.py outputs (WT molecules, >= 10 passes) and writes, per run:
  labeling.tsv          distal m6A fraction histogram (0.01 bins)
  profile.tsv           per element adenine: m6A rate by distal bin, and that rate / the bin's mean distal rate
  clustering.tsv        per A/T position and class: fraction of molecules where the adenine is unmethylated inside
                        an m6A-free span >= L bp, observed vs a within-stratum permutation null, and the excess
  tiles.tsv             the same averaged over 211-bp tiles of the whole plasmid (element = one tile)
and across runs: llr_weight.tsv (HepG2 run-3 emission table) and m_calls.tsv (existing LLR calls in M and C).

Usage: python tools/analysis_v2/a6_proxies.py --run HepG2_run3=<a6 dir>:1234:1445 [--run ...] \
           --naked naked_pEFGFP=<a6 dir> --emissions <emissions.tsv.gz> --llr-dir <llr_vs_fiberhmm dir> \
           --output <new dir>
"""
import argparse
import importlib.util
import json
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('a6x', Path(__file__).with_name('a6_extract.py'))
a6x = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a6x)

STRATUM_WIDTH = 0.05
U_MAX, M_MIN = 0.01, 0.60
SPANS = [20, 80]
N_PERM = 5
TILE = 211
PROFILE_BINS = [0.01, 0.2, 0.4, 0.6, 0.8, 1.0001]
F_BOUND = [0.1, 0.25, 0.5]
CALL_MIN_LLR, CALL_LEN = 5.0, (20, 80)
CHUNK = 20000


def strata(distal):
    return np.minimum((distal / STRATUM_WIDTH).astype(int), int(round(1 / STRATUM_WIDTH)) - 1)


def classes(distal):
    return np.where(distal <= U_MAX, 'U', np.where(distal >= M_MIN, 'M', 'C'))


def span_cover(X, pos, min_len):
    """Bool matrix: the adenine is unmethylated (0) and lies in an m6A-free span >= min_len bp, bounded on both
    sides by methylated adenines of the same molecule (-1 = not observed, skipped)."""
    big = np.int64(1) << 40
    m = X == 1
    prev = np.maximum.accumulate(np.where(m, pos, -1), axis=1)
    nxt = np.minimum.accumulate(np.where(m, pos, big)[:, ::-1], axis=1)[:, ::-1]
    return (X == 0) & (prev >= 0) & (nxt < big) & (nxt - prev - 1 >= min_len)


def permute_within(X, groups, rng):
    """Copy of X with each column's values shuffled among rows of the same group."""
    Y = np.empty_like(X)
    for g in np.unique(groups):
        idx = np.flatnonzero(groups == g)
        Y[idx] = rng.permuted(X[idx], axis=0)
    return Y


def cover_rate(X, pos, rows, min_len):
    """Fraction of the selected rows covered at each position (chunked to bound memory)."""
    tot = np.zeros(X.shape[1])
    for s in range(0, rows.size, CHUNK):
        tot += span_cover(X[rows[s:s + CHUNK]], pos, min_len).sum(axis=0)
    return tot / max(rows.size, 1)


def clustering(pos, X, distal, cls_names, seed):
    """Per position x class x span length: observed and null coverage (null = permutations within strata)."""
    cls = classes(distal) if cls_names is None else np.full(distal.size, cls_names)
    groups = np.char.add(cls.astype(str), strata(distal).astype(str))     # strata never cross a class boundary
    rng = np.random.default_rng(seed)
    perms = [permute_within(X, groups, rng) for _ in range(N_PERM)]
    out = []
    for c in [k for k in ['M', 'C', 'all'] if k in set(cls)]:
        rows = np.flatnonzero(cls == c)
        for L in SPANS:
            obs = cover_rate(X, pos, rows, L)
            null = np.mean([cover_rate(Y, pos, rows, L) for Y in perms], axis=0)
            unm = ((X[rows] == 0).sum(axis=0) / np.maximum((X[rows] >= 0).sum(axis=0), 1))
            out.append(pd.DataFrame(dict(cls=c, span=L, ref_pos0=pos, n_molecules=rows.size, unmethylated=unm,
                                         observed=obs, null=null, excess=obs - null)))
    return pd.concat(out, ignore_index=True)


def tiles(L, element):
    """[start, end) tiles of TILE bp covering the reference, aligned so the element is one tile."""
    s0 = element[0] % TILE if element else 0
    edges = list(range(s0, L - TILE + 1, TILE))
    return [(s, s + TILE) for s in edges]


def tile_table(cl, L, element):
    rows = []
    for (s, e) in tiles(L, element):
        sub = cl[(cl.ref_pos0 >= s) & (cl.ref_pos0 < e)]
        for (c, span), g in sub.groupby(['cls', 'span']):
            rows.append(dict(cls=c, span=span, start0=s, end0=e, n_at=len(g),
                             is_element=bool(element and s == element[0] and e == element[1]),
                             observed=g.observed.mean(), null=g.null.mean(), excess=g.excess.mean(),
                             ratio=g.observed.mean() / g.null.mean() if g.null.mean() > 0 else np.nan))
    return pd.DataFrame(rows)


def profile(pos, X, distal, element):
    el = (pos >= element[0]) & (pos < element[1])
    b = np.digitize(distal, PROFILE_BINS) - 1
    rows = []
    for k in range(len(PROFILE_BINS) - 1):
        sel = b == k
        if not sel.any():
            continue
        Xe = X[sel][:, el]
        rate = (Xe == 1).sum(axis=0) / np.maximum((Xe >= 0).sum(axis=0), 1)
        mean_d = distal[sel].mean()
        lab = f'{PROFILE_BINS[k]:.2f}-{min(PROFILE_BINS[k + 1], 1):.2f}'
        rows.append(pd.DataFrame(dict(distal_bin=lab, n_molecules=int(sel.sum()), mean_distal=mean_d,
                                      ref_pos0=pos[el], rate=rate, rate_rel=rate / mean_d)))
    return pd.concat(rows, ignore_index=True)


def llr_weight(path, element):
    t = pd.read_csv(path, sep='\t')
    t = t[(t.fold == 'A') & (t.ref_pos > element[0]) & (t.ref_pos <= element[1])].copy()
    t['ref_pos0'] = t.ref_pos - 1
    t['weight'] = np.log((1 - t.p_p) / (1 - t.p_a))
    for f in F_BOUND:
        t[f'weight_bound_{f}'] = np.log((1 - t.p_p) / (1 - (1 - f) * t.p_a))
    return t[['ref_pos0', 'strand', 'context', 'p_a', 'p_p', 'weight'] + [f'weight_bound_{f}' for f in F_BOUND]]


def m_calls(llr_dir, element):
    """Per element position: fraction of molecules covered by a 20-80 bp call at min_llr 5 (M: held-out fold)."""
    d = Path(llr_dir)
    rows = []
    for sample, fold_only in [('wt_meth', 1), ('wt_chrom', None)]:
        c = pd.read_parquet(d / f'results/calls/{sample}.llr_global_foldA_pp0_dd25.parquet')
        c = c[(c.min_llr == CALL_MIN_LLR) & c.length.between(*CALL_LEN)]
        ids = [l.strip() for l in open(d / f'config/{sample}.read_ids.txt') if l.strip()]
        if fold_only is not None:
            ids = [i for i in ids if zlib.crc32(i.encode()) % 2 == fold_only]
            c = c[c.fold == fold_only]
        n = len(ids)
        cov = np.zeros(element[1] - element[0])
        for s, e in zip(c.start.to_numpy(), c.end.to_numpy()):     # calls are 1-based inclusive
            a, b = max(s - 1, element[0]), min(e, element[1])
            if a < b:
                cov[a - element[0]:b - element[0]] += 1
        rows.append(pd.DataFrame(dict(sample=sample, n_molecules=n, ref_pos0=np.arange(*element), covered=cov / n)))
    return pd.concat(rows, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', action='append', default=[], help='LABEL=<a6 dir>:start0:end0 (element)')
    ap.add_argument('--naked', help='LABEL=<a6 dir> (no element; all molecules one class)')
    ap.add_argument('--emissions', required=True)
    ap.add_argument('--llr-dir', required=True)
    ap.add_argument('--hepg2-element', default='1234:1445')
    ap.add_argument('--seed', type=int, default=20261029)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    summary = {}
    jobs = [(k, v.split(':')[0], tuple(map(int, v.split(':')[1:]))) for k, v in (r.split('=', 1) for r in a.run)]
    if a.naked:
        k, v = a.naked.split('=', 1)
        jobs.append((k, v, None))
    for label, d, element in jobs:
        pos, X, distal, passes, strand = a6x.load(d)
        L = json.load(open(Path(d) / 'extract.json'))['reference_length']
        o = out / label
        o.mkdir()
        h, edges = np.histogram(distal, bins=np.linspace(0, 1, 101))
        pd.DataFrame(dict(lo=edges[:-1], hi=edges[1:], n=h)).to_csv(o / 'labeling.tsv', sep='\t', index=False)
        cl = clustering(pos, X, distal, None if element else 'all', a.seed)
        cl.to_csv(o / 'clustering.tsv.gz', sep='\t', index=False)
        tt = tile_table(cl, L, element)
        tt.to_csv(o / 'tiles.tsv', sep='\t', index=False)
        s = dict(molecules=int(distal.size), classes={c: int(n) for c, n in zip(*np.unique(classes(distal),
                                                                                       return_counts=True))})
        if element:
            profile(pos, X, distal, element).to_csv(o / 'profile.tsv', sep='\t', index=False)
            for (c, span), g in tt.groupby(['cls', 'span']):
                el = g[g.is_element]
                s[f'{c}_span{span}'] = dict(element_excess=float(el.excess.iloc[0]), element_ratio=float(el.ratio.iloc[0]),
                                            other_tiles_excess_median=float(g[~g.is_element].excess.median()),
                                            other_tiles_excess_max=float(g[~g.is_element].excess.max()),
                                            element_rank_desc=int((g.excess > el.excess.iloc[0]).sum() + 1),
                                            n_tiles=len(g))
        else:
            for (c, span), g in tt.groupby(['cls', 'span']):
                s[f'{c}_span{span}'] = dict(tiles_excess_median=float(g.excess.median()),
                                            tiles_excess_max=float(g.excess.max()), n_tiles=len(g))
        summary[label] = s
        print(label, json.dumps(s), flush=True)
    hep = tuple(map(int, a.hepg2_element.split(':')))
    llr_weight(a.emissions, hep).to_csv(out / 'llr_weight.tsv', sep='\t', index=False)
    m_calls(a.llr_dir, hep).to_csv(out / 'm_calls.tsv', sep='\t', index=False)
    json.dump(dict(summary=summary, settings=dict(stratum_width=STRATUM_WIDTH, spans=SPANS, n_perm=N_PERM, tile=TILE,
                                                  seed=a.seed, call_min_llr=CALL_MIN_LLR, call_len=CALL_LEN)),
              open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

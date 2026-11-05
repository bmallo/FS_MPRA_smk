"""R01 Figure 3 panel b, HEK293T row: plasmid Fiber-seq accessibility change per SNV in HEK293T run 1, computed as
the HepG2 row (repeat3.hepg2) is (settings: fig3_rules.json layout_draft13, written before this was run).

Molecules with labeling >= 0.2, all clones. SNVs with >= 100 such molecules. Endpoint: molecule has a confident
(P_chance < 1e-3) >= 81-bp plasmid-LLR call with midpoint in the element (HEK293T fold-A emissions). The SNV's A7
primary mask (HepG2 catalog mapped to the reversed HEK293T element, as A10 motif mode) is set unobserved in the SNV's
and WT molecules before calling. Writes per SNV: raw variant and WT shares, accessibility change
(-100 x difference, the figure's value) and the labeling-standardised A9 contrast (D, se, p; n, C).
Writes nucleosome_hek.tsv and summary.json to --output.
Usage: python tools/r01_figure/hek_fiberseq.py --root results/analysis_v2_20261028 --procs 3 --output <new dir>
"""
import argparse
import importlib.util
import json
import multiprocessing as mp
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
_spec = importlib.util.spec_from_file_location('fp_features', HERE.parent / 'analysis_v2/fp_features.py')
fpf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fpf)
a10, a9f, a7, a6x = fpf.a10, fpf.a9f, fpf.a7, fpf.a6x

EL = (3236, 3447)                          # contig LDLR, element reversed (DATA_GUIDE coordinates)
EMISSIONS = 'results/hek293t_llr_20261028/results/emissions.tsv.gz'
CATALOG = 'results/rep3_ledger_20261016/catalog/catalog.json'
VARIANTS = 'results/analysis_v2_20261028/a4_informative/visibility/variants_visibility.tsv'
FEATS = [dict(name='nucleosome', mid_lo=0, mid_hi=fpf.WIDTH - 1, len_lo=81, len_hi=10_000)]
MIN_MOLECULES = 100
_G = {}


def calls(X, cols):
    if cols.size:
        X = X.copy()
        X[:, cols] = -1
    return fpf.call_table(X, _G['hit'], _G['miss'], _G['lp'], _G['pos'], *EL, fpf.KMAX, a7.MIN_LLR, a7.MIN_OPPS,
                          a7.DRAWDOWN)


def confident(X, rows, cols):
    return (fpf.occupancy(calls(X[rows], cols), _G['distal'][rows], FEATS)[:, 0] & 2) > 0


def snv_task(v):
    g = _G
    cols = g['masks'][v]['mask'] if v in g['masks'] else np.zeros(0, np.int64)
    rv = g['snv_rows'][v]
    yw = confident(g['X'], g['wt2'], cols) if cols.size else g['yw0']
    yv = confident(g['X'], rv, cols)
    Yw = yw[:, None].astype(float)
    Sy, So = a10.stratum_sums(Yw, g['st'][g['wt2']])
    p = np.divide(Sy, So, out=np.full_like(Sy, np.nan), where=So > 0)
    Y = yv[:, None]
    r = a9f.contrast(Y, np.ones_like(Y), g['st'][rv], g['clone'][rv], p, So, g['ref'])
    clones = pd.Series(g['clone'][rv]).value_counts()
    return dict(snv=v, masked_adenines=int(cols.size), in_catalog=v in g['masks'], n=int(rv.size),
                clones=int(clones.size), top_clone_share=float(clones.iloc[0] / rv.size),
                occ_variant=float(yv.mean()), wt_occ=float(yw.mean()),
                acc_change_pts=float(-100 * (yv.mean() - yw.mean())),
                D=float(r['D'][0]), se=float(r['se'][0]), p=float(r['p'][0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    t0 = time.time()
    ext = Path(a.root) / 'a9_statistics/extract/HEK293T_run1'
    pos, X, distal, _, _ = a6x.load(ext)
    meta = np.load(ext / 'meta.npz')
    label, cls, clone = meta['label'], meta['cls'], meta['clone']
    hit, miss = a7.emission_steps(EMISSIONS, pos)
    t = pd.read_csv(EMISSIONS, sep='\t')
    t = t[t.fold == 'A']
    pa = dict(zip(t.ref_pos - 1, t.p_a))
    lp = np.log(1 - np.clip([pa.get(int(q), 0.0) for q in pos], 0, 0.9999))
    k2 = np.flatnonzero(distal >= a9f.E2_MIN)
    wt2 = k2[cls[k2] == 'WT']
    snv2 = k2[cls[k2] == 'SNV']
    n_snv = pd.Series(label[snv2]).value_counts()
    snvs = sorted(n_snv[n_snv >= MIN_MOLECULES].index)
    masks = a10.snv_masks(CATALOG, VARIANTS, pos, *EL, True)
    st = a9f.strata(distal)
    _G.update(X=X, hit=hit, miss=miss, lp=lp, pos=pos, distal=distal, st=st, clone=clone, wt2=wt2, masks=masks,
              snv_rows={v: snv2[label[snv2] == v] for v in snvs})
    yw0 = confident(X, wt2, np.zeros(0, np.int64))
    _G['yw0'] = yw0
    ref = a9f.Reference(yw0[:, None], np.ones((wt2.size, 1), bool), st[wt2], clone[wt2])
    _G['ref'] = SimpleNamespace(V=ref.V, W2=ref.W2, tau2=ref.tau2)
    print(f'WT {wt2.size} molecules called (confident >= 81 bp share {yw0.mean():.4f}); {len(snvs)} SNVs with '
          f'>= {MIN_MOLECULES} molecules; {time.time() - t0:.0f} s', flush=True)
    rows = []
    with mp.get_context('fork').Pool(a.procs) as pool:
        for i, r in enumerate(pool.imap_unordered(snv_task, snvs)):
            rows.append(r)
            if i % 10 == 0:
                print(f'{i + 1}/{len(snvs)}, {time.time() - t0:.0f} s', flush=True)
    res = pd.DataFrame(rows).sort_values('snv')
    res.to_csv(out / 'nucleosome_hek.tsv', sep='\t', index=False, float_format='%.5g')
    json.dump(dict(wt_molecules=int(wt2.size), wt_confident_nucleosome_share=float(yw0.mean()),
                   snvs_with_molecules=int(n_snv.size), snvs_tested=len(snvs), min_molecules=MIN_MOLECULES,
                   not_in_catalog=[v for v in snvs if v not in masks], seconds=round(time.time() - t0)),
              open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

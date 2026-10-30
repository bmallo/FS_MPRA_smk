"""Step A9 slide figures from the a9_framework tables (no recomputation). Cells are drawn separately.

Usage: python tools/analysis_v2/a9_figures.py --dir results/analysis_v2_20261028/a9_statistics
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
import matplotlib.ticker
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('a7f', Path(__file__).with_name('a7_figures.py'))
a7f = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a7f)
INK, BODY, BLUE, ORANGE, GREY = a7f.INK, a7f.BODY, a7f.BLUE, a7f.ORANGE, a7f.GREY
CELLS = [('HepG2_run3', 'HepG2 run 3'), ('HEK293T_run1', 'HEK293T')]
ENDPOINTS = [('E1', None, 'adenine unmarked (caller-free)'), ('E2', 'TF', '20-40 bp call coverage')]
BIN_NAMES = ['sub_TF', 'TF', 'PIC', 'NUC']
TESTS = [('p', 'random-effects test (primary)', BLUE), ('p_sandwich', 'clone sandwich (prespecified)', '#9A4A12'),
         ('p_naive', 'molecules as independent', ORANGE)]


def thin_qq(p):
    p = np.sort(p[np.isfinite(p)])
    m = p.size
    exp = (np.arange(1, m + 1) - 0.5) / m
    keep = np.unique(np.concatenate([np.arange(min(m, 2000)), np.geomspace(1, m, 3000).astype(int) - 1]))
    return -np.log10(exp[keep]), -np.log10(p[keep])


def qq(d, out, runs):
    fig, ax = plt.subplots(len(runs), 2, figsize=(14, 5.6 * len(runs)), squeeze=False,
                           gridspec_kw=dict(hspace=0.45, wspace=0.25))
    for i, (run, lab) in enumerate(runs):
        cal = pd.read_parquet(d / run / 'calibration.parquet')
        e2_tf = None
        for j, (ep, bn, elab) in enumerate(ENDPOINTS):
            g = cal[cal.endpoint == ep]
            if bn:
                W = 211
                g = g[(g.feature // W) == BIN_NAMES.index(bn)]
            a = ax[i, j]
            a.spines[['top', 'right']].set_visible(False)
            xmax = 0
            for col, tl, c in TESTS:
                x, y = thin_qq(g[col].to_numpy())
                a.scatter(x, y, s=8, color=c, label=tl, rasterized=True)
                xmax = max(xmax, x.max())
            a.plot([0, xmax], [0, xmax], color=GREY, lw=1)
            a.set_xlim(0, xmax * 1.05)
            a.set_ylim(0, 14)
            a.set_xlabel('expected -log10 p (no effect)')
            a.set_ylabel('observed -log10 p')
            a.set_title(f'{lab}: {elab}', loc='left', color=INK, fontsize=14)
            if i == 0 and j == 0:
                a.legend(frameon=False, fontsize=12, loc='upper left')
    fig.savefig(out / 'a9_qq.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def balance(d, out, runs):
    covs = [('smd_distal', 'labeling'), ('smd_passes', 'passes'), ('smd_strand', 'strand')]
    fig, ax = plt.subplots(len(runs), len(covs), figsize=(16, 4.2 * len(runs)), squeeze=False,
                           gridspec_kw=dict(hspace=0.55, wspace=0.25))
    for i, (run, lab) in enumerate(runs):
        b = pd.read_csv(d / run / 'balance.tsv', sep='\t')
        for j, (col, cl) in enumerate(covs):
            a = ax[i, j]
            a.spines[['top', 'right']].set_visible(False)
            real, ps = b[b.draw == -1][col].dropna(), b[b.draw >= 0][col].dropna()
            lim = np.nanpercentile(np.abs(pd.concat([real, ps])), 99.5) * 1.1
            bins = np.linspace(-lim, lim, 41)
            a.hist(ps, bins=bins, density=True, color=GREY, alpha=0.55, label='pseudo-variants (chance)')
            a.hist(real, bins=bins, density=True, histtype='step', lw=2, color=BLUE, label='real SNVs')
            a.axvline(0, color=BODY, lw=0.8)
            a.set_xlabel(f'SMD in {cl} (SNV - WT)')
            a.set_yticks([])
            a.set_title(f'{lab}: {cl}', loc='left', color=INK, fontsize=14)
            if i == 0 and j == 0:
                a.legend(frameon=False, fontsize=11, loc='upper left')
    fig.savefig(out / 'a9_balance.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def mde(d, out, runs):
    fig, ax = plt.subplots(len(runs), 2, figsize=(15, 4.8 * len(runs)), squeeze=False,
                           gridspec_kw=dict(hspace=0.5, wspace=0.25))
    for i, (run, lab) in enumerate(runs):
        m = pd.read_parquet(d / run / 'mde.parquet')
        cl = pd.read_csv(d / run / 'balance.tsv', sep='\t').query('draw == -1').set_index('snv').clones
        e1 = m[m.endpoint == 'E1'].groupby('snv').mde.median()
        tf = m[(m.endpoint == 'E2') & (m['bin'] == 'TF')].groupby('snv')[['mde', 'mde_occupancy']].median()
        a = ax[i, 0]
        a.scatter(cl.reindex(e1.index), 100 * e1, s=12, color=BLUE)
        a.set_ylabel('MDE, fraction unmarked\n(percentage points)')
        a.set_title(f'{lab}: adenine unmarked (median over adenines)', loc='left', color=INK, fontsize=13)
        a = ax[i, 1]
        a.scatter(cl.reindex(tf.index), 100 * tf.mde, s=12, color=BLUE, label='call coverage')
        a.scatter(cl.reindex(tf.index), 100 * tf.mde_occupancy, s=12, color=ORANGE, label='occupancy (÷ spike detection)')
        a.set_ylabel('MDE (percentage points)')
        a.set_title(f'{lab}: 20-40 bp calls (median over positions)', loc='left', color=INK, fontsize=13)
        if i == 0:
            a.legend(frameon=False, fontsize=11)
        for a in ax[i]:
            a.spines[['top', 'right']].set_visible(False)
            a.set_xscale('log')
            a.set_yscale('log')
            a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'))
            a.yaxis.set_minor_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}' if str(v)[0] in '25' else ''))
            a.set_xlabel('clones carrying the SNV')
    fig.savefig(out / 'a9_mde.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    out = d / 'figures'
    out.mkdir(exist_ok=False)
    a7f.a4f.style()
    runs = [(r, lab) for r, lab in CELLS if (d / r / 'calibration.parquet').exists()]
    qq(d, out, runs)
    balance(d, out, runs)
    mde(d, out, runs)


if __name__ == '__main__':
    main()

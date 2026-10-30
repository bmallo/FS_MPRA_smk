"""Step A7 slide figures from the a7_snv_edits tables (no recomputation). Cells are drawn in separate panels.

Usage: python tools/analysis_v2/a7_figures.py --dir results/analysis_v2_20261028/a7_snv_edits \
           --m-calls results/analysis_v2_20261028/a6_proxies/summary/m_calls.tsv
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('a5f', Path(__file__).with_name('a5_figures.py'))
a5f = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a5f)
a4f = a5f.a4f
INK, BODY, BLUE, ORANGE, GREY, MUTED = a5f.INK, a5f.BODY, a5f.BLUE, a5f.ORANGE, a5f.GREY, a5f.MUTED
PURPLE = '#7A4EAB'
BIN_COL = {'sub_TF': ORANGE, 'TF': BLUE, 'PIC': PURPLE, 'NUC': MUTED}
BIN_LAB = {'sub_TF': '10-19 bp', 'TF': '20-40 bp', 'PIC': '41-80 bp', 'NUC': '81+ bp'}
CELLS = [('HepG2_run3', 'HepG2 run 3'), ('HEK293T_run1', 'HEK293T')]
OLIGO0 = 21
REACH = 60


def load(d, run):
    z = np.load(d / run / 'artifacts.npz')
    e = pd.read_csv(d / run / 'edits.tsv', sep='\t')
    return z['delta'], list(z['classes']), list(z['bins']), e


def offsets(d, out, runs):
    """Mean artifact of a single dropped A/T vs distance from it, per size bin."""
    fig, ax = plt.subplots(len(runs), 2, figsize=(16, 4.4 * len(runs)), sharex=True, squeeze=False,
                           gridspec_kw=dict(hspace=0.45, wspace=0.18))
    for i, (run, lab) in enumerate(runs):
        delta, cls, bins, e = load(d, run)
        dr = e[e.kind == 'drop']
        for k, c in enumerate(cls):
            a = ax[i, k]
            a.spines[['top', 'right']].set_visible(False)
            a.axhline(0, color=GREY, lw=1)
            a.axvline(0, color=GREY, lw=1, ls=':')
            for q, b in enumerate(bins):
                acc = np.full((len(dr), 2 * REACH + 1), np.nan)
                for n, r in enumerate(dr.itertuples()):
                    j = r.centre_oligo - OLIGO0
                    for o in range(-REACH, REACH + 1):
                        if 0 <= j + o < delta.shape[-1]:
                            acc[n, o + REACH] = delta[r.edit, k, q, j + o]
                a.plot(np.arange(-REACH, REACH + 1), 100 * np.nanmean(acc, axis=0), color=BIN_COL[b], lw=2,
                       label=BIN_LAB[b])
            a.set_title(f'{lab}, {c} molecules', loc='left', color=INK, fontsize=15)
            a.set_ylabel('coverage change\n(% molecules)')
            if i == len(runs) - 1:
                a.set_xlabel('distance from the dropped A/T (bp)')
        ax[i, 1].legend(frameon=False, fontsize=11, loc='upper right', bbox_to_anchor=(1.0, 1.25), ncol=4)
    fig.savefig(out / 'a7_offsets.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def positions(d, out, runs):
    """Largest |artifact| (any bin) per element position: drop of that A/T, the variant's primary mask, +/-5 mask."""
    fig, ax = plt.subplots(len(runs), 1, figsize=(16.6, 3.9 * len(runs)), sharex=True, squeeze=False,
                           gridspec_kw=dict(hspace=0.4))
    for i, (run, lab) in enumerate(runs):
        e = pd.read_csv(d / run / 'edits.tsv', sep='\t')
        v = pd.read_csv(d / run / 'variants.tsv', sep='\t')
        a = ax[i, 0]
        a5f.shade(a)
        pm = e[e.kind == 'pm5'].sort_values('centre_oligo')
        a.plot(pm.centre_oligo, 100 * pm.C_max_abs, color=GREY, lw=1.5, label='mask +/-5 bp (sensitivity)')
        prim = v[(v.cls == 'SNV') & (v.edit_primary > 0)].merge(e[['edit', 'C_max_abs']], left_on='edit_primary',
                                                                    right_on='edit')
        pp = prim.groupby('index0').C_max_abs.max()
        a.scatter(pp.index + OLIGO0, 100 * pp.values, s=22, color=BLUE, label='primary mask (largest over SNVs)',
                  zorder=3)
        dr = e[e.kind == 'drop']
        a.scatter(dr.centre_oligo, 100 * dr.C_max_abs, s=22, facecolor='none', edgecolor=ORANGE, lw=1.4,
                  label='edited A/T dropped only', zorder=4)
        a.set_ylabel('largest change, C\n(% molecules)')
        a.set_title(lab, loc='left', color=INK, fontsize=15)
        a.set_ylim(bottom=0)
        a.legend(frameon=False, fontsize=11, loc='upper right', bbox_to_anchor=(1.0, 1.22), ncol=3)
    ax[-1, 0].set_xlabel('oligo position of the edit')
    fig.savefig(out / 'a7_positions.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def v1(d, out):
    p = d / 'HepG2_run3' / 'v1_check.tsv'
    if not p.exists():
        return
    t = pd.read_csv(p, sep='\t')
    t = t[t.inside].set_index('caller_class').reindex(['drops its base', 'neither', 'gains an A/T'])
    fig, a = plt.subplots(figsize=(11, 4.6))
    a.spines[['top', 'right']].set_visible(False)
    a.axhline(0, color=GREY, lw=1)
    x = np.arange(len(t))
    for dx, col, c, lab in [(-0.27, INK, 'v1_resid_mean_pp', 'v1 observed (variant molecules)'),
                            (0, ORANGE, 'drop_resid_mean_pp', 'A7: edited base dropped'),
                            (0.27, BLUE, 'context_resid_mean_pp', 'A7: dropped + context model')]:
        a.bar(x + dx, t[c], width=0.25, color=col, label=lab)
    a.set_xticks(x, [f'{k}\n({int(n)} pairs)' for k, n in zip(t.index, t.pairs)])
    a.set_ylabel('mean occupancy residual\n(percentage points)')
    a.set_title('HepG2 run 3: v1 (SNV, feature) pairs with the SNV inside the feature', loc='left', color=INK,
                fontsize=14)
    a.legend(frameon=False, fontsize=11, loc='lower right')
    fig.savefig(out / 'a7_v1_check.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def baseline(d, out, m_calls):
    b = pd.read_csv(d / 'HepG2_run3' / 'baseline.tsv', sep='\t')
    b = b[b['bin'] == '20-80']
    m = pd.read_csv(m_calls, sep='\t')
    m['oligo'] = m.ref_pos0 - 1234 + OLIGO0
    fig, a = plt.subplots(figsize=(16.6, 4.0))
    a5f.shade(a)
    for c, sample, col in [('C', 'wt_chrom', ORANGE), ('M', 'wt_meth', BLUE)]:
        s = b[b.cls == c]
        a.plot(s.oligo_pos1, 100 * s.coverage, color=col, lw=2, label=f'{c}: A7 re-call (A6 molecules)')
        s = m[m['sample'] == sample]
        a.plot(s.oligo, 100 * s.covered, color=col, lw=1.4, ls='--', label=f'{c}: existing calls')
    a.set_ylabel('molecules with a 20-80 bp\ncall (%)')
    a.set_xlabel('oligo position')
    a.set_ylim(bottom=0)
    a.set_title('HepG2 run 3, unedited molecules', loc='left', color=INK, fontsize=15)
    a.legend(frameon=False, fontsize=11, loc='upper right', ncol=2)
    fig.savefig(out / 'a7_baseline.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--m-calls', required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    out = d / 'figures'
    out.mkdir(exist_ok=False)
    a4f.style()
    runs = [(r, lab) for r, lab in CELLS if (d / r / 'artifacts.npz').exists()]
    offsets(d, out, runs)
    positions(d, out, runs)
    v1(d, out)
    baseline(d, out, a.m_calls)


if __name__ == '__main__':
    main()

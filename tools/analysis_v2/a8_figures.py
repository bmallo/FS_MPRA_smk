"""Step A8 slide figures from the a8_spikes and a8_repro tables (no recomputation). Cells are drawn separately.

Usage: python tools/analysis_v2/a8_figures.py --dir results/analysis_v2_20261028/a8_confidence
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('a7f', Path(__file__).with_name('a7_figures.py'))
a7f = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a7f)
a5f, a4f = a7f.a5f, a7f.a4f
INK, BODY, BLUE, ORANGE, GREY, MUTED, PURPLE = a7f.INK, a7f.BODY, a7f.BLUE, a7f.ORANGE, a7f.GREY, a7f.MUTED, a7f.PURPLE
CELLS = [('HepG2_run3', 'HepG2 run 3'), ('HEK293T_run1', 'HEK293T')]
SIZES = [10, 15, 20, 30, 40, 60, 80, 147]
GROUPS = [('C 0.01-0.20', 'C, labeling 0.01-0.2', '#E8B48A'), ('C 0.20-0.40', 'C, 0.2-0.4', ORANGE),
          ('C 0.40-0.60', 'C, 0.4-0.6', '#9A4A12'), ('M', 'M (>= 0.6)', BLUE)]


def spikes(d, run):
    p = d / 'spikes' / run / 'spikes.tsv'
    return pd.read_csv(p, sep='\t') if p.exists() else None


def heatmap(d, out, runs):
    fig, ax = plt.subplots(len(runs), 1, figsize=(16.6, 3.3 * len(runs) + 0.8), squeeze=False,
                           gridspec_kw=dict(hspace=0.45))
    for i, (run, lab) in enumerate(runs):
        t = spikes(d, run)
        t = t[(t.emissions == 'position') & (t['mode'] == 'perfect') & (t.group == 'C')]
        M = t.pivot(index='size', columns='centre_oligo', values='net').reindex(SIZES)
        a = ax[i, 0]
        im = a.imshow(100 * M.to_numpy(), aspect='auto', cmap='Blues', vmin=0, vmax=100, origin='lower',
                      extent=(M.columns.min() - 0.5, M.columns.max() + 0.5, -0.5, len(SIZES) - 0.5))
        a.set_yticks(range(len(SIZES)), [f'{s} bp' for s in SIZES])
        a.set_title(f'{lab}: C molecules, net detection of a spiked footprint', loc='left', color=INK, fontsize=15, pad=24)
        for name, (s, e) in a4f.a4.REGIONS.items():
            a.text((s + e) / 2, 1.02, a4f.SHORT[name], transform=a.get_xaxis_transform(), ha='center', va='bottom',
                   color=BODY, fontsize=12)
    ax[-1, 0].set_xlabel('oligo position of the spike centre')
    cb = fig.colorbar(im, ax=ax[:, 0].tolist(), fraction=0.02, pad=0.01)
    cb.set_label('net detection (% molecules)')
    fig.savefig(out / 'a8_heatmap.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def by_labeling(d, out, runs):
    """Per cell: already called before spiking; net detection (prespecified); detection where nothing was called
    before = net / (1 - unspiked) (post hoc). Position-specific emissions, perfect spikes."""
    fig, ax = plt.subplots(len(runs), 3, figsize=(17, 4.3 * len(runs)), squeeze=False,
                           gridspec_kw=dict(hspace=0.45, wspace=0.28, width_ratios=[0.8, 1, 1]))
    for i, (run, lab) in enumerate(runs):
        t = spikes(d, run)
        t = t[(t.emissions == 'position') & (t['mode'] == 'perfect')].copy()
        t['cond'] = t.net / (1 - t.det_unspiked)
        a = ax[i, 0]
        a.spines[['top', 'right']].set_visible(False)
        u = t[t['size'] == 20].groupby('group').det_unspiked.median()
        a.bar(range(len(GROUPS)), [100 * u[g] for g, _, _ in GROUPS], color=[c for _, _, c in GROUPS], width=0.65)
        a.set_xticks(range(len(GROUPS)), ['C low', 'C mid', 'C high', 'M'])
        a.set_ylim(0, 100)
        a.set_ylabel('molecules already called\nover the spot (%)')
        a.set_title(f'{lab}: before spiking', loc='left', color=INK, fontsize=14)
        for k, (col, ttl, ylab) in enumerate([('net', 'net detection (prespecified)', 'net detection (%)'),
                                              ('cond', 'where nothing was called before (post hoc)',
                                               'detection (%)')]):
            a = ax[i, k + 1]
            a.spines[['top', 'right']].set_visible(False)
            for g, glab, c in GROUPS:
                m = t[t.group == g].groupby('size')[col].median().reindex(SIZES)
                a.plot(range(len(SIZES)), 100 * m.to_numpy(), marker='o', color=c, lw=2, label=glab)
            a.axhline(50, color=GREY, lw=1, ls=':')
            a.set_ylim(0, 102)
            a.set_xticks(range(len(SIZES)), [str(s) for s in SIZES])
            a.set_title(ttl, loc='left', color=INK, fontsize=14)
            a.set_ylabel(ylab)
            if i == len(runs) - 1:
                a.set_xlabel('spiked footprint size (bp)')
    ax[0, 1].legend(frameon=False, fontsize=10.5, loc='center right')
    fig.savefig(out / 'a8_by_labeling.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def size_accuracy(d, out, runs):
    fig, ax = plt.subplots(1, 2, figsize=(16, 4.6), gridspec_kw=dict(wspace=0.25))
    for run, lab in runs:
        t = spikes(d, run)
        s = t[(t.emissions == 'position') & (t['mode'] == 'perfect') & (t.group == 'C')]
        g = s.groupby('size')
        col = BLUE if run.startswith('HepG2') else ORANGE
        frac = (g.net_size_ok.median() / g.net.median()).reindex(SIZES)
        ax[0].plot(range(len(SIZES)), 100 * frac.to_numpy(), marker='o', color=col, lw=2, label=lab)
        err = g.length_error_median.median().reindex(SIZES)
        ax[1].plot(range(len(SIZES)), err.to_numpy(), marker='o', color=col, lw=2, label=lab)
    for a in ax:
        a.spines[['top', 'right']].set_visible(False)
        a.set_xticks(range(len(SIZES)), [str(s) for s in SIZES])
        a.set_xlabel('spiked footprint size (bp)')
    ax[0].set_ylabel('detected calls in the\nright size class (%)')
    ax[1].set_ylabel('called length - true\n(bp, median)')
    ax[1].axhline(0, color=GREY, lw=1)
    ax[0].legend(frameon=False, fontsize=12)
    ax[0].set_title('C molecules, position-specific emissions', loc='left', color=INK, fontsize=14)
    fig.savefig(out / 'a8_size.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def floor(d, out, runs):
    """Post hoc: smallest perfect spike detected in >= 50% of M molecules (the prespecified C floor is never reached,
    because over half of C molecules are already called at any spot), next to A4's prediction (HepG2 table)."""
    fig, ax = plt.subplots(len(runs), 1, figsize=(16.6, 3.4 * len(runs)), sharex=True, squeeze=False,
                           gridspec_kw=dict(hspace=0.45))
    p4 = d / 'spikes' / 'HepG2_run3' / 'floor.tsv'
    a4 = pd.read_csv(p4, sep='\t') if p4.exists() else None
    for i, (run, lab) in enumerate(runs):
        t = spikes(d, run)
        t = t[(t.emissions == 'position') & (t['mode'] == 'perfect') & (t.group == 'M')]
        sm = t[t.net >= 0.5].groupby('centre_oligo')['size'].min()
        a = ax[i, 0]
        a5f.shade(a)
        if a4 is not None and run.startswith('HepG2'):
            a.step(a4.centre_oligo, a4.min_fp_llr_bp, where='mid', color=GREY, lw=1.6,
                   label='A4 prediction (smallest callable, from LLR weights)')
        a.scatter(sm.index, sm.values, s=18, color=BLUE, zorder=3, label='smallest spike found in >= 50% of M')
        a.set_yscale('log')
        a.set_yticks([2, 5, 10, 15, 20], ['2', '5', '10', '15', '20'])
        a.set_ylim(1.5, 30)
        a.set_ylabel('footprint size (bp)')
        a.set_title(lab, loc='left', color=INK, fontsize=15)
        a.legend(frameon=False, fontsize=11, loc='upper right', ncol=2)
    ax[-1, 0].set_xlabel('oligo position (spike centre)')
    fig.savefig(out / 'a8_floor.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def repro(d, out):
    p = d / 'repro' / 'positions.tsv.gz'
    if not p.exists():
        return
    t = pd.read_csv(p, sep='\t')
    s = pd.read_csv(d / 'repro' / 'summary.tsv', sep='\t')
    comps = list(dict.fromkeys(t.comparison))
    eps = [('unmarked', 'adenine unmarked (caller-free)'), ('cov_TF', '20-40 bp call coverage'),
           ('cov_NUC', '81+ bp call coverage')]
    fig, ax = plt.subplots(len(eps), len(comps), figsize=(4.4 * len(comps), 4.0 * len(eps)), squeeze=False,
                           gridspec_kw=dict(hspace=0.5, wspace=0.35))
    for j, c in enumerate(comps):
        for i, (ep, elab) in enumerate(eps):
            a = ax[i, j]
            a.spines[['top', 'right']].set_visible(False)
            g = t[(t.comparison == c) & (t.endpoint == ep)]
            a.scatter(100 * g.est1, 100 * g.est2, s=9, color=BLUE, alpha=0.7)
            lim = [0, 100 * max(g.est1.max(), g.est2.max()) * 1.05]
            a.plot(lim, lim, color=GREY, lw=1)
            r = s[(s.comparison == c) & (s.endpoint == ep)].iloc[0]
            a.text(0.03, 0.97, f'r = {r.r:.3f}\nSD(z) = {r.sd_z:.2f}', transform=a.transAxes, va='top', fontsize=11,
                   color=INK)
            first, second = (c.replace(' halves', ''), 'half 2') if 'halves' in c else \
                [x.replace('_run', ' run ') for x in c.split(' vs ')]
            a.set_xlabel(f'{"half 1" if "halves" in c else first} (%)', fontsize=11)
            a.set_ylabel(f'{second} (%)', fontsize=11)
            a.set_title(f'{c.replace("_run", " run ").replace("HEK293T run 1", "HEK293T")}\n{elab}', loc='left',
                        color=INK, fontsize=11)
    fig.savefig(out / 'a8_repro.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    out = d / 'figures'
    out.mkdir(exist_ok=False)
    a4f.style()
    runs = [(r, lab) for r, lab in CELLS if (d / 'spikes' / r / 'spikes.tsv').exists()]
    heatmap(d, out, runs)
    by_labeling(d, out, runs)
    size_accuracy(d, out, runs)
    floor(d, out, runs)
    repro(d, out)


if __name__ == '__main__':
    main()

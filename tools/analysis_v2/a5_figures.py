"""Step A5 slide figures from the a5_m6a_summary tables (no recomputation). Cells are drawn in separate panels.

Usage: python tools/analysis_v2/a5_figures.py --dir results/analysis_v2_20261028/a5_m6a_qc/summary
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('a4f', Path(__file__).with_name('a4_figures.py'))
a4f = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a4f)
INK, BODY, BLUE, ORANGE, GREY, MUTED, BG = a4f.INK, a4f.BODY, a4f.BLUE, a4f.ORANGE, a4f.GREY, a4f.MUTED, a4f.BG
RUN_STYLE = {'HepG2_run1': (GREY, 'HepG2 run 1 (10 min, 1.5 µL)'), 'HepG2_run2': (ORANGE, 'HepG2 run 2 (10 min, 1.5 µL)'),
             'HepG2_run3': (BLUE, 'HepG2 run 3 (5 min, 3 µL)'), 'HEK293T_run1': (INK, 'HEK293T (10 min)')}
CELLS = [('HepG2', ['HepG2_run1', 'HepG2_run2', 'HepG2_run3']), ('HEK293T', ['HEK293T_run1'])]


def shade(a):
    for name, (s, e) in a4f.a4.REGIONS.items():
        a.axvspan(s - 0.5, e + 0.5, color=a4f.SHADE[name], zorder=0, lw=0)
    a.spines[['top', 'right']].set_visible(False)


def element(d, out, runs_ok):
    t = pd.read_csv(d / 'element_positions.tsv', sep='\t')
    fig, ax = plt.subplots(2, 1, figsize=(16.6, 7.0), sharex=True, gridspec_kw=dict(hspace=0.3))
    for a, (cell, runs) in zip(ax, CELLS):
        shade(a)
        for run in runs:
            if run not in runs_ok:
                continue
            s = t[t.run == run].sort_values('oligo_pos1')
            col, lab = RUN_STYLE[run]
            a.plot(s.oligo_pos1, 100 * s.rate_U, marker='o', ms=4, lw=1.2, color=col, label=lab)
        a.set_ylabel('m6A calls in U\nmolecules (%)')
        a.set_title(cell, loc='left', color=INK, fontsize=16)
        a.set_ylim(0, 1.3 * max(100 * t[t.run.isin(runs)].rate_U.max(), 1))
        a.legend(frameon=False, fontsize=12, loc='upper left', bbox_to_anchor=(0.12, 1.0), ncol=3)
    for name, (s, e) in a4f.a4.REGIONS.items():
        ax[0].text((s + e) / 2, 1.02, a4f.SHORT[name], transform=ax[0].get_xaxis_transform(), ha='center', va='bottom',
                   color=BODY, fontsize=13)
    ax[1].set_xlabel('oligo position (A/T only)')
    fig.savefig(out / 'a5_element_fp.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def passes(d, out, runs_ok):
    t = pd.read_csv(d / 'passes.tsv', sep='\t')
    labels = ['<10', '10-14', '15-19', '20-29', '30+']
    fig, ax = plt.subplots(2, 2, figsize=(15.5, 7.6), sharex=True, gridspec_kw=dict(hspace=0.35, wspace=0.28))
    for r, (cell, runs) in enumerate(CELLS):
        for c, (col_name, ylab) in enumerate([('element_U_rate', 'U molecules: m6A\ncalls in element (%)'),
                                              ('element_M_rate', 'M molecules: m6A\ncalls in element (%)')]):
            a = ax[r, c]
            a.spines[['top', 'right']].set_visible(False)
            a.axvspan(-0.5, 0.5, color='#EDEFF2', lw=0, zorder=0)
            for run in runs:
                if run not in runs_ok:
                    continue
                s = t[t.run == run].set_index('passes').reindex(labels)
                ok = s.element_U_obs >= 2000
                colr, lab = RUN_STYLE[run]
                a.plot(np.arange(5)[ok.values], 100 * s[col_name][ok], marker='o', ms=8, lw=2, color=colr, label=lab)
            a.set_ylabel(ylab)
            a.set_ylim(bottom=0)
            if r == 0 and c == 0:
                a.text(0, 0.03, 'excluded\n(< 10 passes)', transform=a.get_xaxis_transform(), ha='center', va='bottom',
                       color=MUTED, fontsize=11)
            a.set_title(cell, loc='left', color=INK, fontsize=15)
        ax[r, 1].legend(frameon=False, fontsize=11, loc='lower right')
    for a in ax[1]:
        a.set_xticks(range(5), labels)
        a.set_xlabel('CCS passes')
    fig.savefig(out / 'a5_passes.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def gatc(d, out, runs_ok):
    t = pd.read_csv(d / 'gatc_sites.tsv', sep='\t')
    runs = [r for r in ['HepG2_run1', 'HepG2_run3', 'HEK293T_run1'] if r in runs_ok]
    groups = [('GATC', 'GATC'), ('same_3mer', 'GAT / ATC\nnot GATC'), ('isolated_other', 'other isolated\nadenines')]
    fig, ax = plt.subplots(1, len(runs), figsize=(5.4 * len(runs), 5.6), sharey=True, gridspec_kw=dict(wspace=0.12))
    rng = np.random.default_rng(0)
    for a, run in zip(np.atleast_1d(ax), runs):
        a.spines[['top', 'right']].set_visible(False)
        for cls, col, dx in [('U', ORANGE, -0.18), ('M', BLUE, 0.18)]:
            for gi, (g, _) in enumerate(groups):
                v = t[(t.run == run) & (t.cls == cls) & (t.group == g)].rate.dropna().to_numpy()
                a.scatter(gi + dx + rng.uniform(-0.07, 0.07, v.size), 100 * v, s=10, color=col, alpha=0.55, lw=0,
                          label=f'{cls} molecules' if gi == 0 else None)
                a.plot([gi + dx - 0.12, gi + dx + 0.12], [100 * np.median(v)] * 2, color=INK, lw=2.5)
        a.set_xticks(range(3), [g[1] for g in groups], fontsize=12)
        a.set_yscale('symlog', linthresh=1)
        a.set_title(RUN_STYLE[run][1], loc='left', color=INK, fontsize=14)
    np.atleast_1d(ax)[0].set_ylabel('m6A calls per site (%)')
    fig.legend(*np.atleast_1d(ax)[0].get_legend_handles_labels(), frameon=False, fontsize=13, loc='upper center',
               ncol=2, bbox_to_anchor=(0.5, 1.06), markerscale=2.5)
    fig.savefig(out / 'a5_gatc.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--runs', default='HepG2_run1,HepG2_run2,HepG2_run3,HEK293T_run1')
    a = ap.parse_args()
    d = Path(a.dir)
    out = d.parent / 'figures'
    out.mkdir(exist_ok=False)
    a4f.style()
    runs_ok = set(a.runs.split(','))
    element(d, out, runs_ok)
    passes(d, out, runs_ok)
    gatc(d, out, runs_ok)


if __name__ == '__main__':
    main()

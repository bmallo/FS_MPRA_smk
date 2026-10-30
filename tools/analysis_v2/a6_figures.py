"""Step A6 slide figures from the a6_proxies tables (no recomputation). Cells are drawn in separate panels.

Usage: python tools/analysis_v2/a6_figures.py --dir results/analysis_v2_20261028/a6_proxies/summary
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
RUN_STYLE = a5f.RUN_STYLE
ELEMENT = {'HepG2_run1': (1234, 1445), 'HepG2_run2': (1234, 1445), 'HepG2_run3': (1234, 1445),
           'HEK293T_run1': (3236, 3447)}
CELLS = a5f.CELLS


def oligo(run, ref0):
    s, e = ELEMENT[run]
    return ref0 - s + 21 if run.startswith('HepG2') else (e - 1) - ref0 + 21


def labeling(d, out, runs, naked):
    fig, ax = plt.subplots(1, 2, figsize=(15.5, 5.0), sharey=False, gridspec_kw=dict(wspace=0.22))
    for a, (cell, rr) in zip(ax, CELLS):
        a.spines[['top', 'right']].set_visible(False)
        a.axvspan(0, 0.01, color='#FBEFE4', lw=0, zorder=0)
        a.axvspan(0.6, 1.0, color='#E3EDF8', lw=0, zorder=0)
        for run in [r for r in rr if r in runs]:
            t = pd.read_csv(d / run / 'labeling.tsv', sep='\t')
            col, lab = RUN_STYLE[run]
            a.plot((t.lo + t.hi) / 2, 100 * t.n / t.n.sum(), color=col, lw=2, label=lab)
        if naked:
            t = pd.read_csv(d / naked / 'labeling.tsv', sep='\t')
            a.plot((t.lo + t.hi) / 2, 100 * t.n / t.n.sum(), color=GREY, lw=2, ls='--',
                   label='naked pEF-GFP + Hia5 (other plasmid)')
        for x, txt in [(0.005, 'U'), (0.3, 'C'), (0.8, 'M')]:
            a.text(x, 0.97, txt, transform=a.get_xaxis_transform(), ha='center', va='top', color=BODY, fontsize=15)
        a.set_xlabel('distal m6A fraction (backbone)')
        a.set_ylabel('WT molecules (% per 0.01 bin)')
        a.set_title(cell, loc='left', color=INK, fontsize=16)
        a.set_ylim(bottom=0)
        a.legend(frameon=False, fontsize=11, loc='upper left', bbox_to_anchor=(0.08, 0.92))
    fig.savefig(out / 'a6_labeling.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def tiles(d, out, runs, naked):
    """M-molecule excess per 211-bp tile at both span lengths (C is far larger and is reported on the slide)."""
    panels = [r for r in ['HepG2_run3', 'HEK293T_run1'] if r in runs] + ([naked] if naked else [])
    fig, ax = plt.subplots(len(panels), 1, figsize=(16, 2.6 * len(panels) + 0.6), gridspec_kw=dict(hspace=0.7))
    tabs = {p: pd.read_csv(d / p / 'tiles.tsv', sep='\t') for p in panels}
    lim = max(100 * t[t.cls != 'C'].excess.max() for t in tabs.values())
    for a, p in zip(np.atleast_1d(ax), panels):
        t = tabs[p]
        t = t[t.cls == ('all' if p == naked else 'M')]
        a.spines[['top', 'right']].set_visible(False)
        a.axhline(0, color=GREY, lw=1)
        for span, col, dx in [(20, '#9DBBE0', -0.2), (80, BLUE, 0.2)]:
            s = t[t.span == span].sort_values('start0')
            x = np.arange(len(s))
            a.bar(x + dx, 100 * s.excess, width=0.38, color=col, label=f'spans >= {span} bp')
        s = t[t.span == 20].sort_values('start0')
        for i in np.flatnonzero(s.is_element.to_numpy()):
            a.axvspan(i - 0.5, i + 0.5, color='#EDEFF2', zorder=0, lw=0)
            a.text(i, 1.0, 'element', transform=a.get_xaxis_transform(), ha='center', va='bottom', color=INK,
                   fontsize=12)
        a.set_xticks(np.arange(len(s))[::2], [f'{v}' for v in s.start0.to_numpy()[::2]], fontsize=10)
        a.set_ylim(-0.3, 1.1 * lim)
        title = (RUN_STYLE[p][1] + ', M molecules') if p in RUN_STYLE else \
            'naked pEF-GFP + Hia5, all molecules (other plasmid, ft predict-m6a)'
        a.set_title(title, loc='left', color=INK, fontsize=14, pad=14)
        a.set_ylabel('excess\n(% molecules)')
        a.legend(frameon=False, fontsize=11, loc='upper right', ncol=2)
    np.atleast_1d(ax)[-1].set_xlabel('211-bp tile (start, reference bp); excess = observed - permuted')
    fig.savefig(out / 'a6_tiles.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def element(d, out, runs):
    fig, ax = plt.subplots(3, 1, figsize=(16.6, 9.6), sharex=True, gridspec_kw=dict(hspace=0.35))
    for a, (cell, rr) in zip(ax[:2], CELLS):
        a5f.shade(a)
        a.axhline(0, color=GREY, lw=1)
        for run in [r for r in rr if r in runs and r != 'HepG2_run2']:
            t = pd.read_csv(d / run / 'clustering.tsv.gz', sep='\t')
            s, e = ELEMENT[run]
            t = t[(t.span == 20) & (t.cls == 'M') & (t.ref_pos0 >= s) & (t.ref_pos0 < e)].copy()
            t['oligo'] = [oligo(run, p) for p in t.ref_pos0]
            t = t.sort_values('oligo')
            col, lab = RUN_STYLE[run]
            a.plot(t.oligo, 100 * t.excess, marker='o', ms=4, lw=1.2, color=col, label=lab)
        a.set_ylabel('M molecules: excess\n(% molecules)')
        a.set_title(cell, loc='left', color=INK, fontsize=16)
        a.legend(frameon=False, fontsize=12, loc='upper right', ncol=1)
    w = pd.read_csv(d / 'llr_weight.tsv', sep='\t')
    w['oligo'] = w.ref_pos0 - 1234 + 21
    a = ax[2]
    a5f.shade(a)
    for col_name, col, lab in [('weight', INK, 'M emissions as built'), ('weight_bound_0.25', ORANGE, 'if 25% of M bound'),
                               ('weight_bound_0.5', BLUE, 'if 50% of M bound')]:
        a.plot(w.oligo, w[col_name], marker='o', ms=3.5, lw=1.0, color=col, label=lab)
    a.set_ylabel('LLR weight of one\nunmethylated adenine')
    a.set_title('HepG2 run 3 emission table', loc='left', color=INK, fontsize=16)
    a.legend(frameon=False, fontsize=12, loc='lower right', bbox_to_anchor=(1.0, 1.0), ncol=3)
    a.set_ylim(bottom=0)
    for name, (s, e) in a4f.a4.REGIONS.items():
        ax[0].text((s + e) / 2, 1.02, a4f.SHORT[name], transform=ax[0].get_xaxis_transform(), ha='center',
                   va='bottom', color=BODY, fontsize=13)
    ax[2].set_xlabel('oligo position (A/T only)')
    fig.savefig(out / 'a6_element.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def m_calls(d, out):
    t = pd.read_csv(d / 'm_calls.tsv', sep='\t')
    t['oligo'] = t.ref_pos0 - 1234 + 21
    fig, a = plt.subplots(figsize=(16.6, 4.2))
    a5f.shade(a)
    for sample, col, lab in [('wt_chrom', ORANGE, 'C molecules'), ('wt_meth', BLUE, 'M molecules (held-out half)')]:
        s = t[t['sample'] == sample]
        a.plot(s.oligo, 100 * s.covered, color=col, lw=2, label=f'{lab}, n = {s.n_molecules.iloc[0]:,}')
    a.set_ylabel('molecules with a 20-80 bp\nLLR call (%)')
    a.set_xlabel('oligo position')
    a.set_ylim(bottom=0)
    a.set_title('HepG2 run 3, existing LLR calls (min LLR 5)', loc='left', color=INK, fontsize=16)
    a.legend(frameon=False, fontsize=12, loc='upper right')
    fig.savefig(out / 'a6_m_calls.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--runs', default='HepG2_run1,HepG2_run2,HepG2_run3,HEK293T_run1')
    ap.add_argument('--naked', default='naked_pEFGFP')
    a = ap.parse_args()
    d = Path(a.dir)
    out = d.parent / 'figures'
    out.mkdir(exist_ok=False)
    a4f.style()
    runs = [r for r in a.runs.split(',') if (d / r).exists()]
    naked = a.naked if (d / a.naked).exists() else None
    labeling(d, out, runs, naked)
    tiles(d, out, runs, naked)
    element(d, out, runs)
    m_calls(d, out)


if __name__ == '__main__':
    main()

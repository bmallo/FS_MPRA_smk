"""Step A7 explainer figure: two real WT molecules (HepG2 run 3, C) re-called with one element A/T dropped.

For each molecule, the adenines (marked / unmarked), the caller's running score and its calls are drawn before and
after the edit. Molecules are given by their row in the A6 extraction (molecules.npz) (found by scanning for a merge and a lost call).

Usage: python tools/analysis_v2/a7_example.py --extract <a6 HepG2_run3 dir> --emissions <emissions.tsv.gz> \
           --example merge:1090:1350 --example lost:8:1333 --output <png>
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

_spec = importlib.util.spec_from_file_location('a7', Path(__file__).with_name('a7_snv_edits.py'))
a7 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a7)
_spec = importlib.util.spec_from_file_location('a7f', Path(__file__).with_name('a7_figures.py'))
a7f = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a7f)
INK, BODY, BLUE, ORANGE, GREY = a7f.INK, a7f.BODY, a7f.BLUE, a7f.ORANGE, a7f.GREY
EL0 = 1234
HALF = 45
TITLES = {'merge': 'Dropping a marked adenine: two calls merge',
          'lost': 'Dropping an unmarked adenine: a call is lost'}


def running(x, hit, miss):
    """The scan's running score at every A/T column (0 outside a segment), as in llr/kadane.py."""
    out = np.zeros(x.size)
    cur, run, peak = False, 0.0, 0.0
    for j, v in enumerate(x):
        s = hit[j] if v == 1 else (miss[j] if v == 0 else 0.0)
        if not cur:
            if s > 0:
                cur, run, peak = True, s, s
            out[j] = run if cur else 0.0
            continue
        run += s
        peak = max(peak, run)
        if run <= 0 or peak - run >= a7.DRAWDOWN:
            cur, run, peak = False, 0.0, 0.0
        out[j] = run
    return out


def calls(x, pos, hit, miss):
    ii, ih = x >= 0, x == 1
    step = np.where(x == 1, hit, np.where(x == 0, miss, 0.0))
    o = np.flatnonzero(ii)
    s, e, *_ = a7.kadane.scan(step, ih, ii, o[0], o[-1] + 1, a7.MIN_LLR, a7.MIN_OPPS, a7.DRAWDOWN)
    return [(int(pos[a]), int(pos[b - 1]) + 1) for a, b in zip(s, e)]


def panel(a, x, pos, hit, miss, centre, drop, edited):
    w = (pos >= centre - HALF) & (pos <= centre + HALF)
    o = lambda r: r - EL0 + 21
    a.spines[['top', 'right']].set_visible(False)
    for s, e in calls(x, pos, hit, miss):
        if e > centre - HALF and s < centre + HALF:
            a.axvspan(o(max(s, centre - HALF)) - 0.5, o(min(e, centre + HALF + 1)) - 0.5, color='#FBE3CC', lw=0,
                      zorder=0)
            a.text(o((max(s, centre - HALF) + min(e, centre + HALF)) / 2), 1.02, f'call {e - s} bp',
                   transform=a.get_xaxis_transform(), ha='center', va='bottom', color=INK, fontsize=12)
    r = running(x, hit, miss)
    a.step(o(pos[w]), r[w], where='post', color=BLUE, lw=1.8)
    a.axhline(a7.MIN_LLR, color=GREY, lw=1, ls='--')
    a.text(o(centre + HALF), a7.MIN_LLR, ' threshold 5', va='bottom', ha='right', color=BODY, fontsize=11)
    top = max(r[w].max(), a7.MIN_LLR) * 1.15 + 1
    yb = -0.12 * top
    for p, v in zip(pos[w], x[w]):
        if p == drop and edited:
            a.scatter(o(p), yb, marker='x', s=70, color=ORANGE, lw=2, zorder=4)
        elif v == 1:
            a.scatter(o(p), yb, s=38, color=INK, zorder=3)
        elif v == 0:
            a.scatter(o(p), yb, s=38, facecolor='white', edgecolor=INK, lw=1.2, zorder=3)
    if not edited:
        a.annotate('', xy=(o(drop), yb * 0.55), xytext=(o(drop), yb * 0.05 + 0.25 * top),
                   arrowprops=dict(arrowstyle='->', color=ORANGE, lw=1.6))
    a.set_ylim(yb * 1.6, top)
    a.set_ylabel('running score')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--extract', required=True)
    ap.add_argument('--emissions', required=True)
    ap.add_argument('--example', action='append', required=True, help='kind:row:ref_pos0 of the dropped A/T')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    if Path(a.output).exists():
        raise FileExistsError(a.output)
    pos, X, distal, _, _ = a7.a6x.load(a.extract)
    hit, miss = a7.emission_steps(a.emissions, pos)
    a7f.a4f.style()
    ex = [e.split(':') for e in a.example]
    fig, ax = plt.subplots(2, len(ex), figsize=(8.4 * len(ex), 7.6), squeeze=False,
                           gridspec_kw=dict(hspace=0.55, wspace=0.18))
    for k, (kind, row, drop) in enumerate(ex):
        x = X[int(row)]
        if a7.classes(distal[int(row):int(row) + 1])[0] != 'C':
            raise ValueError(f'row {row} is not a C molecule')
        drop = int(drop)
        c = int(np.flatnonzero(pos == drop)[0])
        y = x.copy()
        y[c] = -1
        panel(ax[0, k], x, pos, hit, miss, drop, drop, False)
        panel(ax[1, k], y, pos, hit, miss, drop, drop, True)
        ax[0, k].set_title(f'{TITLES[kind]}\nWT molecule as read', loc='left', color=INK, fontsize=14, pad=26)
        ax[1, k].set_title('same molecule, that adenine dropped', loc='left', color=INK, fontsize=14, pad=26)
        ax[1, k].set_xlabel('oligo position')
    h = [plt.Line2D([], [], marker='o', ls='', color=INK, label='marked adenine (m6A)'),
         plt.Line2D([], [], marker='o', ls='', markerfacecolor='white', markeredgecolor=INK, label='unmarked adenine'),
         plt.Line2D([], [], marker='x', ls='', color=ORANGE, mew=2, label='dropped (as for an SNV there)'),
         plt.Line2D([], [], color=BLUE, lw=1.8, label="caller's running score")]
    fig.legend(handles=h, frameon=False, fontsize=12, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.04))
    fig.savefig(a.output, dpi=100, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()

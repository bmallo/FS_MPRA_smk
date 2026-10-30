"""Step A8 explainer figure: one spike-in on two real WT molecules (HepG2 run 3), one M and one low-labeled C.

A spike makes every observed A/T in [start, start + size) unmarked, as in a8_spikes.py. Each molecule is drawn as
read and spiked, with its adenines, the caller's running score and its calls. Molecules are picked, in a seeded
random order, as the first M molecule with no call over the spike before and a detected call after, and the first
C molecule with labeling 0.01-0.2 already called over the spike before spiking (the typical case for each group).

Usage: python tools/analysis_v2/a8_example.py --extract <a6 HepG2_run3 dir> --emissions <emissions.tsv.gz> \
           --start 1310 --size 20 --output <png>
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

_spec = importlib.util.spec_from_file_location('a7e', Path(__file__).with_name('a7_example.py'))
a7e = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a7e)
a7, a7f = a7e.a7, a7e.a7f
INK, BODY, BLUE, ORANGE, GREY = a7e.INK, a7e.BODY, a7e.BLUE, a7e.ORANGE, a7e.GREY
EL0, HALF = a7e.EL0, 50
SEED = 20261029
LOW_C = (0.01, 0.2)


def spike(x, pos, s0, e0):
    y = x.copy()
    y[(pos >= s0) & (pos < e0) & (y == 1)] = 0
    return y


def detected(cl, s0, e0):
    centre, L = s0 + (e0 - s0 - 1) // 2, e0 - s0
    return any(s <= centre < e and min(e, e0) - max(s, s0) >= 0.5 * L for s, e in cl)


def panel(a, x, pos, hit, miss, s0, e0, centre, spiked):
    w = (pos >= centre - HALF) & (pos <= centre + HALF)
    o = lambda r: r - EL0 + 21
    a.spines[['top', 'right']].set_visible(False)
    a.axvspan(o(s0) - 0.5, o(e0) - 0.5, ymax=0.16, color=ORANGE, alpha=0.18 if spiked else 0.0, lw=0)
    for s, e in a7e.calls(x, pos, hit, miss):
        if e > centre - HALF and s < centre + HALF:
            a.axvspan(o(max(s, centre - HALF)) - 0.5, o(min(e, centre + HALF + 1)) - 0.5, color='#DCE8F5', lw=0,
                      zorder=0)
            a.text(o((max(s, centre - HALF) + min(e, centre + HALF)) / 2), 1.02, f'call {e - s} bp',
                   transform=a.get_xaxis_transform(), ha='center', va='bottom', color=INK, fontsize=12)
    r = a7e.running(x, hit, miss)
    a.step(o(pos[w]), r[w], where='post', color=BLUE, lw=1.8)
    a.axhline(a7.MIN_LLR, color=GREY, lw=1, ls='--')
    a.text(o(centre + HALF), a7.MIN_LLR, ' threshold 5', va='bottom', ha='right', color=BODY, fontsize=11)
    top = max(r[w].max(), a7.MIN_LLR) * 1.15 + 1
    yb = -0.12 * top
    for p, v in zip(pos[w], x[w]):
        if v == 1:
            a.scatter(o(p), yb, s=38, color=INK, zorder=3)
        elif v == 0:
            a.scatter(o(p), yb, s=38, facecolor='white', edgecolor=INK, lw=1.2, zorder=3)
    a.set_ylim(yb * 1.6, top)
    a.set_ylabel('running score')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--extract', required=True)
    ap.add_argument('--emissions', required=True)
    ap.add_argument('--start', type=int, required=True, help='0-based reference start of the spike')
    ap.add_argument('--size', type=int, required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    if Path(a.output).exists():
        raise FileExistsError(a.output)
    pos, X, distal, _, _ = a7.a6x.load(a.extract)
    hit, miss = a7.emission_steps(a.emissions, pos)
    s0, e0 = a.start, a.start + a.size
    centre = s0 + (a.size - 1) // 2
    cls = a7.classes(distal)
    picked = {}
    for row in np.random.default_rng(SEED).permutation(X.shape[0]):
        x = X[row]
        if not ((pos >= s0) & (pos < e0) & (x >= 0)).any():
            continue
        before = detected(a7e.calls(x, pos, hit, miss), s0, e0)
        if 'M' not in picked and cls[row] == 'M' and not before and \
                detected(a7e.calls(spike(x, pos, s0, e0), pos, hit, miss), s0, e0):
            picked['M'] = row
        elif 'C' not in picked and LOW_C[0] < distal[row] < LOW_C[1] and before:
            picked['C'] = row
        if len(picked) == 2:
            break
    a7f.a4f.style()
    fig, ax = plt.subplots(2, 2, figsize=(16.8, 7.6), gridspec_kw=dict(hspace=0.55, wspace=0.18))
    for k, (g, ttl) in enumerate([('M', 'Well-labeled molecule (M)'), ('C', 'Low-labeled molecule (C)')]):
        row = picked[g]
        x = X[row]
        panel(ax[0, k], x, pos, hit, miss, s0, e0, centre, False)
        panel(ax[1, k], spike(x, pos, s0, e0), pos, hit, miss, s0, e0, centre, True)
        ax[0, k].set_title(f'{ttl}, labeling {distal[row]:.2f}\nas read', loc='left', color=INK, fontsize=14, pad=26)
        ax[1, k].set_title(f'{a.size}-bp spike: adenines in the orange band set unmarked', loc='left', color=INK,
                           fontsize=14, pad=26)
        ax[1, k].set_xlabel('oligo position')
        print(f'{g}: row {row}, labeling {distal[row]:.3f}')
    h = [plt.Line2D([], [], marker='o', ls='', color=INK, label='marked adenine (m6A)'),
         plt.Line2D([], [], marker='o', ls='', markerfacecolor='white', markeredgecolor=INK, label='unmarked adenine'),
         plt.Rectangle((0, 0), 1, 1, color=ORANGE, alpha=0.18, label='spike'),
         plt.Line2D([], [], color=BLUE, lw=1.8, label="caller's running score")]
    fig.legend(handles=h, frameon=False, fontsize=12, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.04))
    fig.savefig(a.output, dpi=100, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()

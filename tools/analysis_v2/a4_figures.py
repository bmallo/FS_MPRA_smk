"""Step A4 slide figures from the a4_informative outputs (no recomputation).

Usage: python tools/analysis_v2/a4_figures.py --dir results/analysis_v2_20261028/a4_informative
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('a4', Path(__file__).with_name('a4_informative_bases.py'))
a4 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a4)

INK, BODY, BLUE, ORANGE, GREY, MUTED, BG = '#1B2430', '#4A5568', '#2F6DB5', '#D9782D', '#A7AFB9', '#6B7280', '#F7F6F2'
SHORT = {'Repeat 1 (Sp1)': 'R1 Sp1', 'Repeat 2 (SRE-1)': 'R2 SRE-1', 'Repeat 3 (Sp1)': 'R3 Sp1'}
SHADE = {'Repeat 1 (Sp1)': '#E3EDF8', 'Repeat 2 (SRE-1)': '#FBEFE4', 'Repeat 3 (Sp1)': '#E3EDF8'}


def style():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 15, 'axes.edgecolor': GREY, 'axes.labelcolor': BODY,
                         'xtick.color': BODY, 'ytick.color': BODY, 'axes.facecolor': BG, 'figure.facecolor': BG})


def run_lengths(seq):
    out = np.zeros(len(seq), int)
    i = 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        out[i:j] = j - i
        i = j
    return out


def tracks(d, out):
    pos = pd.read_csv(d / 'HepG2/positions.tsv', sep='\t')
    var = pd.read_csv(d / 'visibility/variants_visibility.tsv', sep='\t')
    wt = ''.join(pos.base)
    run = run_lengths(wt)
    x = pos.oligo_pos1.to_numpy()
    fig, ax = plt.subplots(3, 1, figsize=(16.6, 7.4), sharex=True,
                           gridspec_kw=dict(height_ratios=[1.5, 1.0, 1.6], hspace=0.22))
    for a in ax:
        for name, (s, e) in a4.REGIONS.items():
            a.axvspan(s - 0.5, e + 0.5, color=SHADE[name], zorder=0, lw=0)
        a.spines[['top', 'right']].set_visible(False)
    for name, (s, e) in a4.REGIONS.items():
        ax[0].text((s + e) / 2, 1.02, SHORT[name], transform=ax[0].get_xaxis_transform(), ha='center', va='bottom',
                   color=BODY, fontsize=13)
    a = ax[0]
    inrun = pos.is_at & (run >= 3)
    a.bar(x[pos.is_at & ~inrun], pos.p_a[pos.is_at & ~inrun], width=0.8, color=BLUE, label='adenine, run of 1–2')
    a.bar(x[inrun], pos.p_a[inrun], width=0.8, color=ORANGE, label='adenine inside an A/T run of ≥ 3')
    a.scatter(x[inrun & (pos.p_a < 0.1)], np.full((inrun & (pos.p_a < 0.1)).sum(), 0.04), marker='^', s=40,
              color=ORANGE, zorder=3)
    a.set_ylabel('methylation rate\nin M molecules')
    a.set_ylim(0, 1.05)
    a.legend(loc='upper left', bbox_to_anchor=(0, 1.32), ncol=2, frameon=False, fontsize=13)
    a = ax[1]
    a.step(x, pos.min_fp_llr_bp, where='mid', color=INK, lw=2)
    a.axhline(20, color=MUTED, lw=1, ls='--')
    a.text(x[-1] + 1.5, 20, '20 bp\n(TF)', color=MUTED, va='center', fontsize=12)
    a.set_ylabel('smallest callable\nfootprint (bp)')
    a.set_ylim(0, 25)
    a = ax[2]
    snv = var[var.cls == 'SNV']
    cols = {'A/T loss': ORANGE, 'A<->T swap': INK, 'A/T gain': BLUE, 'G<->C swap': GREY}
    lab = {'A/T loss': 'A/T → G/C', 'A<->T swap': 'A ↔ T', 'A/T gain': 'G/C → A/T', 'G<->C swap': 'G ↔ C'}
    off = {'A/T loss': -0.27, 'A<->T swap': -0.09, 'A/T gain': 0.09, 'G<->C swap': 0.27}
    for c, col in cols.items():
        s = snv[snv.seq_class == c]
        a.scatter(s.oligo_pos1 + off[c], s.delta_visible, s=26, color=col, label=lab[c], zorder=3, edgecolor=BG, lw=0.6)
    a.axhline(0, color=GREY, lw=1)
    a.set_ylabel('predicted change in\nvisible adenines')
    a.set_xlabel('oligo position')
    a.legend(loc='upper left', ncol=4, frameon=False, fontsize=13, bbox_to_anchor=(0, 1.16))
    a.set_ylim(-2, 3.8)
    fig.savefig(out / 'a4_tracks.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def runs(d, out):
    t = pd.read_csv(d / 'context_bias/p_a_by_run_length.tsv', sep='\t')
    t = t[t.scope == 'plasmid']
    fig, a = plt.subplots(figsize=(8.6, 6.0))
    a.spines[['top', 'right']].set_visible(False)
    series = [('naked_pEFGFP', 'naked plasmid (pEF-GFP)', INK), ('HepG2_M', 'HepG2, M molecules', BLUE),
              ('HEK293T_M', 'HEK293T, M molecules', ORANGE)]
    for key, label, col in series:
        s = t[t.table == key].sort_values('run_bin')
        a.plot(s.run_bin, s.p_rel_median, marker='o', ms=9, lw=2, color=col, label=label)
    a.set_xticks([1, 2, 3, 4, 5], ['1', '2', '3', '4', '5+'])
    a.set_xlabel('length of the A (or T) run the adenine sits in')
    a.set_ylabel('methylation rate relative to\nisolated adenines (median)')
    a.set_ylim(0, 1.08)
    a.legend(frameon=False, loc='upper right', fontsize=14)
    a.grid(axis='y', color='#E2E4E8', lw=1)
    a.set_axisbelow(True)
    fig.savefig(out / 'a4_runs.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    out = d / 'figures'
    out.mkdir(exist_ok=False)
    style()
    tracks(d, out)
    runs(d, out)


if __name__ == '__main__':
    main()

"""Slide-sized figures for the HEK293T vs HepG2 WT footprint comparison.

Reads compare_cells.py and stratified.py outputs.
Usage: python tools/hek_compare/deck_figures.py --dir results/hek293t_vs_hepg2_20261028
"""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REGIONS = {'P1 core': (54, 62), 'R1 Sp1': (89, 104), 'R2 SRE1': (124, 139), 'R3 Sp1': (140, 155)}  # oligo, inclusive
BG, INK, MUTED, GRID = '#F7F6F1', '#1B2A41', '#4A5568', '#DCDDD5'
BLUE, ORANGE, GRAY = '#2B6CB0', '#C8641E', '#8A94A3'
plt.rcParams.update({'figure.facecolor': BG, 'axes.facecolor': BG, 'savefig.facecolor': BG, 'font.size': 14,
                     'axes.titlesize': 15, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'axes.labelcolor': MUTED,
                     'axes.edgecolor': GRID, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'text.color': INK, 'axes.titlecolor': INK, 'axes.grid': True, 'grid.color': GRID,
                     'grid.linewidth': 0.7, 'axes.axisbelow': True})


def shade(ax):
    for name, (a, b) in REGIONS.items():
        ax.axvspan(a - 0.5, b + 0.5, color='#ECEAE2', zorder=0)
        ax.text((a + b) / 2, 0.97, name, ha='center', va='top', fontsize=11, color=MUTED, transform=ax.get_xaxis_transform())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', required=True)
    d = Path(ap.parse_args().dir)
    out = d / 'deck'
    out.mkdir(exist_ok=True)
    t = pd.read_csv(d / 'occupancy_oligo.tsv', sep='\t')
    st = pd.read_csv(d / 'stratified_protection.tsv', sep='\t')
    stats = {}

    # 1. Per-site labeling on the M class, HepG2 vs HEK293T.
    fig, ax = plt.subplots(figsize=(7.6, 6.6))
    ok = np.isfinite(t.hepg2_p_a) & np.isfinite(t.hek293t_p_a)
    edge = ok & ((t.oligo <= 26) | (t.oligo >= 226))
    ax.plot([0, 1], [0, 1], color=GRAY, lw=1, zorder=1)
    ax.scatter(t.hepg2_p_a[ok & ~edge], t.hek293t_p_a[ok & ~edge], s=46, color=BLUE, edgecolor=BG, lw=1.2, zorder=3, label='interior A/T')
    ax.scatter(t.hepg2_p_a[edge], t.hek293t_p_a[edge], s=56, color=ORANGE, edgecolor=BG, lw=1.2, zorder=4, label='within 5 bp of an element end')
    for _, r in t[edge & ((t.hepg2_p_a - t.hek293t_p_a).abs() > 0.3)].iterrows():
        ax.annotate(f'oligo {int(r.oligo)}', (r.hepg2_p_a, r.hek293t_p_a), xytext=(10, 4), textcoords='offset points', fontsize=12, color=INK)
    ax.set(xlim=(0, 1.02), ylim=(0, 1.02), xlabel='HepG2  P(m6A | accessible)', ylabel='HEK293T  P(m6A | accessible)')
    ax.legend(loc='upper left', frameon=False, fontsize=12)
    r_all = float(np.corrcoef(t.hepg2_p_a[ok], t.hek293t_p_a[ok])[0, 1])
    r_int = float(np.corrcoef(t.hepg2_p_a[ok & ~edge], t.hek293t_p_a[ok & ~edge])[0, 1])
    stats['labeling'] = dict(sites=int(ok.sum()), edge_sites=int(edge.sum()), pearson_all=r_all, pearson_interior=r_int)
    fig.savefig(out / 'deck_hek_labeling.png', dpi=200, bbox_inches='tight'); plt.close(fig)

    # 2. Occupancy across the element and the difference.
    fig, (a, b) = plt.subplots(2, 1, figsize=(16, 7.6), sharex=True, gridspec_kw={'height_ratios': [1.6, 1], 'hspace': 0.12})
    shade(a)
    a.plot(t.oligo, t.llr_primary_hepg2, color=BLUE, lw=1.2, ls='--', label='HepG2, as observed')
    a.plot(t.oligo, t.llr_primary_hepg2_matched, color=BLUE, lw=2.2, label='HepG2, reweighted to HEK293T labeling')
    a.plot(t.oligo, t.llr_primary_hek293t, color=ORANGE, lw=2.2, label='HEK293T')
    a.legend(loc='lower left', ncol=3, frameon=False, fontsize=12)
    a.set(ylim=(0.3, 0.85), ylabel='fraction protected')
    shade(b)
    b.fill_between(t.oligo, t.llr_primary_diff_matched_lo, t.llr_primary_diff_matched_hi, color=GRID, lw=0)
    b.plot(t.oligo, t.llr_primary_diff_matched, color=INK, lw=1.8)
    b.axhline(0, color=GRAY, lw=1)
    b.axhline(t.llr_primary_diff_matched.mean(), color=ORANGE, lw=1, ls=':')
    b.set(ylim=(-0.16, 0.04), ylabel='HEK293T − HepG2', xlabel='oligo position (shared 211-bp element)', xlim=(20, 232))
    for ax in (a, b):
        for txt in ax.texts[:len(REGIONS)]:
            txt.set_visible(ax is a)
    fig.savefig(out / 'deck_hek_occupancy.png', dpi=200, bbox_inches='tight'); plt.close(fig)
    stats['occupancy'] = dict(mean_diff_matched=float(t.llr_primary_diff_matched.mean()), mean_diff_unmatched=float(t.llr_primary_diff.mean()),
                              resid_sd=float((t.llr_primary_diff_matched - t.llr_primary_diff_matched.mean()).std()),
                              pearson_matched=float(np.corrcoef(t.llr_primary_hepg2_matched, t.llr_primary_hek293t)[0, 1]),
                              diff_153_185=float(t.llr_primary_diff_matched[(t.oligo >= 153) & (t.oligo <= 185)].mean()))

    # 3. Element vs whole plasmid within labeling strata.
    s = st[st.distal_bin != 'all C'].copy()
    mids = {k: i for i, k in enumerate(s[s.cell == 'HepG2'].distal_bin)}
    labels = [k.strip('(]').replace(', ', '–') for k in mids]
    fig, axes = plt.subplots(1, 2, figsize=(16, 5.8), sharey=True, gridspec_kw={'wspace': 0.08})
    for ax, col, title in ((axes[0], 'element_protected', 'Shared 211-bp element'), (axes[1], 'plasmid_protected', 'Whole plasmid')):
        for cell, color in (('HepG2', BLUE), ('HEK293T', ORANGE)):
            c = s[s.cell == cell]
            x = [mids[k] for k in c.distal_bin]
            ax.plot(x, c[col], color=color, lw=2.2, marker='o', ms=9, mec=BG, mew=2, label=cell)
        ax.set(title=title, xticks=range(len(labels)), xlabel='m6A per A/T outside the element (labeling level)')
        ax.set_xticklabels(labels)
    axes[0].set(ylabel='fraction protected', ylim=(0, 1))
    axes[0].legend(frameon=False, loc='upper right')
    fig.savefig(out / 'deck_hek_stratified.png', dpi=200, bbox_inches='tight'); plt.close(fig)
    piv = s.pivot(index='distal_bin', columns='cell', values=['element_protected', 'plasmid_protected'])
    stats['stratified'] = dict(element_gap=(piv.element_protected.HEK293T - piv.element_protected.HepG2).round(3).to_dict(),
                               plasmid_gap=(piv.plasmid_protected.HEK293T - piv.plasmid_protected.HepG2).round(3).to_dict())
    json.dump(stats, open(out / 'deck_stats.json', 'w'), indent=1)
    print(json.dumps(stats, indent=1))


if __name__ == '__main__':
    main()

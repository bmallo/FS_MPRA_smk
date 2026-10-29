"""Slide-sized figures for per-variant effects in HepG2 and HEK293T and their concordance.

Usage: python tools/variant_effects/deck_figures.py --dir results/variant_effects_20261028
"""
import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BG, INK, MUTED, GRID = '#F7F6F1', '#1B2A41', '#4A5568', '#DCDDD5'
BLUE, ORANGE, GRAY = '#2B6CB0', '#C8641E', '#8A94A3'
REGIONS = {'P1 core': (54, 62), 'R1 Sp1': (89, 104), 'R2 SRE1': (124, 139), 'R3 Sp1': (140, 155)}
plt.rcParams.update({'figure.facecolor': BG, 'axes.facecolor': BG, 'savefig.facecolor': BG, 'font.size': 14,
                     'axes.titlesize': 15, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'axes.labelcolor': MUTED,
                     'axes.edgecolor': GRID, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'text.color': INK, 'axes.titlecolor': INK, 'axes.grid': True, 'grid.color': GRID,
                     'grid.linewidth': 0.7, 'axes.axisbelow': True})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', required=True)
    V = Path(ap.parse_args().dir)
    out = V / 'deck'
    out.mkdir(exist_ok=True)
    stats = {}
    hep = pd.read_csv(V / 'hepg2/effects/effects.tsv.gz', sep='\t')
    hep_n = hep.drop_duplicates('source_variant').set_index('source_variant')
    with gzip.open(V / 'hek293t/extract/molecules.tsv.gz', 'rt') as h:
        hek_n = Counter(r['source_variant'] for r in csv.DictReader(h, delimiter='\t'))
    clones = pd.read_csv(V / 'hek293t/extract/clones.tsv', sep='\t')

    # 1. Molecules per shared SNV along the element, both cells; HEK top-clone share.
    fig, ax = plt.subplots(figsize=(16, 5.6))
    pos = hep_n.oligo_pos1.values
    ax.scatter(pos, hep_n.molecules.values, s=12, color=BLUE, label='HepG2 rep3', zorder=3)
    hk = np.array([hek_n.get(v, 0) for v in hep_n.index])
    ax.scatter(pos[hk > 0], hk[hk > 0], s=16, color=ORANGE, label='HEK293T', zorder=4)
    ax.axhline(300, color=INK, lw=1, ls=':')
    ax.text(232, 330, '300-molecule minimum', ha='right', va='bottom', fontsize=12, color=MUTED)
    ax.set(yscale='log', ylim=(0.8, 2e4), xlim=(19, 233), xlabel='oligo position', ylabel='eligible molecules per SNV')
    ax.legend(frameon=False, loc='upper left', ncol=2)
    fig.savefig(out / 'deck_var_coverage.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    stats['coverage'] = dict(hek_observed=int((hk > 0).sum()), hek_ge_300=int((hk >= 300).sum()), hepg2_ge_300=int((hep_n.molecules >= 300).sum()),
                             hek_snv_molecules_median_ge_300=float(np.median(hk[hk >= 300])),
                             hek_wt_clones_ge_300=int(((clones.top_genotype == 'WT') & (clones.molecules >= 300)).sum()))

    # 2. Concordance scatter: caller-free element m6A, MSP, local and flank A_norm.
    con = pd.read_csv(V / 'concordance/concordance.tsv', sep='\t')
    lib = pd.read_csv(V / 'concordance/library_concordance.tsv', sep='\t').set_index('endpoint')
    eps = [('element:C_std', 'Element m6A'), ('msp_C_std', 'Accessibility (MSP)'),
           ('local:A_norm', 'm6A 5–10 bp away'), ('flank:A_norm', 'm6A 11–30 bp away')]
    fig, axes = plt.subplots(1, 4, figsize=(18, 5.2), gridspec_kw={'wspace': 0.28})
    for ax, (ep, title) in zip(axes, eps):
        d = con[con.endpoint == ep]
        rep = d.q_rep_primary_family.lt(.05) if ep in ('local:A_norm', 'flank:A_norm') else d.q_rep_within_endpoint.lt(.05)
        rep = rep.fillna(False)
        ax.axhline(0, color=GRAY, lw=.8)
        ax.axvline(0, color=GRAY, lw=.8)
        ax.scatter(d.hepg2_z[~rep], d.hek293t_z[~rep], s=22, color=GRAY, edgecolor=BG, lw=.6)
        ax.scatter(d.hepg2_z[rep], d.hek293t_z[rep], s=40, color=ORANGE, edgecolor=BG, lw=.8, label='replicated, q < 0.05')
        L = lib.loc[ep]
        ax.set_title(title, fontsize=14)
        ax.text(.03, .97, f"r = {L.pearson_z:.2f}" + ('' if L.usable else '\nHEK gate failed'), transform=ax.transAxes, va='top', fontsize=13)
        ax.set(xlabel='HepG2 z')
    axes[0].set_ylabel('HEK293T z')
    axes[0].legend(frameon=False, loc='lower right', fontsize=11)
    fig.savefig(out / 'deck_var_concordance.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    stats['concordance'] = {ep: dict(r=float(lib.loc[ep].pearson_z), p_rotation=float(lib.loc[ep].p_circular_shift),
                                     usable=bool(lib.loc[ep].usable), sign_agree=f'{int(lib.loc[ep].sign_agree)}/{int(lib.loc[ep].hepg2_significant)}')
                            for ep, _ in eps}

    # 3. Track: element m6A effects along the element, both cells, replicated SNVs marked.
    d = con[con.endpoint == 'element:C_std'].sort_values('oligo_pos1')
    fig, ax = plt.subplots(figsize=(16, 5.4))
    for name, (a, b) in REGIONS.items():
        ax.axvspan(a - .5, b + .5, color='#ECEAE2', zorder=0)
        ax.text((a + b) / 2, .97, name, ha='center', va='top', fontsize=11, color=MUTED, transform=ax.get_xaxis_transform())
    ax.axhline(0, color=GRAY, lw=.8)
    ax.scatter(d.oligo_pos1 - .2, d.hepg2_delta * 100, s=26, color=BLUE, label='HepG2', zorder=3)
    ax.scatter(d.oligo_pos1 + .2, d.hek293t_delta * 100, s=26, color=ORANGE, label='HEK293T', zorder=3)
    for _, r in d[d.q_rep_within_endpoint < .05].iterrows():
        ax.plot([r.oligo_pos1] * 2, [r.hepg2_delta * 100, r.hek293t_delta * 100], color=INK, lw=1, zorder=2)
    ax.set(xlim=(19, 233), xlabel='oligo position', ylabel='Δ element m6A vs WT (pp)')
    ax.legend(frameon=False, loc='lower left', ncol=2)
    fig.savefig(out / 'deck_var_track.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    json.dump(stats, open(out / 'deck_stats.json', 'w'), indent=1)
    print(json.dumps(stats, indent=1))


if __name__ == '__main__':
    main()

"""Track-level layout: variant effects on function, accessibility, nucleosomes and sub-nucleosomal footprints.

Tracks along the shared 211-bp element (oligo 21-231):
  1 published activity (mean log2 of the two sources), functional hits coloured by direction
  2 accessibility change (HepG2 LLR MSP, C_std), HEK293T overlaid for its 105 SNVs
  3 nucleosome-sized protection change (85-250 bp calls)
  4 WT sub-nucleosomal footprint occupancy (TF-sized LLR calls, C_std, per base)
  5 local TF-footprint change at each SNV (observability-scaled residual)
  6 candidate classes per position and the best-ranked motif per region
Plus a footprint-change map (substituted position x element position, TF-sized calls).

Usage: python tools/variant_effects/track_figure.py --dir results/variant_effects_20261028
"""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

BG, INK, MUTED, GRID = '#F7F6F1', '#1B2A41', '#4A5568', '#DCDDD5'
BLUE, ORANGE, GRAY, LIGHT = '#2B6CB0', '#C8641E', '#8A94A3', '#ECEAE2'
REGIONS = {'P1 core': (54, 62), 'R1 Sp1': (89, 104), 'R2 SRE1': (124, 139), 'R3 Sp1': (140, 155)}
CLASS_COLOR = {'activator_accessibility': BLUE, 'activator_footprint': '#1B4F8A', 'repressor_accessibility': ORANGE,
               'repressor_footprint': '#8A3E10', 'created_activator_footprint': '#E9A36B', 'created_repressor_footprint': '#9CC3E6',
               'discordant_accessibility': GRAY}
# Best-ranked WT motif per region (motif_sites_ranked.tsv) and notable created sites (motifs_created.tsv).
MOTIFS = [(52, 62, 'ZNF652-like', BLUE), (82, 102, 'SP4', BLUE), (127, 136, 'SREBF1/2', BLUE), (140, 156, 'SP1/SP3/KLF', BLUE),
          (182, 192, 'CREB1/CREM', BLUE), (157, 157, 'CTCF created\nby 157A>G', ORANGE)]
plt.rcParams.update({'figure.facecolor': BG, 'axes.facecolor': BG, 'savefig.facecolor': BG, 'font.size': 13,
                     'axes.titlesize': 13, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'axes.labelcolor': MUTED,
                     'axes.edgecolor': GRID, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'text.color': INK, 'axes.titlecolor': INK, 'axes.grid': True, 'grid.color': GRID,
                     'grid.linewidth': .6, 'axes.axisbelow': True})


def shade(ax, labels=False):
    for name, (a, b) in REGIONS.items():
        ax.axvspan(a - .5, b + .5, color=LIGHT, zorder=0, lw=0)
        if labels:
            ax.text((a + b) / 2, .97, name, ha='center', va='top', fontsize=11, color=MUTED, transform=ax.get_xaxis_transform())


def stems(ax, x, y, colors, size=16):
    ax.vlines(x, 0, y, colors=colors, lw=.8, alpha=.7)
    ax.scatter(x, y, s=size, c=colors, edgecolor=BG, lw=.4, zorder=3)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', required=True)
    V = Path(ap.parse_args().dir)
    out = V / 'deck'
    out.mkdir(exist_ok=True)
    v = pd.read_csv(V / 'classes/variants.tsv', sep='\t')
    wt = pd.read_csv(V / 'hepg2/footprints/wt_profile.tsv', sep='\t')
    jit = v.alt.map({'A': -.25, 'C': -.08, 'G': .08, 'T': .25}).fillna(0)
    x = v.oligo_pos1 + jit

    def draw(slot):
        # 1 function
        ax = slot.get(0)
        if ax is not None:
            shade(ax, labels=5 not in slot)
            hit = v.functional_hit.astype(bool)
            col = np.where(~hit, GRAY, np.where(v.activity < 0, BLUE, ORANGE))
            stems(ax, x, v.activity, col)
            ax.axhline(0, color=MUTED, lw=.8)
            ax.set(ylabel='activity\n(log2)', title='Function: published MPRA activity per SNV (blue lower, orange higher at q < 0.05 in both sources; grey not significant)')
        # 2 accessibility
        ax = slot.get(1)
        if ax is not None:
            shade(ax)
            sig = v.hepg2_msp_q < .05
            stems(ax, x, v.hepg2_msp_delta * 100, np.where(sig, np.where(v.hepg2_msp_delta < 0, BLUE, ORANGE), GRAY))
            hk = v.hek293t_msp_delta.notna()
            ax.scatter(x[hk] + .35, v.hek293t_msp_delta[hk] * 100, s=14, marker='D', facecolor='none', edgecolor=INK, lw=.8, zorder=4,
                       label='HEK293T (105 SNVs)')
            ax.axhline(0, color=MUTED, lw=.8)
            ax.legend(frameon=False, loc='lower right', fontsize=10)
            ax.set(ylabel='Δ MSP\n(pp)', title='Accessibility: HepG2 change in MSP coverage (blue closes, orange opens at q < 0.05); HEK293T open diamonds')
        # 3 nucleosome-sized
        ax = slot.get(2)
        if ax is not None:
            shade(ax)
            sig = v.hepg2_nuc_q < .05
            stems(ax, x, v.hepg2_nuc_delta * 100, np.where(sig, np.where(v.hepg2_nuc_delta > 0, BLUE, ORANGE), GRAY))
            ax.axhline(0, color=MUTED, lw=.8)
            ax.set(ylabel='Δ nuc\n(pp)', title='Nucleosome-sized protection, 85–250-bp calls (blue more, orange less at q < 0.05)')
        # 4 WT sub-nucleosomal occupancy
        ax = slot.get(3)
        if ax is not None:
            shade(ax)
            ax.fill_between(wt.site_oligo_pos1, wt.wt_tf_C_std * 100, color=INK, alpha=.25, lw=0)
            ax.plot(wt.site_oligo_pos1, wt.wt_tf_C_std * 100, color=INK, lw=1.2)
            ax.set(ylabel='WT TF\n(%)', title='WT sub-nucleosomal footprints: fraction of chromatin molecules with a 10–60-bp call')
        # 5 local footprint change
        ax = slot.get(4)
        if ax is not None:
            shade(ax)
            sig = v.hepg2_tf_local_q < .05
            stems(ax, x, v.hepg2_tf_local_delta * 100, np.where(sig, np.where(v.hepg2_tf_local_delta < 0, BLUE, ORANGE), GRAY))
            ax.axhline(0, color=MUTED, lw=.8)
            ax.set(ylabel='Δ TF local\n(pp)', title='Local sub-nucleosomal footprint change, 5–20 bp from the SNV, net of observability (coloured = q < 0.05)')
        # 6 classes + motifs
        ax = slot.get(5)
        if ax is not None:
            ax.grid(False)
            for _, r in v[v['class'].notna()].iterrows():
                ax.bar(r.oligo_pos1 + jit[_], .8, width=.22, bottom=.95, color=CLASS_COLOR[r['class']], lw=0)
            for a, b, name, c in MOTIFS:
                ax.add_patch(plt.Rectangle((a - .5, .1), b - a + 1, .55, color=c, alpha=.25 if a != b else .9, lw=0))
                if a == b:
                    ax.text(b + 1.5, .38, name, ha='left', va='center', fontsize=9, color=INK)
                else:
                    ax.text((a + b) / 2, .38, name, ha='center', va='center', fontsize=10, color=INK)
            for name, (a, b) in REGIONS.items():
                ax.add_patch(plt.Rectangle((a - .5, 1.95), b - a + 1, .5, color=LIGHT, lw=0))
                ax.text((a + b) / 2, 2.2, name, ha='center', va='center', fontsize=10, color=MUTED)
            ax.set(ylim=(0, 2.5), yticks=[.38, 1.35, 2.2], yticklabels=['motif', 'class', 'region'], xlim=(19, 233),
                   xlabel='oligo position (chr19:11,089,285 = 21)')
            ax.set_title('Candidate classes (blue activator, orange repressor, grey discordant) and best-ranked motifs', loc='left')

    def render(path, keys, figsize, ratios):
        fig, axes = plt.subplots(len(keys), 1, figsize=figsize, sharex=True, gridspec_kw={'height_ratios': ratios, 'hspace': .38})
        draw(dict(zip(keys, axes)))
        for ax in axes:
            ax.tick_params(labelbottom=True)
            ax.set_xlim(19, 233)
        fig.savefig(path, dpi=200, bbox_inches='tight')
        plt.close(fig)

    render(out / 'deck_tracks.png', [0, 1, 2, 3, 4, 5], (17, 15), [1.2, 1.1, 1, .7, 1, .55])
    render(out / 'deck_tracks_a.png', [0, 1, 2], (17, 8), [1.2, 1.1, 1])
    render(out / 'deck_tracks_b.png', [3, 4, 5], (17, 7.4), [.8, 1.1, .75])

    # Footprint-change map: rows = substituted position (mean over alleles), columns = element position.
    z = np.load(V / 'hepg2/footprints/profiles.npz', allow_pickle=True)
    names, wi = z['variants'], int(z['wt_index'])
    b = json.load(open(V / 'classes/summary.json'))['observability_slopes']['tf']
    d = ((z['tf_C_std'] - z['tf_C_std'][wi]) - b * (z['tf_M'] - z['tf_M'][wi])) * 100
    pos = np.array([int(''.join(ch for ch in s.split('>')[0] if ch.isdigit())) if s != 'WT' else 0 for s in names])
    mat = np.full((211, 211), np.nan)
    for p in range(21, 232):
        rows = d[pos == p]
        if len(rows):
            mat[p - 21] = np.nanmean(rows, 0)
    fig, ax = plt.subplots(figsize=(9.5, 8.4))
    lim = np.nanpercentile(np.abs(mat), 99)
    im = ax.imshow(mat, origin='lower', extent=(20.5, 231.5, 20.5, 231.5), cmap='RdBu_r', norm=TwoSlopeNorm(0, -lim, lim), aspect='auto')
    ax.plot([21, 231], [21, 231], color=INK, lw=.6, alpha=.5)
    ax.grid(False)
    ax.set(xlabel='element position', ylabel='substituted position', title=f'TF-sized footprint change (pp), HepG2, net of observability (ΔC − {b:.2f}·ΔM)')
    cb = fig.colorbar(im, ax=ax, shrink=.8)
    cb.set_label('Δ occupancy (pp)', color=MUTED)
    fig.savefig(out / 'deck_footprint_map.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('wrote', out / 'deck_tracks.png', out / 'deck_footprint_map.png')


if __name__ == '__main__':
    main()

"""Slide figures for size-resolved footprints.

deck_size_map.png      WT call density by midpoint x length: HepG2 chromatin, HepG2 near-naked, HEK293T chromatin,
                       frozen feature boxes labelled.
deck_size_effects.png  activity per SNV over a feature x substituted-position heatmap of the occupancy residual
                       (calibrated features only), per-SNV hits marked.
deck_size_f08.png      F08: occupancy change per SNV coloured by activity, and variant - WT size maps for two SNVs.

Usage: python tools/footprint_sizes/figures.py --dir results/footprint_sizes_20261028 \
           --classes results/variant_effects_20261028/classes
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm, TwoSlopeNorm
import numpy as np
import pandas as pd
from scipy import ndimage

BG, INK, MUTED, GRID = '#F7F6F1', '#1B2A41', '#4A5568', '#DCDDD5'
BLUE, ORANGE, GRAY, LIGHT = '#2B6CB0', '#C8641E', '#8A94A3', '#ECEAE2'
REGIONS = {'P1 core': (54, 62), 'R1 Sp1': (89, 104), 'R2 SRE1': (124, 139), 'R3 Sp1': (140, 155)}
SEQ = LinearSegmentedColormap.from_list('seq', [BG, '#9CC3E6', BLUE, INK])
DIV = LinearSegmentedColormap.from_list('div', [BLUE, '#9CC3E6', BG, '#E9A36B', ORANGE])
plt.rcParams.update({'figure.facecolor': BG, 'axes.facecolor': BG, 'savefig.facecolor': BG, 'font.size': 13,
                     'axes.titlesize': 13, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'axes.labelcolor': MUTED,
                     'axes.edgecolor': GRID, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'text.color': INK, 'axes.titlecolor': INK, 'axes.grid': False})
X = (20.5, 231.5)


def regions(ax, y=1.0):
    for name, (a, b) in REGIONS.items():
        ax.axvspan(a - .5, b + .5, ymin=.985, ymax=1, color=GRAY, lw=0, clip_on=False)
        ax.text((a + b) / 2, y, name, ha='center', va='bottom', fontsize=10, color=MUTED, transform=ax.get_xaxis_transform())


def size_map(ax, D, feats, title, vmax, label=True):
    im = ax.imshow(D.T, origin='lower', extent=(20.5, 231.5, 9.5, 90.5), aspect='auto', cmap=SEQ, norm=PowerNorm(.5, 0, vmax))
    for f in feats.itertuples():
        ax.add_patch(plt.Rectangle((f.mid_lo_oligo - .5, f.len_lo - .5), f.mid_hi_oligo - f.mid_lo_oligo + 1, f.len_hi - f.len_lo + 1,
                                   fill=False, ec=ORANGE if f.flag == 'chromatin' else INK, lw=1.1))
        if label:
            ax.text(f.mid_hi_oligo + 1, f.len_hi + .5, f.feature, fontsize=8.5, color=INK, va='bottom')
    ax.set(xlim=X, ylim=(9.5, 90.5), ylabel='call length (bp)')
    ax.set_title(title, loc='left')
    return im


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--classes', required=True)
    a = ap.parse_args()
    R, K = Path(a.dir), Path(a.classes)
    out = R / 'deck'
    out.mkdir(exist_ok=True)
    feats = pd.read_csv(R / 'hepg2_map/features.tsv', sep='\t')
    fs = pd.read_csv(R / 'summary/feature_summary.tsv', sep='\t').set_index('feature')
    hep, hek = np.load(R / 'hepg2_map/map.npz'), np.load(R / 'hek293t_map/map.npz')

    # 1. WT size maps.
    fig, axes = plt.subplots(3, 1, figsize=(16, 12.5), sharex=True, gridspec_kw={'hspace': .3})
    smooth = lambda D: ndimage.gaussian_filter(D, (1.5, 2.0), mode='nearest')   # the peak-calling smoothing (config map)
    vmax = np.percentile(smooth(hep['discovery_C']), 99.9)
    im = size_map(axes[0], smooth(hep['discovery_C']), feats, 'HepG2 chromatin molecules (WT discovery half, 622k molecules); boxes = 18 frozen features', vmax)
    size_map(axes[1], smooth(hep['discovery_M']), feats, 'HepG2 near-naked molecules (M class), same boxes', vmax, label=False)
    size_map(axes[2], smooth(hek['all_C']), feats, 'HEK293T chromatin molecules (all WT, 163k molecules), HepG2 boxes', vmax, label=False)
    axes[2].set_xlabel('call midpoint, oligo position')
    cb = fig.colorbar(im, ax=axes, shrink=.5, pad=.01)
    cb.set_label('calls per molecule per bp² (sqrt scale)', color=MUTED)
    fig.text(.125, .065, 'Smoothed as for peak calling (σ 1.5 bp midpoint, 2 bp length). Orange box: chromatin-enriched (C/M ≥ 2); dark box: also present on near-naked molecules.', fontsize=11, color=MUTED)
    fig.savefig(out / 'deck_size_map.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    # 2. Effects heatmap: calibrated features x substituted position, mean occupancy residual over alleles.
    eff = pd.read_csv(R / 'hepg2_effects/effects.tsv.gz', sep='\t')
    hits = pd.read_csv(R / 'summary/snv_hits.tsv', sep='\t')
    cls = pd.read_csv(K / 'variants.tsv', sep='\t')
    gated = [f for f in feats.feature if fs.loc[f, 'gate_pass']]
    mat = np.full((len(gated), 211), np.nan)
    for r, f in enumerate(gated):
        g = eff[eff.feature == f].groupby('oligo_pos1').occupancy_resid_delta.mean()
        mat[r, g.index - 21] = g.values * 100
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(17, 10), sharex=True, gridspec_kw={'height_ratios': [1, 2.4], 'hspace': .12})
    hit = cls.functional_hit.astype(bool)
    col = np.where(~hit, GRAY, np.where(cls.activity < 0, BLUE, ORANGE))
    jit = cls.alt.map({'A': -.25, 'C': -.08, 'G': .08, 'T': .25}).fillna(0)
    ax0.vlines(cls.oligo_pos1 + jit, 0, cls.activity, colors=col, lw=.8, alpha=.7)
    ax0.scatter(cls.oligo_pos1 + jit, cls.activity, s=14, c=col, edgecolor=BG, lw=.4, zorder=3)
    ax0.axhline(0, color=MUTED, lw=.8)
    ax0.grid(True, color=GRID, lw=.6)
    ax0.set(ylabel='activity (log2)')
    ax0.set_title('Function: published MPRA activity per SNV (blue lower, orange higher at q < 0.05)', loc='left', pad=26)
    regions(ax0, 1.02)
    lim = np.nanpercentile(np.abs(mat), 99)
    im = ax1.imshow(mat, aspect='auto', cmap=DIV, norm=TwoSlopeNorm(0, -lim, lim), extent=(20.5, 231.5, len(gated) - .5, -.5),
                    interpolation='nearest')
    for r, f in enumerate(gated):
        s = fs.loc[f]
        lo, hi = (float(x) for x in s.span_oligo.split('-'))
        ax1.plot([lo, hi], [r + .42, r + .42], color=INK, lw=2.2, solid_capstyle='butt')
        h = hits[hits.feature == f]
        ax1.scatter(h.oligo_pos1, np.full(len(h), r), s=16, facecolor='none', edgecolor=INK, lw=.9, zorder=4)
    ax1.set_yticks(range(len(gated)))
    ax1.set_yticklabels([f'{f}  {int(fs.loc[f, "centre_oligo"])} / {int(fs.loc[f, "length_peak"])} bp' for f in gated], fontsize=11)
    ax1.set(xlim=X, xlabel='substituted position (oligo)')
    ax1.set_title('Feature occupancy change per substituted position (mean over alleles, net of near-naked change); '
                  'bar = footprint span, circle = SNV q < 0.05', loc='left')
    cax = ax1.inset_axes([1.01, .1, .012, .8])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('Δ occupancy (pp)', color=MUTED)
    fig.savefig(out / 'deck_size_effects.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    # 3. F08 focus.
    f8 = feats.set_index('feature').loc['F08']
    e8 = eff[eff.feature == 'F08'].merge(cls[['source_variant', 'activity', 'functional_hit']], on='source_variant')
    sm = np.load(R / 'hepg2_effects/size_maps.npz', allow_pickle=True)
    names, wi = list(sm['variants']), int(sm['wt_index'])
    fig = plt.figure(figsize=(17, 9.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=.38, wspace=.12)
    ax = fig.add_subplot(gs[0, :])
    w = e8[(e8.oligo_pos1 >= 95) & (e8.oligo_pos1 <= 160)].copy()
    w['x'] = w.oligo_pos1 + w.alt.map({'A': -.25, 'C': -.08, 'G': .08, 'T': .25})
    fh = w.functional_hit.fillna(False).astype(bool)
    c = np.where(~fh, GRAY, np.where(w.activity < 0, BLUE, ORANGE))
    top = np.ceil(w.occupancy_resid_delta.max() * 100 + .4) + .6
    ax.axvspan(f8.mid_peak - f8.len_peak / 2 + 21, f8.mid_peak + f8.len_peak / 2 + 21, color=LIGHT, lw=0, zorder=0)
    for nm, (lo, hi) in {'SREBF1/2 125–136': (125, 136), 'SP1/SP3 140–158': (140, 158), '110–121 block': (110, 121)}.items():
        ax.plot([lo - .5, hi + .5], [top - .55, top - .55], color=INK, lw=3, solid_capstyle='butt')
        ax.text((lo + hi) / 2, top - .45, nm, ha='center', va='bottom', fontsize=10, color=MUTED)
    ax.vlines(w.x, 0, w.occupancy_resid_delta * 100, colors=c, lw=.9, alpha=.8)
    ax.scatter(w.x, w.occupancy_resid_delta * 100, s=26, c=c, edgecolor=BG, lw=.5, zorder=3)
    hh = hits[hits.feature == 'F08']
    hx = w[w.source_variant.isin(hh.source_variant)]
    ax.scatter(hx.x, hx.occupancy_resid_delta * 100, s=70, facecolor='none', edgecolor=INK, lw=1, zorder=4)
    ax.axhline(0, color=MUTED, lw=.8)
    ax.axhline(-fs.loc['F08', 'occ_reference'] * 100, color=INK, lw=.8, ls=':')
    ax.text(160, -fs.loc['F08', 'occ_reference'] * 100 + .1, 'complete loss (WT 3.3%)', ha='right', va='bottom', fontsize=10, color=MUTED)
    ax.grid(True, color=GRID, lw=.6)
    ax.set(xlim=(94.5, 160.5), ylim=(-4, top), ylabel='Δ F08 occupancy (pp)', xlabel='substituted position (oligo)')
    ax.set_title('F08 (124 / 12 bp) per SNV: colour = activity (blue lower, orange higher, grey n.s.); ring = q < 0.05; shaded = footprint span', loc='left')
    for k, v in enumerate(['117A>T', '128T>C']):
        a2 = fig.add_subplot(gs[1, k])
        D = ndimage.gaussian_filter((sm['C_std'][names.index(v)] - sm['C_std'][wi]).astype(float), (1.5, 2.0), mode='nearest')
        L = sm['lengths']
        sub = D[95 - 21:161 - 21, :61 - L[0]]
        lim = np.percentile(np.abs(sub), 99.5)
        a2.imshow(sub.T, origin='lower', extent=(94.5, 160.5, 9.5, 60.5), aspect='auto', cmap=DIV, norm=TwoSlopeNorm(0, -lim, lim))
        a2.add_patch(plt.Rectangle((f8.mid_lo + 21 - .5, f8.len_lo - .5), f8.mid_hi - f8.mid_lo + 1, f8.len_hi - f8.len_lo + 1, fill=False, ec=INK, lw=1))
        a2.axvline(int(''.join(ch for ch in v.split('>')[0] if ch.isdigit())), color=INK, lw=.8, ls='--')
        act = cls.set_index('source_variant').loc[v, 'activity']
        a2.set_title(f'{v} minus WT, smoothed (activity {act:+.2f} log2); blue = fewer calls', loc='left', fontsize=12)
        a2.set(xlabel='call midpoint (oligo)', ylabel='call length (bp)' if k == 0 else '')
    fig.savefig(out / 'deck_size_f08.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('wrote', *sorted(p.name for p in out.glob('deck_size_*.png')))


if __name__ == '__main__':
    main()

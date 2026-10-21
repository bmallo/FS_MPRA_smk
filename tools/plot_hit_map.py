"""Static hit map of calibrated replicate-3 chromatin calls across the LDLR element.

Reads only committed evidence (mask-4 site-methylation bundle and the v2 panel):

    python -m tools.plot_hit_map --output docs/figures/ldlr_hit_map_20261021
"""
import argparse
import csv
from pathlib import Path
import numpy as np

EVIDENCE = Path('docs/evidence')
GENOMIC_OFFSET = 11089264
INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
ACCENT, SECOND = '#2a78d6', '#eb6834'
ELEMENTS = [('P1 core', 54, 62, 'data'), ('R1 · Sp1', 89, 104, 'lit'), ('115–117', 115, 117, 'cand'),
            ('R2 · SRE-1', 124, 139, 'lit'), ('R3 · Sp1', 140, 155, 'lit'), ('165–166', 165, 166, 'cand'),
            ('190–191', 190, 191, 'cand')]
TRACKS = [('activity', 'Published activity\n(log2)', None),
          ('msp_C_std', 'Accessibility\n(MSP, pp)', 'accessibility'),
          ('nuc_85_250_C_std', 'Nucleosome-sized\nprotection (pp)', 'nucleosome'),
          ('local:A_norm', 'Local methylation\n5–10 bp (pp)', 'local')]
LABELS = {'msp_C_std': ['146C>G', '149C>G', '98C>G', '117A>T'], 'nuc_85_250_C_std': ['146C>T'],
          'local:A_norm': ['56G>C', '117A>T', '191G>T', '125A>T']}
OFFSETS = {('msp_C_std', '146C>G'): (-46, -16), ('msp_C_std', '149C>G'): (22, -12), ('msp_C_std', '98C>G'): (6, -14),
           ('local:A_norm', '117A>T'): (-48, 6), ('local:A_norm', '125A>T'): (8, 10)}


def load(bundle):
    import gzip
    rows = {}
    with gzip.open(bundle / 'effects' / 'effects.tsv.gz', 'rt') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            if int(r['molecules']) < 300:
                continue
            v = rows.setdefault(r['source_variant'], dict(variant=r['source_variant'], pos=int(r['oligo_pos1']),
                                context=r['context'], activity=float(r['activity_LDLR'])))
            v[r['endpoint']] = (float(r['delta']) * 100 if r['delta'] else np.nan,
                                float(r['q_bh_primary_family']) if r.get('q_bh_primary_family') else np.nan)
    return list(rows.values())


def plot(snvs, kernel, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 9, 'axes.edgecolor': AXIS, 'axes.labelcolor': INK2,
                         'xtick.color': MUTED, 'ytick.color': MUTED, 'text.color': INK})
    fig = plt.figure(figsize=(13, 11.2))
    grid = fig.add_gridspec(7, 2, height_ratios=[.55, 1, 1, 1, 1.15, .62, 1.25], hspace=.32, wspace=.18)
    top = fig.add_subplot(grid[0, :])
    axes = [fig.add_subplot(grid[i + 1, :], sharex=top) for i in range(4)]
    pos = np.array([s['pos'] for s in snvs])
    homo = sorted({s['pos'] for s in snvs if s['context'] == 'homopolymer'})
    for ax in [top] + axes:
        for name, a, b, kind in ELEMENTS:
            ax.axvspan(a - .5, b + .5, color=ACCENT if kind == 'cand' else MUTED, alpha=.07 if kind != 'cand' else .09, lw=0, zorder=0)
    for name, a, b, kind in ELEMENTS:
        if kind == 'lit':                       # wide literature elements: label inside, small gap between neighbours
            top.add_patch(plt.Rectangle((a - .1, .62), b - a + .2, .32, color=INK2, lw=0))
            top.text((a + b) / 2, .78, name, ha='center', va='center', fontsize=7, color='white')
        else:
            y = .5 if kind == 'data' else .22
            top.add_patch(plt.Rectangle((a - .5, y - .1), b - a + 1, .2, color=ACCENT if kind == 'cand' else MUTED, lw=0))
            top.text((a + b) / 2, y + .13, name, ha='center', va='bottom', fontsize=7.5, color=INK2)
    top.vlines(homo, -.02, .06, color=INK2, lw=1.2)
    top.set_ylim(-.05, 1.1)
    top.set_yticks([])
    for side in ('top', 'right', 'left'):
        top.spines[side].set_visible(False)
    top.set_title('Calibrated chromatin calls across the LDLR element, replicate 3 (q < 0.05, prespecified primary family)',
                  loc='left', fontsize=11, color=INK, pad=6)
    for ax, (key, label, _) in zip(axes, TRACKS):
        ax.axhline(0, color=AXIS, lw=1, zorder=1)
        ax.grid(axis='y', color=GRID, lw=.8)
        ax.set_axisbelow(True)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        ax.set_ylabel(label, fontsize=8.5, rotation=0, ha='right', va='center', labelpad=10)
        if key == 'activity':
            y = np.array([s['activity'] for s in snvs])
            ax.scatter(pos, y, s=9, color=MUTED, alpha=.55, lw=0, zorder=2)
            up = sorted(set(pos))
            ax.plot(up, [y[pos == p].mean() for p in up], color=INK2, lw=1.6, zorder=3)
            continue
        d = np.array([s.get(key, (np.nan, np.nan))[0] for s in snvs])
        q = np.array([s.get(key, (np.nan, np.nan))[1] for s in snvs])
        hit = q < .05
        unresolved = (np.array([s['context'] for s in snvs]) == 'homopolymer') & (key == 'local:A_norm')
        base = ~hit & ~unresolved
        ax.scatter(pos[base], d[base], s=9, color=MUTED, alpha=.5, lw=0, zorder=2)
        ax.scatter(pos[unresolved], d[unresolved], s=16, facecolor='none', edgecolor=INK2, lw=.8, alpha=.7, zorder=3)
        ax.scatter(pos[hit & ~unresolved], d[hit & ~unresolved], s=26, color=ACCENT, edgecolor='white', lw=.8, zorder=4)
        lim = np.nanmax(np.abs(d)) * 1.1
        ax.set_ylim(-lim, lim)
        byname = {s['variant']: (s['pos'], s.get(key, (np.nan,))[0]) for s in snvs}
        for v in LABELS.get(key, []):
            if v in byname and np.isfinite(byname[v][1]):
                x, yv = byname[v]
                ax.annotate(v, (x, yv), xytext=OFFSETS.get((key, v), (6, 8 if yv >= 0 else -11)), textcoords='offset points', fontsize=7.5, color=INK,
                            arrowprops=dict(arrowstyle='-', color=MUTED, lw=.6))
        n_hit = int((hit & ~unresolved).sum())
        ax.text(.995, .96, f'{n_hit} calls' + (' outside A runs' if key == 'local:A_norm' else ''), transform=ax.transAxes,
                ha='right', va='top', fontsize=8, color=INK2)
    axes[-1].set_xlim(20, 232)
    ticks = list(range(25, 232, 25))
    axes[-1].set_xticks(ticks)
    axes[-1].set_xticklabels([f'{p}\n{GENOMIC_OFFSET + p:,}' for p in ticks], fontsize=7.5)
    axes[-1].set_xlabel('Oligo position / chr19 GRCh38 position (forward strand)', fontsize=8.5)
    for ax in [top] + axes[:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)
    handles = [plt.Line2D([], [], marker='o', ls='', color=ACCENT, markeredgecolor='white', markersize=6, label='Called (q < 0.05)'),
               plt.Line2D([], [], marker='o', ls='', color=MUTED, alpha=.6, markersize=4, label='Not called'),
               plt.Line2D([], [], marker='o', ls='', markerfacecolor='none', markeredgecolor=INK2, markersize=5, label='Local readout unresolved (A run; ticks on top track)'),
               plt.Line2D([], [], color=INK2, lw=1.6, label='Activity, mean of the 3 substitutions')]
    axes[0].legend(handles=handles, loc='upper left', ncol=4, frameon=False, fontsize=7.5, bbox_to_anchor=(0, 1.32))
    for i, ctx in enumerate(('other', 'homopolymer')):
        ax = fig.add_subplot(grid[6, i])
        ax.axvspan(-4.5, 4.5, color=MUTED, alpha=.1, lw=0)
        ax.axhline(0, color=AXIS, lw=1)
        for sub, color, name in (('gain', SECOND, 'Adds an A/T'), ('loss', ACCENT, 'Removes an A/T')):
            k = sorted((r for r in kernel if r['context'] == ctx and r['at_change'] == sub), key=lambda r: r['offset'])
            x = [r['offset'] for r in k]
            ax.plot([v for v in x if v < 0], [r['mean'] for r in k if r['offset'] < 0], color=color, lw=2, solid_capstyle='round', label=name)
            ax.plot([v for v in x if v > 0], [r['mean'] for r in k if r['offset'] > 0], color=color, lw=2, solid_capstyle='round')
        ax.set_xlim(-14.5, 14.5)
        ax.set_ylim(-50, 50)
        ax.grid(axis='y', color=GRID, lw=.8)
        ax.set_axisbelow(True)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        ax.set_xlabel('Offset of an unchanged A/T site from the substitution (bp)', fontsize=8)
        ax.set_title(('Outside A runs' if ctx == 'other' else 'In or next to A runs') +
                     ': change in calls in heavily methylated molecules (pp)', loc='left', fontsize=9, color=INK)
        ax.text(0, 44, 'masked ±4 bp', ha='center', fontsize=7.5, color=INK2)
        if i == 0:
            ax.legend(frameon=False, fontsize=7.5, loc='lower right')
    fig.text(.01, .045, 'Source: docs/evidence/site_methylation_mask4_20261021 (primary run; SNVs with ≥300 molecules). Activity: published LDLR MPRA. '
             'Single preparation; U/M classes are same-reaction proxies, not prepared controls.', fontsize=7, color=MUTED)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'svg'):
        fig.savefig(f'{output}.{ext}', dpi=170, bbox_inches='tight')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--bundle', default=str(EVIDENCE / 'site_methylation_mask4_20261021'))
    p.add_argument('--output', required=True)
    a = p.parse_args()
    bundle = Path(a.bundle)
    kernel = []
    with open(bundle / 'effects' / 'kernel.tsv') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            if r['measure'] == 'M' and r['at_change'] in ('gain', 'loss') and r['context'] in ('other', 'homopolymer') and abs(int(r['offset'])) <= 14:
                kernel.append(dict(context=r['context'], at_change=r['at_change'], offset=int(r['offset']), mean=100 * float(r['mean_delta'])))
    plot(load(bundle), kernel, a.output)
    print(a.output + '.png')


if __name__ == '__main__':
    main()

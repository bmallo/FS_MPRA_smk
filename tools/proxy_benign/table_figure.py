"""Redraw the user's gene table (mechanism, ClinVar P/LP and B/LB, Fiber-seq accessible lines) with one added column:
proxy benign SNVs as a stacked bar, coloured by the number of supporting proxy benign classes (rules.json).

Table values are transcribed from the user's figure (2026-11-06); the proxy column comes from proxy.py counts.tsv.
User edits 2026-11-06: HEL 92.1.7 removed from the Fiber-seq accessible lines; RBM8A row removed (PAH kept, 15
genes); APC mechanism shown as LOE (was mLOE).
--compact (with --v2): grant-sized version (7 in wide), no footnote.
--v2: AlphaGenome-only proxy benign counts in the 350 bp TSS-side window (proxy_350.py, rules.json v2), one colour.
Usage: python tools/proxy_benign/table_figure.py --root results/proxy_benign_20261106 [--counts | --v2 [--compact]]
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
import pandas as pd  # noqa: E402

plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42, 'svg.fonttype': 'none'})  # editable text
ROWS = [  # gene, mechanisms, P/LP, B/LB, Fiber-seq accessible (user's table)
    ('GATA2', ['LOE'], 2, 85, 'WTC-11'),
    ('SHH', ['LOE', 'GOE'], 17, 27, 'HAP1, WTC-11'),
    ('APC', ['LOE'], 3, 25, 'HEK293T, HAP1, Jurkat, WTC-11'),
    ('MLH1', ['LOE'], 1, 17, 'HEK293T, HAP1, Jurkat, WTC-11'),
    ('MSH2', ['LOE'], 1, 14, 'HEK293T, HAP1, Jurkat, WTC-11'),
    ('GALT', ['LOE'], 0, 12, 'HEK293T, HAP1, Jurkat, WTC-11'),
    ('VHL', ['LOE'], 1, 10, 'HEK293T, HAP1, Jurkat, WTC-11'),
    ('HBB', ['LOE'], 25, 7, 'HUDEP-2'),
    ('LDLR', ['LOE'], 12, 7, 'HEK293T, HAP1, Jurkat, WTC-11'),
    ('TERT', ['GOE'], 1, 7, 'HEK293T, HAP1, Jurkat, WTC-11'),
    ('ANKRD26', ['GOE'], 10, 4, 'HEK293T, HAP1, Jurkat, WTC-11, HUDEP-2'),
    ('HBG1', ['GOE'], 6, 1, 'K562'),
    ('SLC39A4', ['LOE'], 1, 1, 'HEK293T, WTC-11'),
    ('NIPBL', ['LOE'], 1, 0, 'HEK293T, HAP1, Jurkat, WTC-11'),
    ('PAH', ['LOE'], 0, 7, 'HEK293T, Jurkat, WTC-11')]
MECH = {'LOE': ('#009E73', 'white'), 'GOE': ('#CC79A7', 'white'), 'mLOE': ('#F0E442', 'black')}
PLP, BLB = '#D55E00', '#0072B2'
LEVELS = {1: '#c6dbef', 2: '#6baed6', 3: '#08306b'}
# layout: figure size (in), x range, column x, bar widths, font sizes, row band half-height
FULL = dict(fig_w=15.6, row_in=0.45, extra_in=1.5, xmax=156, X=dict(gene=0.6, mech=13.5, plp=27, blb=41, proxy=57,
            fs=101), W=(9, 10.5, 33), hdr=13, gene=13, pill=11.5, num=12, fs=11.5, pill_w=4.0, half=0.3, rule=1.5)
COMPACT = dict(fig_w=7.0, row_in=0.2, extra_in=0.25, xmax=100, X=dict(gene=0.4, mech=11.5, plp=24, blb=35,
               proxy=47.5, fs=66), W=(6, 8, 18), hdr=7.5, gene=7.5, pill=5.8, num=6.5, fs=6.5, pill_w=3.5, half=0.33,
               rule=0.8)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--counts', action='store_true', help='bar in SNV counts instead of percent of possible SNVs')
    ap.add_argument('--v2', action='store_true', help='AlphaGenome-only, 350 bp window (v2_350bp/counts.tsv)')
    ap.add_argument('--cut350', action='store_true', help='with --v2: top-20% cutoff from the 350 bp window')
    ap.add_argument('--compact', action='store_true', help='with --v2: grant-sized, no footnote')
    a = ap.parse_args()
    root = Path(a.root)
    if a.v2:
        root = root / 'v2_350bp'
    c = pd.read_csv(root / ('counts_cut350.tsv' if a.cut350 else 'counts.tsv'), sep='\t').set_index('element')

    n = len(ROWS)
    L = COMPACT if a.compact else FULL
    fig = plt.figure(figsize=(L['fig_w'], L['row_in'] * (n + 1) + L['extra_in']))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, L['xmax'])
    ax.set_ylim(n + (1.0 if a.compact else 2.0), -0.2)
    ax.axis('off')
    X = L['X']
    W_PLP, W_BLB, W_PROXY = L['W']
    hdr = dict(fontsize=L['hdr'], fontweight='bold', va='center')
    ax.text(X['gene'], 0.35, 'Gene', **hdr)
    ax.text(X['mech'], 0.35, 'Mechanism', **hdr)
    ax.text(X['plp'], 0.35, 'P/LP', **hdr)
    ax.text(X['blb'], 0.35, 'B/LB', **hdr)
    ax.text(X['proxy'], 0.35 if a.v2 else 0.15, 'Proxy benign', **hdr)   # v2: no subtitle, footnote explains
    if not a.v2:
        unit = 'SNVs' if a.counts else '% of possible SNVs'
        ax.text(X['proxy'], 0.62, f'{unit}; colour = number of supporting classes', fontsize=8.5, va='center',
                color='#444444')
        lx = X['proxy'] + 15
        for k, col in LEVELS.items():
            ax.add_patch(plt.Rectangle((lx, 0.02), 1.3, 0.28, color=col, lw=0))
            ax.text(lx + 1.7, 0.16, str(k), fontsize=10, va='center')
            lx += 4.2
    ax.text(X['fs'], 0.35, 'Fiber-seq accessible', **hdr)
    ax.plot([0.3, L['xmax'] - 0.5], [0.95, 0.95], color='#333333', lw=L['rule'])

    vmax_plp = max(r[2] for r in ROWS)
    vmax_blb = max(r[3] for r in ROWS)
    if a.v2:
        tot = {r[0]: c.loc[r[0], 'proxy_benign'] for r in ROWS}
        vmax_proxy = 1050                                    # bar length = share of the 350 bp window's SNVs
    else:
        tot = {r[0]: sum(c.loc[r[0], f'level{k}' + ('' if a.counts else '_pct')] for k in LEVELS) for r in ROWS}
        vmax_proxy = max(tot.values()) if a.counts else 100
    for i, (gene, mechs, plp, blb, fs) in enumerate(ROWS):
        y = i + 1.5
        if i % 2 == 0:
            ax.add_patch(plt.Rectangle((0.3, y - 0.5), L['xmax'] - 0.8, 1, color='#f0f0f0', lw=0, zorder=0))
        ax.text(X['gene'], y, gene, fontsize=L['gene'], style='italic', va='center')
        mx = X['mech']
        for m in mechs:
            bg, fg = MECH[m]
            wd = L['pill_w'] * (1.4 if m == 'mLOE' else 1)
            ax.add_patch(FancyBboxPatch((mx, y - L['half']), wd, 2 * L['half'],
                                        boxstyle='round,pad=0,rounding_size=0.25', color=bg, lw=0))
            ax.text(mx + wd / 2, y, m, fontsize=L['pill'], fontweight='bold', color=fg, ha='center', va='center')
            mx += wd + 0.15 * L['pill_w']
        for x0, v, vmax, wmax, col in ((X['plp'], plp, vmax_plp, W_PLP, PLP), (X['blb'], blb, vmax_blb, W_BLB, BLB)):
            wd = wmax * v / vmax
            ax.add_patch(plt.Rectangle((x0, y - L['half']), max(wd, 0.12), 2 * L['half'], color=col, lw=0))
            ax.text(x0 + wd + 0.4, y, str(v), fontsize=L['num'], va='center')
        x0 = X['proxy']
        segs = [(tot[gene], LEVELS[2])] if a.v2 else \
            [(c.loc[gene, f'level{k}' + ('' if a.counts else '_pct')], col) for k, col in LEVELS.items()]
        for v, col in segs:
            wd = W_PROXY * v / vmax_proxy
            ax.add_patch(plt.Rectangle((x0, y - L['half']), wd, 2 * L['half'], color=col, lw=0))
            x0 += wd
        lab = f'{int(tot[gene]):,}' if (a.counts or a.v2) else f'{tot[gene]:.0f}%'
        ax.text(x0 + 0.4, y, lab, fontsize=L['num'] if a.compact else 11, va='center')
        ax.text(X['fs'], y, fs, fontsize=L['fs'], fontweight='bold', va='center')
    foot = ('Proxy benign = SNVs in a 350 bp window containing the TSS whose predicted RNA and DNase changes '
            '(AlphaGenome) are both outside the top 20% of changes across the ' +
            ('350 bp window' if a.cut350 else '1 kb element') + '.\nTSS = MANE Select (APC: promoter 1B); '
            'the window is kept inside the 1 kb element; GATA2 and SHH enhancers: window centred on the element. '
            'SNVs with a ≥ 1-star ClinVar P/LP or B/LB record are excluded.') if a.v2 else None
    if not a.compact:
        ax.text(0.6, n + 1.45, foot or 'Proxy benign classes: (1) in gnomAD at any frequency; (2) tested in GTEx (MAF ≥ 1%) and not '
                'an eQTL for any gene in any tissue; (3) no large AlphaGenome change (RNA and DNase both outside the '
                'element\'s top 20%;\nFiber-seq cell lines available in AlphaGenome). Denominator = all possible SNVs in '
                'the element; SNVs with a ≥ 1-star ClinVar P/LP or B/LB record are not counted.', fontsize=9,
                color='#444444', va='center', linespacing=1.4)
    stem = 'table_proxy_benign' + ('_350bp' + ('_cut350' if a.cut350 else '') + ('_compact' if a.compact else '')
                                   if a.v2 else '_counts' if a.counts else '')
    fig.savefig(root / f'{stem}.png', dpi=300 if a.compact else 200)
    fig.savefig(root / f'{stem}.pdf')
    fig.savefig(root / f'{stem}.svg')
    print('written', root / f'{stem}.png')


if __name__ == '__main__':
    main()

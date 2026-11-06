"""Truth-set figure: confident ClinVar seeds and non-ClinVar candidates per element (definitions: rules.json).

Reads counts.tsv and funnel.tsv from evidence.py (--input). Panels:
  a  per element: possible SNVs; arm A confident ClinVar (>= min stars) P/LP and B/LB SNVs; arm B likely pathogenic and
     likely benign candidates from non-ClinVar evidence
  b  evidence used by each group, with MPRA and plasmid Fiber-seq marked as to be added
  c  evidence steps pooled over elements (benign and pathogenic candidates)
Writes truth_set_figure.png/.pdf to --output.
--no-alphagenome: draft without AlphaGenome (matrix shows it as pending).
Usage: python tools/truth_set_figure/figure.py --input <evidence dir> --output <dir> [--no-alphagenome]
"""
import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

PATH, BEN = '#b2182b', '#2166ac'
PATH_L, BEN_L = '#f4a582', '#92c5de'
GREY = '#666666'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 6.5, 'axes.linewidth': 0.5,
                     'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'pdf.fonttype': 42})


def bar_col(ax, y, v, color, title, sub=None, lim=None):
    lim = lim or max(1, v.max())
    ax.barh(y, v, color=color, height=0.62)
    for yi, vi in zip(y, v):
        ax.text(vi + lim * 0.03, yi, f'{vi:,}', va='center', ha='left', fontsize=6)
    ax.set_xlim(0, lim * 1.35)
    ax.set_ylim(len(y) - 0.5, -0.5)
    ax.set_yticks([])
    ax.set_xticks([])
    for s in ('top', 'right', 'bottom'):
        ax.spines[s].set_visible(False)
    ax.set_title(title, fontsize=6.5, loc='left', pad=2 if sub is None else 9)
    if sub:
        ax.text(0, 1.0, sub, transform=ax.transAxes, fontsize=5.3, color=GREY, va='bottom')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--no-alphagenome', action='store_true')
    a = ap.parse_args()
    inp, out = Path(a.input), Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    c = pd.read_csv(inp / 'counts.tsv', sep='\t')
    f = pd.read_csv(inp / 'funnel.tsv', sep='\t')
    meta = json.load(open(inp / 'summary.json'))
    ms, benign = meta.get('min_stars', 2), meta.get('benign', 'common_not_eqtl')
    star = f'≥ {ms} star' + ('s' if ms > 1 else '')
    conf = 'confident' if 'clinvar_plp_indel_confident' in c else '2star'
    n = len(c)
    y = list(range(n))

    fig = plt.figure(figsize=(7.2, 6.4))
    # panel a: one column of bars per group
    top, bot = 0.87, 0.50
    lab = fig.add_axes([0.02, bot, 0.11, top - bot])
    lab.set_ylim(n - 0.5, -0.5)
    lab.axis('off')
    for yi, r in zip(y, c.itertuples()):
        lab.text(1, yi, r.element, ha='right', va='center', fontsize=6.5, fontweight='bold')
    fig.text(0.02, 0.982, 'a', fontsize=9, fontweight='bold')
    fig.text(0.05, 0.982, 'Truth-set sources per element (SNVs)', fontsize=7.5)
    cols = [(0.14, 0.11, c.possible_snvs, '#bbbbbb', 'All possible\nSNVs', None),
            (0.27, 0.13, c.arm_a_plp_snvs, PATH, 'ClinVar P/LP', star),
            (0.42, 0.13, c.arm_a_blb_snvs, BEN, 'ClinVar B/LB', star),
            (0.59, 0.17, c.lp_candidates, PATH_L, 'Likely pathogenic\ncandidates', None),
            (0.79, 0.17, c.lb_candidates, BEN_L, 'Likely benign\ncandidates', None)]
    lim_a = max(c.arm_a_plp_snvs.max(), c.arm_a_blb_snvs.max(), 1)
    for x0, wd, v, col, title, sub in cols:
        ax = fig.add_axes([x0, bot, wd, top - bot])
        lim = lim_a if col in (PATH, BEN) else None
        bar_col(ax, y, v.values, col, title, sub, lim)
    fig.text(0.41, 0.948, 'A  Confident in ClinVar', fontsize=6.8, ha='center', fontweight='bold')
    fig.text(0.775, 0.948, 'B  Not confident in ClinVar: non-ClinVar evidence', fontsize=6.8, ha='center',
             fontweight='bold')
    fig.add_artist(plt.Line2D([0.27, 0.55], [0.943, 0.943], color='k', lw=0.5))
    fig.add_artist(plt.Line2D([0.59, 0.96], [0.943, 0.943], color='k', lw=0.5))
    notes = []
    below = f'{ms - 1}star'                                     # the star level just under the cutoff
    plp_b, blb_b = c[f'clinvar_plp_snv_{below}'], c[f'clinvar_blb_snv_{below}']
    if (plp_b + blb_b).sum():
        notes.append(f'Not counted in A: {ms - 1}-star P/LP ' + (', '.join(f'{e} {k}' for e, k in zip(c.element, plp_b)
                     if k) or 'none') + f'; {ms - 1}-star B/LB ' +
                     (', '.join(f'{e} {k}' for e, k in zip(c.element, blb_b) if k) or 'none') + '.')
    ind = c[f'clinvar_plp_indel_{conf}'].sum() + c[f'clinvar_blb_indel_{conf}'].sum()
    notes.append(f'Confident ClinVar indels (not in bars): {ind}. B = all SNVs without a {star} P/LP or B/LB record '
                 f'(VUS, conflicting and lower-star records stay in B; ClinVar is not used as evidence there).')
    if c.pai_unscored.sum():
        notes.append('Not scored by PromoterAI (outside TSS ± 500 bp), so cannot be a candidate: ' +
                     ', '.join(f'{e} {k:,}' for e, k in zip(c.element, c.pai_unscored) if k) + ' SNVs.')
    fig.text(0.02, 0.465, '\n'.join(notes), fontsize=5.3, color=GREY, va='top', wrap=True)

    # panel b: evidence matrix
    fig.text(0.02, 0.40, 'b', fontsize=9, fontweight='bold')
    fig.text(0.05, 0.40, 'Evidence per group', fontsize=7.5)
    ax = fig.add_axes([0.03, 0.04, 0.42, 0.33])
    layers = [f'ClinVar {star}', 'gnomAD frequency', 'GTEx eQTL (any gene)', 'PromoterAI', 'AlphaGenome',
              'MPRA', 'Plasmid Fiber-seq']
    groups = ['ClinVar\nP/LP, B/LB', 'Likely\npathogenic', 'Likely\nbenign']
    rules = [['seed', 'context', 'context', 'context', 'context', 'add', 'add'],
             ['', 'absent or\n< 0.01%', '', '|score| ≥ 0.2', 'top 20%', 'add', 'add'],
             ['', '> 1% in a\ngroup', 'not an\neQTL', '|score| < 0.1', 'not top\n20%', 'add', 'add']]
    if a.no_alphagenome:
        rules[1][4] = rules[2][4] = 'pending'
    if benign == 'in_gnomad':                                  # draft 3: no frequency floor, GTEx not used
        rules[2][1] = 'in gnomAD\n(any AF)'
        layers.pop(2)
        rules = [r[:2] + r[3:] for r in rules]
    ax.set_xlim(-0.5, len(groups) - 0.5)
    ax.set_ylim(len(layers) - 0.5, -1.1)
    ax.axis('off')
    for j, gname in enumerate(groups):
        ax.text(j, -0.85, gname, ha='center', va='center', fontsize=6, fontweight='bold',
                color=[PATH, PATH_L, BEN_L][j] if j else 'k')
        for i, (layer, rule) in enumerate(zip(layers, rules[j])):
            if rule == 'add':
                ax.add_patch(plt.Rectangle((j - 0.42, i - 0.38), 0.84, 0.76, fill=False, ls='--', lw=0.6,
                                           ec='#1b7837'))
                ax.text(j, i, 'to add', ha='center', va='center', fontsize=5.5, color='#1b7837')
            elif rule == 'seed':
                ax.add_patch(plt.Rectangle((j - 0.42, i - 0.38), 0.84, 0.76, color='#dddddd', lw=0))
                ax.text(j, i, 'classification', ha='center', va='center', fontsize=5.5)
            elif rule == 'pending':
                ax.add_patch(plt.Rectangle((j - 0.42, i - 0.38), 0.84, 0.76, fill=False, ls=':', lw=0.6, ec=GREY))
                ax.text(j, i, 'pending', ha='center', va='center', fontsize=5.5, color=GREY)
            elif rule == 'context':
                ax.text(j, i, 'shown', ha='center', va='center', fontsize=5.5, color=GREY)
            elif rule:
                ax.add_patch(plt.Rectangle((j - 0.42, i - 0.38), 0.84, 0.76, color='#eeeeee', lw=0))
                ax.text(j, i, rule, ha='center', va='center', fontsize=5.3, linespacing=0.9)
    for i, layer in enumerate(layers):
        ax.text(-0.5, i, layer, ha='right', va='center', fontsize=6)
    ax.add_artist(plt.Line2D([-0.5, len(groups) - 0.5], [len(layers) - 2.5] * 2, color='k', lw=0.4, ls=':'))
    ax.set_position([0.19, 0.04, 0.28, 0.33])

    # panel c: pooled evidence steps
    fig.text(0.53, 0.40, 'c', fontsize=9, fontweight='bold')
    fig.text(0.56, 0.40, 'Evidence steps, all elements pooled (SNVs)', fontsize=7.5)
    p = f.groupby(['side', 'step', 'filter'], sort=False).n.sum().reset_index()
    for k, (side, col, y0) in enumerate((('pathogenic', PATH_L, 0.25), ('benign', BEN_L, 0.035))):
        q = p[p.side == side].sort_values('step')
        ax = fig.add_axes([0.73, y0, 0.20, 0.11 if side == 'pathogenic' else 0.14])
        yy = range(len(q))
        ax.barh(list(yy), q.n, color=col, height=0.65, log=True)
        ax.set_ylim(len(q) - 0.5, -0.5)
        ax.set_yticks(list(yy))
        ax.set_yticklabels(q['filter'], fontsize=5.6)
        for yi, v in zip(yy, q.n):
            ax.text(v * 1.15 if v else 1.2, yi, f'{v:,}', va='center', fontsize=5.6)
        ax.set_xlim(0.8, q.n.max() * 8)
        ax.tick_params(axis='x', labelsize=5.3)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.set_title('Likely pathogenic' if side == 'pathogenic' else 'Likely benign', fontsize=6.3, loc='left',
                     color=PATH if side == 'pathogenic' else BEN)
    fig.savefig(out / 'truth_set_figure.png', dpi=300)
    fig.savefig(out / 'truth_set_figure.pdf')
    print('written', out / 'truth_set_figure.png')


if __name__ == '__main__':
    main()

"""Draft panel b: ClinVar LDLR promoter SNVs by mechanism category (classify.py output; rules classify_rules.json).

One 100%-stacked horizontal bar per ClinVar class. Colour encodes the expression direction (orange family = expression
down, blue family = up, grey = no expression change); shade separates the categories within a direction, and each
segment carries its count, so identity never rests on colour alone.
Usage: python tools/r01_figure/fig_panel_b.py --classify <dir> --output <dir> [--column category]
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

_s = importlib.util.spec_from_file_location('fig_repeat3', Path(__file__).with_name('fig_repeat3.py'))
fr = importlib.util.module_from_spec(_s)
_s.loader.exec_module(fr)
INK, INK2, MUTED = fr.INK, fr.INK2, fr.MUTED
CATS = [('Activator disrupted', '#b8461c', 'white'), ('Repressor gained', '#eb6834', 'white'),
        ('Expression down, mechanism unresolved', '#f6b89a', INK),
        ('Activator gained', '#1c5aa6', 'white'), ('Repressor lost', '#2a78d6', 'white'),
        ('Expression up, mechanism unresolved', '#a9c6ef', INK),
        ('Chromatin change only', '#8a8883', 'white'), ('No detectable effect', '#d6d4ce', INK),
        ('Not assessable by Fiber-seq', '#f1efe9', INK2)]
CLASSES = ['P/LP', 'Conflicting', 'VUS', 'B/LB']
NAMES = {'P/LP': 'Pathogenic /\nlikely pathogenic', 'Conflicting': 'Conflicting', 'VUS': 'VUS',
         'B/LB': 'Benign /\nlikely benign'}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--classify', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--column', default='category')
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    t = pd.read_csv(Path(a.classify) / 'variants.tsv', sep='\t')
    t = t[t.clin.isin(CLASSES)]
    ct = t.groupby(['clin', a.column]).size().unstack(fill_value=0).reindex(CLASSES)
    ct = ct.reindex(columns=[c for c, _, _ in CATS], fill_value=0)
    ct.to_csv(out / 'panel_b_counts.tsv', sep='\t')
    t[['snv', 'c_dot', 'clin', 'clinvar', 'act', 'nuc_call', 'fp_level', 'fp_sign', 'fp_mean_z', 'broken', 'created',
       'at_pm5', a.column, a.column.replace('category', 'evidence')]].to_csv(out / 'panel_b_variants.tsv', sep='\t',
                                                                            index=False)

    fig, ax = plt.subplots(figsize=(7.2, 2.4))
    fig.subplots_adjust(left=0.15, right=0.93, top=0.88, bottom=0.36)
    for i, cls in enumerate(CLASSES):
        y = len(CLASSES) - 1 - i
        n = ct.loc[cls].sum()
        left = 0.0
        for cat, col, tc in CATS:
            k = ct.loc[cls, cat]
            if k == 0:
                continue
            wdt = k / n
            ax.barh(y, wdt, left=left, height=0.72, color=col, edgecolor='white', lw=1.0,
                    hatch='////' if cat == 'Not assessable by Fiber-seq' else None)
            if wdt >= 0.045:
                ax.text(left + wdt / 2, y, str(k), ha='center', va='center', fontsize=6, color=tc)
            left += wdt
        ax.text(1.005, y, f'n={n}', ha='left', va='center', fontsize=6, color=INK2)
    ax.set_yticks(range(len(CLASSES)))
    ax.set_yticklabels([NAMES[c] for c in reversed(CLASSES)], fontsize=6.5)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(['0', '25%', '50%', '75%', '100%'])
    ax.tick_params(axis='y', length=0)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.set_xlabel('Share of ClinVar SNVs in each class (number = SNVs)')
    from matplotlib.patches import Patch
    h = [Patch(facecolor=c, edgecolor='white', hatch='////' if n == 'Not assessable by Fiber-seq' else None)
         for n, c, _ in CATS]
    fig.legend(h, [n for n, _, _ in CATS], loc='lower center', ncol=3, frameon=False, fontsize=5.8,
               bbox_to_anchor=(0.55, 0.0), handletextpad=0.4, columnspacing=1.2)
    fr.panel_letter(ax, 'b')
    fig.savefig(out / 'panel_b.pdf')
    fig.savefig(out / 'panel_b.png', dpi=300)


if __name__ == '__main__':
    main()

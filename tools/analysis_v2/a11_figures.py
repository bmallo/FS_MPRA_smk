"""Step A11 slide figures from the a11_motifs_function tables (no recomputation). HepG2 run 3 only.

Usage: python tools/analysis_v2/a11_figures.py --dir results/analysis_v2_20261028/a11_motifs_function
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('a7f', Path(__file__).with_name('a7_figures.py'))
a7f = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a7f)
INK, BODY, BLUE, ORANGE, GREY = a7f.INK, a7f.BODY, a7f.BLUE, a7f.ORANGE, a7f.GREY


def function(d, out):
    f = pd.read_csv(d / 'function_snvs.tsv', sep='\t').dropna(subset=['act', 'C1'])
    fn = pd.read_csv(d / 'function.tsv', sep='\t')
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.6), gridspec_kw=dict(wspace=0.25))
    a = ax[0]
    o, r = f[~f.region], f[f.region]
    a.scatter(o.C1, o.act, s=12, color='#A7AFBA', label=f'other SNVs ({len(o)})')
    a.scatter(r.C1, r.act, s=26, color=ORANGE, label=f'oligo 143-151, element-wide ({len(r)})')
    a.axhline(0, color=GREY, lw=0.8)
    a.axvline(0, color=GREY, lw=0.8)
    row = fn[(fn.source == 'LDLR') & (fn.model == 'F1') & ~fn.position_fe].iloc[0]
    a.set_title(f'All SNVs · Spearman {row.spearman:.2f}', loc='left', color=INK, fontsize=13)
    a.set_xlabel('element-wide shift: unmarked change > 50 bp away (points)')
    a.set_ylabel('published MPRA effect (log2, LDLR)')
    a.legend(frameon=False, fontsize=11, loc='lower left')
    a = ax[1]
    g = f.copy()
    g['x'] = g.C1 - g.groupby('pos').C1.transform('mean')
    g['y'] = g.act - g.groupby('pos').act.transform('mean')
    a.scatter(g.x, g.y, s=12, color=BLUE, alpha=0.7)
    row = fn[(fn.source == 'LDLR') & (fn.model == 'F1') & fn.position_fe].iloc[0]
    xx = np.linspace(g.x.min(), g.x.max(), 2)
    a.plot(xx, row.slope * xx, color=INK, lw=1.5)
    a.axhline(0, color=GREY, lw=0.8)
    a.axvline(0, color=GREY, lw=0.8)
    a.set_title(f'Alleles at the same position · slope {row.slope:.2f} log2 per point '
                f'({row.ci_lo:.2f} to {row.ci_hi:.2f})', loc='left', color=INK, fontsize=13)
    a.set_xlabel('shift minus the position mean (points)')
    a.set_ylabel('effect minus the position mean (log2)')
    for a in ax:
        a.spines[['top', 'right']].set_visible(False)
    fig.savefig(out / 'a11_function.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def dose(d, out):
    s = pd.read_csv(d / 'dose_snvs.tsv', sep='\t')
    t = pd.read_csv(d / 'dose.tsv', sep='\t')
    names = list(dict.fromkeys(t.motif))
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5), gridspec_kw=dict(hspace=0.62, wspace=0.22))
    for a, name in zip(ax.ravel(), names):
        g = s[(s.motif == name)].dropna(subset=['calls_D'])
        for dr, c, lab in ((0, BLUE, 'gains an A/T or neither'), (1, ORANGE, 'drops its base')):
            h = g[g.drops == dr]
            a.scatter(h.d_pwm, 100 * h.calls_D, s=28, color=c, label=lab)
        r = t[(t.motif == name) & (t.endpoint == 'calls_D')].iloc[0]
        if np.isfinite(r.get('slope', np.nan)):
            xx = np.linspace(min(g.d_pwm.min(), 0), max(g.d_pwm.max(), 0), 2)
            a.plot(xx, 100 * (r.slope * xx + (g.calls_D - r.slope * g.d_pwm).mean()), color=INK, lw=1.2)
            sub = f'slope {100 * r.slope:.2f} · q {r.q:.2g} · n {int(r.snvs)}'
        else:
            sub = f'n {int(r.snvs)}: too few to fit'
        a.axhline(0, color=GREY, lw=0.8)
        a.set_title(f'{name}{" (v1-derived)" if r.set == "secondary" else ""}\n{r.tf} {r.site} · {sub}',
                    loc='left', color=INK, fontsize=11.5)
        a.set_xlabel('motif score change (log2-odds)')
        if a in ax[:, 0]:
            a.set_ylabel('10-40 bp coverage change\n(points)')
        a.spines[['top', 'right']].set_visible(False)
    ax[0, 0].legend(frameon=False, fontsize=11, loc='upper left', bbox_to_anchor=(0, -0.22), ncol=2)
    fig.text(0.99, 0.005, 'slope in points per log2-odds unit', ha='right', color=BODY, fontsize=11)
    fig.savefig(out / 'a11_dose.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    out = d / 'figures'
    out.mkdir(exist_ok=False)
    a7f.a4f.style()
    function(d, out)
    dose(d, out)


if __name__ == '__main__':
    main()

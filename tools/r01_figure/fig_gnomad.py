"""Draft figure for the Aim 1B preview: LDLR promoter SNVs by gnomAD presence vs ClinVar classes (gnomad_check.py output).

a  gnomAD allele count vs HepG2 MPRA expression change for every designed SNV (absent from gnomAD at AC 0, left).
b  expression change and nucleosome-size protection change by class (dots = SNVs, bar = median).
Usage: python tools/r01_figure/fig_gnomad.py --check <gnomad_check dir> --output <dir>
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_s = importlib.util.spec_from_file_location('fig_repeat3', Path(__file__).with_name('fig_repeat3.py'))
fr = importlib.util.module_from_spec(_s)
_s.loader.exec_module(fr)
INK, INK2, MUTED = fr.INK, fr.INK2, fr.MUTED
C_PLP, C_GNOMAD, C_BLB = '#eb6834', '#2a78d6', '#1baf7a'
CLASSES = [('clinvar_plp', 'ClinVar P/LP'), ('clinvar_blb', 'ClinVar B/LB'), ('gnomad_ac2plus', 'gnomAD AC ≥ 2'),
           ('gnomad_singleton', 'gnomAD singleton'), ('not_in_gnomad', 'Not in gnomAD')]
COLS = {'clinvar_plp': C_PLP, 'clinvar_blb': C_BLB, 'gnomad_ac2plus': C_GNOMAD, 'gnomad_singleton': C_GNOMAD,
        'not_in_gnomad': MUTED}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--check', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    t = pd.read_csv(Path(a.check) / 'snvs.tsv', sep='\t')
    fig = plt.figure(figsize=(7.2, 3.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1], left=0.07, right=0.985, top=0.88, bottom=0.2, wspace=0.85)

    ax = fig.add_subplot(gs[0, 0])
    fr.tidy(ax, 'both')
    x = np.where(t.ac > 0, np.log10(t.ac.clip(lower=1)), -0.6) + np.random.default_rng(3).uniform(-0.12, 0.12, len(t))
    base = ~(t.class_clinvar_plp | t.class_clinvar_blb)
    ax.scatter(x[base & (t.ac == 0)], t.act[base & (t.ac == 0)], s=4, color='#d6d4ce', lw=0)
    ax.scatter(x[base & (t.ac > 0)], t.act[base & (t.ac > 0)], s=7, color=C_GNOMAD, lw=0, alpha=0.8)
    ax.scatter(x[t.class_clinvar_plp], t.act[t.class_clinvar_plp], s=20, marker='D', color=C_PLP, edgecolor='white', lw=0.4)
    ax.scatter(x[t.class_clinvar_blb], t.act[t.class_clinvar_blb], s=20, marker='o', facecolor='white',
               edgecolor=C_BLB, lw=1.0)
    for v, dx, dy in (('68C>T', -0.9, -0.45), ('98C>T', 0.35, 0.35), ('146C>G', 0.45, -0.05), ('150C>G', 0.45, -0.2),
                      ('143C>G', -0.15, 0.55)):
        r = t[t.snv == v].index[0]
        ax.annotate(v, (x[r], t.act[r]), xytext=(x[r] + dx, t.act[r] + dy), fontsize=5.5, color=INK2,
                    arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.4))
    ax.axhline(0, color=INK2, lw=0.5)
    ax.set_xticks([-0.6, 0, 1, 2])
    ax.set_xticklabels(['absent', '1', '10', '100'])
    ax.set_xlabel('gnomAD v4.1 allele count (log scale)')
    ax.set_ylabel('Expression change (log2, HepG2 MPRA)')
    ax.text(0.02, 0.98, f'{int((t.ac > 0).sum())} of {len(t)} SNVs in gnomAD;\n'
            f'{int((t.af >= 1e-4).sum())} with AF ≥ 1e-4, {int((t.af >= 0.01).sum())} with AF ≥ 1%',
            transform=ax.transAxes, ha='left', va='top', fontsize=5.5, color=INK2)
    fr.panel_letter(ax, 'a')

    for j, (col, lab, letter) in enumerate((('act', 'Expression change (log2)', 'b'),
                                            ('nuc_D_pts', 'Change in ≥81-bp\nprotection (pts)', 'c'))):
        axx = fig.add_subplot(gs[0, 1 + j])
        fr.tidy(axx, 'x')
        axx.axvline(0, color=INK2, lw=0.5)
        for i, (k, name) in enumerate(CLASSES):
            s = t[t[f'class_{k}']][col].dropna()
            y = len(CLASSES) - i + np.random.default_rng(i).uniform(-0.2, 0.2, len(s))
            axx.scatter(s, y, s=5 if len(s) > 20 else 12, color=COLS[k], lw=0, alpha=0.6 if len(s) > 20 else 0.9)
            axx.scatter(s.median(), len(CLASSES) - i, s=60, marker='|', color=INK, lw=1.4, zorder=5)
        axx.set_yticks([len(CLASSES) - i for i in range(len(CLASSES))])
        axx.set_yticklabels([f'{n} (n={int(t[f"class_{k}"].sum())})' for k, n in CLASSES] if j == 0 else [], fontsize=6)
        axx.set_xlabel(lab)
        fr.panel_letter(axx, letter)
    fig.text(0.985, 0.005, 'draft · bar = median; P/LP and B/LB SNVs may also be in gnomAD', ha='right', fontsize=5,
             color=MUTED)
    fig.savefig(out / 'fig_gnomad.pdf')
    fig.savefig(out / 'fig_gnomad.png', dpi=300)


if __name__ == '__main__':
    main()

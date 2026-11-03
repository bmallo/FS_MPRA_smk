"""Draft R01 panel: the Repeat 3 Sp1 site, clinical classes, nucleosome-size protection and activity.

Reads repeat3.py outputs (+ A11 function table and the step-3 clinical table); writes fig_repeat3.pdf/.png and the
source data behind each panel (panel_*.tsv) to --output. Every number drawn comes from those tables.
Usage: python tools/r01_figure/fig_repeat3.py --root results/analysis_v2_20261028 --repeat3 <dir> --clinical <dir> --output <dir>
"""
import argparse
import subprocess
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#9a9893', '#e4e3df'
COL = {'G1': '#eb6834', 'G2': '#2a78d6', 'G3': '#1baf7a'}          # validated light palette slots 2, 1, 3
MARK = {'G1': 'D', 'G2': 'o', 'G3': '^'}
NAMES = {'G1': 'ClinVar P/LP or conflicting,\nSp1 site', 'G2': 'VUS / unreported,\nbreaks Sp1 site',
         'G3': 'Keeps or strengthens\nSp1 site', 'G4': 'Flanking SNVs\n(138–142, 152–157)',
         'G5': 'ClinVar benign /\nlikely benign (anywhere)', 'G6': 'ClinVar P/LP outside\nSp1 site (SRE-1, other)'}
SITE = (143, 151)
REGIONS = [('R1 Sp1', 89, 104), ('R2 SRE-1', 124, 139), ('R3 Sp1', 140, 155)]

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 8, 'xtick.labelsize': 6.5,
                     'ytick.labelsize': 6.5, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})


def tidy(ax, grid_axis='y'):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.grid(axis=grid_axis, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def panel_letter(ax, s):
    ax.text(-0.02, 1.04, s, transform=ax.transAxes, fontsize=10, fontweight='bold', ha='right', va='bottom', color=INK)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--repeat3', required=True)
    ap.add_argument('--clinical', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, r3, out = Path(a.root), Path(a.repeat3), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(['git', 'status', '--porcelain', 'tools/r01_figure'], capture_output=True, text=True).stdout
    stamp = f'draft {date.today()} · git {commit}{"+uncommitted" if dirty else ""}'

    fn = pd.read_csv(root / 'a11_motifs_function/function_snvs.tsv', sep='\t')
    import json
    gl = json.load(open(r3 / 'groups_membership.json'))
    grp = {v: g for g, vs in gl.items() for v in vs}
    fn['group'] = fn.snv.map(grp)
    clin = pd.read_csv(Path(a.clinical) / 'clinical_variants.tsv', sep='\t')
    clin = clin[clin.type == 'SNV'].set_index('designed').clinvar_class
    fn['clinvar'] = fn.snv.map(clin).fillna('')
    se = pd.read_csv(r3 / 'snv_endpoints.tsv', sep='\t')
    nuc = se[(se.endpoint == 'nucleosome') & (se.calls == 'confident') & (se['mask'] == 'A7')].copy()
    nuc['D_pts'], nuc['se_pts'] = 100 * nuc.D, 100 * nuc.se
    g = pd.read_csv(r3 / 'groups.tsv', sep='\t')
    hc = pd.read_csv(r3 / 'hek_clones.tsv', sep='\t')

    fig = plt.figure(figsize=(7.2, 8.2))
    gs_a = fig.add_gridspec(1, 1, left=0.085, right=0.985, top=0.955, bottom=0.775)
    gs_b = fig.add_gridspec(1, 2, left=0.225, right=0.985, top=0.63, bottom=0.425, wspace=0.10)
    gs_c = fig.add_gridspec(1, 2, left=0.085, right=0.985, top=0.33, bottom=0.075, wspace=0.30)

    # A: activity along the element, clinical classes
    ax = fig.add_subplot(gs_a[0, 0])
    tidy(ax)
    for name, lo, hi in REGIONS:
        ax.axvspan(lo - 0.5, hi + 0.5, color='#f1efe9', zorder=0, lw=0)
        ax.text((lo + hi) / 2, 1.02, name, transform=ax.get_xaxis_transform(), ha='center', va='bottom',
                fontsize=6, color=INK2)
    ax.axvspan(SITE[0] - 0.5, SITE[1] + 0.5, color='#fbe3d9', zorder=0, lw=0)
    rest = fn[fn.group.isna() | fn.group.isin(['G4'])]
    ax.scatter(rest.pos, rest.act, s=5, color='#c9c7c1', lw=0, zorder=2, label='Other SNVs')
    for k in ('G2', 'G3', 'G1'):
        t = fn[fn.group == k]
        ax.scatter(t.pos, t.act, s=16 if k == 'G1' else 10, marker=MARK[k], color=COL[k], lw=0.4,
                   edgecolor='white', zorder=4, label=NAMES[k].replace('\n', ' '))
    t = fn[fn.group == 'G6']
    ax.scatter(t.pos, t.act, s=16, marker='D', facecolor='white', edgecolor=INK, lw=0.7, zorder=4,
               label=NAMES['G6'].replace('\n', ' '))
    t = fn[fn.group == 'G5']
    ax.scatter(t.pos, t.act, s=14, marker='o', facecolor='white', edgecolor=MUTED, lw=0.7, zorder=4,
               label=NAMES['G5'].replace('\n', ' '))
    ax.axhline(0, color=INK2, lw=0.5)
    ax.set_xlim(20, 232)
    ax.set_xlabel('Promoter position (oligo; c. = oligo − 285)')
    ax.set_ylabel('Expression change\n(log2 variant / WT, HepG2 MPRA)')
    ax.legend(loc='upper center', ncol=3, frameon=False, fontsize=5.8, handletextpad=0.3, columnspacing=1.2,
              bbox_to_anchor=(0.5, -0.30))
    panel_letter(ax, 'a')
    fn[['snv', 'pos', 'act', 'act2', 'group', 'clinvar']].to_csv(out / 'panel_a_activity.tsv', sep='\t', index=False)

    # B: pooled nucleosome-size occupancy change and activity per group
    order = ['G1', 'G2', 'G3', 'G4', 'G5', 'G6']
    yy = {k: len(order) - i for i, k in enumerate(order)}
    axb = fig.add_subplot(gs_b[0, 0])
    axc = fig.add_subplot(gs_b[0, 1], sharey=axb)
    rows = []
    for axx, measure, col, xl in ((axb, 'nucleosome_confident_A7', 'D_pts', 'Change in molecules with a ≥81-bp\n'
                                   'protected stretch (pts, variant − WT)'),
                                  (axc, 'act', 'act', 'Expression change (log2, HepG2 MPRA)')):
        tidy(axx, 'x')
        axx.axvline(0, color=INK2, lw=0.5)
        for k in order:
            pts = nuc[nuc.group == k][col] if col == 'D_pts' else fn[fn.group == k].act
            jit = np.random.default_rng(1).uniform(-0.18, 0.18, len(pts))
            c = COL.get(k, MUTED)
            axx.scatter(pts, yy[k] + jit, s=6, color=c, alpha=0.55, lw=0, zorder=2)
            r = g[(g.group == k) & (g.measure == measure)].iloc[0]
            if col == 'D_pts':
                axx.errorbar(r.est, yy[k], xerr=[[r.est - r.lo], [r.hi - r.est]], fmt='none', ecolor=INK, lw=1.0,
                             capsize=0, zorder=4)
                axx.scatter(r.est, yy[k], s=26, marker=MARK.get(k, 's'), color=c, edgecolor=INK, lw=0.6, zorder=5)
                rows.append(dict(group=k, measure='nucleosome-size occupancy change, pooled (DL)', k=r.k,
                                 est=r.est, lo=r.lo, hi=r.hi))
            else:
                axx.scatter(r.est, yy[k], s=26, marker='|', color=INK, lw=1.2, zorder=5)
                rows.append(dict(group=k, measure='activity mean', k=r.k, est=r.est, lo=r.lo, hi=r.hi))
        axx.set_xlabel(xl)
    axb.set_yticks([yy[k] for k in order])
    axb.set_yticklabels([f'{NAMES[k]} (n={int(g[(g.group == k) & (g.measure == "act")].k.iloc[0])})' for k in order],
                        fontsize=6)
    plt.setp(axc.get_yticklabels(), visible=False)
    axb.text(0.02, 0.98, 'dot = SNV; black = random-effects\npooled estimate, bar = 95% CI', transform=axb.transAxes,
             ha='left', va='top', fontsize=5.5, color=INK2)
    axc.text(0.02, 0.02, 'dot = SNV; bar = group mean', transform=axc.transAxes, ha='left', va='bottom',
             fontsize=5.5, color=INK2)
    panel_letter(axb, 'b')
    pd.DataFrame(rows).to_csv(out / 'panel_b_groups.tsv', sep='\t', index=False)
    nuc[['snv', 'group', 'n', 'C', 'D_pts', 'se_pts', 'act']].to_csv(out / 'panel_b_snvs.tsv', sep='\t', index=False)

    # C: within the site, nucleosome-size change vs activity
    axd = fig.add_subplot(gs_c[0, 0])
    tidy(axd, 'both')
    s = nuc[nuc.group.isin(['G1', 'G2', 'G3'])].copy()
    for k in ('G2', 'G3', 'G1'):
        t = s[s.group == k]
        axd.errorbar(t.D_pts, t.act, xerr=1.96 * t.se_pts, fmt='none', ecolor=COL[k], alpha=0.5, lw=0.6, zorder=2)
        axd.scatter(t.D_pts, t.act, s=22 if k == 'G1' else 14, marker=MARK[k], color=COL[k], edgecolor='white',
                    lw=0.4, zorder=3)
    rho = stats.spearmanr(s.D_pts, s.act)
    for v, dx, dy in (('147T>G', 0.5, 0.12), ('144T>C', 0.6, 0.15), ('147T>A', 0.6, 0.18), ('151C>A', 0.6, -0.25),
                      ('149C>T', -3.2, 0.25), ('150C>G', -2.2, 0.55), ('145C>T', 1.2, 0.35), ('146C>G', 1.2, -0.1)):
        r = s[s.snv == v].iloc[0]
        axd.annotate(v, (r.D_pts, r.act), xytext=(r.D_pts + dx, r.act + dy), fontsize=5.5, color=INK2,
                     arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.4))
    axd.axhline(0, color=INK2, lw=0.5)
    axd.axvline(0, color=INK2, lw=0.5)
    axd.set_xlim(-5, 12.5)
    axd.set_xlabel('Change in molecules with a ≥81-bp protected stretch\n(pts, variant − WT; bar = 95% CI)')
    axd.set_ylabel('Expression change (log2)')
    axd.text(0.98, 0.97, f'27 SNVs in the Sp1 site (143–151)\nSpearman ρ = {rho.correlation:.2f}, '
             f'P = {rho.pvalue:.1g}', transform=axd.transAxes, ha='right', va='top', fontsize=6, color=INK)
    panel_letter(axd, 'c')
    s[['snv', 'group', 'D_pts', 'se_pts', 'act']].to_csv(out / 'panel_c_site.tsv', sep='\t', index=False)

    # D: HEK293T clones, element-wide inaccessibility far from the site
    axe = fig.add_subplot(gs_c[0, 1])
    tidy(axe, 'x')
    w = hc[hc.group == 'WT']
    q = w.far_resid_pts.quantile([0.05, 0.95])
    axe.axvspan(q.iloc[0], q.iloc[1], color='#efeeea', zorder=0, lw=0)
    axe.scatter(w.far_resid_pts, np.random.default_rng(2).uniform(-0.3, 0.3, len(w)), s=5, color='#c9c7c1', lw=0,
                zorder=2)
    site = hc[hc.group.isin(['G1', 'G2', 'G3'])].sort_values('far_resid_pts')
    ys = np.arange(1, len(site) + 1)
    for (i, r), y in zip(site.iterrows(), ys):
        axe.scatter(r.far_resid_pts, y, s=22 if r.group == 'G1' else 14, marker=MARK[r.group], color=COL[r.group],
                    edgecolor='white', lw=0.4, zorder=3)
        axe.text(r.far_resid_pts + 0.25, y, f'{r.label} ({r.molecules} mol.)', va='center', fontsize=5.5, color=INK2)
    axe.axvline(0, color=INK2, lw=0.5)
    axe.set_yticks([])
    axe.spines['left'].set_visible(False)
    axe.text(q.iloc[1] + 0.2, -0.05, f'WT clones (n={len(w)})', va='center', fontsize=5.5, color=INK2)
    axe.set_xlabel('HEK293T clones: inaccessible fraction >50 bp\nfrom the site (pts vs WT, labeling-matched)')
    axe.text(0.02, 0.98, 'grey band = WT clones 5–95%', transform=axe.transAxes, ha='left', va='top', fontsize=5.5,
             color=INK2)
    axe.set_xlim(min(-4, site.far_resid_pts.min() - 0.5), max(7.5, site.far_resid_pts.max() + 3.5))
    panel_letter(axe, 'd')
    hc.to_csv(out / 'panel_d_hek_clones.tsv', sep='\t', index=False)

    fig.text(0.985, 0.003, stamp, ha='right', va='bottom', fontsize=5, color=MUTED)
    meta = {'Title': 'LDLR Repeat 3 Sp1 site: clinical variants, chromatin and activity', 'Subject': stamp}
    fig.savefig(out / 'fig_repeat3.pdf', metadata=meta)
    fig.savefig(out / 'fig_repeat3.png', dpi=300)


if __name__ == '__main__':
    main()

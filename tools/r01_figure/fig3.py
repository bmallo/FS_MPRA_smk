"""Grant Figure 3 (draft 5): LDLR promoter truth-set validation across readouts (fig3_rules.json: layout_draft5).
One Illustrator-editable SVG (live Arial text) plus PDF and PNG renders.

a  HepG2 MPRA expression change per SNV along the promoter (HGVS c.), ClinVar class by marker, five shaded regions
b  which readout flags which ClinVar control: P/LP (8) and B/LB (6) SNVs x MPRA, plasmid Fiber-seq, AlphaGenome,
   PromoterAI and 'any readout'. Flag = |score| in that readout's top 20% across all designed SNVs (same budget for
   every readout); glyph = direction of the readout (higher / lower than WT). Margin: P/LP and B/LB flagged, share of
   gnomAD-observed and of all SNVs flagged.
c  MPRA x plasmid Fiber-seq for every SNV, ClinVar class by marker, dashed lines at the two flag cutoffs, counts of
   SNVs flagged by both, one or neither readout
d  HepG2 vs HEK293T local inaccessibility change (separate cells, separate axes)
Readout scores as fig_A.py: MPRA |log2|; Fiber-seq max(|nucleosome-size z|, max |footprint z| over A7 confident
fp_features); AlphaGenome max(|DNase|/SD, |LDLR RNA|/SD), HepG2; PromoterAI |score| (MANE ENST00000558518.6).
Usage: python tools/r01_figure/fig3.py --work results/r01_figure_20261102 --root results/analysis_v2_20261028
           --output <dir>
"""
import argparse
import json
import os
import subprocess
import sys
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib import patheffects as pe_  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HPRC = Path('/gscratch/stergachislab/bmallo/large_home/projects/HPRC_driver_variants')
RSVG = '/mmfs1/home/bmallo/miniconda3/envs/fiberseq-qc/bin/rsvg-convert'
FONTS_CONF = Path(__file__).with_name('fonts.conf')
sys.path.insert(0, str(HPRC / 'scripts/plasmid_background'))
from svg_fonts import expand_fonts  # noqa: E402
for f in ('Arimo-Regular.ttf', 'Arimo-Bold.ttf'):
    font_manager.fontManager.addfont(str(HPRC / 'scripts/fonts' / f))
plt.rcParams.update({'font.family': 'Arimo', 'svg.fonttype': 'none', 'font.size': 6.5, 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7, 'xtick.major.size': 2.5,
                     'ytick.major.size': 2.5, 'xtick.labelsize': 6, 'ytick.labelsize': 6})

INK, INK2, MUTED, GRID = '#2E2E2E', '#666666', '#9a9893', '#E6E4DF'
ORANGE, BLUE, AQUA, FLAG = '#eb6834', '#2a78d6', '#1baf7a', '#3B3B3B'
C0, TSS_C = 285, -86                      # c. = oligo - 285; NM_000527.5 TSS at c.-86
# regions (fig3_rules.json regions_draft4): the three motif sites (P < 1e-5), the CRE half-site and the F04 footprint;
# (label in a, short name in b, oligo lo, oligo hi). Shaded in a, named under the panel-b columns.
REGIONS = [('F04 footprint\n(no motif)', 'F04', 66, 82), ('SP-like\n(ZFP3 match)', 'SP-like', 89, 100),
           ('SRE-1\n(SREBP)', 'SRE-1', 127, 136), ('GC box\n(SP1)', 'GC box', 143, 151),
           ('CRE\n(CREB/ATF)', 'CRE', 182, 192)]
DATA_DEFINED = {'F04'}                     # dashed box in a (footprint from the data, not a motif)
CLS = {'none': ('Not in ClinVar', dict(s=5, color='#C9C6BF', lw=0, zorder=2)),
       'vus': ('VUS', dict(s=11, color=BLUE, lw=0, zorder=3)),
       'conf': ('Conflicting', dict(s=16, marker='D', facecolor='white', edgecolor=ORANGE, lw=0.9, zorder=4)),
       'blb': ('Benign / likely benign', dict(s=15, marker='o', facecolor='white', edgecolor=AQUA, lw=1.0, zorder=4)),
       'plp': ('Pathogenic / likely pathogenic', dict(s=18, marker='D', color=ORANGE, edgecolor='white', lw=0.4,
                                                      zorder=5))}
READOUTS = ['MPRA (HepG2)', 'Plasmid Fiber-seq (HepG2)', 'AlphaGenome (HepG2)', 'PromoterAI']
TOP = 0.8                                  # flag = |score| >= the readout's 0.8 quantile across all designed SNVs
W_IN, H_IN = 7.2, 7.3


def readout_scores(w, root):
    """Per SNV and readout: |score| (fig_A.py definitions), signed value of the component that set it, and flag."""
    t = pd.read_csv(w / 'evidence_v3/evidence.tsv', sep='\t')
    fe = pd.read_parquet(root / 'fp_features/effects.parquet').query('mask == "A7" and calls == "confident"')
    fe = fe[fe.p.notna()].assign(zs=lambda d: d.D / d.se)
    fp = fe.loc[fe.zs.abs().groupby(fe.snv).idxmax(), ['snv', 'zs']].set_index('snv').zs    # signed, largest |z|
    fpz = t.snv.map(fp)
    nuc_dom = fpz.isna() | (t.nuc_z.abs() >= fpz.abs())
    ag_d, ag_r = t.ag_dnase_hepg2 / t.ag_dnase_hepg2.std(), t.ag_rna_seq_hepg2 / t.ag_rna_seq_hepg2.std()
    signed = {'MPRA (HepG2)': t.act,
              'Plasmid Fiber-seq (HepG2)': t.nuc_z.where(nuc_dom, fpz),       # + = more protection than WT
              'AlphaGenome (HepG2)': ag_d.where(ag_d.abs() >= ag_r.abs(), ag_r),
              'PromoterAI': t.pai_score}
    S = t[['snv', 'oligo', 'c_dot', 'act', 'clinvar', 'gnomad_ac']].copy()
    for k in READOUTS:
        score = signed[k].abs()
        S[f'{k}|signed'] = signed[k]
        S[f'{k}|score'] = score
        S[f'{k}|flag'] = score >= score.quantile(TOP)
        S[f'{k}|cutoff'] = score.quantile(TOP)
    S['Any readout|flag'] = S[[f'{k}|flag' for k in READOUTS]].any(axis=1)
    S['gnomad'] = S.gnomad_ac.fillna(0) > 0
    return S


def region_of(oligo):
    return next((short for _, short, o0, o1 in REGIONS if o0 <= oligo <= o1), '–')


def cname(snv, c):
    return f'c.{c}{snv.lstrip("0123456789")}'


def scatter_cls(ax, x, y, cls, keys=('none', 'vus', 'conf', 'blb', 'plp')):
    for k in keys:
        m = cls == k
        if m.any():
            ax.scatter(x[m], y[m], **CLS[k][1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--work', required=True)
    ap.add_argument('--root', required=True, help='results/analysis_v2_20261028 (fp_features for the Fiber-seq readout)')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    w, out = Path(a.work), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    pa = pd.read_csv(w / 'fig_A_v9/panel_a.tsv', sep='\t')
    pe = pd.read_csv(w / 'fig_A_v9/panel_d_hek.tsv', sep='\t')
    corr = pd.read_csv(w / 'correlations_v2/pairs.tsv', sep='\t')
    re_ = corr[(corr.x == 'local_inacc_pts') & (corr.y == 'hek_local_pts') & (corr.subset == 'all')].iloc[0]
    S = readout_scores(w, Path(a.root)).merge(pa[['snv', 'cls']], on='snv', how='left')
    assert len(pa) == 633 and len(S) == 633 and len(pe) == int(re_.n) and S.cls.notna().all()
    ctrl = pd.concat([S[S.cls == 'plp'].sort_values('oligo'), S[S.cls == 'blb'].sort_values('oligo')])
    ctrl['region'] = ctrl.oligo.map(region_of)
    ctrl['label'] = [cname(s, c) for s, c in zip(ctrl.snv, ctrl.c_dot)]
    long = []
    for _, r in ctrl.iterrows():
        for k in READOUTS + ['Any readout']:
            long.append(dict(snv=r.snv, label=r.label, clinvar_class=r.cls, region=r.region, readout=k,
                             signed=r.get(f'{k}|signed', np.nan), score=r.get(f'{k}|score', np.nan),
                             flag=bool(r[f'{k}|flag'])))
    pd.DataFrame(long).to_csv(out / 'panel_b_matrix.tsv', sep='\t', index=False, float_format='%.4g')
    rows = READOUTS + ['Any readout']
    margin = {k: dict(plp=int(ctrl[ctrl.cls == 'plp'][f'{k}|flag'].sum()), blb=int(ctrl[ctrl.cls == 'blb'][f'{k}|flag'].sum()),
                      gnomad=float(S[S.gnomad][f'{k}|flag'].mean()), all=float(S[f'{k}|flag'].mean())) for k in rows}
    n_plp, n_blb, n_gn = int((ctrl.cls == 'plp').sum()), int((ctrl.cls == 'blb').sum()), int(S.gnomad.sum())
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()

    fig = plt.figure(figsize=(W_IN, H_IN))
    Wpt, Hpt = W_IN * 72, H_IN * 72

    def axpt(x0, y0, x1, y1):
        return fig.add_axes([x0 / Wpt, 1 - y1 / Hpt, (x1 - x0) / Wpt, (y1 - y0) / Hpt])

    def lab(ch, x, y):
        fig.text(x / Wpt, 1 - y / Hpt, ch, fontsize=10, fontweight='bold', va='top', ha='left')

    ov = fig.add_axes([0, 0, 1, 1], facecolor='none')            # point-coordinate overlay (y down) for panel b
    ov.set_xlim(0, Wpt)
    ov.set_ylim(Hpt, 0)
    ov.axis('off')

    # ---- a: MPRA along the promoter
    lab('a', 4, 6)
    xlo, xhi = pa.c_dot.min() - 2, pa.c_dot.max() + 2
    axt = axpt(50, 12, 508, 34)
    axt.set_xlim(xlo, xhi)
    axt.set_ylim(0, 1)
    axt.axis('off')
    ax = axpt(50, 34, 508, 132)
    for name, short, o0, o1 in REGIONS:
        c0, c1 = o0 - C0 - 0.5, o1 - C0 + 0.5
        ax.axvspan(c0, c1, color='#ECE9F7', lw=0, zorder=0)
        axt.add_patch(Rectangle((c0, 0.08), c1 - c0, 0.32, fc='white', ec='#6A4FC9', lw=0.8,
                                ls=(0, (2, 1.2)) if short in DATA_DEFINED else '-'))
        axt.text((c0 + c1) / 2, 0.48, name, ha='center', va='bottom', fontsize=5.8, color=INK, linespacing=1.05)
    axt.plot([TSS_C, TSS_C, TSS_C + 9], [0.08, 0.45, 0.45], color=INK, lw=0.9)
    axt.annotate('', xy=(TSS_C + 13, 0.45), xytext=(TSS_C + 8, 0.45),
                 arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.9, mutation_scale=6, shrinkA=0, shrinkB=0))
    axt.text(TSS_C + 15, 0.45, 'TSS', ha='left', va='center', fontsize=5.8, color=INK, fontweight='bold')
    ax.axvline(TSS_C, color=INK2, lw=0.6, ls=(0, (2, 2)), zorder=1)
    ax.axhline(0, color=INK2, lw=0.5, zorder=1)
    scatter_cls(ax, pa.c_dot.to_numpy(), pa.act.to_numpy(), pa.cls.to_numpy())
    ax.set_xlim(xlo, xhi)
    ax.set_xticks(range(-250, -54, 25), [f'c.{v}' for v in range(-250, -54, 25)])
    ax.set_xlabel('Position relative to the start codon (HGVS c., NM_000527.5)', fontsize=6.5)
    ax.set_ylabel('Expression change\n(log2, HepG2 MPRA)', fontsize=6.5)
    ax.spines[['top', 'right']].set_visible(False)
    n = pa.cls.value_counts()
    hs = [Line2D([], [], ls='none', marker=CLS[k][1].get('marker', 'o'), ms=4 if k != 'none' else 3,
                 mfc=CLS[k][1].get('facecolor', CLS[k][1].get('color')),
                 mec=CLS[k][1].get('edgecolor', CLS[k][1].get('color')), mew=0.9) for k in CLS]
    fig.legend(hs, [f'{CLS[k][0]} ({n.get(k, 0)})' for k in CLS], loc='upper center', ncol=5, frameon=False,
               fontsize=6, bbox_to_anchor=(279 / Wpt, 1 - 154 / Hpt), handletextpad=0.3, columnspacing=1.4)

    # ---- b: control-variant detection matrix (drawn in points on the overlay)
    lab('b', 4, 174)
    MX = (424, 450, 480, 512)                                     # right edges of the margin columns (pt)
    cw, rh, gap, x0, ytop = 20.0, 16.0, 10.0, 112.0, 252.0
    xs = [x0 + i * cw + (gap if i >= n_plp else 0) for i in range(len(ctrl))]
    for (lo, hi, txt, col) in ((0, n_plp, f'Pathogenic / likely pathogenic ({n_plp})', ORANGE),
                               (n_plp, len(ctrl), f'Benign / likely benign ({n_blb})', AQUA)):
        xa, xb = xs[lo] + 1, xs[hi - 1] + cw - 1
        ov.plot([xa, xb], [190, 190], color=col, lw=1.4, solid_capstyle='butt')
        ov.text((xa + xb) / 2, 187, txt, ha='center', va='bottom', fontsize=6.2, fontweight='bold', color=INK)
    for xi, (_, r) in zip(xs, ctrl.iterrows()):
        ov.text(xi + cw / 2, 236, r.label, rotation=90, ha='center', va='bottom', fontsize=5.6, color=INK)
        ov.text(xi + cw / 2, 245, r.region, ha='center', va='center', fontsize=5.0, color=INK2)
    for j, k in enumerate(rows):
        yc = ytop + j * rh + rh / 2 + (4 if k == 'Any readout' else 0)
        ov.text(x0 - 6, yc, k, ha='right', va='center', fontsize=6.2, color=INK,
                fontweight='bold' if k == 'Any readout' else 'normal')
        for xi, (_, r) in zip(xs, ctrl.iterrows()):
            fl = bool(r[f'{k}|flag'])
            s = 11.0
            ov.add_patch(Rectangle((xi + (cw - s) / 2, yc - s / 2), s, s, fc=FLAG if fl else 'white',
                                   ec=FLAG if fl else '#BDBAB3', lw=0.6))
            if fl and k != 'Any readout':
                ov.text(xi + cw / 2, yc + 0.3, '▲' if r[f'{k}|signed'] > 0 else '▼', ha='center', va='center',
                        fontsize=5.0, color='white')
        mg = margin[k]
        for xm, txt in zip(MX, (f'{mg["plp"]}/{n_plp}', f'{mg["blb"]}/{n_blb}',
                                                   f'{100 * mg["gnomad"]:.0f}%', f'{100 * mg["all"]:.0f}%')):
            ov.text(xm, yc, txt, ha='right', va='center', fontsize=6.0, color=INK,
                    fontweight='bold' if k == 'Any readout' else 'normal')
    ov.plot([x0 - 80, MX[-1]], [ytop + 4 * rh + 2, ytop + 4 * rh + 2], color=GRID, lw=0.6)
    for xm, txt in zip(MX, ('P/LP', 'B/LB', 'gnomAD', 'All SNVs')):
        ov.text(xm, 245, txt, ha='right', va='center', fontsize=5.3, color=INK2, fontweight='bold')
    ov.text((MX[0] + MX[-1]) / 2 - 8, 232, 'Flagged', ha='center', va='bottom', fontsize=5.6, color=INK2,
            fontweight='bold')
    ov.text(x0, ytop + 5 * rh + 12, f'■ flagged = |effect| in that readout\'s top 20% of all {len(S)} SNVs '
            f'(▲ / ▼ readout higher / lower than WT);  □ not flagged', ha='left', va='center', fontsize=5.4,
            color=INK2)
    ov.text(x0, ytop + 5 * rh + 20, f'gnomAD = the {n_gn} designed SNVs observed in gnomAD (all rare, the only '
            f'proxy-benign class available here)', ha='left', va='center', fontsize=5.4, color=INK2)

    # ---- c: MPRA x Fiber-seq for every SNV
    lab('c', 4, 356)
    axc = axpt(44, 366, 300, 486)
    xm_, ym_ = 'MPRA (HepG2)', 'Plasmid Fiber-seq (HepG2)'
    tx, ty = float(S[f'{xm_}|cutoff'].iloc[0]), float(S[f'{ym_}|cutoff'].iloc[0])
    X, Y = S[f'{xm_}|signed'].to_numpy(), S[f'{ym_}|signed'].to_numpy()
    for v in (-tx, tx):
        axc.axvline(v, color=MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)
    for v in (-ty, ty):
        axc.axhline(v, color=MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)
    axc.add_patch(Rectangle((-tx, -ty), 2 * tx, 2 * ty, fc='#F3F1EC', lw=0, zorder=0))
    scatter_cls(axc, X, Y, S.cls.to_numpy())
    fm, ff = S[f'{xm_}|flag'], S[f'{ym_}|flag']
    quad = {'both': int((fm & ff).sum()), 'mpra': int((fm & ~ff).sum()), 'fs': int((~fm & ff).sum()),
            'neither': int((~fm & ~ff).sum())}
    axc.text(0.99, 0.98, f'Flagged by both: {quad["both"]}\nMPRA only: {quad["mpra"]}\n'
             f'Fiber-seq only: {quad["fs"]}\nNeither: {quad["neither"]}', transform=axc.transAxes, ha='right',
             va='top', fontsize=5.6, color=INK, linespacing=1.25)
    # region tag next to each P/LP point (one tag for the GC box pair); offsets in pt, chosen to avoid overlaps
    halo = [pe_.withStroke(linewidth=1.6, foreground='white')]
    tags = {'150C>G': (5, 3, 'left'), '136C>A': (5, 0, 'left'), '133C>T': (5, 0, 'left'),
            '134C>G': (-5, -4, 'right'), '186A>G': (-5, 3, 'right'), '187C>T': (5, 2, 'left'),
            '117A>G': (5, 0, 'left')}
    for _, r in ctrl[ctrl.cls == 'plp'].iterrows():
        if r.snv in tags:
            dx, dy, ha = tags[r.snv]
            axc.annotate(r.region if r.region != '–' else r.label, (r[f'{xm_}|signed'], r[f'{ym_}|signed']),
                         xytext=(dx, dy), textcoords='offset points', fontsize=5.4, color=INK, ha=ha, va='center',
                         zorder=7, path_effects=halo)
    axc.set_xlabel('MPRA: expression change (log2, HepG2)', fontsize=6.3)
    axc.set_ylabel('Plasmid Fiber-seq: protection change\n(signed z, HepG2)', fontsize=6.3)
    axc.spines[['top', 'right']].set_visible(False)

    # ---- d: HepG2 vs HEK293T
    lab('d', 352, 356)
    axe = axpt(392, 380, 508, 486)
    site = pe.oligo.between(143, 151)
    axe.scatter(pe.local_inacc_pts[~site], pe.hek_local_pts[~site], s=6, color=MUTED, lw=0, alpha=0.85)
    axe.scatter(pe.local_inacc_pts[site], pe.hek_local_pts[site], s=12, color=ORANGE, edgecolor='white', lw=0.3)
    axe.axhline(0, color=INK2, lw=0.5)
    axe.axvline(0, color=INK2, lw=0.5)
    axe.set_xlabel('HepG2 (pts)', fontsize=6.3)
    axe.set_ylabel('HEK293T (pts)', fontsize=6.3)
    axe.spines[['top', 'right']].set_visible(False)
    axe.set_title(f'Local inaccessibility change\nSpearman ρ = {re_.rho:.2f} (95% CI {re_.lo:.2f}–{re_.hi:.2f}), '
                  f'{int(re_.n)} SNVs', fontsize=5.8, loc='right', pad=3, linespacing=1.2)
    axe.legend([Line2D([], [], ls='none', marker='o', ms=3, color=MUTED),
                Line2D([], [], ls='none', marker='o', ms=3.5, color=ORANGE)], ['other SNVs', 'GC box (SP1)'],
               frameon=False, fontsize=5.6, loc='lower right', handletextpad=0.2, borderaxespad=0.2)

    S.drop(columns=['clinvar']).to_csv(out / 'panel_c_points.tsv', sep='\t', index=False, float_format='%.4g')
    svg = out / 'fig3.svg'
    fig.savefig(svg, facecolor='white')
    plt.close(fig)
    svg.write_text(expand_fonts(svg.read_text().replace("'Arimo'", 'Arial, Helvetica, sans-serif')
                                .replace('Arimo', 'Arial'))[0])
    env = dict(os.environ, FONTCONFIG_FILE=str(FONTS_CONF))
    subprocess.run([RSVG, '-f', 'pdf', '-o', str(out / 'fig3.pdf'), str(svg)], check=True, env=env)
    subprocess.run([RSVG, '-d', '300', '-p', '300', '-f', 'png', '-o', str(out / 'fig3.png'), str(svg)], check=True,
                   env=env)
    json.dump(dict(snvs=len(S), clinvar=n.to_dict(), flag_quantile=TOP,
                   cutoffs={k: float(S[f'{k}|cutoff'].iloc[0]) for k in READOUTS}, margin=margin,
                   controls=dict(plp=n_plp, blb=n_blb, gnomad_observed=n_gn), panel_c_quadrants=quad,
                   hek=dict(n=int(re_.n), rho=re_.rho, lo=re_.lo, hi=re_.hi), date=str(date.today()), commit=commit),
              open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

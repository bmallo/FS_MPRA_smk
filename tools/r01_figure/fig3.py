"""Grant Figure 3 (draft 7): LDLR promoter truth-set validation across readouts (fig3_rules.json: layout_draft5-7).
One Illustrator-editable SVG (live Arial text) plus PDF and PNG renders.

a  HepG2 MPRA expression change per SNV along the promoter (HGVS c.), ClinVar class by colour and shape, five
   shaded regions
b  which readout shows a large change for which ClinVar control: P/LP (8) and B/LB (6) SNVs x readouts grouped by the
   quantity measured (variant - WT): expression = MPRA (measured), AlphaGenome LDLR RNA and PromoterAI (predicted);
   accessibility = plasmid Fiber-seq (measured), AlphaGenome DNase (predicted); then 'any readout'. Large change =
   |variant - WT| in that row's top 20% across all designed SNVs (same budget for every row; not a pathogenicity
   call); triangle = direction. Margin: share of P/LP, B/LB and rare gnomAD-observed SNVs with a large change.
c  MPRA expression change x plasmid Fiber-seq accessibility change for every SNV, dashed lines at the two cutoffs,
   counts of SNVs with a large change in both, one or neither
Readouts (variant - WT): MPRA log2 (Kircher 2019); Fiber-seq % of molecules without a confident >= 81-bp call in the
element (A7 mask, >= 100 molecules; plain difference, no z); AlphaGenome HepG2 LDLR RNA-seq LFC and DNase 501-bp
DIFF_LOG2_SUM; PromoterAI score (MANE ENST00000558518.6).
Usage: python tools/r01_figure/fig3.py --work results/r01_figure_20261102 --output <dir>
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
PLP_C, BLB_C, CONF_C, VUS_C, FLAG = '#D55E00', '#009E73', '#E69F00', '#56B4E9', '#3B3B3B'   # Okabe-Ito
C0, TSS_C = 285, -86                      # c. = oligo - 285; NM_000527.5 TSS at c.-86
# regions (fig3_rules.json regions_draft4, labels layout_draft8): the three motif sites (P < 1e-5), the CRE half-site
# and the data-defined protected site (F04);
# (label in a, short name in b, oligo lo, oligo hi). Shaded in a, named under the panel-b columns.
REGIONS = [('Protected site\n(no motif match)', 'Protected site', 66, 82), ('SP-like\n(ZFP3 match)', 'SP-like', 89, 100),
           ('SRE-1\n(SREBP)', 'SRE-1', 127, 136), ('GC box\n(SP1)', 'GC box', 143, 151),
           ('CRE\n(CREB/ATF)', 'CRE', 182, 192)]
DATA_DEFINED = {'Protected site'}                     # dashed box in a (footprint from the data, not a motif)
CLS = {'none': ('Not in ClinVar', dict(s=4, color='#CFCCC5', lw=0, zorder=2)),
       'vus': ('VUS', dict(s=8, color=VUS_C, lw=0, zorder=3)),
       'conf': ('Conflicting', dict(s=20, marker='^', color=CONF_C, edgecolor='white', lw=0.4, zorder=4)),
       'blb': ('Benign / likely benign', dict(s=20, marker='o', color=BLB_C, edgecolor='white', lw=0.4, zorder=5)),
       'plp': ('Pathogenic / likely pathogenic', dict(s=24, marker='D', color=PLP_C, edgecolor='white', lw=0.4,
                                                      zorder=6))}
# panel b rows: (key, group, name, detail); every value is variant - WT
ROWS = [('MPRA', 'Expression', 'MPRA', 'measured, HepG2'),
        ('AG RNA', 'Expression', 'AlphaGenome', 'predicted, HepG2 RNA'),
        ('PromoterAI', 'Expression', 'PromoterAI', 'predicted'),
        ('Fiber-seq', 'Accessibility', 'Plasmid Fiber-seq', 'measured, HepG2'),
        ('AG DNase', 'Accessibility', 'AlphaGenome', 'predicted, HepG2 DNase')]
READOUTS = [r[0] for r in ROWS]
TOP = 0.8                                  # large change = |value| >= the row's 0.8 quantile across all designed SNVs
MIN_MOLECULES = 100                        # Fiber-seq readout needs this many variant molecules (post hoc)
W_IN, H_IN = 7.2, 7.4


def readout_scores(w):
    """Per SNV and readout: signed variant - WT value, |value| and large-change flag (top 20%)."""
    t = pd.read_csv(w / 'evidence_v3/evidence.tsv', sep='\t')
    nuc = pd.read_csv(w / 'gnomad_check/nucleosome_all.tsv', sep='\t').query(
        'mask == "A7" and calls == "confident" and endpoint == "nucleosome"').set_index('snv')
    nuc = nuc[nuc.n >= MIN_MOLECULES]
    acc_pts = t.snv.map(-100 * (nuc.occ_variant - nuc.wt_occ))    # % of molecules without nucleosome-sized protection
    signed = {'MPRA': t.act, 'AG RNA': t.ag_rna_seq_hepg2, 'PromoterAI': t.pai_score,
              'Fiber-seq': acc_pts, 'AG DNase': t.ag_dnase_hepg2}
    S = t[['snv', 'oligo', 'c_dot', 'act', 'clinvar', 'gnomad_ac', 'gnomad_af']].copy()
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
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    w, out = Path(a.work), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    pa = pd.read_csv(w / 'fig_A_v9/panel_a.tsv', sep='\t')
    S = readout_scores(w).merge(pa[['snv', 'cls']], on='snv', how='left')
    assert len(pa) == 633 and len(S) == 633 and S.cls.notna().all()
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
    hs = [Line2D([], [], ls='none', marker=CLS[k][1].get('marker', 'o'), ms={'none': 2.5, 'vus': 3}.get(k, 4.5),
                 color=CLS[k][1]['color'], mec='white', mew=0.3) for k in CLS]
    fig.legend(hs, [f'{CLS[k][0]} ({n.get(k, 0)})' for k in CLS], loc='upper center', ncol=5, frameon=False,
               fontsize=6, bbox_to_anchor=(279 / Wpt, 1 - 154 / Hpt), handletextpad=0.3, columnspacing=1.4)

    # ---- b: which readout shows a large change for which ClinVar control (drawn in points on the overlay)
    lab('b', 4, 174)
    MX = (460, 486, 514)                                     # right edges of the margin columns (pt)
    cw, rh, gap, x0, ytop = 19.0, 16.0, 10.0, 150.0, 252.0
    xs = [x0 + i * cw + (gap if i >= n_plp else 0) for i in range(len(ctrl))]
    for (lo, hi, txt, col) in ((0, n_plp, f'Pathogenic / likely pathogenic ({n_plp})', PLP_C),
                               (n_plp, len(ctrl), f'Benign / likely benign ({n_blb})', BLB_C)):
        xa, xb = xs[lo] + 1, xs[hi - 1] + cw - 1
        ov.plot([xa, xb], [190, 190], color=col, lw=1.4, solid_capstyle='butt')
        ov.text((xa + xb) / 2, 187, txt, ha='center', va='bottom', fontsize=6.2, fontweight='bold', color=INK)
    for xi, (_, r) in zip(xs, ctrl.iterrows()):
        ov.text(xi + cw / 2, 236, r.label, rotation=90, ha='center', va='bottom', fontsize=5.6, color=INK)
    i = 0                                                         # one region name per run of consecutive columns
    while i < len(ctrl):
        j = i
        while j + 1 < len(ctrl) and j + 1 != n_plp and ctrl.region.iloc[j + 1] == ctrl.region.iloc[i]:
            j += 1
        xa, xb = xs[i] + 2, xs[j] + cw - 2
        if j > i and ctrl.region.iloc[i] != '–':
            ov.plot([xa, xb], [240, 240], color=MUTED, lw=0.5)
        ov.text((xa + xb) / 2, 245.5, ctrl.region.iloc[i], ha='center', va='center', fontsize=4.8, color=INK2)
        i = j + 1
    for xm, txt in zip(MX, (f'P/LP\n({n_plp})', f'B/LB\n({n_blb})', f'Rare in\ngnomAD\n({n_gn})')):
        ov.text(xm, 248.5, txt, ha='right', va='bottom', fontsize=5.0, color=INK2, fontweight='bold', linespacing=1.05)
    ov.text((MX[0] + MX[-1]) / 2 - 6, 229, 'Share with a\nlarge change', ha='center', va='bottom', fontsize=5.6,
            color=INK2, fontweight='bold', linespacing=1.1)

    def margin_row(k, yc, bold=False):
        mg = margin[k]
        for xm, txt in zip(MX, (f'{mg["plp"]}/{n_plp}', f'{mg["blb"]}/{n_blb}', f'{100 * mg["gnomad"]:.0f}%')):
            ov.text(xm, yc, txt, ha='right', va='center', fontsize=6.0, fontweight='bold' if bold else 'normal',
                    color=INK)

    y = ytop
    for grp in ('Expression', 'Accessibility'):
        ov.text(8, y + 5, f'{grp} change (variant − WT)', ha='left', va='center', fontsize=6.2, fontweight='bold',
                color=INK)
        ov.plot([8, MX[-1]], [y + 10, y + 10], color=GRID, lw=0.5)
        y += 10
        for k, g, name, detail in ROWS:
            if g != grp:
                continue
            yc = y + rh / 2
            ov.text(x0 - 8, yc - 2.6, name, ha='right', va='center', fontsize=6.2, color=INK)
            ov.text(x0 - 8, yc + 4.2, detail, ha='right', va='center', fontsize=4.9, color=INK2)
            for xi, (_, r) in zip(xs, ctrl.iterrows()):
                if r[f'{k}|flag']:
                    ov.scatter(xi + cw / 2, yc, marker='^' if r[f'{k}|signed'] > 0 else 'v', s=30, color=FLAG, lw=0)
                else:
                    ov.scatter(xi + cw / 2, yc, s=3, color='#BDBAB3', lw=0)
            margin_row(k, yc)
            y += rh
        y += 3
    ov.plot([8, MX[-1]], [y, y], color=INK2, lw=0.6)
    yc = y + 4 + rh / 2
    ov.text(x0 - 8, yc, 'Any readout', ha='right', va='center', fontsize=6.2, color=INK, fontweight='bold')
    for xi, (_, r) in zip(xs, ctrl.iterrows()):
        if r['Any readout|flag']:
            ov.scatter(xi + cw / 2, yc, marker='s', s=34, color=FLAG, lw=0)
        else:
            ov.scatter(xi + cw / 2, yc, s=3, color='#BDBAB3', lw=0)
    margin_row('Any readout', yc, bold=True)
    yk = yc + rh / 2 + 9
    ov.text(30, yk, f'▲ / ▼ large increase / decrease (|variant − WT| in that row\'s top 20% of {len(S)} SNVs; any readout: '
            f'{100 * margin["Any readout"]["all"]:.0f}% of SNVs);  · smaller change', ha='left', va='center', fontsize=5.4, color=INK2)
    ov.text(30, yk + 8, f'Fiber-seq accessibility = % of molecules without a nucleosome-sized (≥ 81 bp) protected '
            f'stretch;  rare in gnomAD: allele frequency ≤ {100 * S.gnomad_af.max():.2g}%', ha='left', va='center',
            fontsize=5.4, color=INK2)

    # ---- c: MPRA x Fiber-seq for every SNV
    y0c = yk + 20
    lab('c', 4, y0c)
    axc = axpt(50, y0c + 8, 508, y0c + 96)
    xm_, ym_ = 'MPRA', 'Fiber-seq'
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
    axc.text(0.995, 0.98, f'Large change in both: {quad["both"]}\nMPRA only: {quad["mpra"]}\n'
             f'Fiber-seq only: {quad["fs"]}\nNeither: {quad["neither"]}', transform=axc.transAxes, ha='right',
             va='top', fontsize=5.6, color=INK, linespacing=1.25)
    # region tag next to each P/LP point (one tag each for the GC box and CRE pairs); offsets in pt, chosen to avoid overlaps
    halo = [pe_.withStroke(linewidth=1.6, foreground='white')]
    tags = {'150C>G': (6, 0, 'left'), '136C>A': (5, 3, 'left'), '133C>T': (5, 2, 'left'),
            '134C>G': (-5, 4, 'right'), '186A>G': (-5, -5, 'right'), '117A>G': (5, 2, 'left')}
    for _, r in ctrl[ctrl.cls == 'plp'].iterrows():
        if r.snv in tags:
            dx, dy, ha = tags[r.snv]
            axc.annotate(r.region if r.region != '–' else r.label, (r[f'{xm_}|signed'], r[f'{ym_}|signed']),
                         xytext=(dx, dy), textcoords='offset points', fontsize=5.4, color=INK, ha=ha, va='center',
                         zorder=7, path_effects=halo)
    axc.set_xlabel('MPRA: expression change (log2, HepG2)', fontsize=6.3)
    axc.set_ylabel('Plasmid Fiber-seq: accessibility\nchange (% of molecules)', fontsize=6.3)
    axc.set_yticks(range(-10, 11, 5))
    axc.spines[['top', 'right']].set_visible(False)

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
                   date=str(date.today()), commit=commit),
              open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

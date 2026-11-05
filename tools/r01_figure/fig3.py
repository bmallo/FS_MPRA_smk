"""Grant Figure 3 (draft 14): LDLR promoter truth-set validation across readouts and truth-set potential of the new
promoters (fig3_rules.json: layout_draft5-14). One Illustrator-editable SVG (live Arial text) plus PDF and PNG renders.

a  HepG2 MPRA expression change per SNV along the promoter (HGVS c.), ClinVar class by colour and shape, five
   shaded regions
b  evidence per LDLR ClinVar control: P/LP (8) and B/LB (6) SNVs x rows grouped by evidence type: population
   (gnomAD v4.1 allele frequency, genomes + exomes; GTEx V10: eQTL / tested / in the variant table but not tested /
   absent), predicted (AlphaGenome HepG2 RNA and DNase, PromoterAI) and measured (MPRA, plasmid Fiber-seq). HEK293T
   plasmid Fiber-seq: same accessibility readout from tools/r01_figure/hek_fiberseq.py (SNVs with >= 100 molecules,
   any number of clones; top 20% within HEK293T), its own row.
   Large change = |variant - WT| in that row's top 20% across all designed SNVs (same
   budget for every row; not a pathogenicity call); triangle = direction. Margin: P/LP and B/LB with a large change
   (population rows: in gnomAD / eQTL).
c  truth-set potential per element (LDLR, HBB, GATA2 +9.5 enhancer, SLC39A4): ClinVar P/LP, B/LB, VUS/conflicting and the
   candidate-benign funnel (possible SNVs -> gnomAD -> common -> not an eQTL -> predicted no effect) from
   tools/promoter_variant_census/truth_sets.py (truth_set_rules.json)
d  MPRA expression change x plasmid Fiber-seq accessibility change for every SNV, dashed lines at the two cutoffs,
   counts of SNVs with a large change in both, one or neither
Readouts (variant - WT): MPRA log2 (Kircher 2019); Fiber-seq % of molecules without a confident >= 81-bp call in the
element (A7 mask, >= 100 molecules; plain difference, no z); AlphaGenome HepG2 LDLR RNA-seq LFC and DNase 501-bp
DIFF_LOG2_SUM; PromoterAI score (MANE ENST00000558518.6).
Usage: python tools/r01_figure/fig3.py --work results/r01_figure_20261102
           --census results/promoter_variant_census_20261104 --hek results/r01_figure_20261102/hek_fiberseq_20261104
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
# panel b functional rows: (key, group, name, detail); every value is variant - WT
ROWS = [('AG RNA', 'Predicted', 'AlphaGenome', 'expression, HepG2 RNA'),
        ('AG DNase', 'Predicted', 'AlphaGenome', 'accessibility, HepG2 DNase'),
        ('PromoterAI', 'Predicted', 'PromoterAI', 'expression'),
        ('MPRA', 'Measured', 'MPRA', 'expression, HepG2'),
        ('Fiber-seq', 'Measured', 'Plasmid Fiber-seq', 'accessibility, HepG2')]
GROUPS = {'Population': 'Population evidence', 'Predicted': 'Predicted change (variant − WT)',
          'Measured': 'Measured change (variant − WT)'}
READOUTS = [r[0] for r in ROWS]
TOP = 0.8                                  # large change = |value| >= the row's 0.8 quantile across all designed SNVs
MIN_MOLECULES = 100                        # Fiber-seq readout needs this many variant molecules (post hoc)
OLIGO_OFFSET = 11089264                    # hg38 = 11,089,264 + oligo
# panel c: (truth_sets.tsv element, label, note)
ELEMENTS = [('LDLR', 'LDLR', 'promoter; b = 211-bp measured element'), ('HBB', 'HBB', 'promoter'),
            ('GATA2_E95', 'GATA2', '+9.5 intronic enhancer, 2 windows, 503 bp*'),
            ('SLC39A4', 'SLC39A4', 'promoter, 2 windows, 542 bp')]
W_IN, H_IN = 7.2, 10.51


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


def population(census, snvs):
    """Per control SNV: gnomAD v4.1 genome status (common = AF > 1% in any group with AN >= 2000 / rare / absent), joint
    genome + exome allele frequency, and GTEx V10 status (eQTL for any gene / tested, MAF >= 1% in GTEx / in the
    variant table, not tested / absent), from the census variants.tsv."""
    v = pd.read_csv(census / 'variants.tsv', sep='\t').query('window == "LDLR"')
    v = v.set_index(v.pos.astype(str) + v.ref + v.alt)
    rows = []
    for s in snvs:
        k = f'{OLIGO_OFFSET + int(s[:-3])}{s[-3]}{s[-1]}'
        r = v.loc[k] if k in v.index else None
        ing = r is not None and bool(r.in_gnomad)
        gn = 'common' if ing and bool(r.af_gt_1pct_any_group) else 'rare' if ing else 'absent'
        sig = r is not None and isinstance(r.gtex_v10_sig_any_genes, str) and r.gtex_v10_sig_any_genes != ''
        gt = ('eqtl' if sig else 'tested' if r is not None and bool(r.gtex_v10_maf01) else
              'untested' if r is not None and bool(r.gtex_v10_present) else 'absent')
        ac, an = (r.g_ac + r.e_ac, r.g_an + r.e_an) if r is not None else (0, 0)
        rows.append(dict(snv=s, pop_gnomad=gn, pop_gtex=gt, joint_ac=int(ac), joint_an=int(an),
                         joint_af=ac / an if ac else 0.0, g_af=r.g_af if ing else 0.0,
                         g_grpmax_af=r.g_grpmax_af if ing else 0.0,
                         g_grpmax_group=r.g_grpmax_group if ing else ''))
    return pd.DataFrame(rows)


def hek_scores(hek):
    """HEK293T run 1 per SNV (hek_fiberseq.py, >= 100 molecules): accessibility change and top-20% flag within HEK293T."""
    h = pd.read_csv(hek / 'nucleosome_hek.tsv', sep='\t')
    score = h.acc_change_pts.abs()
    return pd.DataFrame(dict(snv=h.snv, hek_signed=h.acc_change_pts, hek_flag=score >= score.quantile(TOP),
                             hek_cutoff=score.quantile(TOP), hek_n=h.n, hek_clones=h.clones,
                             hek_top_clone_share=h.top_clone_share, hek_D_pts=-100 * h.D))


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
    ap.add_argument('--census', required=True)
    ap.add_argument('--hek', required=True, help='hek_fiberseq.py output dir')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    w, out, census = Path(a.work), Path(a.output), Path(a.census)
    out.mkdir(parents=True, exist_ok=False)
    pa = pd.read_csv(w / 'fig_A_v9/panel_a.tsv', sep='\t')
    S = readout_scores(w).merge(pa[['snv', 'cls']], on='snv', how='left')
    assert len(pa) == 633 and len(S) == 633 and S.cls.notna().all()
    ctrl = pd.concat([S[S.cls == 'plp'].sort_values('oligo'), S[S.cls == 'blb'].sort_values('oligo')])
    ctrl['region'] = ctrl.oligo.map(region_of)
    ctrl['label'] = [cname(s, c) for s, c in zip(ctrl.snv, ctrl.c_dot)]
    ctrl = ctrl.merge(population(census, ctrl.snv), on='snv', how='left', validate='1:1')
    H = hek_scores(Path(a.hek))
    ctrl = ctrl.merge(H, on='snv', how='left', validate='1:1')
    ctrl['hek_covered'] = ctrl.hek_signed.notna()
    long = []
    for _, r in ctrl.iterrows():
        for k in READOUTS + ['Any readout']:
            long.append(dict(snv=r.snv, label=r.label, clinvar_class=r.cls, region=r.region, readout=k,
                             signed=r.get(f'{k}|signed', np.nan), score=r.get(f'{k}|score', np.nan),
                             flag=bool(r[f'{k}|flag'])))
    pd.DataFrame(long).to_csv(out / 'panel_b_matrix.tsv', sep='\t', index=False, float_format='%.4g')
    ctrl[['snv', 'label', 'cls', 'joint_ac', 'joint_an', 'joint_af', 'pop_gnomad', 'g_af', 'g_grpmax_af', 'g_grpmax_group',
          'pop_gtex']].to_csv(
        out / 'panel_b_population.tsv', sep='\t', index=False, float_format='%.4g')
    ctrl[['snv', 'label', 'cls', 'hek_covered', 'hek_n', 'hek_clones', 'hek_top_clone_share', 'hek_signed', 'hek_D_pts',
          'hek_flag']].to_csv(out / 'panel_b_hek293t.tsv', sep='\t', index=False, float_format='%.4g')
    hg = H.merge(S[['snv', 'gnomad']], on='snv', how='left')
    hek_margin = dict(plp=int(ctrl[ctrl.cls == 'plp'].hek_flag.fillna(False).sum()),
                      plp_n=int(ctrl[ctrl.cls == 'plp'].hek_covered.sum()),
                      blb=int(ctrl[ctrl.cls == 'blb'].hek_flag.fillna(False).sum()),
                      blb_n=int(ctrl[ctrl.cls == 'blb'].hek_covered.sum()),
                      gnomad=int(hg[hg.gnomad.fillna(False)].hek_flag.sum()), gnomad_n=int(hg.gnomad.fillna(False).sum()),
                      snvs=len(H), cutoff=float(H.hek_cutoff.iloc[0]))
    ts = pd.read_csv(census / 'truth_sets/truth_sets.tsv', sep='\t').set_index('element')
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
    MX = (486, 514)                                          # right edges of the margin columns (pt)
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
    for xm, txt in zip(MX, (f'P/LP\n({n_plp})', f'B/LB\n({n_blb})')):
        ov.text(xm, 248.5, txt, ha='right', va='bottom', fontsize=5.0, color=INK2, fontweight='bold', linespacing=1.05)
    ov.text((MX[0] + MX[-1]) / 2 - 8, 229, 'Large\nchange', ha='center', va='bottom', fontsize=5.6,
            color=INK2, fontweight='bold', linespacing=1.1)

    def margin_row(k, yc, bold=False):
        mg = margin[k]
        for xm, txt in zip(MX, (f'{mg["plp"]}/{n_plp}', f'{mg["blb"]}/{n_blb}')):
            ov.text(xm, yc, txt, ha='right', va='center', fontsize=6.0, fontweight='bold' if bold else 'normal',
                    color=INK)

    def group_head(grp, y):
        ov.text(8, y + 5, GROUPS[grp], ha='left', va='center', fontsize=6.2, fontweight='bold', color=INK)
        ov.plot([8, MX[-1]], [y + 10, y + 10], color=GRID, lw=0.5)
        return y + 10

    def pop_margin(yc, n_p, n_b):
        for xm, txt in zip(MX, (f'{n_p}/{n_plp}', f'{n_b}/{n_blb}')):
            ov.text(xm, yc, txt, ha='right', va='center', fontsize=6.0, color=INK)

    y = group_head('Population', ytop)
    yc = y + rh / 2                                              # gnomAD: allele frequency, genomes + exomes
    ov.text(x0 - 8, yc - 2.6, 'gnomAD v4.1', ha='right', va='center', fontsize=6.2, color=INK)
    ov.text(x0 - 8, yc + 4.2, 'allele frequency', ha='right', va='center', fontsize=4.9, color=INK2)
    for xi, (_, r) in zip(xs, ctrl.iterrows()):
        if r.joint_ac:
            m_, e_ = f'{r.joint_af:.1e}'.split('e')
            ov.text(xi + cw / 2, yc, f'{m_}e{int(e_)}', ha='center', va='center', fontsize=4.7, color=INK)
        else:
            ov.scatter(xi + cw / 2, yc, s=3, color='#BDBAB3', lw=0)
    pop_margin(yc, int(((ctrl.cls == 'plp') & (ctrl.joint_ac > 0)).sum()),
               int(((ctrl.cls == 'blb') & (ctrl.joint_ac > 0)).sum()))
    y += rh
    for name, detail, col, hit in [('GTEx V10', 'eQTL for any gene', 'pop_gtex', 'eqtl')]:
        yc = y + rh / 2
        ov.text(x0 - 8, yc - 2.6, name, ha='right', va='center', fontsize=6.2, color=INK)
        ov.text(x0 - 8, yc + 4.2, detail, ha='right', va='center', fontsize=4.9, color=INK2)
        for xi, (_, r) in zip(xs, ctrl.iterrows()):
            st = r[col]
            if st == hit:
                ov.scatter(xi + cw / 2, yc, marker='o', s=26, color=FLAG, lw=0)
            elif st in ('rare', 'tested'):
                ov.scatter(xi + cw / 2, yc, marker='o', s=16, facecolor='white', edgecolor=FLAG, lw=0.6)
            elif st == 'untested':
                ov.text(xi + cw / 2, yc, 'n.t.', ha='center', va='center', fontsize=4.6, color=MUTED)
            else:
                ov.scatter(xi + cw / 2, yc, s=3, color='#BDBAB3', lw=0)
        pop_margin(yc, int(((ctrl.cls == 'plp') & (ctrl[col] == hit)).sum()),
                   int(((ctrl.cls == 'blb') & (ctrl[col] == hit)).sum()))
        y += rh
    y += 3
    for grp in ('Predicted', 'Measured'):
        y = group_head(grp, y)
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
        if grp == 'Measured':                                   # HEK293T: own row, top 20% within HEK293T
            yc = y + rh / 2
            ov.text(x0 - 8, yc - 2.6, 'Plasmid Fiber-seq', ha='right', va='center', fontsize=6.2, color=INK)
            ov.text(x0 - 8, yc + 4.2, 'accessibility, HEK293T', ha='right', va='center', fontsize=4.9, color=INK2)
            for xi, (_, r) in zip(xs, ctrl.iterrows()):
                if not r.hek_covered:
                    ov.text(xi + cw / 2, yc, '–', ha='center', va='center', fontsize=5.6, color=MUTED)
                elif r.hek_flag:
                    ov.scatter(xi + cw / 2, yc, marker='^' if r.hek_signed > 0 else 'v', s=30, color=FLAG, lw=0)
                else:
                    ov.scatter(xi + cw / 2, yc, s=3, color='#BDBAB3', lw=0)
            hm = hek_margin
            for xm, txt in zip(MX, (f'{hm["plp"]}/{hm["plp_n"]}', f'{hm["blb"]}/{hm["blb_n"]}')):
                ov.text(xm, yc, txt, ha='right', va='center', fontsize=6.0, color=INK)
            y += rh
        y += 3
    ov.plot([8, MX[-1]], [y, y], color=GRID, lw=0.5)
    yk = y + 9
    ov.text(30, yk, f'▲ / ▼ large increase / decrease (|variant − WT| in that row\'s top 20% of {len(S)} SNVs);  · smaller '
            f'change;  margins: number with a large change (population rows: in gnomAD / eQTL)', ha='left', va='center',
            fontsize=5.4, color=INK2)
    ov.text(30, yk + 8, 'Fiber-seq accessibility = % of molecules without a nucleosome-sized (≥ 81 bp) protected stretch',
            ha='left', va='center', fontsize=5.4, color=INK2)
    top = ctrl.loc[ctrl.g_grpmax_af.idxmax()]
    ov.text(30, yk + 16, f'gnomAD: allele frequency, genomes + exomes ({top.label}: {100 * top.g_grpmax_af:.2f}% in '
            f'African/African American genomes);  GTEx: n.t. genotyped, not tested (MAF < 1% in donors);  · absent',
            ha='left', va='center', fontsize=5.4, color=INK2)
    ov.text(30, yk + 24, f'HEK293T: same measure, top 20% of its {hek_margin["snvs"]} SNVs with ≥ 100 molecules; – fewer '
            f'molecules; most molecules from 1–2 barcode clones (clone and variant effects not separable)',
            ha='left', va='center', fontsize=5.4, color=INK2)

    # ---- c: truth-set potential per element (census truth_sets.tsv)
    y0t = yk + 40
    lab('c', 4, y0t)
    ov.text(30, y0t + 4, 'Candidate truth-set variants in each promoter or enhancer (300-bp windows)', ha='left', va='top',
            fontsize=6.4, fontweight='bold', color=INK)
    # column right edges (pt) and headers; two groups: ClinVar, candidate benign (each step adds one criterion)
    CX = [176, 212, 254, 304, 346, 390, 440, 500]
    heads = ['P/LP', 'B/LB', 'VUS /\nconflicting', 'Possible\nSNVs', 'In gnomAD\ngenomes', 'Common\n(> 1%)',
             '+ not an\neQTL', '+ predicted\nno effect']
    yh = y0t + 50
    for (lo, hi, txt, col) in ((0, 3, 'ClinVar (SNVs + indels)', INK2), (3, 8, 'Candidate benign SNVs', BLB_C)):
        xa, xb = (CX[lo - 1] + 6 if lo else 128), CX[hi - 1]
        ov.plot([xa, xb], [y0t + 29, y0t + 29], color=col, lw=1.2, solid_capstyle='butt')
        ov.text((xa + xb) / 2, y0t + 26, txt, ha='center', va='bottom', fontsize=6.0, fontweight='bold', color=INK)
    for xm, h in zip(CX, heads):
        ov.text(xm, yh, h, ha='right', va='bottom', fontsize=5.2, color=INK2, fontweight='bold', linespacing=1.05)
    ov.plot([8, CX[-1]], [yh + 3, yh + 3], color=GRID, lw=0.5)
    rt = 13.0
    for i, (el, name, note) in enumerate(ELEMENTS):
        r = ts.loc[el]
        yc = yh + 3 + rt * (i + 0.5)
        ov.text(12, yc, name, ha='left', va='center', fontsize=6.2, color=INK, style='italic')
        if note:
            ov.text(12 + 5.6 * len(name) + 4, yc, note, ha='left', va='center', fontsize=4.9, color=INK2)
        vals = [r.plp_snv + r.plp_indel, r.blb_snv + r.blb_indel, r.vus + r.conflicting, r.possible_snvs,
                r.gnomad_snvs, r.b1_snvs, r.b1_b2_snvs, r.b1_b2_b3_snvs]
        cols = [PLP_C, BLB_C, INK, INK2, INK, INK, INK, INK]
        for j, (xm, val, c) in enumerate(zip(CX, vals, cols)):
            ov.text(xm, yc, f'{int(val):,}', ha='right', va='center', fontsize=6.2, color=c,
                    fontweight='bold' if j in (0, 1, 7) else 'normal')
    yb = yh + 3 + rt * len(ELEMENTS)
    ov.plot([8, CX[-1]], [yb, yb], color=GRID, lw=0.5)
    tot = ts.loc[[e for e, _, _ in ELEMENTS]]
    for k, line in enumerate((
            f'Common = AF > 1% in ≥ 1 ancestry group, gnomAD v4.1 genomes (> 1% overall: '
            f'{int(tot.b1_overall_snvs.sum())} SNVs across the four elements)',
            f'Not an eQTL = not significant for any gene in GTEx V10; {int(tot.b1_untested.sum())} of '
            f'{int(tot.b1_snvs.sum())} common SNVs were not testable (MAF < 1% in GTEx donors)',
            'Predicted no effect = no top-20% change in AlphaGenome RNA or DNase (gene-relevant cell types) and '
            'PromoterAI |score| < 0.1 (* enhancer: AlphaGenome only; PromoterAI is a promoter model)',
            'ClinVar germline classifications, accessed 2026-11-04;  windows cover the FIRE peaks (HBB: no peak; '
            'Kircher 2019 element)')):
        ov.text(30, yb + 8 + 8 * k, line, ha='left', va='center', fontsize=5.4, color=INK2)
    yb += 16

    # ---- d: MPRA x Fiber-seq for every SNV
    y0c = yb + 34
    lab('d', 4, y0c)
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

    S.drop(columns=['clinvar']).to_csv(out / 'panel_d_points.tsv', sep='\t', index=False, float_format='%.4g')
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
                   controls=dict(plp=n_plp, blb=n_blb, gnomad_observed=n_gn), panel_d_quadrants=quad,
                   panel_b_population=ctrl.set_index('label')[['pop_gnomad', 'pop_gtex']].to_dict('index'),
                   panel_b_hek293t=dict(margin=hek_margin, controls=ctrl.set_index('label')[
                       ['hek_covered', 'hek_n', 'hek_clones', 'hek_signed', 'hek_flag']].to_dict('index')),
                   panel_c=ts.loc[[e for e, _, _ in ELEMENTS]].to_dict('index'),
                   date=str(date.today()), commit=commit),
              open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

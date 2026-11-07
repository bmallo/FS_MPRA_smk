"""Standalone Fig 3b (draft 15; fig3_rules.json panel_b_isolated_draft15): which evidence flags each LDLR ClinVar
control across the whole Kircher 2019 element (chr19:11,089,231-11,089,548, c.-318 to c.-1; 954 SNVs).

Columns: ClinVar P/LP and B/LB SNVs in the element, 5'->3' within class. Rows:
  population (yes/no)  gnomAD v4.1 (genomes or exomes), GTEx V10 tested for eQTL (drafts 15-16 also: alt allele in a primate
                       reference genome, Ensembl EPO; removed in draft 17)
  predicted change     AlphaGenome HepG2 RNA (LDLR) and DNase, PromoterAI
  measured change      MPRA (Kircher 2019, HepG2), plasmid Fiber-seq accessibility HepG2 and HEK293T (library SNVs)
Large change = |variant - WT| in that row's top 20% over every SNV of the element the row scores (954; Fiber-seq rows:
the plasmid-library SNVs with >= 100 molecules, as draft 14). Margins: P/LP and B/LB with yes / a large change, over
the controls the row covers.
Inputs: tools/r01_figure/fig3b_data.py outputs (scores.tsv, controls_full.tsv).
Draft 18 (compact_draft18): grant-size layout, ~3.9 in tall, one-line labels, section labels in the bands.
Draft 17 (panel_b_draft17): primate row removed; each row group coloured (heading, symbols, tint band).
Draft 16 (panel_ab_draft16): panel a (MPRA along c.-318..c.-1, ClinVar class by marker, Fig 3 a regions) added
above b; b's explanatory text replaced by one symbol key.
Usage: python tools/r01_figure/fig3b.py --data results/r01_figure_20261102/fig3b_data --output <dir>
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
import pandas as pd  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HPRC = Path('/gscratch/stergachislab/bmallo/large_home/projects/HPRC_driver_variants')
RSVG = '/mmfs1/home/bmallo/miniconda3/envs/fiberseq-qc/bin/rsvg-convert'
FONTS_CONF = Path(__file__).with_name('fonts.conf')
sys.path.insert(0, str(HPRC / 'scripts/plasmid_background'))
from svg_fonts import expand_fonts  # noqa: E402
for f in ('Arimo-Regular.ttf', 'Arimo-Bold.ttf'):
    font_manager.fontManager.addfont(str(HPRC / 'scripts/fonts' / f))
plt.rcParams.update({'font.family': 'Arimo', 'svg.fonttype': 'none'})

INK, INK2, MUTED, GRID, DOT = '#2E2E2E', '#666666', '#9a9893', '#E6E4DF', '#BDBAB3'
PLP_C, BLB_C, FLAG = '#D55E00', '#009E73', '#3B3B3B'
CONF_C, VUS_C = '#E69F00', '#56B4E9'
TSS_C = -86                               # NM_000527.5 TSS at c.-86
# panel a markers (as Fig 3 a); 'other' ClinVar records drawn as not in ClinVar
CLS = {'none': ('Not in ClinVar', dict(s=4, color='#CFCCC5', lw=0, zorder=2)),
       'vus': ('VUS', dict(s=8, color=VUS_C, lw=0, zorder=3)),
       'conflicting': ('Conflicting', dict(s=20, marker='^', color=CONF_C, edgecolor='white', lw=0.4, zorder=4)),
       'blb': ('Benign / likely benign', dict(s=20, marker='o', color=BLB_C, edgecolor='white', lw=0.4, zorder=5)),
       'plp': ('Pathogenic / likely pathogenic', dict(s=24, marker='D', color=PLP_C, edgecolor='white', lw=0.4,
                                                      zorder=6))}
A_LABELS = {k: k for k in ('Protected site', 'SP-like', 'SRE-1', 'GC box', 'CRE')}   # one line (draft 18)
A_ALIGN = {'Protected site': 'right', 'SP-like': 'left', 'SRE-1': 'right', 'GC box': 'left', 'CRE': 'center'}
FS, FS2 = 7.0, 6.0                       # main text, secondary text (row details, region names, key)
TOP = 0.8
C0 = 285
# regions named under the columns (oligo coordinates, as panel a of Fig 3)
REGIONS = [('Protected site', 66, 82), ('SP-like', 89, 100), ('SRE-1', 127, 136), ('GC box', 143, 151),
           ('CRE', 182, 192)]
# (key, group, name, detail, column); population rows are yes/no, the rest top-20% flags
POP = [('gnomad', 'gnomAD v4.1', 'observed'),
       ('gtex_tested', 'GTEx V10', 'eQTL tested')]           # primate row removed in draft 17 (user)
FUNC = [('ag_rna', 'Predicted', 'AlphaGenome RNA', 'HepG2'),
        ('ag_dnase', 'Predicted', 'AlphaGenome DNase', 'HepG2'),
        ('pai', 'Predicted', 'PromoterAI', ''),
        ('mpra', 'Measured', 'MPRA', 'HepG2'),
        ('fs_hepg2', 'Measured', 'Fiber-seq', 'HepG2'),
        ('fs_hek293t', 'Measured', 'Fiber-seq', 'HEK293T')]
# section colours (panel_b_draft17): heading, symbols, and a light tint band behind the rows
SEC = {'Population': '#9A6A12', 'Predicted': '#6A55B5', 'Measured': '#00809A'}
TINT = 0.10                               # band = this fraction of the hue over white
GROUPS = {'Population': ('Population', 'yes / no'), 'Predicted': ('Predicted', 'variant − WT'),
          'Measured': ('Measured', 'variant − WT')}
W_IN = 7.2


def region_of(oligo):
    return next((r for r, o0, o1 in REGIONS if o0 <= oligo <= o1), '–')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--data', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    d, out = Path(a.data), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    S = pd.read_csv(d / 'scores.tsv', sep='\t')
    C = pd.read_csv(d / 'controls_full.tsv', sep='\t')
    assert len(S) == 954 and len(C) == len(C.snv.unique())
    cut = {}
    for k, *_ in FUNC:
        v = S[k].dropna().abs()
        cut[k] = float(v.quantile(TOP))
        S[f'{k}|flag'] = S[k].abs() >= cut[k]
    C = C.merge(S[['snv'] + [f'{k}|flag' for k, *_ in FUNC]], on='snv', validate='1:1')
    C = pd.concat([C[C.cls == 'plp'].sort_values('pos'), C[C.cls == 'blb'].sort_values('pos')]).reset_index(drop=True)
    C['label'] = [f'c.{c}{s[-3]}>{s[-1]}' for s, c in zip(C.snv, C.c_dot)]
    C['region'] = C.oligo.map(region_of)
    n_plp, n_blb = int((C.cls == 'plp').sum()), int((C.cls == 'blb').sum())
    cv = pd.read_csv(d / 'clinvar_all.tsv', sep='\t')
    A = S[['snv', 'c_dot', 'mpra', 'in_library']].merge(cv[['snv', 'cls']], on='snv', how='left', validate='1:1')
    A['cls'] = A.cls.where(A.cls.isin(list(CLS)), 'none')
    assert set(A[A.cls.isin(['plp', 'blb'])].snv) == set(C.snv)

    # ---- layout (pt, y down)
    Wpt = W_IN * 72
    AH = 114.0                                                # panel a block height; panel b starts below
    x0, cw, gap, rh, sg = 150.0, 15.5, 8.0, 10.5, 2.5         # sg = gap between sections
    xs = [x0 + i * cw + (gap if i >= n_plp else 0) for i in range(len(C))]
    MX = (xs[-1] + cw + 30, xs[-1] + cw + 60)               # right edges of the margin columns
    assert MX[-1] <= Wpt - 4, MX
    y_bar, y_lab, y_reg = AH + 9.0, AH + 51.0, AH + 55.5
    ytop = AH + 64.0
    n_rows = len(POP) + len(FUNC)
    y_end = ytop + n_rows * rh + 2 * sg
    Hpt = y_end + 13
    fig = plt.figure(figsize=(W_IN, Hpt / 72))

    def axpt(x0_, y0_, x1_, y1_):
        return fig.add_axes([x0_ / Wpt, 1 - y1_ / Hpt, (x1_ - x0_) / Wpt, (y1_ - y0_) / Hpt])

    def letter(ch, y):
        fig.text(2 / Wpt, 1 - y / Hpt, ch, fontsize=9, fontweight='bold', va='top', ha='left')

    # ---- a: MPRA along the whole measured element (as Fig 3 a, extended to c.-318..c.-1)
    letter('a', 2)
    xlo, xhi = -321, 2
    axt = axpt(54, 2, Wpt - 6, 16)
    axt.set_xlim(xlo, xhi)
    axt.set_ylim(0, 1)
    axt.axis('off')
    ax = axpt(54, 16, Wpt - 6, 74)
    for short, o0, o1 in REGIONS:
        c0, c1 = o0 - C0 - 0.5, o1 - C0 + 0.5
        ax.axvspan(c0, c1, color='#ECE9F7', lw=0, zorder=0)
        axt.add_patch(Rectangle((c0, 0.06), c1 - c0, 0.36, fc='white', ec='#6A4FC9', lw=0.8,
                                ls=(0, (2, 1.2)) if short == 'Protected site' else '-'))
        al = A_ALIGN[short]                                   # neighbours anchored at their outer box edges
        axt.text({'right': c1, 'left': c0, 'center': (c0 + c1) / 2}[al], 0.5, A_LABELS[short], ha=al, va='bottom',
                 fontsize=FS2, color=INK, linespacing=1.0, multialignment=al)
    axt.plot([TSS_C, TSS_C, TSS_C + 7], [0.06, 0.6, 0.6], color=INK, lw=0.9)
    axt.annotate('', xy=(TSS_C + 11, 0.6), xytext=(TSS_C + 6, 0.6),
                 arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.9, mutation_scale=6, shrinkA=0, shrinkB=0))
    axt.text(TSS_C + 13, 0.6, 'TSS', ha='left', va='center', fontsize=FS2, color=INK, fontweight='bold')
    ax.axvline(TSS_C, color=INK2, lw=0.6, ls=(0, (2, 2)), zorder=1)
    ax.axhline(0, color=INK2, lw=0.5, zorder=1)
    for k in CLS:
        m = A.cls == k
        ax.scatter(A.c_dot[m], A.mpra[m], **CLS[k][1])
    ax.set_xlim(xlo, xhi)
    ax.set_xticks(range(-300, 0, 50), [f'c.{v}' for v in range(-300, 0, 50)], fontsize=FS)
    ax.tick_params(axis='y', labelsize=FS, length=2.5, width=0.7)
    ax.tick_params(axis='x', length=2.5, width=0.7)
    ax.set_xlabel('Position relative to the start codon (HGVS c., NM_000527.5)', fontsize=FS, labelpad=1)
    ax.set_ylabel('MPRA, HepG2\n(log2 change)', fontsize=FS, labelpad=2, linespacing=1.1)
    ax.set_yticks([-3, -2, -1, 0, 1])
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_linewidth(0.7)
    n = A.cls.value_counts()
    hs = [Line2D([], [], ls='none', marker=CLS[k][1].get('marker', 'o'), ms={'none': 2.5, 'vus': 3}.get(k, 4.5),
                 color=CLS[k][1]['color'], mec='white', mew=0.3) for k in CLS]
    fig.legend(hs, [f'{CLS[k][0]} ({n.get(k, 0)})' for k in CLS], loc='upper center', ncol=5, frameon=False,
               fontsize=FS, bbox_to_anchor=(0.5, 1 - 97 / Hpt), handletextpad=0.3, columnspacing=1.2)

    ov = fig.add_axes([0, 0, 1, 1], facecolor='none')
    ov.set_xlim(0, Wpt)
    ov.set_ylim(Hpt, 0)
    ov.axis('off')
    letter('b', AH)

    for (lo, hi, txt, col) in ((0, n_plp, f'Pathogenic / likely pathogenic ({n_plp})', PLP_C),
                               (n_plp, len(C), f'Benign / likely benign ({n_blb})', BLB_C)):
        xa, xb = xs[lo] + 1, xs[hi - 1] + cw - 1
        ov.plot([xa, xb], [y_bar + 3, y_bar + 3], color=col, lw=1.4, solid_capstyle='butt')
        ov.text((xa + xb) / 2, y_bar, txt, ha='center', va='bottom', fontsize=FS, fontweight='bold', color=INK)
    for xi, r in zip(xs, C.itertuples()):
        ov.text(xi + cw / 2, y_lab, r.label, rotation=90, ha='center', va='bottom', fontsize=FS, color=INK)
    i = 0                                                     # one region name per run of consecutive columns
    while i < len(C):
        j = i
        while j + 1 < len(C) and j + 1 != n_plp and C.region[j + 1] == C.region[i]:
            j += 1
        if C.region[i] != '–':
            xa, xb = xs[i] + 1.5, xs[j] + cw - 1.5
            if j > i:
                ov.plot([xa, xb], [y_reg - 3.5, y_reg - 3.5], color=MUTED, lw=0.5)
            ov.text((xa + xb) / 2, y_reg, C.region[i], ha='center', va='top', fontsize=FS2, color=INK2)
        i = j + 1
    for xm, txt in zip(MX, ('P/LP', 'B/LB')):
        ov.text(xm, ytop - 1, txt, ha='right', va='bottom', fontsize=FS, color=INK, fontweight='bold')

    def head(grp, y, n):                                      # tinted band, accent edge and section label (left)
        col = SEC[grp]
        rgb = matplotlib.colors.to_rgb(col)
        tint = tuple(1 - TINT * (1 - c) for c in rgb)
        ov.add_patch(Rectangle((2, y), MX[-1] + 2 - 2, n * rh, fc=tint, ec='none', zorder=0))
        ov.add_patch(Rectangle((2, y), 2.0, n * rh, fc=col, ec='none', zorder=0.5))
        name, sub = GROUPS[grp]
        ym = y + n * rh / 2
        ov.text(8, ym - 3.2, name, ha='left', va='center', fontsize=FS, fontweight='bold', color=col)
        ov.text(8, ym + 3.8, sub, ha='left', va='center', fontsize=FS2, color=col)
        return y

    def row_label(name, detail, yc):
        ov.text(x0 - 6, yc, f'{name} ({detail})' if detail else name, ha='right', va='center', fontsize=FS,
                color=INK)

    def margin(yc, mask, hit):
        for xm, cls in zip(MX, ('plp', 'blb')):
            m = (C.cls == cls) & mask
            ov.text(xm, yc, f'{int((m & hit).sum())}/{int(m.sum())}', ha='right', va='center', fontsize=FS, color=INK)

    rows_out, y = [], head('Population', ytop, len(POP))
    for k, name, detail in POP:
        yc = y + rh / 2
        row_label(name, detail, yc)
        for xi, r in zip(xs, C.itertuples()):
            if getattr(r, k):
                ov.scatter(xi + cw / 2, yc, marker='o', s=24, color=SEC['Population'], lw=0)
            else:
                ov.scatter(xi + cw / 2, yc, s=3, color=DOT, lw=0)
        margin(yc, pd.Series(True, index=C.index), C[k].astype(bool))
        rows_out.append(dict(row=name, kind='population', plp=int((C[k] & (C.cls == 'plp')).sum()),
                             blb=int((C[k] & (C.cls == 'blb')).sum())))
        y += rh
    y += sg
    for grp in ('Predicted', 'Measured'):
        y = head(grp, y, sum(g == grp for _, g, _, _ in FUNC))
        for k, g, name, detail in FUNC:
            if g != grp:
                continue
            yc = y + rh / 2
            row_label(name, detail, yc)
            cov = C[k].notna()
            for xi, r, c in zip(xs, C.itertuples(), cov):
                if not c:
                    ov.text(xi + cw / 2, yc, '–', ha='center', va='center', fontsize=FS, color=MUTED)
                elif C.loc[r.Index, f'{k}|flag']:
                    ov.scatter(xi + cw / 2, yc, marker='^' if C.loc[r.Index, k] > 0 else 'v', s=28, color=SEC[grp],
                               lw=0)
                else:
                    ov.scatter(xi + cw / 2, yc, s=3, color=DOT, lw=0)
            hit = C[f'{k}|flag'] & cov
            margin(yc, cov, hit)
            rows_out.append(dict(row=f'{name} {detail}'.strip(), kind='change', key=k, cutoff=cut[k],
                                 n_scored=int(S[k].notna().sum()), plp=int((hit & (C.cls == 'plp')).sum()),
                                 plp_n=int((cov & (C.cls == 'plp')).sum()), blb=int((hit & (C.cls == 'blb')).sum()),
                                 blb_n=int((cov & (C.cls == 'blb')).sum())))
            y += rh
        y += sg
    y -= sg
    sp = '\u2003\u2003'                                       # em spaces: SVG collapses ordinary runs
    ov.text(4, y + 7.5, f'●\u2002yes{sp}▲ / ▼\u2002large increase / decrease (top 20% of |variant − WT|){sp}·\u2002no / '
            f'smaller change{sp}–\u2002not measured', ha='left', va='center', fontsize=FS2, color=INK2)
    assert y + 11 <= Hpt

    svg = out / 'fig3.svg'
    fig.savefig(svg, facecolor='white')
    plt.close(fig)
    svg.write_text(expand_fonts(svg.read_text().replace("'Arimo'", 'Arial, Helvetica, sans-serif')
                                .replace('Arimo', 'Arial'))[0])
    env = dict(os.environ, FONTCONFIG_FILE=str(FONTS_CONF))
    subprocess.run([RSVG, '-f', 'pdf', '-o', str(out / 'fig3.pdf'), str(svg)], check=True, env=env)
    subprocess.run([RSVG, '-d', '300', '-p', '300', '-f', 'png', '-o', str(out / 'fig3.png'), str(svg)], check=True,
                   env=env)
    cols = ['label', 'cls', 'clinvar_desc', 'stars', 'accession', 'region', 'in_library', 'gnomad', 'gtex_present',
            'gtex_tested', 'primate_alt', 'primate_alt_species'] + [c for k, *_ in FUNC for c in (k, f'{k}|flag')]
    C[cols].to_csv(out / 'matrix.tsv', sep='\t', index=False, float_format='%.4g')
    pd.DataFrame(rows_out).to_csv(out / 'margins.tsv', sep='\t', index=False, float_format='%.4g')
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()
    json.dump(dict(snvs=len(S), controls=dict(plp=n_plp, blb=n_blb), cutoffs=cut, rows=rows_out,
                   size_in=[W_IN, round(Hpt / 72, 3)], date=str(date.today()), commit=commit),
              open(out / 'summary.json', 'w'), indent=1)
    print(pd.DataFrame(rows_out).to_string(index=False))


if __name__ == '__main__':
    main()

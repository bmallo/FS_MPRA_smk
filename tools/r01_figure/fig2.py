"""Grant Figure 2 (draft 3): SLC39A4 c.-169A>G disrupts the NF-Y footprint and promoter activity (fig2_rules.json).
Output is one Illustrator-editable SVG (live Arial text), plus PDF and PNG renders (rsvg-convert with fonts.conf).

Inputs are the SLC39A4_variant_slide project's outputs (HEK293T plasmid Fiber-seq, WT and c.-169A>G co-transfected,
3 replicates): analysis/occ_<rep>_<allele>.npy (reads with a footprint of size s covering plasmid position p),
n_<rep>_<allele>.txt (kept reads), reads.tsv.gz (per read: allele, kept, site_fp), motifs/MA0060.4.jaspar, and
Fig. S15A (HEK293T half) as vector: fig2_s15_svg.py exports Supplemental.pdf page 27 to S15_page_runs.svg, which is
embedded here clipped to S15_CROP (the raw luciferase values are not on klone). TF-sized footprint = 20-60 bp.
a  construct map (SnapGene style, to scale, plasmid MAP_P in gene orientation) with c.-169A>G and the NF-Y logo
b  luciferase (Fig. S15A, HEK293T)
c  TF-sized occupancy around the NF-Y motif, WT vs variant, reps pooled, Gaussian sigma SIGMA bp
d  per replicate: % of kept reads with a 20-60 bp footprint over the variant base; Fisher exact test per replicate
e  whole plasmid: WT - variant TF-sized occupancy, mean of the 3 replicates, Gaussian sigma SIGMA bp
Usage: python tools/r01_figure/fig2.py --project <SLC39A4_variant_slide> --output <dir>
"""
import argparse
import json
import os
import re
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
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.patches import PathPatch, Polygon, Rectangle  # noqa: E402
from matplotlib.textpath import TextPath  # noqa: E402
from matplotlib.transforms import Affine2D  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

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

INK, INK2 = '#2E2E2E', '#666666'
WT_C, VAR_C, TF_C = '#2E3192', '#E3191F', '#338B8C'          # WT / variant as Fig. S15A; TF teal as Figure 1
PROM_C, LUC_C = '#B9BCE6', '#F2C14E'
LETTER_C = {'A': '#109648', 'C': '#255C99', 'G': '#F7B32B', 'T': '#D62839'}
REPS = ['rep_1', 'rep_2', 'rep_3']
L = 6466                                  # plasmid length (bp)
PROM, LUC, VPOS1 = (5215, 6273), (2916, 4568), 6120          # 1-based plasmid features; variant c.-169
ATG_A = 5951                              # 1-based plasmid position of the ATG A; c.-x = ATG_A + x
MOTIF_P = ATG_A + 170                     # NF-Y CCAAT core centre (c.-170)
TF = (20, 60)
SIGMA, HALF = 3, 150
W_IN, H_IN = 7.2, 5.5
S15_CROP = (56.0, 68.0, 140.0, 216.0)        # x, y, w, h (PDF pt): Fig. S15A HEK293T half, left of the divider
S15_FRAME_Y = (80.07, 252.83)                 # S15A axes frame (rect x 100.09-322.07, stroke 1 pt): close it at the crop
MAP_P = (6400, 2800)                          # panel a map window, plasmid coordinates (gene orientation: left = 6400)


def letter(ax, ch, x, y, w, h, color, prop):
    tp = TextPath((0, 0), ch, size=1, prop=prop)
    bb = tp.get_extents()
    t = Affine2D().translate(-bb.x0, -bb.y0).scale(w / bb.width, h / bb.height).translate(x, y)
    ax.add_patch(PathPatch(t.transform_path(tp), fc=color, ec='none', zorder=2))


def feature_arrow(ax, x0, x1, y, h, color, label, tc, gid, lx=None):
    """SnapGene-style feature: thin box with a pointed 3' end, label inside."""
    tip = min(6.0, (x1 - x0) * 0.3)
    ax.add_patch(Polygon([[x0, y - h / 2], [x1 - tip, y - h / 2], [x1, y], [x1 - tip, y + h / 2], [x0, y + h / 2]],
                         fc=color, ec='#4A4A4A', lw=0.4, zorder=2, gid=gid))
    ax.text((x0 + x1 - tip) / 2 if lx is None else lx, y + 0.3, label, fontsize=6.3, fontweight='bold', color=tc,
            ha='center' if lx is None else 'left', va='center', zorder=3)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--project', required=True)
    ap.add_argument('--s15', required=True, help='S15_page_runs.svg from fig2_s15_svg.py')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    P, out = Path(a.project), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    A = P / 'analysis'
    occ = {al: {r: np.load(A / f'occ_{r}_{al}.npy') for r in REPS} for al in ('WT', 'Var')}
    n = {al: {r: int(open(A / f'n_{r}_{al}.txt').read()) for r in REPS} for al in ('WT', 'Var')}
    tf = {al: {r: occ[al][r][TF[0]:TF[1] + 1].sum(0) / n[al][r] for r in REPS} for al in occ}   # share per position
    pooled = {al: sum(occ[al][r][TF[0]:TF[1] + 1].sum(0) for r in REPS) / sum(n[al].values()) for al in occ}
    sm = lambda v: gaussian_filter1d(v.astype(float), SIGMA, mode='wrap')   # noqa: E731  (circular plasmid)

    # d: per-replicate site test from per-read calls
    reads = pd.read_csv(A / 'reads.tsv.gz', sep='\t')
    reads = reads[reads.kept == 1]
    reads['tf_site'] = reads.site_fp.between(*TF)
    stats = []
    for r in REPS:
        g = reads[reads.rep == r]
        t = pd.crosstab(g.allele, g.tf_site).reindex(index=['WT', 'Var'], columns=[True, False], fill_value=0)
        orr, p = fisher_exact(t.values)
        stats.append(dict(rep=r, n_wt=int(t.loc['WT'].sum()), n_var=int(t.loc['Var'].sum()),
                          fp_wt=int(t.loc['WT', True]), fp_var=int(t.loc['Var', True]),
                          pct_wt=100 * t.loc['WT', True] / t.loc['WT'].sum(),
                          pct_var=100 * t.loc['Var', True] / t.loc['Var'].sum(), odds_ratio=orr, fisher_p=p))
    S = pd.DataFrame(stats)
    S.to_csv(out / 'site_stats.tsv', sep='\t', index=False, float_format='%.4g')
    assert all(S.n_wt == [n['WT'][r] for r in REPS]) and all(S.n_var == [n['Var'][r] for r in REPS])

    # c: source table
    p1 = np.arange(MOTIF_P - HALF, MOTIF_P + HALF + 1)
    xc = MOTIF_P - p1                                              # gene orientation: upstream (5') on the left
    csm = {al: sm(pooled[al])[p1 - 1] for al in pooled}
    pd.DataFrame({'plasmid_pos': p1, 'bp_from_nfy': xc, 'wt_share': pooled['WT'][p1 - 1],
                  'var_share': pooled['Var'][p1 - 1], 'wt_smooth': csm['WT'], 'var_smooth': csm['Var']}).to_csv(
        out / 'panel_c_occupancy.tsv', sep='\t', index=False, float_format='%.5g')
    # e: source table
    dlt = {r: sm(tf['WT'][r] - tf['Var'][r]) * 100 for r in REPS}
    dmean = np.mean([dlt[r] for r in REPS], axis=0)                # mean of the 3 replicates (each weighted equally)
    pd.DataFrame({'plasmid_pos': np.arange(1, L + 1), **{f'delta_pp_{r}': dlt[r] for r in REPS},
                  'delta_pp_mean': dmean}).to_csv(
        out / 'panel_e_delta.tsv', sep='\t', index=False, float_format='%.4g')
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()

    fig = plt.figure(figsize=(W_IN, H_IN))
    Wpt, Hpt = W_IN * 72, H_IN * 72

    def axpt(x0, y0, x1, y1, **kw):
        return fig.add_axes([x0 / Wpt, 1 - y1 / Hpt, (x1 - x0) / Wpt, (y1 - y0) / Hpt], **kw)

    def letter_lab(ch, x, y):
        fig.text(x / Wpt, 1 - y / Hpt, ch, fontsize=10, fontweight='bold', va='top', ha='left')

    # ---- a: schematic (points, y down)
    letter_lab('a', 4, 6)
    ax = axpt(14, 14, 318, 150)
    ax.set_xlim(14, 318)
    ax.set_ylim(150, 14)
    ax.axis('off')
    yb, bh = 40, 10                                                # to scale, gene orientation (plasmid MAP_P0 -> MAP_P1)
    mx0, mx1 = 24.0, 306.0
    xm = lambda p: mx0 + (MAP_P[0] - p) / (MAP_P[0] - MAP_P[1]) * (mx1 - mx0)   # noqa: E731
    ax.plot([mx0 - 6, mx1 + 6], [yb, yb], color='#4A4A4A', lw=0.8, solid_capstyle='butt', zorder=1, gid='backbone')
    for xe in (mx0 - 6, mx1 + 6):                                  # backbone continues (break marks)
        ax.plot([xe - 1.5, xe + 1.5], [yb + 3, yb - 3], color='#4A4A4A', lw=0.6)
    feature_arrow(ax, xm(PROM[1]), xm(PROM[0]), yb, bh, PROM_C, 'SLC39A4 promoter', WT_C, 'feature-promoter',
                  lx=xm(VPOS1) + 3)                             # label starts right of the variant tick
    feature_arrow(ax, xm(LUC[1]), xm(LUC[0]), yb, bh, LUC_C, 'Luciferase', INK, 'feature-luciferase')
    xv = xm(VPOS1)
    ax.plot([xv, xv], [yb - bh / 2 - 1, yb + bh / 2 + 3], color=VAR_C, lw=0.9, zorder=4)
    ax.add_patch(Polygon([[xv, yb + bh / 2 + 3], [xv - 3.2, yb + bh / 2 + 8.5], [xv + 3.2, yb + bh / 2 + 8.5]],
                         fc=VAR_C, ec='none', zorder=4))
    ax.text(xv + 5, yb + bh / 2 + 7, 'c.-169A>G', fontsize=6.5, fontweight='bold', color=VAR_C, va='center')
    ax.plot([mx1 - xm(0) + xm(500), mx1], [yb - 13, yb - 13], color=INK, lw=0.7)
    ax.text(mx1 - (xm(0) - xm(500)) / 2, yb - 15, '500 bp', fontsize=5.5, color=INK, ha='center', va='bottom')
    cnt = np.array([[float(x) for x in ln.split('[')[1].split(']')[0].split()]
                    for ln in open(P / 'motifs/MA0060.4.jaspar').read().splitlines()[1:5]])
    pr = (cnt + 0.8) / (cnt.sum(0) + 3.2)
    ic = 2 + (pr * np.log2(pr)).sum(0)
    cw, bit, ly0 = 13.0, 15.0, 128.0                               # column width, pt per bit, logo baseline
    lx0 = xv + 40
    ax.add_patch(Polygon([[xv - 1.5, yb + bh / 2 + 11], [xv + 1.5, yb + bh / 2 + 11],
                          [lx0 + 8 * cw, ly0 - 2 * bit - 6], [lx0, ly0 - 2 * bit - 6]], fc='#E9E9F4', ec='none',
                         zorder=0))
    bold = FontProperties(family='Arimo', weight='bold')
    lax = axpt(14, 14, 318, 150)                                   # logo in a y-up copy of the same frame
    lax.set_xlim(14, 318)
    lax.set_ylim(150, 14)
    lax.axis('off')
    for j in range(8):
        y = ly0
        for k in np.argsort(pr[:, j] * ic[j]):
            h = pr[k, j] * ic[j] * bit
            if h > 0.7:
                tp = TextPath((0, 0), 'ACGT'[k], size=1, prop=bold)
                bb = tp.get_extents()
                tr = Affine2D().translate(-bb.x0, -bb.y0).scale((cw - 2) / bb.width, -h / bb.height).translate(
                    lx0 + j * cw + 1, y)
                lax.add_patch(PathPatch(tr.transform_path(tp), fc=LETTER_C['ACGT'[k]], ec='none', zorder=2))
            y -= h
    ax.plot([lx0 - 2, lx0 + 8 * cw + 2], [ly0 + 1.5, ly0 + 1.5], color=INK2, lw=0.6)
    vj = 3                                                         # CCA[A]TCAC: variant = motif position 4
    ax.add_patch(Rectangle((lx0 + vj * cw, ly0 - 2 * bit - 3), cw, 2 * bit + 7, fc='none', ec=VAR_C, lw=1.0,
                           ls=(0, (3, 1.5)), zorder=3))
    ax.text(lx0 + (vj + 0.5) * cw, ly0 + 6, 'A>G', fontsize=6.5, fontweight='bold', color=VAR_C, ha='center',
            va='top')
    ax.text(lx0 + 8 * cw + 8, ly0 - bit, 'NF-Y motif', fontsize=6.5, fontweight='bold', color=INK, va='center')

    # ---- b: luciferase (Fig. S15A, HEK293T)
    letter_lab('b', 330, 6)
    b_slot, b_h = (352.0, 12.0), 140.0                              # S15 vector embedded here after saving
    b_sc = b_h / S15_CROP[3]

    # ---- c: TF occupancy around the motif
    letter_lab('c', 4, 160)
    cax = axpt(40, 172, 300, 290)
    for al, col, lab in (('WT', WT_C, 'WT'), ('Var', VAR_C, 'c.-169A>G')):
        cax.fill_between(xc, 0, csm[al] * 100, color=col, alpha=0.12, lw=0)
        cax.plot(xc, csm[al] * 100, color=col, lw=1.2, label=lab)
    cax.axvspan(-4, 4, color='#E9E9F4', lw=0, zorder=0)
    cax.set_xlim(-HALF, HALF)
    cax.set_ylim(0, None)
    cax.set_xticks([-150, -100, -50, 0, 50, 100, 150], ['−150', '−100', '−50', 'NF-Y', '+50', '+100', '+150'])
    cax.set_xlabel('Distance from NF-Y motif (bp; 5′ ← → 3′ of SLC39A4)', fontsize=6.5)
    cax.set_ylabel('Reads with a 20–60 bp\nfootprint (%)', fontsize=6.5)
    cax.spines[['top', 'right']].set_visible(False)
    cax.legend(frameon=False, fontsize=6, loc='upper right', handlelength=1.5)

    # ---- d: per-replicate site occupancy
    letter_lab('d', 330, 160)
    dax = axpt(392, 172, 500, 290)
    x = np.arange(2)
    for i, row in S.iterrows():
        dax.plot(x, [row.pct_wt, row.pct_var], color='#9A9A9A', lw=0.8, zorder=1)
        dax.scatter(x, [row.pct_wt, row.pct_var], c=[WT_C, VAR_C], s=18, zorder=2, edgecolor='white', lw=0.4)
        dax.text(-0.09, row.pct_wt, f'Rep {i + 1}', fontsize=5.5, color=INK2, va='center', ha='right')
    dax.set_xticks(x, ['WT', 'c.-169A>G'])
    dax.get_xticklabels()[0].set_color(WT_C)
    dax.get_xticklabels()[1].set_color(VAR_C)
    dax.set_xlim(-0.5, 1.4)
    dax.set_ylim(0, max(S.pct_wt) * 1.25)
    dax.set_ylabel('Reads with a 20–60 bp footprint\nover c.-169 (%)', fontsize=6.5)
    dax.spines[['top', 'right']].set_visible(False)
    pmax = S.fisher_p.max()
    ptxt = f'Fisher p < 10$^{{{int(np.floor(np.log10(pmax))) + 1}}}$' if pmax > 0 else 'Fisher p < 1e-300'
    dax.text(0.5, max(S.pct_wt) * 1.16, f'each replicate\n{ptxt}', ha='center', va='center', fontsize=5.8,
             color=INK)

    # ---- e: whole-plasmid delta
    letter_lab('e', 4, 314)
    eax = axpt(40, 324, 500, 380)
    pos = np.arange(1, L + 1)
    eax.plot(pos, dmean, color=TF_C, lw=0.8)
    eax.axhline(0, color=INK, lw=0.5)
    eax.set_xlim(1, L)
    eax.set_ylabel('Δ occupancy\n(WT − var., %)', fontsize=6.5)
    eax.spines[['top', 'right']].set_visible(False)
    eax.set_xticks([])
    yl = eax.get_ylim()
    eax.annotate('c.-169 (NF-Y)', xy=(VPOS1, dmean[VPOS1 - 1]), xytext=(VPOS1 - 700, yl[1] * 0.85),
                 fontsize=6, ha='right', va='center', arrowprops=dict(arrowstyle='-|>', lw=0.6, color=INK,
                                                                      mutation_scale=6))
    fax = axpt(40, 382, 500, 392)                                  # plasmid feature bar
    fax.set_xlim(1, L)
    fax.set_ylim(0, 1)
    fax.axis('off')
    fax.add_patch(Rectangle((1, 0.42), L - 1, 0.16, fc='#8B9094', lw=0))
    fax.add_patch(Rectangle((PROM[0], 0.15), PROM[1] - PROM[0], 0.7, fc=PROM_C, lw=0))
    fax.add_patch(Rectangle((LUC[0], 0.15), LUC[1] - LUC[0], 0.7, fc=LUC_C, lw=0))
    fax.text((PROM[0] + PROM[1]) / 2, 0.5, 'promoter', fontsize=5.5, ha='center', va='center', color=WT_C)
    fax.text((LUC[0] + LUC[1]) / 2, 0.5, 'luciferase', fontsize=5.5, ha='center', va='center', color=INK)

    svg = out / 'fig2.svg'
    fig.savefig(svg, facecolor='white', dpi=300)
    plt.close(fig)
    s = expand_fonts(svg.read_text().replace("'Arimo'", 'Arial, Helvetica, sans-serif').replace('Arimo', 'Arial'))[0]
    body = re.sub(r'^.*?<svg[^>]*>', '', Path(a.s15).read_text(), count=1, flags=re.S).rsplit('</svg>', 1)[0]
    body = re.sub(r'id="([^"]+)"', r'id="s15-\1"', body)
    body = re.sub(r'url\(#([^)]+)\)', r'url(#s15-\1)', body)
    body = re.sub(r'href="#([^"]+)"', r'href="#s15-\1"', body)
    cx, cy, cw, ch = S15_CROP
    group = (f'<g id="panel-b-luciferase-S15A" transform="translate({b_slot[0]:.3f} {b_slot[1]:.3f}) '
             f'scale({b_sc:.5f}) translate({-cx:.3f} {-cy:.3f})"><defs><clipPath id="s15-crop"><rect x="{cx}" '
             f'y="{cy}" width="{cw}" height="{ch}"/></clipPath></defs><g clip-path="url(#s15-crop)">{body}</g>'
             f'<line id="s15-frame-right" x1="{cx + cw - 0.5}" y1="{S15_FRAME_Y[0]}" x2="{cx + cw - 0.5}" '
             f'y2="{S15_FRAME_Y[1]}" stroke="#231f20" stroke-width="1"/></g>')
    s = s.replace('</svg>', group + '\n</svg>')
    if 'xmlns:inkscape' not in s.split('>', 2)[1]:
        s = re.sub(r'<svg ', '<svg xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" ', s, count=1)
    svg.write_text(s)
    env = dict(os.environ, FONTCONFIG_FILE=str(FONTS_CONF))
    subprocess.run([RSVG, '-f', 'pdf', '-o', str(out / 'fig2.pdf'), str(svg)], check=True, env=env)
    subprocess.run([RSVG, '-d', '300', '-p', '300', '-f', 'png', '-o', str(out / 'fig2.png'), str(svg)], check=True,
                   env=env)
    peak = {al: (int(xc[np.argmax(csm[al])]), float(csm[al].max() * 100)) for al in csm}
    json.dump(dict(reads=n, tf_bin=TF, sigma=SIGMA, peak_bp_from_nfy_and_pct=peak,
                   at_variant_pct={al: float(sm(pooled[al])[VPOS1 - 1] * 100) for al in pooled},
                   delta_at_variant_pp={r: float(dlt[r][VPOS1 - 1]) for r in REPS},
                   max_abs_delta_outside_promoter_pp={r: float(np.abs(np.r_[dlt[r][:PROM[0] - 1], dlt[r][PROM[1]:]]).max())
                                                      for r in REPS},
                   commit=commit), open(out / 'summary.json', 'w'), indent=1)
    print(S.to_string(index=False))


if __name__ == '__main__':
    main()

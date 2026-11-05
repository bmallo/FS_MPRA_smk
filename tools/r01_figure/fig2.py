"""Grant Figure 2 (draft 8): SLC39A4 c.-169A>G disrupts the NF-Y footprint and promoter activity (fig2_rules.json,
revision_draft5-8). Output is one Illustrator-editable SVG (live Arial text), plus PDF and PNG renders (rsvg-convert with
fonts.conf).

Inputs: the SLC39A4_variant_slide project's outputs (HEK293T plasmid Fiber-seq, WT and c.-169A>G co-transfected,
3 replicates) and the same outputs for HepG2 (4 replicates; fig2_hepg2.sh, identical processing):
analysis/occ_<rep>_<allele>.npy (reads with a footprint of size s covering plasmid position p), n_<rep>_<allele>.txt
(kept reads), reads.tsv.gz (per read: allele, kept, site_fp); motifs/MA0060.4.jaspar; and Fig. S15A (HEK293T half) as vector: fig2_s15_svg.py exports Supplemental.pdf page 27 to
S15_page_runs.svg, embedded here clipped to S15_CROP (the raw luciferase values are not on klone).
TF-sized footprint = 20-60 bp. Cell lines are never pooled or overlaid.
a  whole plasmid: WT - variant TF-sized occupancy (HEK293T, mean of the per-replicate tracks, Gaussian sigma SIGMA bp)
   above the construct map (to scale, gene orientation from the promoter's 5' end, wrapping round), variant tick and
   the NF-Y logo as a zoom from the c.-169 peak
b  luciferase (Fig. S15A, HEK293T)
Per-replicate site tests (site_stats_<cell>.tsv) and the HepG2 whole-plasmid track are written for the legend only.
Usage: python tools/r01_figure/fig2.py --project <SLC39A4_variant_slide> --hepg2 <fig2/hepg2> --s15 <svg> --output <dir>
"""
import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker  # noqa: E402
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
REPS = {'HEK293T': ['rep_1', 'rep_2', 'rep_3'], 'HepG2': ['rep_1', 'rep_2', 'rep_3', 'rep_4']}
L = 6466                                  # plasmid length (bp)
PROM, LUC, VPOS1 = (5215, 6273), (2916, 4568), 6120          # 1-based plasmid features; variant c.-169
ATG_A = 5951                              # 1-based plasmid position of the ATG A; c.-x = ATG_A + x
MOTIF_P = ATG_A + 170                     # NF-Y CCAAT core centre (c.-170)
TF = (20, 60)
SIGMA = 3
W_IN, H_IN = 7.2, 2.1
S15_CROP = (56.0, 68.0, 140.0, 216.0)        # x, y, w, h (PDF pt): Fig. S15A HEK293T half, left of the divider
S15_FRAME_Y = (80.07, 252.83)                 # S15A axes frame (rect x 100.09-322.07, stroke 1 pt): close it at the crop
B_H = 147.0                                   # panel b height (pt): S15 text 8.9-12 pt -> >= 6 pt


def gx(p):
    """Plasmid position -> bp from the promoter fragment's 5' end, gene orientation, wrapping round the circle."""
    return (PROM[1] - np.asarray(p)) % L


def feature_arrow(ax, x0, x1, y, h, color, label, tc, gid, lx=None, fs=6.3):
    """SnapGene-style feature: thin box with a pointed 3' end, label inside."""
    tip = min(6.0, (x1 - x0) * 0.3)
    ax.add_patch(Polygon([[x0, y - h / 2], [x1 - tip, y - h / 2], [x1, y], [x1 - tip, y + h / 2], [x0, y + h / 2]],
                         fc=color, ec='#4A4A4A', lw=0.4, zorder=2, gid=gid))
    ax.text((x0 + x1 - tip) / 2 if lx is None else lx, y + 0.3, label, fontsize=fs, fontweight='bold', color=tc,
            ha='center' if lx is None else 'left', va='center', zorder=3)


def load(adir, reps):
    occ = {al: {r: np.load(adir / f'occ_{r}_{al}.npy') for r in reps} for al in ('WT', 'Var')}
    n = {al: {r: int(open(adir / f'n_{r}_{al}.txt').read()) for r in reps} for al in ('WT', 'Var')}
    return occ, n


def site_stats(adir, reps, n):
    """Per replicate: kept reads with a 20-60 bp footprint covering the variant base, WT vs variant; Fisher exact."""
    reads = pd.read_csv(adir / 'reads.tsv.gz', sep='\t')
    reads = reads[reads.kept == 1]
    reads['tf_site'] = reads.site_fp.between(*TF)
    stats = []
    for r in reps:
        g = reads[reads.rep == r]
        t = pd.crosstab(g.allele, g.tf_site).reindex(index=['WT', 'Var'], columns=[True, False], fill_value=0)
        orr, p = fisher_exact(t.values)
        stats.append(dict(rep=r, n_wt=int(t.loc['WT'].sum()), n_var=int(t.loc['Var'].sum()),
                          fp_wt=int(t.loc['WT', True]), fp_var=int(t.loc['Var', True]),
                          pct_wt=100 * t.loc['WT', True] / t.loc['WT'].sum(),
                          pct_var=100 * t.loc['Var', True] / t.loc['Var'].sum(), odds_ratio=orr, fisher_p=p))
    S = pd.DataFrame(stats)
    assert all(S.n_wt == [n['WT'][r] for r in reps]) and all(S.n_var == [n['Var'][r] for r in reps])
    return S


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--project', required=True)
    ap.add_argument('--hepg2', required=True, help='fig2_hepg2.sh work dir (holds analysis/)')
    ap.add_argument('--s15', required=True, help='S15_page_runs.svg from fig2_s15_svg.py')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    P, out = Path(a.project), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    adirs = {'HEK293T': P / 'analysis', 'HepG2': Path(a.hepg2) / 'analysis'}
    sm = lambda v: gaussian_filter1d(v.astype(float), SIGMA, mode='wrap')   # noqa: E731  (circular plasmid)
    occ, n, S, dlt, spec = {}, {}, {}, {}, {}
    for cell, reps in REPS.items():
        occ[cell], n[cell] = load(adirs[cell], reps)
        S[cell] = site_stats(adirs[cell], reps, n[cell])
        S[cell].to_csv(out / f'site_stats_{cell}.tsv', sep='\t', index=False, float_format='%.4g')
        tf = {al: {r: occ[cell][al][r][TF[0]:TF[1] + 1].sum(0) / n[cell][al][r] for r in reps} for al in ('WT', 'Var')}
        dlt[cell] = {r: sm(tf['WT'][r] - tf['Var'][r]) * 100 for r in reps}           # panel a (HEK293T)
        dmean = np.mean([dlt[cell][r] for r in reps], axis=0)                       # replicates weighted equally
        pd.DataFrame({'plasmid_pos': np.arange(1, L + 1), **{f'delta_pp_{r}': dlt[cell][r] for r in reps},
                      'delta_pp_mean': dmean}).to_csv(out / f'whole_plasmid_delta_{cell}.tsv', sep='\t', index=False,
                                                      float_format='%.4g')
        outside = np.r_[dmean[:PROM[0] - 1], dmean[PROM[1]:]]
        spec[cell] = dict(mean_delta_at_variant_pp=float(dmean[VPOS1 - 1]),
                          max_abs_mean_delta_outside_promoter_pp=float(np.abs(outside).max()),
                          delta_at_variant_pp={r: float(dlt[cell][r][VPOS1 - 1]) for r in reps})

    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()
    dmean = pd.read_csv(out / 'whole_plasmid_delta_HEK293T.tsv', sep='\t').delta_pp_mean.to_numpy()   # as drawn
    order = np.argsort(gx(np.arange(1, L + 1)))
    xg, yg = gx(np.arange(1, L + 1))[order], dmean[order]

    fig = plt.figure(figsize=(W_IN, H_IN))
    Wpt, Hpt = W_IN * 72, H_IN * 72

    def axpt(x0, y0, x1, y1, **kw):
        return fig.add_axes([x0 / Wpt, 1 - y1 / Hpt, (x1 - x0) / Wpt, (y1 - y0) / Hpt], **kw)

    def letter_lab(ch, x, y):
        fig.text(x / Wpt, 1 - y / Hpt, ch, fontsize=10, fontweight='bold', va='top', ha='left')

    # ---- a: whole-plasmid WT - variant track over the construct map (shared x, to scale)
    letter_lab('a', 2, 2)
    X0, X1, Y0, Y1 = 50.0, 398.0, 8.0, 104.0                      # track axes (pt)
    ylo, yhi = -1.5, 8.5
    tax = axpt(X0, Y0, X1, Y1)
    tax.axhline(0, color=INK2, lw=0.5, zorder=1)
    tax.plot(xg, yg, color=TF_C, lw=0.8, zorder=2)
    tax.set_xlim(0, L)
    tax.set_ylim(ylo, yhi)
    tax.set_xticks([])
    tax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5, integer=True))
    tax.spines[['top', 'right', 'bottom']].set_visible(False)
    tax.set_ylabel('Δ 20–60 bp footprints\n(WT − variant, % reads)', fontsize=6.3, labelpad=2)
    tax.text(1.0, 0.98, 'HEK293T', transform=tax.transAxes, fontsize=6.3, fontweight='bold', color=INK, ha='right',
             va='top')
    ov = axpt(0, 0, Wpt, Hpt)                                       # overlay in pt (map, logo, labels)
    ov.set_xlim(0, Wpt)
    ov.set_ylim(Hpt, 0)
    ov.axis('off')
    xp = lambda g: X0 + np.asarray(g) / L * (X1 - X0)              # noqa: E731
    yp = lambda v: Y0 + (yhi - v) / (yhi - ylo) * (Y1 - Y0)         # noqa: E731
    yb, bh = Y1 + 10, 9                                             # map backbone y, feature height
    ov.plot([X0, X1], [yb, yb], color='#4A4A4A', lw=0.8, zorder=1, solid_capstyle='butt')
    feature_arrow(ov, xp(0), xp(gx(PROM[0])), yb, bh, PROM_C, '', WT_C, 'feature-promoter')
    feature_arrow(ov, xp(gx(LUC[1])), xp(gx(LUC[0])), yb, bh, LUC_C, 'Luciferase', INK, 'feature-luciferase')
    ov.text(xp(gx(PROM[0])) / 2 + X0 / 2, yb + bh / 2 + 2.5, 'SLC39A4 promoter', fontsize=6.3, fontweight='bold',
            color=WT_C, ha='center', va='top')
    xv = xp(gx(VPOS1))
    ov.plot([xv, xv], [yb - bh / 2 - 1, yb + bh / 2 + 1], color=VAR_C, lw=1.0, zorder=4)
    kb = xp(1000) - xp(0)
    ov.plot([X1 - kb, X1], [yb + bh / 2 + 5, yb + bh / 2 + 5], color=INK, lw=0.7)
    ov.text(X1 - kb / 2, yb + bh / 2 + 6.5, '1 kb', fontsize=6, color=INK, ha='center', va='top')
    # NF-Y logo, zoomed from the c.-169 peak
    cnt = np.array([[float(x) for x in ln.split('[')[1].split(']')[0].split()]
                    for ln in open(P / 'motifs/MA0060.4.jaspar').read().splitlines()[1:5]])
    pr = (cnt + 0.8) / (cnt.sum(0) + 3.2)
    ic = 2 + (pr * np.log2(pr)).sum(0)
    cw, bit = 9.0, 11.0                                             # column width, pt per bit
    apex = (xv, yp(dmean[VPOS1 - 1]))
    lx0, ly0 = xv + 34, Y0 + 44                                      # logo left edge, baseline
    ov.add_patch(Polygon([apex, [lx0 - 2, ly0 - 2 * bit - 3], [lx0 - 2, ly0 + 2]], fc='#E9E9F4', ec='none', zorder=0))
    bold = FontProperties(family='Arimo', weight='bold')
    for j in range(8):
        y = ly0
        for k in np.argsort(pr[:, j] * ic[j]):
            h = pr[k, j] * ic[j] * bit
            if h > 0.5:
                tp = TextPath((0, 0), 'ACGT'[k], size=1, prop=bold)
                bb = tp.get_extents()
                tr = Affine2D().translate(-bb.x0, -bb.y0).scale((cw - 1.4) / bb.width, -h / bb.height).translate(
                    lx0 + j * cw + 0.7, y)
                ov.add_patch(PathPatch(tr.transform_path(tp), fc=LETTER_C['ACGT'[k]], ec='none', zorder=2))
            y -= h
    ov.plot([lx0 - 1, lx0 + 8 * cw + 1], [ly0 + 1, ly0 + 1], color=INK2, lw=0.5)
    vj = 3                                                          # CCA[A]TCAC: variant = motif position 4
    ov.add_patch(Rectangle((lx0 + vj * cw, ly0 - 2 * bit - 2), cw, 2 * bit + 5, fc='none', ec=VAR_C, lw=0.8,
                           ls=(0, (2.5, 1.2)), zorder=3))
    ov.text(lx0 + (vj + 0.5) * cw, ly0 + 4.5, 'c.-169A>G', fontsize=6.3, fontweight='bold', color=VAR_C,
            ha='center', va='top')
    ov.text(lx0 + 8 * cw + 5, ly0 - bit, 'NF-Y motif', fontsize=6.3, fontweight='bold', color=INK, va='center')

    # ---- b: luciferase (Fig. S15A, HEK293T)
    b_sc = B_H / S15_CROP[3]
    b_slot = (Wpt - S15_CROP[2] * b_sc - 2, 2.0)                     # S15 vector embedded here after saving
    letter_lab('b', b_slot[0] - 14, 2)

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
    json.dump(dict(reads=n, tf_bin=TF, sigma=SIGMA, whole_plasmid_specificity=spec, commit=commit),
              open(out / 'summary.json', 'w'), indent=1)
    for cell in REPS:
        print(cell)
        print(S[cell].to_string(index=False))
    print(json.dumps(spec, indent=1))


if __name__ == '__main__':
    main()

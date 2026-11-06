"""Grant Figure 2 (draft 14, half page): SLC39A4 c.-169A>G disrupts the NF-Y footprint and promoter activity
(fig2_rules.json, revision_draft5-14). Laid out at its final size (3.5 in wide, text >= 7 pt). Output is one
Illustrator-editable SVG (live Arial text), plus PDF and PNG renders (rsvg-convert with fonts.conf).

Inputs: the SLC39A4_variant_slide project's outputs (HEK293T plasmid Fiber-seq, WT and c.-169A>G co-transfected,
3 replicates) and the same outputs for HepG2 (4 replicates; fig2_hepg2.sh, identical processing):
analysis/occ_<rep>_<allele>.npy (reads with a footprint of size s covering plasmid position p), n_<rep>_<allele>.txt
(kept reads), reads.tsv.gz (per read: allele, kept, site_fp); motifs/MA0060.4.jaspar; and the HEK293T luciferase
points and box geometry read from the vector Fig. S15A by fig2_luc_digitize.py (points.tsv, box.tsv; the raw values
are not on klone).
TF-sized footprint = 20-60 bp. Cell lines are never pooled or overlaid.
a  WT - variant TF-sized occupancy (HEK293T, mean of the per-replicate tracks, Gaussian sigma SIGMA bp) over the
   promoter and luciferase (WIN), above the plasmid map (to scale, gene orientation from the promoter's 5' end); a
   dotted line joins the c.-169 peak to the variant tick; NF-Y logo as a zoom from the peak
b  luciferase (HEK293T; box and points as in Fig. S15A, p as stated there)
Per-replicate site tests (site_stats_<cell>.tsv) and the HepG2 whole-plasmid track are written for the legend only.
Usage: python tools/r01_figure/fig2.py --project <SLC39A4_variant_slide> --hepg2 <fig2/hepg2> --luc <dir> --output <dir>
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
FS = 7                                    # every label (pt); panel letters 9 pt
plt.rcParams.update({'font.family': 'Arimo', 'svg.fonttype': 'none', 'font.size': FS, 'axes.linewidth': 0.7,
                     'xtick.major.width': 0.7, 'ytick.major.width': 0.7, 'xtick.major.size': 2.5,
                     'ytick.major.size': 2.5, 'xtick.labelsize': FS, 'ytick.labelsize': FS,
                     'xtick.major.pad': 1.5, 'ytick.major.pad': 1.5})

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
W_IN, H_IN = 3.5, 1.62
WIN = (-100, 3460)                            # panel a window, bp from the promoter's 5' end (promoter + luciferase)


def gx(p):
    """Plasmid position -> bp from the promoter fragment's 5' end, gene orientation, wrapping round the circle."""
    return (PROM[1] - np.asarray(p)) % L


def feature_arrow(ax, x0, x1, y, h, color, label, tc, gid, lx=None, fs=FS):
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
    ap.add_argument('--luc', required=True, help='fig2_luc_digitize.py output dir (points.tsv, box.tsv)')
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
    pts = pd.read_csv(Path(a.luc) / 'points.tsv', sep='\t')
    box = pd.read_csv(Path(a.luc) / 'box.tsv', sep='\t').set_index('group')
    order = np.argsort(gx(np.arange(1, L + 1)))
    xg, yg = gx(np.arange(1, L + 1))[order], dmean[order]
    xg, yg = np.r_[xg - L, xg], np.r_[yg, yg]                      # circular: the window may start before 0

    fig = plt.figure(figsize=(W_IN, H_IN))
    Wpt, Hpt = W_IN * 72, H_IN * 72

    def axpt(x0, y0, x1, y1, **kw):
        return fig.add_axes([x0 / Wpt, 1 - y1 / Hpt, (x1 - x0) / Wpt, (y1 - y0) / Hpt], **kw)

    def letter_lab(ch, x, y):
        fig.text(x / Wpt, 1 - y / Hpt, ch, fontsize=9, fontweight='bold', va='top', ha='left')

    # ---- a: WT - variant track over the plasmid map, promoter + luciferase window (shared x, to scale)
    letter_lab('a', 0, 1)
    X0, X1, Y0, Y1 = 31.0, 168.0, 6.0, 78.0                        # track axes (pt)
    ylo, yhi = -1.5, 8.5
    tax = axpt(X0, Y0, X1, Y1)
    tax.axhline(0, color=INK2, lw=0.5, zorder=1)
    tax.plot(xg, yg, color=TF_C, lw=1.0, zorder=2)
    tax.set_xlim(*WIN)
    tax.set_ylim(ylo, yhi)
    tax.set_xticks([])
    tax.set_yticks([0, 4, 8])
    tax.spines[['top', 'right', 'bottom']].set_visible(False)
    tax.set_ylabel('Footprint loss (%)\n(WT − variant)', fontsize=FS, labelpad=1.5, linespacing=1.1)
    ov = axpt(0, 0, Wpt, Hpt)                                       # overlay in pt (map, logo, labels)
    ov.set_xlim(0, Wpt)
    ov.set_ylim(Hpt, 0)
    ov.axis('off')
    xp = lambda g: X0 + (np.asarray(g) - WIN[0]) / (WIN[1] - WIN[0]) * (X1 - X0)   # noqa: E731
    yp = lambda v: Y0 + (yhi - v) / (yhi - ylo) * (Y1 - Y0)         # noqa: E731
    yb, bh = Y1 + 10, 9                                             # map backbone y, feature height
    loop_rx, loop_ry = 5.0, 9.0                                     # plasmid loop ends (half-ellipses), not to scale
    t = np.linspace(0, np.pi, 60)
    lx = np.r_[X1, X0, X0 - loop_rx * np.sin(t), X1, X1 + loop_rx * np.sin(t)]
    ly = np.r_[yb, yb, yb + loop_ry - loop_ry * np.cos(t), yb + 2 * loop_ry, yb + loop_ry + loop_ry * np.cos(t)]
    ov.plot(lx, ly, color='#4A4A4A', lw=0.8, solid_joinstyle='round', zorder=1, gid='backbone')
    xv = xp(gx(VPOS1))
    feature_arrow(ov, xp(0), xp(gx(PROM[0])), yb, bh, PROM_C, 'SLC39A4', WT_C, 'feature-promoter')
    feature_arrow(ov, xp(gx(LUC[1])), xp(gx(LUC[0])), yb, bh, LUC_C, 'Luciferase', INK, 'feature-luciferase')
    ov.text((X0 + X1) / 2, yb + loop_ry + bh / 4, 'Reporter plasmid', fontsize=FS, color=INK2, ha='center',
            va='center')
    dots = dict(color=VAR_C, lw=0.5, ls=(0, (1.0, 1.2)))           # peak -> variant: behind the track, then to the map
    tax.axvline(gx(VPOS1), zorder=1.5, **dots)
    ov.plot([xv, xv], [Y1, yb - bh / 2], zorder=0.5, **dots)
    # NF-Y logo, zoomed from the c.-169 peak
    cnt = np.array([[float(x) for x in ln.split('[')[1].split(']')[0].split()]
                    for ln in open(P / 'motifs/MA0060.4.jaspar').read().splitlines()[1:5]])
    pr = (cnt + 0.8) / (cnt.sum(0) + 3.2)
    ic = 2 + (pr * np.log2(pr)).sum(0)
    cw, bit = 7.5, 8.5                                              # column width, pt per bit
    apex = (xv, yp(dmean[VPOS1 - 1]))
    lx0, ly0 = xv + 20, Y0 + 40                                      # logo left edge, baseline
    ov.add_patch(Polygon([apex, [lx0 - 2, ly0 - 2 * bit - 2], [lx0 - 2, ly0 + 1.5]], fc='#E9E9F4', ec='none',
                         zorder=0))
    bold = FontProperties(family='Arimo', weight='bold')
    for j in range(8):
        y = ly0
        for k in np.argsort(pr[:, j] * ic[j]):
            h = pr[k, j] * ic[j] * bit
            if h > 0.5:
                tp = TextPath((0, 0), 'ACGT'[k], size=1, prop=bold)
                bb = tp.get_extents()
                tr = Affine2D().translate(-bb.x0, -bb.y0).scale((cw - 1.2) / bb.width, -h / bb.height).translate(
                    lx0 + j * cw + 0.6, y)
                ov.add_patch(PathPatch(tr.transform_path(tp), fc=LETTER_C['ACGT'[k]], ec='none', zorder=2))
            y -= h
    ov.plot([lx0 - 1, lx0 + 8 * cw + 1], [ly0 + 1, ly0 + 1], color=INK2, lw=0.5)
    vj = 3                                                          # CCA[A]TCAC: variant = motif position 4
    ov.add_patch(Rectangle((lx0 + vj * cw, ly0 - 2 * bit - 1.5), cw, 2 * bit + 4, fc='none', ec=VAR_C, lw=0.7,
                           ls=(0, (2.2, 1.1)), zorder=3))
    ov.text(lx0 + (vj + 0.5) * cw, ly0 + 3.5, 'c.-169A>G', fontsize=FS, fontweight='bold', color=VAR_C,
            ha='center', va='top')
    ov.text(lx0 + 4 * cw, ly0 - 2 * bit - 3, 'NF-Y motif', fontsize=FS, fontweight='bold', color=INK, ha='center',
            va='bottom')

    # ---- b: luciferase (HEK293T; points and box geometry from the vector Fig. S15A)
    letter_lab('b', 176, 1)
    bax = axpt(206, Y0 + 4, 250, 88)
    for i, (g, col) in enumerate((('WT', WT_C), ('Var', VAR_C))):
        b = box.loc[g]
        bw = 0.56
        bax.add_patch(Rectangle((i - bw / 2, b.q1), bw, b.q3 - b.q1, fc=col, alpha=0.35, ec=col, lw=0.7, zorder=1))
        bax.add_patch(Rectangle((i - bw / 2, b.q1), bw, b.q3 - b.q1, fc='none', ec=col, lw=0.7, zorder=2))
        bax.plot([i - bw / 2, i + bw / 2], [b['median']] * 2, color=INK, lw=1.2, zorder=3)
        for lo, hi in ((b.whisker_lo, b.q1), (b.q3, b.whisker_hi)):
            bax.plot([i, i], [lo, hi], color=col, lw=0.7, zorder=1)
        for yw in (b.whisker_lo, b.whisker_hi):
            bax.plot([i - bw / 4, i + bw / 4], [yw, yw], color=col, lw=0.7, zorder=1)
        g_pts = pts[pts.group == g]
        bax.scatter(i + g_pts.x_offset * bw, g_pts.value, s=7, color=col, edgecolor='white', lw=0.3, zorder=4)
    yt = box.whisker_hi.max() + 0.25                                # p bracket (as in S15A: 'p < 0.01')
    bax.plot([0, 0, 1, 1], [box.loc['WT', 'whisker_hi'] + 0.08, yt, yt, box.loc['Var', 'whisker_hi'] + 0.08],
             color=INK2, lw=0.6)
    bax.text(0.5, yt + 0.04, 'p < 0.01', fontsize=FS, color=INK2, ha='center', va='bottom')
    bax.set_xlim(-0.6, 1.6)
    bax.set_ylim(0, 2.15)
    bax.set_yticks([0, 1, 2])
    bax.set_xticks([0, 1], ['WT', 'c.-169A>G'], rotation=35, ha='right', rotation_mode='anchor')
    bax.get_xticklabels()[0].set_color(WT_C)
    bax.get_xticklabels()[1].set_color(VAR_C)
    bax.spines[['top', 'right']].set_visible(False)
    bax.set_ylabel('Luciferase activity\n(firefly/renilla)', fontsize=FS, labelpad=1.5, linespacing=1.1)

    svg = out / 'fig2.svg'
    fig.savefig(svg, facecolor='white', dpi=300)
    plt.close(fig)
    s = expand_fonts(svg.read_text().replace("'Arimo'", 'Arial, Helvetica, sans-serif').replace('Arimo', 'Arial'))[0]
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

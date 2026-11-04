"""Grant Figure 1 (DPAGT1 draft 7; WASF1 draft 1): plasmid Fiber-seq workflow, motif-footprinted 300-bp construct
molecules along the whole plasmid, and occupancy tracks (fig1_rules.json: fp_bins_v2, panel_a_ctcf; wasf1.draft1 for
--gene WASF1: WASF12 construct, SP/KLF GC box, FIRE from fig1_fire_profiles.py). Output is one Illustrator-editable
SVG (live Arial text; workflow icons are the HPRC main-figure vectors), plus PDF and PNG renders of it (rsvg-convert
with fonts.conf: Arial -> Arimo, metric-compatible, so the renders match the SVG layout).

Footprint bins (fp_bins_v2, chosen from fig1_fpsizes.py): 15-30, 31-60, 61-89 bp, nucleosome >= 90 bp; calls < 15 bp
are not drawn or counted.
a  workflow: fig1_schematic.py's schematic.svg embedded as vector, with the shared colour legend beside it; then the
   linearized 300-bp construct (4,774 bp: vector backbone, FIRE peak insert, motif; DPAGT11 + CTCF or WASF12 + GC box)
   and N_SHOW of its molecules that carry a 31-60 bp call over the motif (>= MIN_OVERLAP bp), drawn at
   random (seed) and ordered by Ward clustering of 31-60 bp occupancy over the motif +/- CLUSTER_FLANK (leaf order).
   Nucleosome and 31-60 bp calls only, along the whole plasmid (no m6A marks, no FIRE).
b  gene model(s) (with the motif) and occupancy tracks over X0-X1: endogenous Fiber-seq shows nucleosomes and FIRE
   only; the plasmid 300-bp track shows all bins, sub-nucleosome bins x FP_SCALE. Footprint shares from this project's
   FiberHMM calls (molecules with a call of the bin covering a base / molecules covering it), FIRE share from the LRSS
   slide 5/6 table (DPAGT1) or fig1_fire_profiles.py (WASF1); Gaussian sigma SIGMA bp within covered runs; solid fills
   drawn tallest first.
Usage: python tools/r01_figure/fig1.py --work results/r01_figure_20261102 --output <dir> [--gene DPAGT1|WASF1]
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
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from scipy.cluster.hierarchy import leaves_list, linkage  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402

HPRC = Path('/gscratch/stergachislab/bmallo/large_home/projects/HPRC_driver_variants')
LRSS = Path('/gscratch/stergachislab/bmallo/large_home/resources/presentations/lrss2026_slides56')
RSVG = '/mmfs1/home/bmallo/miniconda3/envs/fiberseq-qc/bin/rsvg-convert'
FONTS_CONF = Path(__file__).with_name('fonts.conf')
sys.path.insert(0, str(HPRC / 'scripts/plasmid_background'))
from svg_fonts import expand_fonts  # noqa: E402
for f in ('Arimo-Regular.ttf', 'Arimo-Bold.ttf', 'Arimo-Italic.ttf', 'Arimo-BoldItalic.ttf'):
    p_ = HPRC / 'scripts/fonts' / f
    if p_.exists():
        font_manager.fontManager.addfont(str(p_))
plt.rcParams.update({'font.family': 'Arimo', 'svg.fonttype': 'none', 'font.size': 6})

INK, INK2 = '#2E2E2E', '#555555'
COL = {'nuc': '#C7C8CA', 'FIRE': '#BA7C69', 'fp_61_89': '#E77C34', 'fp_31_60': '#338B8C', 'fp_15_30': '#755B9D'}
NAMES = {'nuc': 'Nucleosome (≥90 bp)', 'FIRE': 'FIRE', 'fp_61_89': '61–89 bp footprint',
         'fp_31_60': '31–60 bp footprint ({})', 'fp_15_30': '15–30 bp footprint'}
BINS = {'fp_15_30': (15, 30), 'fp_31_60': (31, 60), 'fp_61_89': (61, 89), 'nuc': (90, 10 ** 6)}
CODE = {0: 'accessible', 1: 'fp_15_30', 2: 'fp_31_60', 3: 'fp_61_89', 4: 'nuc'}
C_GENE, HILITE, C_MOTIF, C_BACKBONE = '#1f6fb4', '#FBE3D3', '#B8860B', '#8B9094'
GENES = {
    # DPAGT1 gene body in the 4.1-kb construct reads (insert start .. TSS); FIRE peak 1 300-bp element (DPAGT11);
    # HOCOMOCO v14 CTCF.H14CORE.0.P.B, P 3.8e-6 (fig1/motifs.tsv)
    'DPAGT1': dict(work='fig1', window=(119100224, 119101853), peak=(119100271, 119100570),
                   motif=(119100436, 119100455, '-'), motif_name='CTCF motif', fp_name='CTCF-sized', read_src='DPAGT11',
                   genes=[('NM_001382.4', 'DPAGT1')], refseq='raw/ucsc_ncbiRefSeq_dpagt1.json',
                   fire=LRSS / 'tables/profiles_DPAGT1.tsv.gz'),
    # WASF1 exons 1-2 .. CDC40 exon 1; FIRE peak 4 300-bp element (WASF12); SP/KLF GC box, union of the best SP1/SP4/KLF
    # matches (HOCOMOCO v14 KLF10 P 7e-10, SP1 P 1.3e-8; fig1_wasf1/motifs/motifs.tsv)
    'WASF1': dict(work='fig1_wasf1', window=(110178550, 110180700), peak=(110179911, 110180211),
                  motif=(110180054, 110180097, '.'), motif_name='GC box (SP/KLF)', fp_name='TF-sized',
                  read_src='WASF12',
                  genes=[('NM_003931.3', 'WASF1'), ('NM_015891.3', 'CDC40')], refseq='raw/ucsc_ncbiRefSeq_wasf1.json',
                  fire='fire/fire_profiles.tsv.gz'),
}
X0 = X1 = PEAK1 = MOTIF = None          # set from GENES in main(): window, 300-bp element, motif
SEED, N_SHOW, SPAN_SLACK, MIN_OVERLAP, CLUSTER_FLANK = 20261102, 15, 20, 10, 60
SIGMA, FP_SCALE = 3, 5
W_IN, H_IN = 7.2, 5.45
CROP = (58.63 + 1, 46.63 + 82 + 1.5, 185.0, 54)   # schematic crop (fig1_schematic.CROP)


def fp_code(ln):
    return next((k for k, c in CODE.items() if k and BINS[c][0] <= ln <= BINS[c][1]), None)


def select_ctcf_reads(mol):
    """Molecules spanning the motif (CTCF for DPAGT1) with a motif-sized (31-60 bp) call overlapping it by
    >= MIN_OVERLAP bp (fig1_ctcf_reads.py); returns their indices, all footprint lists, and per-base 31-60 bp
    occupancy over the motif
    +/- CLUSTER_FLANK (clustering features)."""
    lo, hi = MOTIF[0] - CLUSTER_FLANK, MOTIF[1] + CLUSTER_FLANK
    keep, fps_all, feats = [], [], []
    for i, r in enumerate(mol.itertuples()):
        f = json.loads(r.fp)
        fps_all.append(f)
        at = np.array(json.loads(r.at), int)
        if not at.size or at.min() > MOTIF[0] - SPAN_SLACK or at.max() < MOTIF[1] + SPAN_SLACK:
            continue
        hit = [(s, e) for s, e, ln in f if fp_code(ln) == 2 and min(e, MOTIF[1]) - max(s, MOTIF[0]) + 1 >= MIN_OVERLAP]
        if hit:
            v = np.zeros(hi - lo + 1)
            for s, e, ln in f:
                if fp_code(ln) == 2 and e >= lo and s <= hi:
                    v[max(s, lo) - lo:min(e, hi) - lo + 1] = 1
            keep.append(i)
            feats.append(v)
    return keep, fps_all, np.array(feats)


def track_shares(mol, lo, hi):
    """Per base X0..X1: molecules covering it (span of observed A/T, limited to [lo, hi]) and, per bin, molecules with
    a call of that bin covering it; shares smoothed with a Gaussian inside covered runs."""
    L = X1 - X0 + 1
    cov = np.zeros(L)
    cnt = {c: np.zeros(L) for c in BINS}
    for r in mol.itertuples():
        at = np.array(json.loads(r.at), int)
        if not at.size:
            continue
        a, b = max(at.min(), lo, X0), min(at.max(), hi, X1)
        if b < a:
            continue
        cov[a - X0:b - X0 + 1] += 1
        for c in BINS:
            m = np.zeros(L, bool)
            for s, e, ln in json.loads(r.fp):
                if BINS[c][0] <= ln <= BINS[c][1]:
                    s_, e_ = max(s, a), min(e, b)
                    if e_ >= s_:
                        m[s_ - X0:e_ - X0 + 1] = True
            cnt[c] += m
    out = {'pos': np.arange(X0, X1 + 1), 'molecules': cov}
    for c in BINS:
        out[f'{c}_share'] = np.where(cov > 0, cnt[c] / np.maximum(cov, 1), np.nan)
        out[f'{c}_smooth'] = smooth(out[f'{c}_share'])
    return pd.DataFrame(out)


def smooth(raw):
    out = np.full(raw.size, np.nan)
    ok = ~np.isnan(raw)
    edges = np.flatnonzero(np.diff(np.r_[0, ok.astype(int), 0]))
    for a, b in zip(edges[::2], edges[1::2]):
        out[a:b] = gaussian_filter1d(raw[a:b], SIGMA, mode='nearest')
    return out


def gene_model(ax, ref, genes, y=0.4, h=0.44, size=7, tss_h=0.48, arrow=6):
    """Each (transcript, label) clipped to X0-X1, intron chevrons by strand, a bent TSS arrow; the label sits on the
    side where the gene runs out of the window."""
    for name, label in genes:
        r = next(r for r in ref if r['name'] == name)
        a, b = max(r['txStart'], X0), min(r['txEnd'], X1)
        d = -1 if r['strand'] == '-' else 1
        ax.plot([a, b], [y, y], color=C_GENE, lw=0.8)
        for xx in np.arange(a + 60, b - 40, 120):
            ax.text(xx, y, '‹' if d < 0 else '›', ha='center', va='center', fontsize=size - 1.5, color=C_GENE)
        for es, ee in zip(map(int, r['exonStarts'].strip(',').split(',')),
                          map(int, r['exonEnds'].strip(',').split(','))):
            if ee > X0 and es < X1:
                ax.add_patch(Rectangle((max(es, X0), y - h / 2), min(ee, X1) - max(es, X0), h, fc=C_GENE, lw=0))
        tss = r['txEnd'] if d < 0 else r['txStart']
        ax.plot([tss, tss, tss + d * 60], [y, y + tss_h, y + tss_h], color=C_GENE, lw=1.0)
        ax.annotate('', xy=(tss + d * 95, y + tss_h), xytext=(tss + d * 55, y + tss_h),
                    arrowprops=dict(arrowstyle='-|>', color=C_GENE, lw=1.0, shrinkA=0, shrinkB=0, mutation_scale=arrow))
        left = a == X0
        ax.text(X0 - 20 if left else X1 + 20, y, label, ha='right' if left else 'left', va='center', fontsize=size,
                color=C_GENE, fontweight='bold', style='italic', clip_on=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--work', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--gene', default='DPAGT1', choices=sorted(GENES))
    a = ap.parse_args()
    global X0, X1, PEAK1, MOTIF
    G = GENES[a.gene]
    (X0, X1), PEAK1, MOTIF, READ_SRC = G['window'], G['peak'], G['motif'], G['read_src']
    NAMES['fp_31_60'] = NAMES['fp_31_60'].format(G['fp_name'])
    w, out = Path(a.work), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    f1 = w / G['work']
    lrss = pd.read_csv(f1 / G['fire'], sep='\t')
    ref = json.load(open(f1 / G['refseq']))['ncbiRefSeq']
    ins = pd.read_csv(f1 / 'inserts.tsv', sep='\t').set_index('construct').loc[READ_SRC]
    assert ins.strand == '+' and ins.hg38_start == PEAK1[0]
    off = int(ins.hg38_start - ins.construct_start)                 # plasmid coordinate (1-based) = hg38 - off
    plen = int(ins.length)
    mol = pd.read_parquet(f1 / f'profiles/{READ_SRC}.molecules.parquet')
    endo = pd.read_parquet(f1 / 'profiles/endogenous.molecules.parquet')
    keep, fps_all, feats = select_ctcf_reads(mol)
    rng = np.random.default_rng(SEED)
    pick = np.sort(rng.choice(len(keep), N_SHOW, replace=False))
    order = [keep[pick[j]] for j in leaves_list(linkage(feats[pick], 'ward'))]
    pd.DataFrame({'row': range(1, N_SHOW + 1), 'read': mol.read.iloc[order].to_numpy()}).to_csv(
        out / 'panel_a_reads.tsv', sep='\t', index=False)
    tracks = {'genome': track_shares(endo, X0, X1), '300bp': track_shares(mol, *PEAK1)}
    for k, t in tracks.items():
        lr = lrss[lrss.track == k].set_index('pos')
        t['FIRE_share'] = lr.FIRE_share.reindex(t.pos).to_numpy()
        t['FIRE_smooth'] = smooth(t.FIRE_share.to_numpy())
    pd.concat([t.assign(track=k) for k, t in tracks.items()]).to_csv(out / 'tracks.tsv', sep='\t', index=False,
                                                                      float_format='%.5g')
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()

    fig = plt.figure(figsize=(W_IN, H_IN))
    Wpt, Hpt = W_IN * 72, H_IN * 72

    def axpt(x0, y0, x1, y1):            # axes from a box in points, origin top left
        return fig.add_axes([x0 / Wpt, 1 - y1 / Hpt, (x1 - x0) / Wpt, (y1 - y0) / Hpt])

    ov = fig.add_axes([0, 0, 1, 1], facecolor='none')
    ov.set_xlim(0, Wpt)
    ov.set_ylim(Hpt, 0)
    ov.axis('off')

    # ---- a: schematic slot + legend, then the linear plasmid map and the clustered molecules
    sch_h = 86.0
    sc = sch_h / CROP[3]
    sch_w = CROP[2] * sc
    slot = (16.0, 12.0)
    fig.text(4 / Wpt, 1 - 6 / Hpt, 'a', fontsize=10, fontweight='bold', va='top')
    fig.legend([Patch(color=COL[c]) for c in COL], [NAMES[c] for c in COL], loc='upper left',
               bbox_to_anchor=((slot[0] + sch_w + 28) / Wpt, 1 - 26 / Hpt), ncol=1, frameon=False, fontsize=6,
               handlelength=1.3, handleheight=0.9, labelspacing=0.55, title='Fiber-seq calls',
               title_fontproperties={'weight': 'bold', 'size': 6.3}, alignment='left')
    rx0, rx1 = 40.0, Wpt - 14
    my0, my1, ry0, ry1 = 104.0, 124.0, 128.0, 214.0
    axm = axpt(rx0, my0, rx1, my1)
    axm.set_xlim(0.5, plen + 0.5)
    axm.set_ylim(0, 1)
    axm.axis('off')
    i0, i1 = int(ins.construct_start), int(ins.construct_end)
    c0, c1 = MOTIF[0] - off, MOTIF[1] - off
    axm.add_patch(Rectangle((0.5, 0.22), plen, 0.16, fc=C_BACKBONE, lw=0, gid='backbone'))
    axm.add_patch(Rectangle((i0 - 0.5, 0.12), i1 - i0 + 1, 0.36, fc=C_GENE, lw=0, gid='insert'))
    axm.add_patch(Rectangle((c0 - 0.5, 0.56), c1 - c0 + 1, 0.3, fc=C_MOTIF, lw=0, gid='ctcf-motif-map'))
    axm.text(c1 + 25, 0.71, G['motif_name'], ha='left', va='center', fontsize=5.6, color=INK, fontweight='bold')
    axm.text(i0 - 30, 0.71, f'{a.gene} 300-bp element', ha='right', va='center', fontsize=5.6, color=C_GENE,
             fontweight='bold')
    axm.text((i1 + plen) / 2, 0.71, f'vector backbone ({plen:,} bp plasmid, linearized)', ha='center', va='center',
             fontsize=5.6, color=INK2)
    for xx, ha in ((0.5, 'left'), (plen + 0.5, 'right')):
        axm.text(xx, -0.05, 'NotI', ha=ha, va='top', fontsize=5, color=INK2)
    axr = axpt(rx0, ry0, rx1, ry1)
    axr.axvspan(i0 - 0.5, i1 + 0.5, color=HILITE, lw=0, zorder=0)
    axr.axvspan(c0 - 0.5, c1 + 0.5, color=C_MOTIF, alpha=0.25, lw=0, zorder=0)
    for y, i in enumerate(reversed(order)):
        at = np.array(json.loads(mol['at'].iloc[i]), int) - off
        axr.plot([max(at.min(), 1), min(at.max(), plen)], [y, y], color='#BDBDBD', lw=0.4, zorder=1)
        for code in (4, 2):                                         # nucleosomes and CTCF-sized (31-60 bp) calls only
            for s, e, ln in fps_all[i]:
                if fp_code(ln) == code:
                    axr.add_patch(Rectangle((s - off - 0.5, y - 0.36), e - s + 1, 0.72, fc=COL[CODE[code]], lw=0,
                                            zorder=2))
    axr.set_xlim(0.5, plen + 0.5)
    axr.set_ylim(-0.6, N_SHOW - 0.4)
    axr.axis('off')
    axr.plot([plen - 500, plen], [-1.4, -1.4], color=INK, lw=0.8, clip_on=False)
    axr.text(plen - 250, -1.9, '500 bp', ha='center', va='top', fontsize=5.5, color=INK)
    axr.text(0.5, -1.4, f'300-bp construct molecules with a 31–60 bp footprint on the {G["motif_name"]} '
             f'({N_SHOW} of {len(keep):,} such molecules; {len(mol):,} total), clustered on that footprint',
             ha='left', va='top', fontsize=5.5, color=INK2)

    # ---- b: gene model (with the CTCF motif) + occupancy tracks
    bx0, bx1 = 112.0, Wpt - 52
    t1, t2 = (270, 316), (322, 368)
    fig.text(4 / Wpt, 1 - 232 / Hpt, 'b', fontsize=10, fontweight='bold', va='top')
    axg = axpt(bx0, 236, bx1, 264)
    axg.set_xlim(X0, X1)
    axg.set_ylim(0, 1.05)
    axg.axis('off')
    gene_model(axg, ref, G['genes'])
    axg.add_patch(Rectangle((MOTIF[0] - 0.5, 0.74), MOTIF[1] - MOTIF[0] + 1, 0.2, fc=C_MOTIF, lw=0, gid='ctcf-motif-b'))
    axg.text(MOTIF[1] + 12, 0.84, G['motif_name'], ha='left', va='center', fontsize=5.6, color=INK, fontweight='bold')
    axg.plot([X0 + 700, X0 + 900], [0.9, 0.9], color=INK, lw=0.8)
    axg.text(X0 + 800, 0.96, '200 bp', ha='center', va='bottom', fontsize=5.5, color=INK)
    axes = []
    for track, cats, scale, (y0, y1) in (('genome', ['nuc', 'FIRE'], 1, t1),
                                         ('300bp', list(COL), FP_SCALE, t2)):
        ax = axpt(bx0, y0, bx1, y1)
        axes.append(ax)
        t = tracks[track]
        x = t.pos.to_numpy() + 0.5
        layers = {c: np.minimum(t[f'{c}_smooth'].to_numpy() * (scale if c.startswith('fp_') else 1), 1.0)
                  for c in cats}
        for cat in sorted(layers, key=lambda c: -np.nanmax(layers[c])):   # tallest behind, shortest in front
            yv = layers[cat]
            ax.fill_between(x, 0, np.nan_to_num(yv), where=~np.isnan(yv), color=COL[cat], lw=0,
                            gid=f'{track}-{cat}-fill')
        ax.axvspan(MOTIF[0] - 0.5, MOTIF[1] + 0.5, color=C_MOTIF, alpha=0.3, lw=0, zorder=3, gid=f'{track}-ctcf')
        if track == '300bp':
            for xx in PEAK1:
                ax.axvline(xx, color=INK2, lw=0.6, ls=(0, (3, 2)))
        ax.set_xlim(X0, X1)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.yaxis.tick_right()
        ax.set_yticks([0, 1])
        ax.tick_params(axis='y', labelsize=6.5, length=2, pad=1.5)
        for lab in ax.get_yticklabels():
            lab.set_fontweight('bold')
        for s in ax.spines.values():
            s.set_linewidth(0.9)
    axes[0].text(-0.015, 0.6, 'Fiber-seq', transform=axes[0].transAxes, ha='right', va='center', fontsize=7,
                 fontweight='bold', color=INK)
    axes[0].text(-0.015, 0.36, '(genome-wide)', transform=axes[0].transAxes, ha='right', va='center', fontsize=6.3,
                 style='italic', color=INK)
    axes[1].text(-0.015, 0.5, '300bp\nregion', transform=axes[1].transAxes, ha='right', va='center', fontsize=7,
                 fontweight='bold', color=INK, multialignment='right', linespacing=1.15)
    bk = bx0 - 40
    ov.plot([bk, bk], list(t2), color=INK, lw=0.9)
    for yy in t2:
        ov.plot([bk, bk + 4], [yy, yy], color=INK, lw=0.9)
    ov.text(bk - 8, sum(t2) / 2, 'Plasmid\nFiber-seq', rotation=90, ha='center', va='center', fontsize=6.5,
            fontweight='bold', color=INK, linespacing=1.3)
    ov.text(bx1 + 22, (t1[0] + t2[1]) / 2, 'Proportion of Reads Bound', rotation=270, ha='center', va='center',
            fontsize=7, fontweight='bold', color=INK, rotation_mode='anchor')
    ov.text((bx0 + bx1) / 2, t2[1] + 6, f'Endogenous track: nucleosomes and FIRE only. Plasmid track: sub-nucleosome '
            f'footprint shares ×{FP_SCALE}. Calls <15 bp not shown.', ha='center', va='top', fontsize=5.4, color=INK2)
    fig.text(0.995, 0.997, f'draft {date.today()} · git {commit}+uncommitted', ha='right', va='top', fontsize=4.5,
             color='#9a9893', gid='draft-stamp')

    # ---- write SVG, then embed the schematic vector in its slot
    svg = out / 'fig1.svg'
    fig.savefig(svg, facecolor='white')
    plt.close(fig)
    s = expand_fonts(svg.read_text().replace("'Arimo'", 'Arial, Helvetica, sans-serif').replace('Arimo', 'Arial'))[0]
    sch = (w / 'fig1/schematic2/schematic.svg').read_text()
    body = re.sub(r'^.*?<svg[^>]*>', '', sch, count=1, flags=re.S).rsplit('</svg>', 1)[0]
    body = re.sub(r'<metadata>.*?</metadata>', '', body, flags=re.S)
    body = re.sub(r'id="([^"]+)"', r'id="sch-\1"', body)
    body = re.sub(r'url\(#([^)]+)\)', r'url(#sch-\1)', body)
    body = re.sub(r'href="#([^"]+)"', r'href="#sch-\1"', body)
    cx, cy, cw, ch = CROP
    group = (f'<g id="panel-a-schematic" transform="translate({slot[0]:.3f} {slot[1]:.3f}) scale({sc:.5f}) '
             f'translate({-cx:.3f} {-cy:.3f})"><defs><clipPath id="sch-crop"><rect x="{cx}" y="{cy}" width="{cw}" '
             f'height="{ch}"/></clipPath></defs><g clip-path="url(#sch-crop)">{body}</g></g>')
    s = s.replace('</svg>', group + '\n</svg>')
    if 'xmlns:inkscape' not in s.split('>', 2)[1]:      # the HPRC layers carry inkscape:label attributes
        s = re.sub(r'<svg ', '<svg xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" ', s, count=1)
    svg.write_text(s)
    env = dict(os.environ, FONTCONFIG_FILE=str(FONTS_CONF))
    subprocess.run([RSVG, '-f', 'pdf', '-o', str(out / 'fig1.pdf'), str(svg)], check=True, env=env)
    subprocess.run([RSVG, '-d', '300', '-p', '300', '-f', 'png', '-o', str(out / 'fig1.png'), str(svg)], check=True,
                   env=env)
    json.dump(dict(gene=a.gene, motif_name=G['motif_name'], read_source=READ_SRC, molecules=len(mol),
                   ctcf_footprint_molecules=len(keep), shown=N_SHOW,
                   seed=SEED, plasmid_length=plen, insert=[i0, i1], ctcf_plasmid=[c0, c1],
                   endogenous_molecules=len(endo), endogenous_max_coverage=int(tracks['genome'].molecules.max()),
                   bins=BINS, sigma=SIGMA, fp_scale_300bp=FP_SCALE, window=[X0, X1], peak1=list(PEAK1),
                   ctcf=list(MOTIF), commit=commit), open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

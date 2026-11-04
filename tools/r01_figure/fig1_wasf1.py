"""Grant Figure 1, WASF1 draft 5 (fig1_rules.json: wasf1.draft3, draft4, draft5): plasmid Fiber-seq workflow with 4.5-kb construct
molecules beside it (footprints on the GC box and at the 20-30 bp site), and occupancy tracks over the WASF1 gene body,
FIRE peak 4 and its downstream nucleosome. Same style as fig1.py (DPAGT1 draft 7); one Illustrator-editable SVG (live Arial text; workflow icons are
the HPRC main-figure vectors) plus PDF and PNG renders (rsvg-convert with fonts.conf).

Footprint bins (wasf1.draft3, user-chosen after fig1_fpsizes.py --gene WASF1): 20-30, 31-50, 51-80 bp, nucleosome
>= 81 bp (as Figure 3); calls < 20 bp are not drawn or counted.
a  workflow schematic; to its right a map (300-bp element, GC box, 20-30 bp site) and N_SHOW 4.5-kb construct molecules
   with a 20-80 bp call over the GC box (>= MIN_OVERLAP bp) and a separate 20-30 bp call with midpoint within
   SMALL_SITE +/- SMALL_WIN, drawn at random (seed), ordered by Ward clustering of one-hot per-base call class over
   the 300-bp element; every bin and nucleosomes drawn over X0-READ_X1 (hg38, as b; cropped before the 20-30 bp
   calls past the downstream nucleosome), with FIRE elements (fig1_fire.sh on the 4.5-kb reads) drawn between the
   nucleosomes and the footprints, as in b; arrow from the sequencer to the reads; legend row below. No draft stamp (commit and date in summary.json).
b  WASF1 gene model with the GC box and occupancy tracks over X0-X1: endogenous Fiber-seq nucleosomes and FIRE; plasmid
   300-bp all bins x --fp-scale. Shares = molecules with a call of the bin covering a base / molecules covering it,
   Gaussian sigma SIGMA bp within covered runs; FIRE from fig1_fire_profiles.py. Draw order: nucleosome, FIRE (alpha
   FIRE_ALPHA), then footprint bins tallest first.
Usage: python tools/r01_figure/fig1_wasf1.py --work results/r01_figure_20261102 --output <dir> [--fp-scale 1]
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
BINS = {'fp_20_30': (20, 30), 'fp_31_50': (31, 50), 'fp_51_80': (51, 80), 'nuc': (81, 10 ** 6)}
COL = {'nuc': '#C7C8CA', 'FIRE': '#BA7C69', 'fp_51_80': '#E77C34', 'fp_31_50': '#338B8C', 'fp_20_30': '#755B9D'}
NAMES = {'nuc': 'Nucleosome (≥81 bp)', 'FIRE': 'FIRE', 'fp_51_80': '51–80 bp footprint',
         'fp_31_50': '31–50 bp footprint', 'fp_20_30': '20–30 bp footprint'}
CODE = {0: 'none', 1: 'fp_20_30', 2: 'fp_31_50', 3: 'fp_51_80', 4: 'nuc'}
C_GENE, HILITE, C_MOTIF, C_BACKBONE = '#1f6fb4', '#FBE3D3', '#B8860B', '#8B9094'
X0, X1 = 110178550, 110180380          # WASF1 exons 1-2 .. peak 4 downstream nucleosome (wasf1.draft2 window_b)
READ_X1 = 110180358                    # reads in a end here (wasf1.draft4 read_window)
PEAK = (110179911, 110180211)          # FIRE peak 4 300-bp element (WASF12 insert)
MOTIF = (110180054, 110180097)         # SP/KLF GC box (wasf1.draft1 motif)
SMALL_SITE, SMALL_WIN = 110180148, 10   # 20-30 bp call-midpoint density peak in WASF12 (wasf1.draft3 small_peak)
MOTIF_NAME, GENE_TX = 'GC box (SP/KLF)', 'NM_003931.3'
READ_SRC, PLASMID_SRC = 'WASF1_full', 'WASF12'   # reads in a: 4.5-kb construct; plasmid track in b: 300-bp construct
SEED, N_SHOW, SPAN_SLACK, MIN_OVERLAP = 20261102, 15, 20, 10
SIGMA, FIRE_ALPHA = 3, 0.65
W_IN, H_IN = 7.2, 4.15
CROP = (58.63 + 1, 46.63 + 82 + 1.5, 185.0, 54)   # schematic crop (fig1_schematic.CROP)


def fp_code(ln):
    return next((k for k, c in CODE.items() if k and BINS[c][0] <= ln <= BINS[c][1]), None)


def select_reads(mol, lo, hi):
    """Molecules spanning the motif with a 20-80 bp call overlapping it by >= MIN_OVERLAP bp and a separate 20-30 bp
    call with midpoint within SMALL_SITE +/- SMALL_WIN; returns their indices, all footprint lists, and one-hot
    per-base call class over [lo, hi] (clustering features)."""
    keep, fps_all, feats = [], [], []
    for i, r in enumerate(mol.itertuples()):
        f = json.loads(r.fp)
        fps_all.append(f)
        at = np.array(json.loads(r.at), int)
        if not at.size or at.min() > MOTIF[0] - SPAN_SLACK or at.max() < MOTIF[1] + SPAN_SLACK:
            continue
        on_motif = [min(e, MOTIF[1]) - max(s, MOTIF[0]) + 1 >= MIN_OVERLAP for s, e, ln in f]
        if any(fp_code(ln) in (1, 2, 3) and om for (s, e, ln), om in zip(f, on_motif)) and \
                any(fp_code(ln) == 1 and not om and abs((s + e) / 2 - SMALL_SITE) <= SMALL_WIN
                    for (s, e, ln), om in zip(f, on_motif)):
            v = np.zeros(hi - lo + 1, int)
            for s, e, ln in f:
                c = fp_code(ln)
                if c and e >= lo and s <= hi:
                    v[max(s, lo) - lo:min(e, hi) - lo + 1] = c
            keep.append(i)
            feats.append(np.eye(5)[v].ravel())
    return keep, fps_all, np.array(feats)


def fire_elements(path, contig, off, reads):
    """{read: [(hg38 start, end), ...]} FIRE elements of the given reads from ft extract --all (0-based, half-open
    reference coordinates on the construct contig; hg38 = position + off + 1)."""
    x = pd.read_csv(path, sep='\t', usecols=['#ct', 'fiber', 'ref_fire_starts', 'ref_fire_lengths'])
    x = x[(x['#ct'] == contig) & x.fiber.isin(reads)].drop_duplicates('fiber')
    out = {}
    for r in x.itertuples():
        st = [int(v) for v in str(r.ref_fire_starts).split(',') if v not in ('', '.', 'nan')]
        ln = [int(v) for v in str(r.ref_fire_lengths).split(',') if v not in ('', '.', 'nan')]
        out[r.fiber] = [(s + off + 1, s + off + l_) for s, l_ in zip(st, ln) if s >= 0 and l_ > 0]
    return out


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


def gene_model(ax, ref, y=0.4, h=0.44, size=7, tss_h=0.48, arrow=6):
    r = next(r for r in ref if r['name'] == GENE_TX)
    tss = r['txEnd']
    ax.plot([X0, tss], [y, y], color=C_GENE, lw=0.8)
    for xx in np.arange(X0 + 60, tss - 40, 120):
        ax.text(xx, y, '‹', ha='center', va='center', fontsize=size - 1.5, color=C_GENE)
    for es, ee in zip(map(int, r['exonStarts'].strip(',').split(',')), map(int, r['exonEnds'].strip(',').split(','))):
        if ee > X0 and es < X1:
            ax.add_patch(Rectangle((max(es, X0), y - h / 2), min(ee, X1) - max(es, X0), h, fc=C_GENE, lw=0))
    ax.plot([tss, tss, tss - 60], [y, y + tss_h, y + tss_h], color=C_GENE, lw=1.0)
    ax.annotate('', xy=(tss - 95, y + tss_h), xytext=(tss - 55, y + tss_h),
                arrowprops=dict(arrowstyle='-|>', color=C_GENE, lw=1.0, shrinkA=0, shrinkB=0, mutation_scale=arrow))
    ax.text(X0 - 20, y, 'WASF1', ha='right', va='center', fontsize=size, color=C_GENE, fontweight='bold',
            style='italic')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--work', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--fp-scale', type=float, default=1.0)
    a = ap.parse_args()
    w, out = Path(a.work), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    f1 = w / 'fig1_wasf1'
    fire = pd.read_csv(f1 / 'fire/fire_profiles.tsv.gz', sep='\t')
    ref = json.load(open(f1 / 'raw/ucsc_ncbiRefSeq_wasf1.json'))['ncbiRefSeq']
    ins = pd.read_csv(f1 / 'inserts.tsv', sep='\t').set_index('construct')
    assert ins.loc[PLASMID_SRC, 'hg38_start'] == PEAK[0]
    assert ins.loc['WASF1', 'hg38_start'] <= X0 and ins.loc['WASF1', 'hg38_end'] >= X1   # 4.5-kb insert holds X0-X1
    mol = pd.read_parquet(f1 / f'profiles/{READ_SRC}.molecules.parquet')
    p300 = pd.read_parquet(f1 / f'profiles/{PLASMID_SRC}.molecules.parquet')
    endo = pd.read_parquet(f1 / 'profiles/endogenous.molecules.parquet')
    keep, fps_all, feats = select_reads(mol, *PEAK)
    rng = np.random.default_rng(SEED)
    pick = np.sort(rng.choice(len(keep), N_SHOW, replace=False))
    order = [keep[pick[j]] for j in leaves_list(linkage(feats[pick], 'ward'))]
    pd.DataFrame({'row': range(1, N_SHOW + 1), 'read': mol.read.iloc[order].to_numpy()}).to_csv(
        out / 'panel_a_reads.tsv', sep='\t', index=False)
    r45 = ins.loc['WASF1']
    fire_rd = fire_elements(f1 / 'fire/WASF1_full.extract.tsv.gz', 'WASF1', int(r45.hg38_start - r45.construct_start),
                            set(mol.read.iloc[order]))
    tracks = {'genome': track_shares(endo, X0, X1), '300bp': track_shares(p300, *PEAK)}
    for k, t in tracks.items():
        fr = fire[fire.track == k].set_index('pos')
        t['FIRE_share'] = fr.FIRE_share.reindex(t.pos).to_numpy()
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

    # ---- a: schematic slot, then the construct map and molecules beside the PacBio icon
    sch_h = 86.0
    sc = sch_h / CROP[3]
    sch_w = CROP[2] * sc
    slot = (12.0, 12.0)
    fig.text(4 / Wpt, 1 - 6 / Hpt, 'a', fontsize=10, fontweight='bold', va='top')
    rx0, rx1 = slot[0] + sch_w + 34, Wpt - 10
    my0, my1, ry0, ry1 = 10.0, 26.0, 30.0, 96.0
    axm = axpt(rx0, my0, rx1, my1)
    axm.set_xlim(X0 - 0.5, READ_X1 + 0.5)
    axm.set_ylim(0, 1)
    axm.axis('off')
    axm.add_patch(Rectangle((X0 - 0.5, 0.2), READ_X1 - X0 + 1, 0.16, fc=C_BACKBONE, lw=0, gid='construct'))
    axm.add_patch(Rectangle((PEAK[0] - 0.5, 0.1), PEAK[1] - PEAK[0] + 1, 0.36, fc=C_GENE, lw=0, gid='peak4'))
    axm.add_patch(Rectangle((MOTIF[0] - 0.5, 0.56), MOTIF[1] - MOTIF[0] + 1, 0.3, fc=C_MOTIF, lw=0, gid='motif-map'))
    axm.text(PEAK[0] - 15, 0.71, '300-bp element', ha='right', va='center', fontsize=5.4, color=C_GENE,
             fontweight='bold')
    axm.text(X0, 0.71, '4.5-kb construct', ha='left', va='center', fontsize=5.0, color=INK2)
    axm.text(MOTIF[1] + 15, 0.71, 'GC box', ha='left', va='center', fontsize=5.4, color=INK, fontweight='bold',
             clip_on=False)
    axr = axpt(rx0, ry0, rx1, ry1)
    axr.axvspan(PEAK[0] - 0.5, PEAK[1] + 0.5, color=HILITE, lw=0, zorder=0)
    axr.axvspan(MOTIF[0] - 0.5, MOTIF[1] + 0.5, color=C_MOTIF, alpha=0.25, lw=0, zorder=0)
    axr.axvspan(SMALL_SITE - SMALL_WIN - 0.5, SMALL_SITE + SMALL_WIN + 0.5, color=COL['fp_20_30'], alpha=0.15, lw=0,
                zorder=0)
    n_fire_drawn = 0
    for y, i in enumerate(reversed(order)):
        at = np.array(json.loads(mol['at'].iloc[i]), int)
        axr.plot([max(at.min(), X0), min(at.max(), READ_X1)], [y, y], color='#BDBDBD', lw=0.4, zorder=1)
        for s, e, ln in fps_all[i]:                                  # nucleosomes
            if fp_code(ln) == 4 and e >= X0 and s <= READ_X1:
                s_, e_ = max(s, X0), min(e, READ_X1)
                axr.add_patch(Rectangle((s_ - 0.5, y - 0.36), e_ - s_ + 1, 0.72, fc=COL['nuc'], lw=0, zorder=2))
        for s, e in fire_rd.get(mol.read.iloc[i], []):              # FIRE elements
            if e >= X0 and s <= READ_X1:
                s_, e_ = max(s, X0), min(e, READ_X1)
                axr.add_patch(Rectangle((s_ - 0.5, y - 0.36), e_ - s_ + 1, 0.72, fc=COL['FIRE'], alpha=FIRE_ALPHA,
                                        lw=0, zorder=3))
                n_fire_drawn += 1
        for code in (3, 2, 1):                                       # footprints
            for s, e, ln in fps_all[i]:
                if fp_code(ln) == code and e >= X0 and s <= READ_X1:
                    s_, e_ = max(s, X0), min(e, READ_X1)
                    axr.add_patch(Rectangle((s_ - 0.5, y - 0.36), e_ - s_ + 1, 0.72, fc=COL[CODE[code]], lw=0,
                                            zorder=4))
    axr.set_xlim(X0 - 0.5, READ_X1 + 0.5)
    axr.set_ylim(-0.6, N_SHOW - 0.4)
    axr.axis('off')
    axr.plot([READ_X1 - 200, READ_X1], [-1.3, -1.3], color=INK, lw=0.8, clip_on=False)
    axr.text(READ_X1 - 100, -1.7, '200 bp', ha='center', va='top', fontsize=5.3, color=INK)
    ov.text(rx0, ry1 + 4, f'4.5-kb construct molecules with footprints on the\nGC box and the 20–30 bp site '
            f'({N_SHOW} of {len(keep):,}; {len(mol):,} total)', ha='left', va='top', fontsize=5.2, color=INK2,
            linespacing=1.15)
    ov.annotate('', xy=(rx0 - 4, 56), xytext=(slot[0] + sch_w - 8, 56), gid='arrow-to-reads',
                arrowprops=dict(arrowstyle='-|>', color='#3B3B3B', lw=1.3, shrinkA=0, shrinkB=0, mutation_scale=8))
    fig.legend([Patch(color=COL[c], alpha=FIRE_ALPHA if c == 'FIRE' else 1) for c in COL],
               [NAMES[c] for c in COL], loc='center', bbox_to_anchor=(0.5, 1 - 122 / Hpt), ncol=5, frameon=False,
               fontsize=6, handlelength=1.3, handleheight=0.9, columnspacing=1.6, handletextpad=0.5)

    # ---- b: gene model (with the motif) + occupancy tracks
    bx0, bx1 = 112.0, Wpt - 52
    yb = 134
    t1, t2 = (yb + 40, yb + 86), (yb + 92, yb + 138)
    fig.text(4 / Wpt, 1 - yb / Hpt, 'b', fontsize=10, fontweight='bold', va='top')
    axg = axpt(bx0, yb + 6, bx1, yb + 34)
    axg.set_xlim(X0, X1)
    axg.set_ylim(0, 1.05)
    axg.axis('off')
    gene_model(axg, ref)
    axg.add_patch(Rectangle((MOTIF[0] - 0.5, 0.74), MOTIF[1] - MOTIF[0] + 1, 0.2, fc=C_MOTIF, lw=0, gid='motif-b'))
    axg.text(MOTIF[1] + 12, 0.84, MOTIF_NAME, ha='left', va='center', fontsize=5.6, color=INK, fontweight='bold',
             clip_on=False)
    axg.plot([X0 + 150, X0 + 350], [0.9, 0.9], color=INK, lw=0.8)
    axg.text(X0 + 250, 0.96, '200 bp', ha='center', va='bottom', fontsize=5.5, color=INK)
    axes = []
    for track, cats, (y0, y1) in (('genome', ['nuc', 'FIRE'], t1), ('300bp', list(COL), t2)):
        ax = axpt(bx0, y0, bx1, y1)
        axes.append(ax)
        t = tracks[track]
        x = t.pos.to_numpy() + 0.5
        layers = {c: np.minimum(t[f'{c}_smooth'].to_numpy() * (a.fp_scale if c.startswith('fp_') else 1), 1.0)
                  for c in cats}
        fps = sorted([c for c in cats if c.startswith('fp_')], key=lambda c: -np.nanmax(layers[c]))
        for cat in [c for c in ('nuc', 'FIRE') if c in cats] + fps:      # nucleosome, FIRE, footprints (tallest first)
            yv = layers[cat]
            ax.fill_between(x, 0, np.nan_to_num(yv), where=~np.isnan(yv), color=COL[cat], lw=0,
                            alpha=FIRE_ALPHA if cat == 'FIRE' else 1, gid=f'{track}-{cat}-fill')
        ax.axvspan(MOTIF[0] - 0.5, MOTIF[1] + 0.5, color=C_MOTIF, alpha=0.3, lw=0, zorder=3, gid=f'{track}-motif')
        if track == '300bp':
            for xx in PEAK:
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
    scale_txt = 'unscaled' if a.fp_scale == 1 else f'×{a.fp_scale:g}'
    ov.text((bx0 + bx1) / 2, t2[1] + 6, f'Endogenous track: nucleosomes and FIRE only. Plasmid track: sub-nucleosome '
            f'footprint shares {scale_txt}. Calls <20 bp not shown.', ha='center', va='top', fontsize=5.4, color=INK2)

    # ---- write SVG, then embed the schematic vector in its slot
    svg = out / 'fig1.svg'
    fig.savefig(svg, facecolor='white')
    plt.close(fig)
    s = expand_fonts(svg.read_text().replace("'Arimo'", 'Arial, Helvetica, sans-serif').replace('Arimo', 'Arial'))[0]
    sch = (w / 'fig1/schematic2/schematic.svg').read_text()
    body = re.sub(r'^.*?<svg[^>]*>', '', sch, count=1, flags=re.S).rsplit('</svg>', 1)[0]
    body = re.sub(r'<metadata>.*?</metadata>', '', body, flags=re.S)
    assert body.count('>PacBio</text>') == 1
    body = body.replace('>PacBio</text>', '>PacBio HiFi</text>')
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
    json.dump(dict(gene='WASF1', motif_name=MOTIF_NAME, read_source=READ_SRC, plasmid_track=PLASMID_SRC,
                   molecules=len(mol), selected_molecules=len(keep), shown=N_SHOW, seed=SEED,
                   small_site=[SMALL_SITE - SMALL_WIN, SMALL_SITE + SMALL_WIN], read_window=[X0, READ_X1],
                   reads_with_fire=sum(bool(v) for v in fire_rd.values()), fire_elements_drawn=n_fire_drawn, plasmid_300bp_molecules=len(p300),
                   endogenous_molecules=len(endo), endogenous_max_coverage=int(tracks['genome'].molecules.max()),
                   bins=BINS, sigma=SIGMA, fp_scale_300bp=a.fp_scale, fire_alpha=FIRE_ALPHA, window=[X0, X1],
                   peak=list(PEAK), motif=list(MOTIF), commit=commit, date=str(date.today())),
              open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

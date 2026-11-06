"""Grant Figure 1, WASF1 draft 13 (half page; fig1_rules.json: wasf1.draft3-5, draft7-13_halfpage): plasmid Fiber-seq
workflow above one shared WASF1 axis holding the gene model, occupancy tracks and single 300-bp construct molecules
(the element aligned with the tracks).
Laid out at its final size, 3.5 in wide, text >= 7 pt (as Figure 2). One Illustrator-editable SVG (live Arial text;
workflow icons are the HPRC main-figure vectors) plus PDF and PNG renders (rsvg-convert with fonts.conf).

Footprint bins (wasf1.draft3, user-chosen after fig1_fpsizes.py --gene WASF1): 20-30, 31-50, 51-80 bp, nucleosome
>= 81 bp (as Figure 3); calls < 20 bp are not drawn or counted.
a  workflow schematic (labels raised to 7 pt in the embedded copy); arrow from the sequencer into b
b  shared hg38 axis X0-X1, top to bottom: WASF1 gene model with the GC box; endogenous Fiber-seq nucleosomes and FIRE;
   plasmid 300-bp all bins x --fp-scale. Shares = molecules with a call of the bin covering a base / molecules covering
   it, Gaussian sigma SIGMA bp within covered runs; FIRE from fig1_fire_profiles.py. Draw order: nucleosome, FIRE
   (alpha FIRE_ALPHA), then footprint bins tallest first. Below them, on the same axis (outside the element the read
   is plasmid sequence in construct coordinates + offset): N_SHOW 300-bp construct molecules (all kept reads span the
   element) with a 20-80 bp call over the GC box (>= MIN_OVERLAP bp) and a separate 20-30 bp call with midpoint within
   SMALL_SITE +/- SMALL_WIN, then N_NUC molecules whose nucleosome calls cover >= NUC_COVER of the element. Both groups
   are the molecules whose nucleosome occupancy outside the element (X0-X1) is nearest (Hamming) the medoid of the
   footprinted candidates, with swaps so the footprinted set holds a 20-30, 31-50 and 51-80 bp call inside the element
   a footprint at the left peak LEFT_SITE and a 31-50 bp call covering it (pick_typical); each group ordered by Ward clustering of one-hot per-base call class over the element; every bin and
   nucleosomes drawn,
   with FIRE elements (fig1_fire.sh) between the nucleosomes and the footprints. Legend below. Commit and date in
   summary.json.
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
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.textpath import TextPath  # noqa: E402
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
FS = 7                                    # every label (pt); panel letters 9 pt
plt.rcParams.update({'font.family': 'Arimo', 'svg.fonttype': 'none', 'font.size': FS})

INK, INK2 = '#2E2E2E', '#555555'
BINS = {'fp_20_30': (20, 30), 'fp_31_50': (31, 50), 'fp_51_80': (51, 80), 'nuc': (81, 10 ** 6)}
COL = {'nuc': '#C7C8CA', 'FIRE': '#BA7C69', 'fp_51_80': '#E77C34', 'fp_31_50': '#338B8C', 'fp_20_30': '#755B9D'}
NAMES = {'nuc': 'Nucleosome (≥81 bp)', 'FIRE': 'FIRE', 'fp_51_80': '51–80 bp footprint',
         'fp_31_50': '31–50 bp footprint', 'fp_20_30': '20–30 bp footprint'}
CODE = {0: 'none', 1: 'fp_20_30', 2: 'fp_31_50', 3: 'fp_51_80', 4: 'nuc'}
C_GENE, HILITE, C_MOTIF, C_BACKBONE = '#1f6fb4', '#FBE3D3', '#B8860B', '#8B9094'
X0, X1 = 110178550, 110180380          # WASF1 exons 1-2 .. peak 4 downstream nucleosome (wasf1.draft2 window_b)
PEAK = (110179911, 110180211)          # FIRE peak 4 300-bp element (WASF12 insert)
MOTIF = (110180054, 110180097)         # SP/KLF GC box (wasf1.draft1 motif)
SMALL_SITE, SMALL_WIN = 110180148, 10   # 20-30 bp call-midpoint density peak in WASF12 (wasf1.draft3 small_peak)
LEFT_SITE = 110179959                   # left sub-nucleosome footprint peak of the 300-bp track (wasf1.draft11 peaks)
MOTIF_NAME, GENE_TX = 'GC box (SP/KLF)', 'NM_003931.3'
READ_SRC, PLASMID_SRC = 'WASF12', 'WASF12'       # molecules and plasmid track in b: 300-bp construct
SEED, N_SHOW, SPAN_SLACK, MIN_OVERLAP = 20261102, 10, 20, 10
N_NUC, NUC_COVER = 3, 0.8
SIGMA, FIRE_ALPHA = 3, 0.65
W_IN = 3.5
CROP = (58.63 + 1, 46.63 + 82 + 1.5, 185.0, 54)   # schematic crop (fig1_schematic.CROP)
SCH_FS = 5                                # schematic label size (pt) in schematic.svg
PACBIO_LAB = (180.5 + 46.63 + 2.0, 133.6 + 46.63)   # "sequencing" label (centre x, baseline y), schematic page pt


def fp_code(ln):
    return next((k for k, c in CODE.items() if k and BINS[c][0] <= ln <= BINS[c][1]), None)


def onehot(f, lo, hi):
    """One-hot per-base call class (none, 20-30, 31-50, 51-80, nucleosome) over [lo, hi]: clustering features."""
    v = np.zeros(hi - lo + 1, int)
    for s, e, ln in f:
        c = fp_code(ln)
        if c and e >= lo and s <= hi:
            v[max(s, lo) - lo:min(e, hi) - lo + 1] = c
    return np.eye(5)[v].ravel()


def nuc_cover(f, lo, hi):
    """Share of [lo, hi] covered by nucleosome calls (>= 81 bp)."""
    cov = np.zeros(hi - lo + 1, bool)
    for s, e, ln in f:
        if fp_code(ln) == 4 and e >= lo and s <= hi:
            cov[max(s, lo) - lo:min(e, hi) - lo + 1] = True
    return cov.mean()


def flank_nuc(f):
    """Nucleosome occupancy (calls >= 81 bp) per base over X0-X1 outside the 300-bp element (the flank features)."""
    v = np.zeros(X1 - X0 + 1, bool)
    for s, e, ln in f:
        if fp_code(ln) == 4 and e >= X0 and s <= X1:
            v[max(s, X0) - X0:min(e, X1) - X0 + 1] = True
    return np.r_[v[:PEAK[0] - X0], v[PEAK[1] - X0 + 1:]]


def hamming(A, B):
    """Share of flank bases where each row of A disagrees with each row of B (boolean matrices)."""
    A, B = A.astype(np.float32), B.astype(np.float32)
    return (A @ (1 - B).T + (1 - A) @ B.T) / A.shape[1]


def pick_typical(fps_all, keep, nuc_ok, in_el):
    """Footprinted and nucleosome-covered molecules nearest the medoid flank pattern of the footprinted candidates.
    Set-level constraints on the footprinted picks: a 20-30, a 31-50 and a 51-80 bp call inside the element, and a
    20-80 bp call with midpoint within LEFT_SITE +/- SMALL_WIN, and a 31-50 bp call covering LEFT_SITE. For each unmet constraint the farthest pick that is not
    the only carrier of another met constraint is swapped for the nearest candidate that meets it."""
    cons = [lambda f, b=b: b in in_el(f) for b in (1, 2, 3)] + \
           [lambda f: any(fp_code(ln) in (1, 2, 3) and abs((s + e) / 2 - LEFT_SITE) <= SMALL_WIN for s, e, ln in f),
            lambda f: any(fp_code(ln) == 2 and s <= LEFT_SITE <= e for s, e, ln in f)]
    Fk = np.array([flank_nuc(fps_all[i]) for i in keep])
    med = int(np.argmin(hamming(Fk, Fk).mean(1)))
    dk = hamming(Fk, Fk[[med]])[:, 0]
    rank = list(np.argsort(dk, kind='stable'))
    pick, swaps = rank[:N_SHOW], 0
    for c in cons:
        if any(c(fps_all[keep[j]]) for j in pick):
            continue
        sole = {j for c2 in cons for j in pick if c2(fps_all[keep[j]])
                and sum(c2(fps_all[keep[k]]) for k in pick) == 1}
        drop = max((j for j in pick if j not in sole), key=lambda j: dk[j])
        add = next(j for j in rank if j not in pick and c(fps_all[keep[j]]))
        pick = [j for j in pick if j != drop] + [add]
        swaps += 1
    assert all(any(c(fps_all[keep[j]]) for j in pick) for c in cons)
    Fn = np.array([flank_nuc(fps_all[i]) for i in nuc_ok])
    dn = hamming(Fn, Fk[[med]])[:, 0]
    npick = list(np.argsort(dn, kind='stable')[:N_NUC])
    return sorted(pick), npick, dict(medoid=int(keep[med]), swaps=swaps, dist_footprinted=[float(dk[j]) for j in pick],
                                     dist_nucleosome=[float(dn[j]) for j in npick])


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
            keep.append(i)
            feats.append(onehot(f, lo, hi))
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
    in_el = lambda f: {fp_code(ln) for s, e, ln in f if e >= PEAK[0] and s <= PEAK[1]}   # noqa: E731
    nuc_ok = [i for i, f in enumerate(fps_all) if nuc_cover(f, *PEAK) >= NUC_COVER]
    pick, npick, typ = pick_typical(fps_all, keep, nuc_ok, in_el)
    pick = np.array(pick)
    order = [keep[pick[j]] for j in leaves_list(linkage(feats[pick], 'ward'))]
    Fs = np.array([flank_nuc(fps_all[i]) for i in order + [nuc_ok[j] for j in npick]])
    Fr = np.array([flank_nuc(fps_all[i]) for i in keep])
    rng = np.random.default_rng(SEED)                              # reference: random sets of N_SHOW candidates
    rand = [hamming(Fr[r], Fr[r])[np.triu_indices(N_SHOW, 1)].mean()
            for r in (rng.choice(len(keep), N_SHOW, replace=False) for _ in range(1000))]
    typ.update(mean_pairwise_shown_footprinted=float(hamming(Fs[:N_SHOW], Fs[:N_SHOW])[np.triu_indices(N_SHOW, 1)].mean()),
               mean_pairwise_shown_all=float(hamming(Fs, Fs)[np.triu_indices(len(Fs), 1)].mean()),
               mean_pairwise_random10_median=float(np.median(rand)), mean_pairwise_random10_min=float(np.min(rand)))
    order_nuc = [nuc_ok[npick[j]] for j in leaves_list(linkage(np.array([onehot(fps_all[nuc_ok[j]], *PEAK)
                                                                           for j in npick]), 'ward'))]
    shown = order + order_nuc
    pd.DataFrame({'row': range(1, len(shown) + 1), 'read': mol.read.iloc[shown].to_numpy(),
                  'group': ['footprinted'] * N_SHOW + ['nucleosome'] * N_NUC,
                  'nuc_cover_element': [round(nuc_cover(fps_all[i], *PEAK), 3) for i in shown]}).to_csv(
        out / 'molecules_shown.tsv', sep='\t', index=False)
    rc = ins.loc[READ_SRC]
    fire_rd = fire_elements(f1 / f'fire/{READ_SRC}.extract.tsv.gz', READ_SRC, int(rc.hg38_start - rc.construct_start),
                            set(mol.read.iloc[shown]))
    tracks = {'genome': track_shares(endo, X0, X1), '300bp': track_shares(p300, *PEAK)}
    for k, t in tracks.items():
        fr = fire[fire.track == k].set_index('pos')
        t['FIRE_share'] = fr.FIRE_share.reindex(t.pos).to_numpy()
        t['FIRE_smooth'] = smooth(t.FIRE_share.to_numpy())
    pd.concat([t.assign(track=k) for k, t in tracks.items()]).to_csv(out / 'tracks.tsv', sep='\t', index=False,
                                                                      float_format='%.5g')
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()

    # layout (pt, origin top left), computed before the figure so the height fits
    Wpt = W_IN * 72
    sc = (Wpt - 10) / CROP[2]                                      # a: schematic scale
    slot = (8.0, 2.0)
    sch_fs = FS / sc                                               # embedded label size that prints at FS pt
    sch_bottom = slot[1] + CROP[3] * sc
    yb0 = sch_bottom + 8
    bx0, bx1 = 62.0, Wpt - 32                                      # b: data x range
    g0, g1 = yb0 + 2, yb0 + 16                                     # gene model
    t1 = (g1 + 7, g1 + 37)                                         # endogenous track
    t2 = (t1[1] + 4, t1[1] + 34)                                   # 300-bp plasmid track
    GAP = 0                                                        # rows between footprinted and nucleosome reads
    r0, r1 = t2[1] + 5, t2[1] + 5 + 4.0 * (N_SHOW + GAP + N_NUC)   # molecules (same axis as the tracks)
    ly = r1 + 13                                                   # legend (200-bp bar above it)
    Hpt = ly + 18
    fig = plt.figure(figsize=(W_IN, Hpt / 72))

    def axpt(x0, y0, x1, y1):            # axes from a box in points, origin top left
        return fig.add_axes([x0 / Wpt, 1 - y1 / Hpt, (x1 - x0) / Wpt, (y1 - y0) / Hpt])

    def letter_lab(ch, x, y):
        fig.text(x / Wpt, 1 - y / Hpt, ch, fontsize=9, fontweight='bold', va='top', ha='left')

    ov = fig.add_axes([0, 0, 1, 1], facecolor='none')
    ov.set_xlim(0, Wpt)
    ov.set_ylim(Hpt, 0)
    ov.axis('off')

    # ---- a: workflow schematic (embedded after saving); arrow from the sequencer into b
    letter_lab('a', 7, 1)
    ax_ = slot[0] + (PACBIO_LAB[0] - CROP[0]) * sc
    ay_ = slot[1] + (PACBIO_LAB[1] - CROP[1]) * sc
    xe, yh = bx1 - 4, t1[0] - 5.5                                  # bent: down, left (under the GC box label), down
    ov.plot([ax_, ax_, xe], [ay_ + 3, yh, yh], color='#3B3B3B', lw=1.0, solid_joinstyle='round', gid='arrow-to-b')
    ov.annotate('', xy=(xe, t1[0] - 0.4), xytext=(xe, yh), gid='arrow-to-b-head',
                arrowprops=dict(arrowstyle='-|>', color='#3B3B3B', lw=1.0, shrinkA=0, shrinkB=0, mutation_scale=7))

    # ---- b: shared WASF1 axis: gene model, tracks, then the 300-bp molecules
    letter_lab('b', 7, yb0)
    axg = axpt(bx0, g0, bx1, g1)
    axg.set_xlim(X0, X1)
    axg.set_ylim(0, 1.05)
    axg.axis('off')
    gene_model(axg, ref, y=0.38, h=0.6, size=FS, tss_h=0.6, arrow=5)
    axg.add_patch(Rectangle((MOTIF[0] - 0.5, 0.16), MOTIF[1] - MOTIF[0] + 1, 0.44, fc=C_MOTIF, lw=0, gid='motif-b'))
    axg.text(MOTIF[1] + 25, 0.38, 'GC box', ha='left', va='center', fontsize=FS, color=INK, fontweight='bold')
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
                ax.axvline(xx, color=INK2, lw=0.5, ls=(0, (2.5, 1.5)))
        ax.set_xlim(X0, X1)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.yaxis.tick_right()
        ax.set_yticks([0, 1])
        ax.tick_params(axis='y', labelsize=FS, length=2, pad=1.5, width=0.6)
        lab0, lab1 = ax.get_yticklabels()                          # keep 0 / 1 inside the track: no overlap
        lab0.set_va('bottom')
        lab1.set_va('top')
        for sp in ax.spines.values():
            sp.set_linewidth(0.6)
    scale_txt = '' if a.fp_scale == 1 else f'\n(footprints ×{a.fp_scale:g})'
    for ax, lab in ((axes[0], 'Endogenous\nlocus'), (axes[1], '300-bp\nplasmid' + scale_txt)):
        ax.text(-0.03, 0.5, lab, transform=ax.transAxes, ha='right', va='center', fontsize=FS, color=INK,
                multialignment='right', linespacing=1.1, fontweight='bold')
    ov.text(bx1 + 18, (t1[0] + t2[1]) / 2, 'Fraction of molecules\nwith each feature', rotation=270, ha='center',
            va='center', fontsize=FS, color=INK, rotation_mode='anchor', linespacing=1.1)
    axr = axpt(bx0, r0, bx1, r1)
    axr.axvspan(MOTIF[0] - 0.5, MOTIF[1] + 0.5, color=C_MOTIF, alpha=0.25, lw=0, zorder=0)
    for xx in PEAK:
        axr.axvline(xx, color=INK2, lw=0.5, ls=(0, (2.5, 1.5)), zorder=5)
    ys = [N_NUC + GAP + N_SHOW - 1 - k for k in range(N_SHOW)] + [N_NUC - 1 - k for k in range(N_NUC)]
    n_fire_drawn = 0
    for y, i in zip(ys, shown):
        at = np.array(json.loads(mol['at'].iloc[i]), int)
        axr.plot([max(at.min(), X0), min(at.max(), X1)], [y, y], color='#BDBDBD', lw=0.3, zorder=1)
        for s_e_ln in fps_all[i]:                                    # nucleosomes
            s, e, ln = s_e_ln
            if fp_code(ln) == 4 and e >= X0 and s <= X1:
                s_, e_ = max(s, X0), min(e, X1)
                axr.add_patch(Rectangle((s_ - 0.5, y - 0.38), e_ - s_ + 1, 0.76, fc=COL['nuc'], lw=0, zorder=2))
        for s, e in fire_rd.get(mol.read.iloc[i], []):              # FIRE elements
            if e >= X0 and s <= X1:
                s_, e_ = max(s, X0), min(e, X1)
                axr.add_patch(Rectangle((s_ - 0.5, y - 0.38), e_ - s_ + 1, 0.76, fc=COL['FIRE'], alpha=FIRE_ALPHA,
                                        lw=0, zorder=3))
                n_fire_drawn += 1
        for code in (3, 2, 1):                                       # footprints
            for s, e, ln in fps_all[i]:
                if fp_code(ln) == code and e >= X0 and s <= X1:
                    s_, e_ = max(s, X0), min(e, X1)
                    axr.add_patch(Rectangle((s_ - 0.5, y - 0.38), e_ - s_ + 1, 0.76, fc=COL[CODE[code]], lw=0,
                                            zorder=4))
    axr.set_xlim(X0, X1)
    axr.set_ylim(-0.6, N_NUC + GAP + N_SHOW - 0.4)
    axr.axis('off')
    xb2 = lambda g: bx0 + (g - X0) / (X1 - X0) * (bx1 - bx0)      # noqa: E731  (shared axis, pt)
    ov.plot([xb2(X1 - 200), xb2(X1)], [r1 + 5, r1 + 5], color=INK, lw=0.7)   # 200-bp bar under the molecules
    ov.text(xb2(X1 - 200) - 3, r1 + 5, '200 bp', ha='right', va='center', fontsize=FS, color=INK)
    ov.text(bx0 - 4.5, (r0 + r1) / 2, '300-bp\nplasmid\nmolecules', ha='right', va='center', fontsize=FS, color=INK,
            fontweight='bold', multialignment='right', linespacing=1.1)
    leg = [(('nuc', 'Nucleosome'), ('FIRE', 'FIRE')),
           (('fp_20_30', '20–30 bp footprint'), ('fp_31_50', '31–50 bp footprint')),
           (('fp_51_80', '51–80 bp footprint'), ('motif', 'SP/KLF GC box'))]
    tw = lambda t: TextPath((0, 0), t, size=FS, prop=FontProperties(family='Arimo')).get_extents().width   # noqa: E731
    cw_ = [9 + max(tw(lab) for _, lab in col) for col in leg]      # swatch + gap + widest label per column
    xs = (Wpt - sum(cw_) - 10 * (len(cw_) - 1)) / 2 + np.r_[0, np.cumsum(np.array(cw_[:-1]) + 10)]
    for (col_items, xx) in zip(leg, xs):
        for k, (c, lab) in enumerate(col_items):
            yy = ly + 9 * k
            fc, al = (C_MOTIF, 0.5) if c == 'motif' else (COL[c], FIRE_ALPHA if c == 'FIRE' else 1)
            ov.add_patch(Rectangle((xx, yy), 7, 6, fc=fc, alpha=al, lw=0))
            ov.text(xx + 9, yy + 3, lab, fontsize=FS, color=INK, va='center', ha='left')

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
    body = re.sub(r'font-size="5"', f'font-size="{sch_fs:.3f}"', body)       # schematic labels print at FS pt
    body, nbg = re.subn(r'<g id="patch_1">\s*<path[^>]*style="fill: #ffffff"\s*/>\s*</g>', '', body)
    assert nbg == 1                                                # drop the page background: keeps letter a visible
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
                   molecules=len(mol), selected_molecules=len(keep), shown=N_SHOW, seed=SEED, typical_selection=typ,
                   nucleosome_candidates=len(nuc_ok), nuc_cover_min=NUC_COVER, shown_nucleosome=N_NUC,
                   small_site=[SMALL_SITE - SMALL_WIN, SMALL_SITE + SMALL_WIN], read_window=list(PEAK),
                   reads_with_fire=sum(bool(v) for v in fire_rd.values()), fire_elements_drawn=n_fire_drawn, plasmid_300bp_molecules=len(p300),
                   endogenous_molecules=len(endo), endogenous_max_coverage=int(tracks['genome'].molecules.max()),
                   bins=BINS, sigma=SIGMA, fp_scale_300bp=a.fp_scale, fire_alpha=FIRE_ALPHA, window=[X0, X1],
                   peak=list(PEAK), motif=list(MOTIF), commit=commit, date=str(date.today())),
              open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

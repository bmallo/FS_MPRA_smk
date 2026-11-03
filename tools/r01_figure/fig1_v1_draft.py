"""Draft grant Figure 1: what plasmid Fiber-seq is and what it detects (DPAGT1 locus, HEK293T).

a  workflow schematic                      b  scale: molecules and clones per SNV in the LDLR library (HepG2 run 3)
c  DPAGT1 locus: gene model, plasmid inserts, motifs; accessibility (m6A fraction, 51-bp smoothed) for the endogenous
   genome, the 4.1-kb construct and the two ~300-bp constructs, each on its own track (no pooling)
d  FIRE peak 1 zoom: relative accessibility (m6A fraction / its window mean), endogenous vs ~300-bp construct
e  single molecules over peak 1: endogenous and ~300-bp construct; m6A marks and footprints by size class
f  footprint sizes on the ~300-bp constructs (all calls inside the inserts): TF-size, PIC-size, nucleosome
Inputs: fig1/profiles (fig1_profiles.py), fig1/inserts.tsv, fig1/motifs.tsv, UCSC RefSeq json, LDLR clones.tsv.
Usage: python tools/r01_figure/fig1.py --work results/r01_figure_20261102 --root results/analysis_v2_20261028 --output <dir>
"""
import argparse
import importlib.util
import json
import subprocess
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

_s = importlib.util.spec_from_file_location('fig_repeat3', Path(__file__).with_name('fig_repeat3.py'))
fr = importlib.util.module_from_spec(_s)
_s.loader.exec_module(fr)
INK, INK2, MUTED, GRID = fr.INK, fr.INK2, fr.MUTED, fr.GRID
C_ENDO, C_FULL, C_300 = INK, '#2a78d6', '#eb6834'
C_TF, C_PIC, C_NUC = '#1baf7a', '#6250d6', '#9a9893'
WIN = (119100000, 119104500)
PEAK1, PEAK3 = (119100271, 119100570), (119102160, 119102459)
CTCF1 = (119100436, 119100455)
SEED, NPILE = 20261102, 30


def smooth(s, w=51):
    return s.rolling(w, center=True, min_periods=w // 3).mean()


def smooth_ratio(p, w=51):
    """m6A fraction smoothed as (sum of m6A calls) / (sum of observed A/T) over a w-bp window: no gaps at G/C runs."""
    meth = (p.m6a_frac.fillna(0) * p.at_observed)
    num = meth.rolling(w, center=True, min_periods=1).sum()
    den = p.at_observed.rolling(w, center=True, min_periods=1).sum()
    return (num / den).where(den > 0)


def footprint_class(ln):
    return 'tf' if 10 <= ln <= 40 else 'pic' if 41 <= ln <= 80 else 'nuc' if ln >= 85 else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--work', required=True)
    ap.add_argument('--root', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    w, out = Path(a.work), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    f1 = w / 'fig1'
    P = {d: pd.read_csv(f1 / f'profiles/{d}.profile.tsv', sep='\t').set_index('pos')
         for d in ('endogenous', 'DPAGT1_full', 'DPAGT11', 'DPAGT13')}
    M = {d: pd.read_parquet(f1 / f'profiles/{d}.molecules.parquet') for d in ('endogenous', 'DPAGT11', 'DPAGT13')}
    ins = pd.read_csv(f1 / 'inserts.tsv', sep='\t').set_index('construct')
    mot = pd.read_csv(f1 / 'motifs.tsv', sep='\t')
    ref = json.load(open(f1 / 'raw/ucsc_ncbiRefSeq_dpagt1.json'))['ncbiRefSeq']
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()

    fig = plt.figure(figsize=(7.2, 7.4))
    gA = fig.add_gridspec(1, 2, left=0.06, right=0.985, top=0.965, bottom=0.80, width_ratios=[2.2, 1], wspace=0.3)
    gC = fig.add_gridspec(1, 1, left=0.13, right=0.985, top=0.735, bottom=0.43)
    gD = fig.add_gridspec(1, 3, left=0.08, right=0.985, top=0.35, bottom=0.07, wspace=0.42, width_ratios=[1, 1.15, 1])

    # a schematic
    ax = fig.add_subplot(gA[0, 0])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3)
    ax.axis('off')
    steps = [('Reporter plasmids\n(WT or variant\nlibraries)', 0.0), ('Transfect cells;\nHia5 m6A-marks\naccessible DNA', 2.55),
             ('PacBio HiFi reads\nwhole plasmids', 5.1), ('Per molecule:\ngenotype + footprints', 7.65)]
    for i, (txt, x0) in enumerate(steps):
        ax.add_patch(FancyBboxPatch((x0, 1.25), 2.2, 1.55, boxstyle='round,pad=0.02,rounding_size=0.15',
                                    fc='#f4f3ef', ec='#d6d4ce', lw=0.6))
        ax.text(x0 + 1.1, 2.03, txt, ha='center', va='center', fontsize=5.8, color=INK, linespacing=1.15)
        if i < 3:
            ax.add_patch(FancyArrowPatch((x0 + 2.22, 2.03), (x0 + 2.53, 2.03), arrowstyle='-|>', mutation_scale=7,
                                         color=INK2, lw=0.8))
    # cartoon molecule
    y = 0.55
    ax.plot([0.3, 9.7], [y, y], color=INK2, lw=0.8)
    for xc in (1.4, 2.6, 6.6, 7.8):
        ax.add_patch(Ellipse((xc, y), 0.95, 0.42, fc=C_NUC, ec='none', alpha=0.85))
    ax.add_patch(Rectangle((4.0, y - 0.13), 0.45, 0.26, fc=C_TF, ec='none'))
    ax.add_patch(Rectangle((5.0, y - 0.15), 0.9, 0.30, fc=C_PIC, ec='none'))
    for xm in np.r_[np.linspace(3.25, 3.85, 4), np.linspace(4.6, 4.9, 2), np.linspace(6.0, 6.1, 1), np.linspace(8.5, 9.5, 5)]:
        ax.plot(xm, y + 0.27, marker='v', ms=2.6, color=INK, lw=0)
    for x_, lab, c in ((2.0, 'nucleosomes', C_NUC), (4.22, 'TF', C_TF), (5.45, 'Pol II / PIC', C_PIC),
                       (9.0, 'm6A = accessible', INK)):
        ax.text(x_, y - 0.42, lab, ha='center', va='top', fontsize=5.4, color=c if c != C_NUC else INK2)
    fr.panel_letter(ax, 'a')

    # b scale (LDLR library)
    axs = fig.add_subplot(gA[0, 1])
    fr.tidy(axs, 'y')
    cl = pd.read_csv(Path(a.root) / 'a9_statistics/extract/HepG2_run3/clones.tsv', sep='\t')
    snv = cl[cl.label != 'WT'].groupby('label').agg(molecules=('molecules', 'sum'), clones=('clone', 'size'))
    axs.hist(np.log10(snv.molecules), bins=30, color=C_FULL, alpha=0.9, lw=0)
    axs.set_xticks([2, 3, 4])
    axs.set_xticklabels(['100', '1,000', '10,000'])
    axs.set_xlabel('Molecules per SNV', fontsize=6)
    axs.set_ylabel('SNVs', fontsize=6)
    axs.text(0.02, 0.98, f'LDLR promoter library (HepG2)\n{len(snv)} of 633 SNVs\nmedian {int(snv.molecules.median()):,} molecules,\n'
             f'{int(snv.clones.median())} independent clones', transform=axs.transAxes, ha='left', va='top', fontsize=5.4,
             color=INK)
    fr.panel_letter(axs, 'b')
    snv.to_csv(out / 'panel_b_scale.tsv', sep='\t')

    # c locus
    gc = gC[0, 0].subgridspec(5, 1, height_ratios=[0.9, 1, 1, 0.6, 0.6], hspace=0.12)
    axg = fig.add_subplot(gc[0])
    axg.set_xlim(*WIN)
    axg.set_ylim(0, 3)
    axg.axis('off')
    for r in ref:
        if r['name2'] in ('DPAGT1', 'H2AX') and r['name'].startswith('NM_'):
            s0, e0 = max(r['txStart'], WIN[0]), min(r['txEnd'], WIN[1])
            if s0 >= e0:
                continue
            axg.plot([s0, e0], [2.3, 2.3], color=INK2, lw=0.7)
            for es, ee in zip(map(int, r['exonStarts'].strip(',').split(',')), map(int, r['exonEnds'].strip(',').split(','))):
                if ee > WIN[0] and es < WIN[1]:
                    axg.add_patch(Rectangle((max(es, WIN[0]), 2.05), min(ee, WIN[1]) - max(es, WIN[0]), 0.5, fc=INK2, lw=0))
            tss = r['txEnd'] if r['strand'] == '-' else r['txStart']
            if WIN[0] < tss < WIN[1]:
                axg.annotate('', xy=(tss - 250, 2.85), xytext=(tss, 2.85),
                             arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.8, shrinkA=0, shrinkB=0))
                axg.plot([tss, tss], [2.3, 2.85], color=INK, lw=0.8)
            axg.text((s0 + e0) / 2 if r['name2'] == 'H2AX' else 119101300, 2.7, r['name2'], fontsize=6, style='italic',
                     ha='center', va='bottom', color=INK)
            break
    rowy = {'DPAGT1': 1.35, 'DPAGT11': 0.55, 'DPAGT13': 0.55}
    lab = {'DPAGT1': '4.1-kb construct', 'DPAGT11': '~300-bp construct (FIRE peak 1)', 'DPAGT13': '~300-bp (FIRE peak 3)'}
    for c, rr in ins.iterrows():
        col = C_FULL if c == 'DPAGT1' else C_300
        axg.add_patch(Rectangle((rr.hg38_start, rowy[c] - 0.18), rr.hg38_end - rr.hg38_start, 0.36, fc=col, lw=0))
        if c == 'DPAGT1':
            axg.text(rr.hg38_start, rowy[c] + 0.3, lab[c], fontsize=5.4, va='bottom', ha='left', color=col)
        else:
            axg.text(rr.hg38_end + 40, rowy[c], lab[c], fontsize=5.4, va='center', color=col)
    for (s0, e0), name in ((CTCF1, 'CTCF'), ((119102288, 119102317), 'GC box')):
        axg.plot([s0, e0], [0.05, 0.05], color=INK, lw=3, solid_capstyle='butt')
        axg.text(e0 + 30, 0.05, name, fontsize=5.2, va='center', color=INK)
    fr.panel_letter(axg, 'c')
    tracks = [('endogenous', 'Endogenous\ngenome', C_ENDO), ('DPAGT1_full', '4.1-kb\nconstruct', C_FULL)]
    for i, (d, name, col) in enumerate(tracks):
        axt = fig.add_subplot(gc[1 + i], sharex=axg)
        fr.tidy(axt, 'y')
        p = P[d][(P[d].index >= WIN[0]) & (P[d].index <= WIN[1])]
        y_ = smooth_ratio(p).where(p.molecules > 0)
        axt.fill_between(p.index, 0, y_, color=col, alpha=0.75, lw=0)
        axt.set_ylabel(name, fontsize=5.6, rotation=0, ha='right', va='center', labelpad=4)
        axt.set_ylim(0, None)
        axt.tick_params(labelbottom=False, labelsize=5)
        axt.text(1.0, 0.95, f'n={int(p.molecules.max())}', transform=axt.transAxes, ha='right', va='top', fontsize=5,
                 color=INK2)
    for i, (d, name) in enumerate((('DPAGT11', 'peak 1'), ('DPAGT13', 'peak 3'))):
        pass
    axt = fig.add_subplot(gc[3], sharex=axg)
    fr.tidy(axt, 'y')
    for d in ('DPAGT11', 'DPAGT13'):
        p = P[d][(P[d].molecules > 0)]
        axt.fill_between(p.index, 0, smooth_ratio(p), color=C_300, alpha=0.8, lw=0)
    axt.set_ylabel('~300-bp\nconstructs', fontsize=5.6, rotation=0, ha='right', va='center', labelpad=4)
    axt.tick_params(labelbottom=False, labelsize=5)
    axt.text(1.0, 0.95, f'n={len(M["DPAGT11"]):,} / {len(M["DPAGT13"]):,}', transform=axt.transAxes, ha='right',
             va='top', fontsize=5, color=INK2)
    axx = fig.add_subplot(gc[4], sharex=axg)
    axx.axis('off')
    ticks = np.arange(119100500, WIN[1] + 1, 1000)
    for t in ticks:
        axx.text(t, 0.8, f'{t:,}', ha='center', va='top', fontsize=5, color=INK2)
    axx.text(WIN[0], 0.8, 'chr11', ha='left', va='top', fontsize=5, color=INK2)
    axx.text(0.5, -0.1, 'm6A fraction at A/T bases (accessibility), 51-bp smoothed; each track on its own scale',
             transform=axx.transAxes, ha='center', va='top', fontsize=5.3, color=INK2)

    # d peak 1 zoom
    axz = fig.add_subplot(gD[0, 0])
    fr.tidy(axz, 'both')
    zz = []
    for d, col, lab_ in (('endogenous', C_ENDO, 'endogenous genome'), ('DPAGT11', C_300, '~300-bp construct')):
        p = P[d].loc[PEAK1[0]:PEAK1[1]]
        y_ = smooth_ratio(p).where(p.molecules > 0)
        y_ = y_ / y_.mean()
        axz.plot((p.index - PEAK1[0]), y_, color=col, lw=1.2, label=lab_)
        zz.append(y_.rename(d))
    axz.axvspan(CTCF1[0] - PEAK1[0], CTCF1[1] - PEAK1[0], color='#eceae4', lw=0, zorder=0)
    axz.text((CTCF1[0] + CTCF1[1]) / 2 - PEAK1[0], 0.02, 'CTCF', transform=axz.get_xaxis_transform(), ha='center',
             va='bottom', fontsize=5, color=INK2)
    zz = pd.concat(zz, axis=1).dropna()
    from scipy import stats
    rho = stats.spearmanr(zz.endogenous, zz.DPAGT11).correlation
    axz.text(0.02, 0.98, f'FIRE peak 1 (300 bp)\nSpearman ρ = {rho:.2f}', transform=axz.transAxes, ha='left', va='top',
             fontsize=5.6, color=INK)
    axz.set_xlabel('Position in insert (bp)', fontsize=6)
    axz.set_ylabel('Relative accessibility\n(m6A / window mean)', fontsize=6)
    axz.legend(frameon=False, fontsize=5.4, loc='upper right', bbox_to_anchor=(1.0, 0.86))
    fr.panel_letter(axz, 'd')
    zz.to_csv(out / 'panel_d_peak1.tsv', sep='\t')

    # e single molecules
    sub = gD[0, 1].subgridspec(2, 1, height_ratios=[0.8, 1.2], hspace=0.12)
    rng = np.random.default_rng(SEED)
    for i, (d, name) in enumerate((('endogenous', 'Endogenous'), ('DPAGT11', '~300-bp construct'))):
        axp = fig.add_subplot(sub[i])
        m = M[d]
        if d == 'endogenous':
            rows = [r for r in m.itertuples() if min(json.loads(r.at), default=1e12) <= PEAK1[0]
                    and max(json.loads(r.at), default=0) >= PEAK1[1]]
        else:
            rows = list(m.itertuples())
            rows = [rows[j] for j in rng.choice(len(rows), NPILE, replace=False)]
        for yy, r in enumerate(rows):
            axp.plot([0, PEAK1[1] - PEAK1[0]], [yy, yy], color='#ecebe7', lw=0.5, zorder=0)
            for s0, e0, ln in json.loads(r.fp):
                c = footprint_class(ln)
                if c is None or e0 < PEAK1[0] or s0 > PEAK1[1]:
                    continue
                axp.add_patch(Rectangle((max(s0, PEAK1[0]) - PEAK1[0], yy - 0.4), min(e0, PEAK1[1]) - max(s0, PEAK1[0]),
                                        0.8, fc={'tf': C_TF, 'pic': C_PIC, 'nuc': C_NUC}[c], lw=0, zorder=1))
            mm = np.array(json.loads(r.m6a))
            mm = mm[(mm >= PEAK1[0]) & (mm <= PEAK1[1])] - PEAK1[0]
            axp.scatter(mm, np.full(mm.size, yy), s=1.2, color=INK, lw=0, zorder=2)
        axp.axvspan(CTCF1[0] - PEAK1[0], CTCF1[1] - PEAK1[0], color='#fbe3d9', lw=0, zorder=-1)
        axp.set_xlim(0, PEAK1[1] - PEAK1[0])
        axp.set_ylim(-0.8, len(rows) - 0.2)
        axp.set_yticks([])
        axp.set_ylabel(f'{name}\n{len(rows)} molecules', fontsize=5.6)
        for s_ in ('top', 'right', 'left'):
            axp.spines[s_].set_visible(False)
        if i == 0:
            axp.tick_params(labelbottom=False)
            fr.panel_letter(axp, 'e')
        else:
            axp.set_xlabel('Position in FIRE peak 1 (bp)', fontsize=6)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    axp.legend([Line2D([], [], marker='o', ls='none', color=INK, ms=2), Patch(color=C_TF), Patch(color=C_PIC),
                Patch(color=C_NUC)], ['m6A', 'TF-size', 'PIC-size', 'nucleosome'], loc='upper center',
               bbox_to_anchor=(0.5, -0.32), ncol=4, frameon=False, fontsize=5.2, handletextpad=0.2, columnspacing=0.6)

    # f footprint sizes
    axf = fig.add_subplot(gD[0, 2])
    fr.tidy(axf, 'y')
    lens = []
    for d in ('DPAGT11', 'DPAGT13'):
        lo, hi = ins.loc[d, 'hg38_start'], ins.loc[d, 'hg38_end']
        for r in M[d].itertuples():
            lens += [ln for s0, e0, ln in json.loads(r.fp) if s0 >= lo and e0 <= hi or (s0 < hi and e0 > lo and ln >= 85)]
    lens = np.array(lens)
    n_short = int((lens < 10).sum())
    lens = lens[lens >= 10]                     # < 10-bp pieces are not footprints of any class (fiberhmm fragments)
    bins = np.arange(10, 260, 5)
    hist, edges = np.histogram(lens, bins=bins)
    cols = [C_TF if 10 <= e < 41 else C_PIC if 41 <= e < 81 else C_NUC if e >= 85 else '#d6d4ce' for e in edges[:-1]]
    axf.bar(edges[:-1], hist / hist.sum() * 100, width=4.6, align='edge', color=cols, lw=0)
    axf.set_xlabel('Footprint length (bp)', fontsize=6)
    axf.set_ylabel('% of footprints', fontsize=6)
    for x_, t_, c_ in ((25, 'TF', C_TF), (60, 'PIC', C_PIC), (150, 'nucleosome', C_NUC)):
        axf.text(x_, 0.97, t_, transform=axf.get_xaxis_transform(), ha='center', va='top', fontsize=5.4,
                 color=c_ if c_ != C_NUC else INK2)
    axf.text(0.98, 0.80, f'{len(lens):,} footprints ≥10 bp,\n~300-bp constructs\n({n_short:,} shorter pieces\nnot shown)',
             transform=axf.transAxes, ha='right', va='top', fontsize=5.2, color=INK2)
    fr.panel_letter(axf, 'f')
    pd.Series(lens, name='length').to_csv(out / 'panel_f_lengths.tsv', sep='\t', index=False)

    fig.text(0.985, 0.995, f'draft {date.today()} · git {commit}+uncommitted · HEK293T (plasmid cell line to confirm)',
             ha='right', va='top', fontsize=5, color=MUTED)
    fig.savefig(out / 'fig1.pdf')
    fig.savefig(out / 'fig1.png', dpi=300)
    json.dump(dict(peak1_rho=float(rho), footprints=int(len(lens)), footprints_under_10bp=n_short,
                   endo_peak1_molecules=sum(1 for r in M['endogenous'].itertuples()
                                            if min(json.loads(r.at), default=1e12) <= PEAK1[0]
                                            and max(json.loads(r.at), default=0) >= PEAK1[1])),
              open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

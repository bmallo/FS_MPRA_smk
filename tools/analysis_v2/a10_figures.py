"""Step A10 slide figures from the a10_effects tables (no recomputation). Cells are drawn separately.

Usage: python tools/analysis_v2/a10_figures.py --dir results/analysis_v2_20261028/a10_effects [--example 61A>T]
The example SNV defaults to the largest max |z| in the +/-5 bp run, E1, group 'rest' (Simes p ties at 0).
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('a7f', Path(__file__).with_name('a7_figures.py'))
a7f = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a7f)
INK, BODY, BLUE, ORANGE, GREY = a7f.INK, a7f.BODY, a7f.BLUE, a7f.ORANGE, a7f.GREY
RUNS = [('HepG2_run3', 'A7 mask'), ('HepG2_run3_pm5', '±5 bp mask')]
GROUPS = [('rest', 'gains an A/T or neither'), ('drops its base', 'drops its base')]
BANDS = [(0, 5), (6, 10), (11, 20), (21, 50), (51, 300)]


def load(d, run):
    e = pd.read_parquet(d / run / 'effects.parquet')
    e['dist'] = (e.oligo - e.snv_oligo).abs()
    return e, pd.read_csv(d / run / 'snvs.tsv', sep='\t')


def hitmap(d, out):
    fig, ax = plt.subplots(2, 2, figsize=(15, 13), gridspec_kw=dict(hspace=0.3, wspace=0.2))
    for j, (run, rl) in enumerate(RUNS):
        e, _ = load(d, run)
        for i, (grp, gl) in enumerate(GROUPS):
            a = ax[i, j]
            h = e[(e.endpoint == 'E1') & e.hit & (e.group == grp)]
            for sign, c in ((1, ORANGE), (-1, BLUE)):
                s = h[np.sign(h.D) == sign]
                a.scatter(s.snv_oligo, s.oligo, s=np.clip(np.abs(s.D) * 300, 2, 40), color=c, alpha=0.6, lw=0,
                          label='more unmarked' if sign > 0 else 'less unmarked', rasterized=True)
            a.plot([21, 231], [21, 231], color=GREY, lw=0.8)
            a.set_xlim(20, 232)
            a.set_ylim(20, 232)
            a.set_aspect('equal')
            a.spines[['top', 'right']].set_visible(False)
            n = h.snv.nunique()
            a.set_title(f'{rl} · {gl}: {len(h)} hits in {n} SNVs', loc='left', color=INK, fontsize=13)
            a.set_xlabel('SNV position (oligo)')
            a.set_ylabel('adenine with a significant change')
            if i == 0 and j == 0:
                a.legend(frameon=False, fontsize=11, loc='upper left', markerscale=1.5)
    fig.suptitle('HepG2 run 3 · adenine unmarked (E1) · hierarchical FDR 5%', x=0.06, ha='left', color=INK, fontsize=15)
    fig.savefig(out / 'a10_hitmap.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def distance(d, out):
    fig, ax = plt.subplots(1, 2, figsize=(15, 4.6), gridspec_kw=dict(wspace=0.25))
    labels = [f'{lo}-{hi}' if hi < 300 else f'>{lo - 1}' for lo, hi in BANDS]
    for k, (ep, el) in enumerate([('E1', 'adenine unmarked'), ('E2', '20-40 bp call coverage')]):
        a = ax[k]
        for j, (run, rl) in enumerate(RUNS):
            e, _ = load(d, run)
            g = e[(e.endpoint == ep) & e.testable]
            if ep == 'E2':
                g = g[g['bin'] == 'TF']
            share = [g[(g.dist >= lo) & (g.dist <= hi)].hit.mean() * 100 for lo, hi in BANDS]
            a.bar(np.arange(len(BANDS)) + (j - 0.5) * 0.38, share, width=0.36, color=[BLUE, ORANGE][j], label=rl)
        a.set_xticks(range(len(BANDS)), labels)
        a.set_xlabel('distance from the SNV (bp)')
        a.set_ylabel('tested features that are hits (%)')
        a.set_title(f'HepG2 run 3 · {el}', loc='left', color=INK, fontsize=13)
        a.spines[['top', 'right']].set_visible(False)
        if k == 0:
            a.legend(frameon=False, fontsize=11)
    fig.savefig(out / 'a10_distance.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def manhattan(d, out):
    fig, ax = plt.subplots(2, 1, figsize=(16, 7.6), sharex=True, gridspec_kw=dict(hspace=0.35))
    for i, (run, rl) in enumerate(RUNS):
        _, s = load(d, run)
        s = s[s.endpoint == 'E1']
        a = ax[i]
        for (grp, gl), c in zip(GROUPS, (INK, '#A7AFBA')):
            g = s[s.group == grp]
            a.scatter(g.snv_oligo, -np.log10(g.simes_p.clip(lower=1e-30)), s=12, color=c, label=gl, alpha=0.85)
        a.set_ylabel('-log10 Simes p')
        a.set_title(f'HepG2 run 3 · E1 · {rl}: {int(s.selected.sum())} of {len(s)} SNVs selected', loc='right',
                    color=INK, fontsize=13)
        a.spines[['top', 'right']].set_visible(False)
        a.set_ylim(0, 31)
        if i == 0:
            a.legend(frameon=False, fontsize=11, ncol=2, loc='lower left', bbox_to_anchor=(0, 1.08))
    ax[-1].set_xlabel('SNV position (oligo) · values above 30 drawn at 30')
    fig.savefig(out / 'a10_snvs.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def example(d, out, snv):
    e, s = load(d, 'HepG2_run3_pm5')
    if snv is None:
        r = s[(s.endpoint == 'E1') & (s.group == 'rest')].sort_values('max_abs_z', ascending=False).iloc[0]
        snv = r.snv
    g = e[e.snv == snv]
    fig, ax = plt.subplots(2, 1, figsize=(16, 7.8), sharex=True, gridspec_kw=dict(hspace=0.35))
    for a, (ep, bn, yl) in zip(ax, [('E1', 'adenine', 'fraction unmarked\n(SNV - WT, points)'),
                                    ('E2', 'TF', '20-40 bp call coverage\n(SNV - WT, points)')]):
        h = g[(g.endpoint == ep) & (g['bin'] == bn)].sort_values('oligo')
        t = h[h.testable]
        a.axhline(0, color=GREY, lw=0.8)
        a.errorbar(t.oligo, 100 * t.D, yerr=196 * t.se, fmt='o', ms=3.5, color=BODY, ecolor='#C9CED6', lw=1)
        hh = t[t.hit]
        a.scatter(hh.oligo, 100 * hh.D, s=40, color=ORANGE, zorder=4, label='hit (FDR 5%)')
        m = h[h.masked]
        a.scatter(m.oligo, np.zeros(len(m)), marker='x', color=GREY, s=30, label='masked')
        a.axvline(g.snv_oligo.iloc[0], color=BLUE, lw=1.2, ls='--')
        a.set_ylabel(yl)
        a.spines[['top', 'right']].set_visible(False)
        a.set_title(f'{snv} · {"adenine unmarked" if ep == "E1" else "20-40 bp calls"} · ±5 bp mask', loc='left',
                    color=INK, fontsize=13)
    ax[0].legend(frameon=False, fontsize=11, loc='upper right')
    ax[-1].set_xlabel('oligo position (dashed: the SNV)')
    fig.savefig(out / 'a10_example.png', dpi=100, bbox_inches='tight')
    plt.close(fig)
    return snv


def hek(d, out):
    t = pd.read_csv(d / 'HEK293T_run1' / 'motifs.tsv', sep='\t')
    t = t[t.cm_p.notna()]
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.4), sharey=True, gridspec_kw=dict(wspace=0.08))
    names = list(dict.fromkeys(t.motif))
    for a, (ep, el) in zip(ax, [('E1', 'adenine unmarked'), ('E2', '10-80 bp call coverage')]):
        for k, (grp, gl) in enumerate(GROUPS):
            g = t[(t.endpoint == ep) & (t.group == grp)].set_index('motif')
            y = np.array([names.index(m) for m in g.index]) + (k - 0.5) * 0.3
            c = [BLUE, ORANGE][k]
            a.errorbar(100 * g.cm_effect, y, xerr=196 * g.cm_null_sd, fmt='o', color=c, ecolor=c, ms=7, lw=1.5,
                       label=gl)
            for yy, (m, r) in zip(y, g.iterrows()):
                if r.cm_bh:
                    a.text(100 * r.cm_effect, yy - 0.06, '*', ha='center', va='bottom', color=INK, fontsize=16)
        a.axvline(0, color=GREY, lw=0.8)
        a.set_yticks(range(len(names)), names)
        a.set_ylim(len(names) - 0.4, -0.75)
        a.set_xlabel('SNV clones - WT clones (points)')
        a.set_title(f'HEK293T · {el}', loc='left', color=INK, fontsize=13, pad=14)
        a.spines[['top', 'right']].set_visible(False)
    ax[0].legend(frameon=False, fontsize=11, loc='upper left', bbox_to_anchor=(0, -0.16), ncol=2)
    fig.text(0.9, -0.06, 'bars: ±1.96 SD of the WT-clone null · * BH 5%', ha='right', color=BODY, fontsize=11)
    fig.savefig(out / 'a10_hek_motifs.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--example')
    a = ap.parse_args()
    d = Path(a.dir)
    out = d / 'figures'
    out.mkdir(exist_ok=False)
    a7f.a4f.style()
    hitmap(d, out)
    distance(d, out)
    manhattan(d, out)
    print('example', example(d, out, a.example))
    hek(d, out)


if __name__ == '__main__':
    main()

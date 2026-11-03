"""Single-variant breakdown figure (HepG2 run 3 unless stated). Descriptive; reuses stored tests.

a  sequence around the SNV: A/T bases (visible to Hia5) vs G/C, motif matches covering the SNV (P < 1e-4, best per
   family), masked adenines (A7), WT -> variant motif P.
b  HepG2 MPRA expression change of every SNV in the window, ClinVar class by marker, the SNV highlighted.
c  single molecules: WT and variant, 40 each, labeling 0.2-0.6, random sample (seed), accessible adenines and LLR calls
   (unmasked calls from fp_features/calls.npz; 10-40 bp orange, darker if confident; >= 81 bp grey).
d  change in 10-19 and 20-40 bp call coverage per position (A10 E2, A7 mask; 95% CI; filled = A10 hit).
e  change in inaccessible fraction per adenine (A10 E1, A7 mask filled, +/-5 bp mask open; 95% CI).
f  share of molecules with a 10-40 bp call covering the motif core, WT vs variant, by labeling (all and confident
   calls; unmasked).  A real protein footprint should not shrink with labeling; a sequence artefact often does.
g  summary of every evidence type (evidence table, classify, HEK293T, AlphaGenome, ClinVar/LOVD, gnomAD).
Usage: python tools/r01_figure/fig_variant.py --snv 129C>T --root results/analysis_v2_20261028 --work results/r01_figure_20261102 --output <dir>
"""
import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).parent


def _load(name):
    s = importlib.util.spec_from_file_location(name, HERE / f'{name}.py')
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m
    s.loader.exec_module(m)
    return m


fr = _load('fig_repeat3')
fa = _load('fig_A')
INK, INK2, MUTED, GRID = fr.INK, fr.INK2, fr.MUTED, fr.GRID
ORANGE, BLUE, VIOLET = '#eb6834', '#2a78d6', '#6250d6'
FAMILY = re.compile(r'^(SP[1-4]|KLF\d*|SREBF[12]|SRBP[12]|CREB\d|ATF\d|NFY[ABC])$')
SEED, NPILE, WIN = 20261102, 40, 30
C0 = 285


def motif_p(ref, snv, db):
    ms, cl = _load('motif_scan'), _load('clinical')
    k, alt = int(snv[:-3]), snv[-1]
    seq = cl.edit(ref, k, k, alt)
    bg = np.array([ref.count(b) for b in ms.BASES], float)
    bg = (bg + bg[::-1]) / (2 * len(ref))
    fam = {json.loads(x)['name']: json.loads(x)['masterlist_info']['tf'] for x in open(db / ms.HOCO_ANN)}
    rows = []
    for m in ms.read_motifs(db / ms.JASPAR, 'JASPAR2026') + ms.read_motifs(db / ms.HOCO, 'HOCOMOCOv14'):
        tf = fam.get(m['id'], m['name']).upper()
        mn, sf, I = ms.pvalue_table(ms.log_odds(m['counts'], bg), bg)
        L = I.shape[0]
        best = []
        for sq in (ref, seq):
            pf, pr = ms.scan(sq, I, mn, sf)
            st = np.arange(len(pf)) + cl.LO
            cov = (st <= k) & (st + L - 1 >= k)
            i = np.argmin(np.minimum(pf, pr) + (~cov) * 9)
            best.append((min(pf[i], pr[i]), st[i], st[i] + L - 1))
        rows.append(dict(tf=tf, motif=m['id'], p_wt=best[0][0], lo=best[0][1], hi=best[0][2], p_var=best[1][0]))
    t = pd.DataFrame(rows)
    return t[(t.p_wt < 1e-4) | (t.p_var < 1e-4)].sort_values('p_wt')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--snv', required=True)
    ap.add_argument('--root', required=True)
    ap.add_argument('--work', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, w, out = Path(a.root), Path(a.work), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    from tools.analysis_v2 import access
    v = a.snv
    k = int(v[:-3])
    lo, hi = k - WIN, k + WIN
    cl = _load('clinical')
    ref, _ = cl.reference(root)
    ev = pd.read_csv(w / 'classify/variants.tsv', sep='\t').set_index('snv')
    r = ev.loc[v]
    cname = f'c.{int(r.c_dot)}{v[-3:]}'
    mot = motif_p(ref, v, w / 'motif_db')
    mot.to_csv(out / 'motifs.tsv', sep='\t', index=False, float_format='%.3g')
    cache = w / 'cache/a7_masks.json'                 # a7_masks rebuilds every SNV's mask (~minutes); cache once
    if not cache.exists():
        cache.parent.mkdir(exist_ok=True)
        json.dump({m: {s_: list(map(int, x)) for s_, x in access.a7_masks('HepG2', pm5=(m == 'pm5')).items()}
                   for m in ('A7', 'pm5')}, open(cache, 'w'))
    masks = {m: np.array(json.load(open(cache))[m][v]) for m in ('A7', 'pm5')}
    e = pd.read_parquet(root / 'a10_effects/HepG2_run3/effects.parquet',
                        columns=['snv', 'endpoint', 'oligo', 'bin', 'D', 'se', 'p', 'masked', 'hit'],
                        filters=[('snv', '==', v)])
    e5 = pd.read_parquet(root / 'a10_effects/HepG2_run3_pm5/effects.parquet',
                         columns=['snv', 'endpoint', 'oligo', 'bin', 'D', 'se', 'p', 'masked', 'hit'],
                         filters=[('snv', '==', v), ('endpoint', '==', 'E1')])

    d = access.load_element('HepG2')
    z = np.load(root / 'fp_features/calls.npz')
    rows = z['rows']
    lab = d.label[rows]
    L = d.distal[rows]
    conf = (np.log(1 - np.clip(L, .01, .99)) / np.log(1 - 0.70))[:, None] * z['S'] / np.log(10) < -3
    rng = np.random.default_rng(SEED)

    fig = plt.figure(figsize=(7.2, 9.6))
    gs = fig.add_gridspec(6, 2, height_ratios=[0.55, 0.9, 2.6, 0.95, 0.95, 1.15], width_ratios=[1.0, 0.62],
                          left=0.09, right=0.985, top=0.945, bottom=0.05, hspace=0.55, wspace=0.28)
    xl = (lo - C0 - 0.5, hi - C0 + 0.5)

    # a sequence + motifs
    ax = fig.add_subplot(gs[0, 0])
    for j in range(lo, hi + 1):
        b = ref[j - 21]
        col = INK if b in 'AT' else '#b9b7b1'
        ax.text(j - C0, 0.2, b, ha='center', va='center', fontsize=5.2, family='monospace', color=col,
                fontweight='bold' if b in 'AT' else 'normal')
    ax.text(k - C0, -0.55, v[-1], ha='center', va='center', fontsize=6, family='monospace', color=ORANGE,
            fontweight='bold')
    ax.add_patch(plt.Rectangle((k - C0 - 0.5, -0.35), 1, 1.1, fill=False, edgecolor=ORANGE, lw=0.9))
    for mk in masks['A7']:
        ax.plot(mk - C0, 0.85, marker='v', color=MUTED, ms=3)
    fams = mot[mot.tf.str.match(FAMILY)].drop_duplicates('tf').head(4)
    grp = fams.groupby(['lo', 'hi']).agg(tf=('tf', '/'.join), p_wt=('p_wt', 'min'), p_var=('p_var', 'min')).reset_index()
    grp = grp.sort_values('p_wt')
    for i, m in enumerate(grp.head(2).itertuples()):
        y = 1.3 + 0.7 * i
        ax.plot([m.lo - C0 - 0.5, m.hi - C0 + 0.5], [y, y], color=VIOLET, lw=2.2, solid_capstyle='butt')
        ax.text(m.hi - C0 + 1.5, y, f'{m.tf} motif: P {m.p_wt:.0e} (WT) → {m.p_var:.0e} (variant)', va='center',
                fontsize=5.8, color=INK)
    ax.set_xlim(*xl)
    ax.set_ylim(-0.9, 2.5)
    ax.axis('off')
    fig.text(0.09, 0.985, f'{cname} ({v}; oligo {k}) · ClinVar: {r.clinvar if isinstance(r.clinvar, str) else "not reported"}',
             fontsize=8, fontweight='bold', color=INK, va='top')
    ax.text(xl[1], -0.85, 'bold = A/T (Hia5-markable); ▾ = masked adenine (A7)', ha='right', va='bottom', fontsize=5.2,
            color=INK2)
    fr.panel_letter(ax, 'a')

    # b MPRA neighbourhood
    ax = fig.add_subplot(gs[1, 0])
    fr.tidy(ax)
    nb = ev[(ev.oligo >= lo) & (ev.oligo <= hi)].copy()
    nb['cls'] = nb.clinvar.map(fa.clin_class).fillna('none')
    fa.scatter_clin(ax, (nb.oligo - C0).values, nb.act.values, nb.cls.values, small=6)
    ax.scatter(k - C0, r.act, s=70, facecolor='none', edgecolor=INK, lw=1.0, zorder=6)
    ax.axhline(0, color=INK2, lw=0.5)
    ax.set_xlim(*xl)
    ax.set_ylabel('Expression\n(log2, MPRA)', fontsize=6.5)
    ax.tick_params(labelbottom=False)
    fr.panel_letter(ax, 'b')

    # c pileups
    sub = gs[2, 0].subgridspec(2, 1, hspace=0.12)
    pile = {}
    for i, (name, sel) in enumerate((('WT', 'WT'), (cname, v))):
        axp = fig.add_subplot(sub[i, 0])
        cand = np.flatnonzero((lab == sel) & (L >= 0.2) & (L < 0.6))
        pick = np.sort(rng.choice(cand, min(NPILE, cand.size), replace=False))
        pick = pick[np.argsort(-L[pick])]
        pile[name] = pick
        cols = (d.oligo >= lo) & (d.oligo <= hi)
        for y, ix in enumerate(pick):
            x = d.X[rows[ix]]
            axp.plot([lo - C0, hi - C0], [y, y], color='#ecebe7', lw=0.6, zorder=0)
            acc = d.oligo[cols & (x == 1)]
            axp.scatter(acc - C0, np.full(acc.size, y), s=3, color=INK, lw=0, zorder=3)
            for c in range(z['mid'].shape[1]):
                if z['mid'][ix, c] < 0:
                    continue
                s, en, ln = z['start'][ix, c], z['end'][ix, c], z['length'][ix, c]
                if en < lo or s > hi:
                    continue
                col = ('#b8461c' if conf[ix, c] else '#f6b89a') if 10 <= ln <= 40 else (
                    '#8a8883' if ln >= 81 else '#d6d4ce')
                axp.add_patch(plt.Rectangle((max(s, lo) - C0 - 0.5, y - 0.38), min(en, hi) - max(s, lo) + 1, 0.76,
                                            color=col, lw=0, zorder=2))
        axp.set_xlim(*xl)
        axp.set_ylim(-0.8, len(pick) - 0.2)
        axp.set_yticks([])
        axp.set_ylabel(f'{name}\n{len(pick)} molecules', fontsize=6)
        for s_ in ('top', 'right', 'left'):
            axp.spines[s_].set_visible(False)
        axp.axvline(k - C0, color=ORANGE, lw=0.5, ls=(0, (2, 2)), zorder=1)
        axp.tick_params(labelbottom=False)
        if i == 0:
            fr.panel_letter(axp, 'c')
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    axp.legend([Line2D([], [], marker='o', ls='none', color=INK, ms=2.5), Patch(color='#b8461c'),
                Patch(color='#f6b89a'), Patch(color='#8a8883')],
               ['accessible (marked) adenine', '10-40 bp call, confident', '10-40 bp call, other', '≥81 bp call'],
               loc='upper center', bbox_to_anchor=(0.5, -0.04), ncol=4, frameon=False, fontsize=5.5, handletextpad=0.3)

    # d TF-size coverage change
    ax = fig.add_subplot(gs[3, 0])
    fr.tidy(ax)
    e2 = e[(e.endpoint == 'E2') & e['bin'].isin(['sub_TF', 'TF']) & e.oligo.between(lo, hi) & ~e.masked]
    for b, col in (('sub_TF', '#f6b89a'), ('TF', '#b8461c')):
        g = e2[e2['bin'] == b].sort_values('oligo')
        x = g.oligo - C0
        ax.fill_between(x, 100 * (g.D - 1.96 * g.se), 100 * (g.D + 1.96 * g.se), color=col, alpha=0.25, lw=0)
        ax.plot(x, 100 * g.D, color=col, lw=1.1, label='10-19 bp' if b == 'sub_TF' else '20-40 bp')
        hh = g[g.hit]
        ax.scatter(hh.oligo - C0, 100 * hh.D, s=6, color=col, zorder=3)
    ax.axhline(0, color=INK2, lw=0.5)
    ax.axvline(k - C0, color=ORANGE, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(*xl)
    ax.set_ylabel('Footprint-size call\ncoverage change (pts)', fontsize=6.5)
    ax.legend(frameon=False, fontsize=5.6, loc='upper left', ncol=2)
    ax.tick_params(labelbottom=False)
    fr.panel_letter(ax, 'd')

    # e inaccessible fraction per adenine
    ax = fig.add_subplot(gs[4, 0])
    fr.tidy(ax)
    for frame, filled, lab_ in ((e[(e.endpoint == 'E1')], True, 'A7 mask'), (e5, False, '±5 bp mask')):
        g = frame[frame.oligo.between(lo, hi) & frame.p.notna()]
        x = g.oligo - C0 + (0 if filled else 0.25)
        ax.errorbar(x, 100 * g.D, yerr=196 * g.se, fmt='none', ecolor='#c9c7c1', lw=0.6)
        ax.scatter(x, 100 * g.D, s=12, facecolor=np.where(g.hit, BLUE, '#a9c6ef') if filled else 'white',
                   edgecolor=BLUE, lw=0.7, zorder=3, label=lab_)
    ax.axhline(0, color=INK2, lw=0.5)
    ax.axvline(k - C0, color=ORANGE, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(*xl)
    ax.set_ylabel('Inaccessible fraction\nchange (pts)', fontsize=6.5)
    ax.legend(frameon=False, fontsize=5.6, loc='upper left', ncol=2)
    ticks = [t_ for t_ in range(((lo - C0) // 10) * 10, hi - C0 + 1, 10) if xl[0] <= t_ <= xl[1]]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'c.{t_}' for t_ in ticks])
    ax.set_xlabel('Position (HGVS c., NM_000527.5)')
    fr.panel_letter(ax, 'e')

    # f labeling dependence (right column, top)
    axf = fig.add_subplot(gs[1:3, 1])
    fr.tidy(axf, 'y')
    core = (int(grp.lo.iloc[0]), int(grp.hi.iloc[0])) if len(grp) else (k - 5, k + 5)   # best motif match
    cover = ((z['start'] <= core[1]) & (z['end'] >= core[0]) & (z['length'] >= 10) & (z['length'] <= 40)
             & (z['mid'] >= 0))
    bins = [(0.2, 0.4), (0.4, 0.6), (0.6, 1.01)]
    ft = []
    for j, (b0, b1) in enumerate(bins):
        for name, sel, col, dx in (('WT', 'WT', '#9a9893', -0.18), (cname, v, ORANGE, 0.18)):
            m = (lab == sel) & (L >= b0) & (L < b1)
            for kind, cv, mk in (('all', cover, 'o'), ('confident', cover & conf, 's')):
                p = cv[m].any(1).mean() if m.any() else np.nan
                ft.append(dict(bin=f'{b0}-{min(b1, 1)}', group=name, calls=kind, n=int(m.sum()), share=p))
                axf.scatter(j + dx + (0.06 if kind == 'confident' else -0.06), 100 * p, s=18, marker=mk,
                            facecolor=col if kind == 'all' else 'white', edgecolor=col, lw=0.9, zorder=3)
    ft = pd.DataFrame(ft)
    ft.to_csv(out / 'panel_f_labeling.tsv', sep='\t', index=False, float_format='%.5g')
    axf.set_xticks(range(len(bins)))
    axf.set_xticklabels(['0.2–0.4', '0.4–0.6', '≥0.6'])
    axf.set_xlabel('Molecule labeling (m6A fraction)')
    axf.set_ylabel(f'Molecules with a 10–40 bp call over\nc.{core[0] - C0}..c.{core[1] - C0} (%)', fontsize=6.5)
    axf.legend([Line2D([], [], marker='o', ls='none', color='#9a9893'), Line2D([], [], marker='o', ls='none', color=ORANGE),
                Line2D([], [], marker='o', ls='none', color=INK, mfc=INK), Line2D([], [], marker='s', ls='none', color=INK, mfc='white')],
               ['WT', cname, 'all calls', 'confident calls'], frameon=False, fontsize=5.6, loc='upper right')
    fr.panel_letter(axf, 'f')

    # g summary card
    axg = fig.add_subplot(gs[3:6, 1])
    axg.axis('off')
    hep = pd.read_csv(root / 'a9_statistics/extract/HepG2_run3/clones.tsv', sep='\t')
    hep = hep[hep.label == v]
    best = mot[mot.tf.str.match(FAMILY)].head(1)
    lines = [
        ('Database', None),
        ('ClinVar', r.clinvar if isinstance(r.clinvar, str) else 'not reported'),
        ('LOVD', str(r.lovd_effect).split(';')[0] if isinstance(r.lovd_effect, str) else '–'),
        ('gnomAD v4.1', f'AC {int(r.gnomad_ac)}' if r.gnomad_ac else 'absent'),
        ('GTEx', 'no eQTL in element'),
        ('Prediction', None),
        ('AlphaGenome', f'LDLR RNA {r.ag_rna_seq_hepg2:+.3f}, DNase {r.ag_dnase_hepg2:+.2f}'),
        ('Motif', f'{best.tf.iloc[0]} P {best.p_wt.iloc[0]:.0e} → {best.p_var.iloc[0]:.0e}' if len(best) else 'none'),
        ('Function (HepG2)', None),
        ('MPRA', f'{r.act:+.2f} log2 (rep. {r.act2:+.2f})'),
        ('Molecules', f'{int(hep.molecules.sum())} ({len(hep)} clones)'),
        ('Footprint ≤10 bp', f'{"gain" if r.fp_sign > 0 else "loss" if r.fp_sign < 0 else "none"} ({r.fp_level})'),
        ('', f'z {r.fp_mean_z:+.1f}; {int(r.fp_hits)} sig. positions'),
        ('Inaccessible ≤10 bp', f'{r.local_inacc_pts:+.1f} pts'),
        ('≥81-bp protection', f'{r.nuc_D_pts:+.2f} pts (q {r.nuc_q:.2f})'),
        ('Element-wide shift', f'{r.C1:+.2f} pts'),
        ('HEK293T', None),
        ('Molecules', f'{int(r.hek_molecules)} ({int(r.hek_clones)} clones)'),
        ('Inaccessible ≤10 bp', f'{r.hek_local_pts:+.1f} pts' if np.isfinite(r.hek_local_pts) else '–'),
        ('Far (>50 bp)', f'{r.hek_far_pts:+.1f} pts' if np.isfinite(r.hek_far_pts) else '–'),
        ('Call', None),
        ('Category', r.category),
        ('Deciding evidence', r.evidence),
    ]
    y = 1.0
    for kk, val in lines:
        if val is None:
            y -= 0.012
            axg.text(0, y, kk, fontsize=6.6, fontweight='bold', color=INK, va='top', transform=axg.transAxes)
            y -= 0.048
            continue
        axg.text(0.0, y, kk, fontsize=5.9, color=INK2, va='top', transform=axg.transAxes)
        axg.text(0.42, y, val, fontsize=5.9, color=INK, va='top', transform=axg.transAxes, wrap=True)
        y -= 0.043
    fr.panel_letter(axg, 'g')

    fig.savefig(out / f'variant_{v.replace(">", "_")}.pdf')
    fig.savefig(out / f'variant_{v.replace(">", "_")}.png', dpi=300)
    json.dump(dict(snv=v, c=cname, pileup_wt=len(pile['WT']), pileup_var=len(pile[cname]), core=list(map(int, core)),
                   masks={m: list(map(int, x)) for m, x in masks.items()}), open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

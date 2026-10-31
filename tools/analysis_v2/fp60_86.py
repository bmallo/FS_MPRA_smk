"""Exploratory follow-up (after Phase A): is the 20-40 bp footprint at oligo ~60-86 (HepG2 run 3) real?

Parts (all from existing outputs and the A9 molecules; no settings change):
  adenines.tsv    each A/T at oligo 50-95: run class, naked relative marking (A4 lookup), run-3 emissions p_a/p_p,
                  evidence per unmarked adenine w = ln((1-p_p)/(1-p_a)), WT unmarked fraction
  anatomy.tsv     WT calls (labeling >= 0.2, sample) covering the core 70-78: start, end, length, LLR, opportunities,
                  as read and with the A run 59-65 (and also the TT 85-86) set unobserved
  labeling.tsv    20-40 bp coverage at oligo 73 by labeling bin: WT observed, and naked DNA simulated from this
                  sequence (marking = q * naked relative rate, q fitted per bin to WT marking elsewhere in the element);
                  the same for a comparison footprint position
  variants.tsv    every SNV at oligo 55-91: WT - variant change over the footprint (66-84) in 10-19/20-40/41-80/81+
                  coverage (A7-mask and +/-5 bp-mask tests), core unmarked change, hits, activity, caller class,
                  A7-predicted caller artifact of dropping its base
  summary.json, figures/

Usage: python tools/analysis_v2/fp60_86.py --root results/analysis_v2_20261028 --catalog <catalog.json> --output <new dir>
"""
import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).parent


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


a7 = _load('a7', 'a7_snv_edits.py')
a4v = _load('a4v', 'a4_variant_visibility.py')
a7f = _load('a7f', 'a7_figures.py')
a6x, kadane = a7.a6x, a7.kadane
INK, BODY, BLUE, ORANGE, GREY, PURPLE = a7f.INK, a7f.BODY, a7f.BLUE, a7f.ORANGE, a7f.GREY, a7f.PURPLE
TEAL = '#23877C'

EL0, OLIGO0 = 1234, 21
WIN = (50, 95)
CORE = (70, 78)
FOOT = (66, 84)
ARUN = (59, 65)
TT = (85, 86)
PROBE = 73
TF = (20, 40)
LAB_BINS = [(0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 0.6), (0.6, 0.7), (0.7, 0.8), (0.8, 1.01)]
N_WT, N_SIM = 40000, 20000
SEED = 20261030
BINS = {'sub_TF': (10, 19), 'TF': (20, 40), 'PIC': (41, 80), 'NUC': (81, 100000)}


def ref(o):
    return o - OLIGO0 + EL0


def olig(r):
    return r - EL0 + OLIGO0


def calls(x, hit, miss, pos):
    """Calls of one molecule: (start ref, end ref exclusive, llr, opportunities)."""
    ii = x >= 0
    if not ii.any():
        return []
    step = np.where(x == 1, hit, np.where(x == 0, miss, 0.0))
    o = np.flatnonzero(ii)
    s, e, llr, nop, *_ = kadane.scan(step, x == 1, ii, o[0], o[-1] + 1, a7.MIN_LLR, a7.MIN_OPPS, a7.DRAWDOWN)
    return [(int(pos[a]), int(pos[b - 1]) + 1, float(l), int(n)) for a, b, l, n in zip(s, e, llr, nop)]


def covering(cl, r):
    return [c for c in cl if c[0] <= r < c[1]]


def anatomy(Xs, dist, hit, miss, pos):
    """Per WT molecule: the call covering the probe, as read / A run hidden / A run + TT hidden."""
    cols_arun = np.flatnonzero((pos >= ref(ARUN[0])) & (pos <= ref(ARUN[1])))
    cols_tt = np.flatnonzero((pos >= ref(TT[0])) & (pos <= ref(TT[1])))
    rows = []
    for x, d in zip(Xs, dist):
        out = dict(labeling=d)
        for tag, hide in (('read', None), ('no_arun', cols_arun), ('no_arun_tt', np.r_[cols_arun, cols_tt])):
            y = x
            if hide is not None:
                y = x.copy()
                y[hide] = -1
            c = covering(calls(y, hit, miss, pos), ref(PROBE))
            if c:
                s, e, l, n = c[0]
                out.update({f'{tag}_start': olig(s), f'{tag}_end': olig(e - 1), f'{tag}_len': e - s, f'{tag}_llr': l,
                            f'{tag}_opps': n})
        rows.append(out)
    return pd.DataFrame(rows)


def naked_rates(seq_ref, pos, lookup):
    """Relative naked marking per reference A/T column from its (run length, position in run) class, either strand."""
    cls = a4v.adenine_classes(seq_ref)
    rel = np.ones(pos.size)
    for c, p in enumerate(pos):
        for strand in ('top', 'bottom'):
            k = cls.get((strand, int(p)))
            if k is not None:
                rel[c] = lookup.get(f'{k[0]}|{k[1]}', 1.0)
    return rel


def simulate(rel, hit, miss, pos, target, n, rng, probes):
    """Naked molecules over the columns given: marking min(1, q * rel), q matched to the target mean marking.
    Returns 20-40 bp coverage at each probe position."""
    lo, hi = 0.0, 5.0
    for _ in range(50):
        q = (lo + hi) / 2
        lo, hi = (q, hi) if np.minimum(1, q * rel).mean() < target else (lo, q)
    p = np.minimum(1, q * rel)
    cov = np.zeros(len(probes))
    for _ in range(n):
        x = (rng.random(p.size) < p).astype(np.int8)
        cl = calls(x, hit, miss, pos)
        for k, pr in enumerate(probes):
            c = covering(cl, ref(pr))
            if c and TF[0] <= c[0][1] - c[0][0] <= TF[1]:
                cov[k] += 1
    return cov / n, q


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'figures').mkdir()
    rng = np.random.default_rng(SEED)
    ext = root / 'a9_statistics/extract/HepG2_run3'
    cfg = json.load(open(root / 'a6_proxies/extract/HepG2_run3/extract.json'))
    seq_ref = a6x.a5x.ref_seq(cfg['reference'])
    pos, X, distal, _, _ = a6x.load(ext)
    meta = np.load(ext / 'meta.npz')
    cls, label = meta['cls'], meta['label']
    em = 'results/llr_vs_fiberhmm_20261026/results/emissions.tsv.gz'
    hit, miss = a7.emission_steps(em, pos)
    t = pd.read_csv(em, sep='\t')
    t = t[t.fold == 'A']
    pa, pp = dict(zip(t.ref_pos - 1, t.p_a)), dict(zip(t.ref_pos - 1, t.p_p))
    lookup = pd.read_csv(root / 'a4_informative/visibility/naked_lookup.tsv', sep='\t', index_col=0).p_rel.to_dict()
    rel = naked_rates(seq_ref, pos, lookup)
    wt = np.flatnonzero(cls == 'WT')
    summ = {}

    # 1. adenines
    rows = []
    for c in np.flatnonzero((pos >= ref(WIN[0])) & (pos <= ref(WIN[1]))):
        x = X[wt, c]
        k1 = distal[wt] > 0.01
        km = distal[wt] >= 0.6
        obs = (x >= 0)
        rows.append(dict(oligo=olig(int(pos[c])), base=seq_ref[pos[c]], naked_rel=round(rel[c], 3), p_a=round(pa[pos[c]], 3),
                         p_p=round(pp[pos[c]], 4), w_unmarked=round(miss[c], 2), w_marked=round(hit[c], 2),
                         wt_unmarked=round(((x == 0) & k1).sum() / max((obs & k1).sum(), 1), 4),
                         wt_unmarked_M=round(((x == 0) & km).sum() / max((obs & km).sum(), 1), 4)))
    ad = pd.DataFrame(rows)
    ad.to_csv(out / 'adenines.tsv', sep='\t', index=False)
    core = ad[ad.oligo.between(*FOOT) | ad.oligo.between(*TT)]
    summ['adenines'] = dict(window=WIN, at_in_60_86=int(ad.oligo.between(60, 86).sum()),
                            informative_in_60_86=int((ad.oligo.between(60, 86) & (ad.w_unmarked >= 0.5)).sum()),
                            evidence_60_86=round(float(ad[ad.oligo.between(60, 86)].w_unmarked.sum()), 2),
                            evidence_65_86=round(float(ad[ad.oligo.between(65, 86)].w_unmarked.sum()), 2),
                            at_free_gaps=[[66, 69], [79, 84]])
    print('adenines', summ['adenines'], flush=True)

    # 2. anatomy of WT calls covering the probe
    k2 = wt[distal[wt] >= 0.2]
    samp = rng.choice(k2, min(N_WT, k2.size), replace=False)
    an = anatomy(X[samp], distal[samp], hit, miss, pos)
    an.to_csv(out / 'anatomy.tsv.gz', sep='\t', index=False, compression='gzip')
    an_s = {}
    for tag in ('read', 'no_arun', 'no_arun_tt'):
        L = an[f'{tag}_len'].dropna()
        an_s[tag] = dict(molecules=len(an), covered=round(len(L) / len(an), 4),
                         tf_size=round(float(L.between(*TF).sum() / len(an)), 4),
                         sub_tf=round(float(L.between(10, 19).sum() / len(an)), 4),
                         pic=round(float(L.between(41, 80).sum() / len(an)), 4), nuc=round(float((L >= 81).sum() / len(an)), 4),
                         median_len=float(L.median()) if len(L) else None)
        if tag == 'read':
            tfm = an[an.read_len.between(*TF)]
            an_s[tag].update(start_in_arun=round(float(tfm.read_start.between(*ARUN).mean()), 3),
                             start_at_or_after_70=round(float((tfm.read_start >= 70).mean()), 3),
                             end_at_85_86=round(float(tfm.read_end.between(*TT).mean()), 3),
                             end_at_or_before_78=round(float((tfm.read_end <= 78).mean()), 3),
                             start_mode=int(tfm.read_start.mode()[0]), end_mode=int(tfm.read_end.mode()[0]),
                             llr_median=round(float(tfm.read_llr.median()), 2), opps_median=float(tfm.read_opps.median()))
    summ['anatomy'] = an_s
    print('anatomy', json.dumps(an_s), flush=True)

    # 3. labeling dependence, observed vs simulated naked
    elcols = np.flatnonzero((pos >= EL0) & (pos < EL0 + 211) & ~((pos >= ref(WIN[0])) & (pos <= ref(WIN[1]))))
    wcols = np.flatnonzero((pos >= EL0 - 150) & (pos < EL0 + 211 + 150))
    base = json.load(open(root / 'explorer_v3/base.json'))['HepG2']['coverage']['TF']
    cand = [(v, OLIGO0 + i) for i, v in enumerate(base) if v is not None and not (WIN[0] - 10 <= OLIGO0 + i <= WIN[1] + 10)]
    comp = max(cand)[1]
    probes = [PROBE, comp]
    lrows = []
    for lo, hi in LAB_BINS:
        m = k2[(distal[k2] >= lo) & (distal[k2] < hi)]
        if m.size < 200:
            continue
        mm = rng.choice(m, min(m.size, 6000), replace=False)
        obs = np.zeros(len(probes))
        for x in X[mm]:
            cl = calls(x, hit, miss, pos)
            for k, pr in enumerate(probes):
                c = covering(cl, ref(pr))
                if c and TF[0] <= c[0][1] - c[0][0] <= TF[1]:
                    obs[k] += 1
        obs /= mm.size
        xe = X[mm][:, elcols]
        target = float((xe == 1).sum() / max((xe >= 0).sum(), 1))
        sim, q = simulate(rel[wcols], hit[wcols], miss[wcols], pos[wcols], target, N_SIM // len(LAB_BINS), rng, probes)
        for k, pr in enumerate(probes):
            lrows.append(dict(labeling=f'{lo:.1f}-{min(hi, 1):.1f}', probe=pr, wt_molecules=int(m.size), element_marking=round(target, 3),
                              q=round(q, 3), observed=round(float(obs[k]), 4), simulated_naked=round(float(sim[k]), 4)))
        print('labeling', lo, hi, [r for r in lrows[-2:]], flush=True)
    lab = pd.DataFrame(lrows)
    lab.to_csv(out / 'labeling.tsv', sep='\t', index=False)
    summ['labeling'] = dict(comparison_position=comp, rows=lab.to_dict('records'))

    # 4. variants across the footprint
    e7 = pd.read_parquet(root / 'a10_effects/HepG2_run3/effects.parquet')
    e5 = pd.read_parquet(root / 'a10_effects/HepG2_run3_pm5/effects.parquet')
    fn = pd.read_csv(root / 'a11_motifs_function/function_snvs.tsv', sep='\t').set_index('snv')
    v7 = pd.read_csv(root / 'a7_snv_edits/HepG2_run3/variants.tsv', sep='\t').set_index('source_variant')
    art = np.load(root / 'a7_snv_edits/HepG2_run3/artifacts.npz')
    snv = e7[(e7.endpoint == 'E1') & e7.snv_oligo.between(55, 91)].drop_duplicates('snv')[['snv', 'snv_oligo', 'caller_class']]
    vrows = []
    foot = slice(FOOT[0] - OLIGO0, FOOT[1] - OLIGO0 + 1)
    for v, o, cc in snv.itertuples(index=False):
        row = dict(snv=v, oligo=o, caller_class=cc, act=fn.act.get(v), act2=fn.act2.get(v),
                   zone=('A run 59-65' if ARUN[0] <= o <= ARUN[1] else 'gap 66-69' if 66 <= o <= 69 else 'core 70-78' if 70 <= o <= 78
                         else 'gap 79-84' if 79 <= o <= 84 else 'TT 85-86' if 85 <= o <= 86 else 'flank'))
        for tag, e in (('A7', e7), ('pm5', e5)):
            g = e[(e.snv == v) & (e.endpoint == 'E2') & e.oligo.between(*FOOT)]
            for b in BINS:
                h = g[(g['bin'] == b) & g.testable]
                row[f'{b}_{tag}'] = -100 * h.D.mean() if len(h) else np.nan      # WT - variant, points
                row[f'{b}_{tag}_hits'] = int(h.hit.sum())
            g = e[(e.snv == v) & (e.endpoint == 'E1') & e.oligo.between(*CORE) & e.testable]
            row[f'unmarked_core_{tag}'] = -100 * g.D.mean() if len(g) else np.nan
            row[f'unmarked_core_{tag}_hits'] = int(g.hit.sum())
            row[f'unmarked_core_{tag}_n'] = int(len(g))
        if v in v7.index and int(v7.loc[v, 'edit_drop']) >= 0:
            ed = int(v7.loc[v, 'edit_drop'])
            for ci, cn in enumerate(art['classes']):
                row[f'predicted_drop_TF_{cn}'] = float(-100 * art['delta'][ed, ci, 1, foot].mean())
        vrows.append(row)
    var = pd.DataFrame(vrows).sort_values(['oligo', 'snv'])
    var.to_csv(out / 'variants.tsv', sep='\t', index=False)
    vs = {}
    for zone, g in var.groupby('zone'):
        vs[zone] = dict(snvs=len(g), TF_A7_mean=round(g.TF_A7.mean(), 2), TF_A7_sd=round(g.TF_A7.std(), 2),
                        PIC_A7_mean=round(g.PIC_A7.mean(), 2), unmarked_core_A7_mean=round(g.unmarked_core_A7.mean(), 2),
                        TF_hits_A7=int(g.TF_A7_hits.sum()))
    cg = var[var.caller_class == 'neither']
    at = var[var.caller_class != 'neither']
    multi = var.groupby('oligo').filter(lambda g: len(g) >= 2)
    icc = None
    if len(multi):
        between = multi.groupby('oligo').TF_A7.mean().var()
        icc = float(between / multi.TF_A7.var())
    ok = var.dropna(subset=['TF_A7', 'act'])
    fe = ok.assign(x=ok.TF_A7 - ok.groupby('oligo').TF_A7.transform('mean'), y=ok.act - ok.groupby('oligo').act.transform('mean'))
    summ['variants'] = dict(
        zones=vs,
        neither_vs_at_changing=dict(neither_n=len(cg), neither_TF_A7=round(cg.TF_A7.mean(), 2), neither_TF_hits=int(cg.TF_A7_hits.sum()),
                                    at_n=len(at), at_TF_A7=round(at.TF_A7.mean(), 2)),
        allele_share_of_variance=None if icc is None else round(icc, 2),
        edge_vs_core=dict(edge_abs=round(var[var.zone.isin(['A run 59-65', 'gap 66-69', 'gap 79-84', 'TT 85-86'])].TF_A7.abs().mean(), 2),
                          core_abs=round(var[var.zone == 'core 70-78'].TF_A7.abs().mean(), 2)),
        activity=dict(spearman=round(float(stats.spearmanr(ok.TF_A7, ok.act).correlation), 3), n=len(ok),
                      within_position_spearman=round(float(stats.spearmanr(fe.x, fe.y).correlation), 3)),
        merging=dict(corr_TF_PIC=round(float(var[['TF_A7', 'PIC_A7']].corr().iloc[0, 1]), 3),
                     tf_lost_pic_gained=int(((var.TF_A7 > 0.5) & (var.PIC_A7 < -0.5)).sum())))
    print('variants', json.dumps(summ['variants'], default=float), flush=True)
    json.dump(summ, open(out / 'summary.json', 'w'), indent=1, default=float)
    figures(out, ad, an, lab, var, comp)


def figures(out, ad, an, lab, var, comp):
    a7f.a4f.style()
    shade = lambda a: [a.axvspan(lo - .5, hi + .5, color=c, alpha=.12, lw=0) for lo, hi, c in
                       ((52, 62, ORANGE), (89, 104, BLUE), (*ARUN, GREY), (*TT, GREY))]
    # anatomy
    fig, ax = plt.subplots(3, 1, figsize=(15, 9.5), sharex=True, gridspec_kw=dict(hspace=0.35))
    a = ax[0]
    shade(a)
    a.bar(ad.oligo, ad.w_unmarked, color=[INK if b in 'A' else BLUE for b in ad.base], width=0.7)
    a.set_ylabel('evidence per\nunmarked A/T')
    a.set_title('Oligo 50-95 · bars: ln((1-p_p)/(1-p_a)) per A/T (top-strand A dark, T = bottom-strand A blue); grey: A run 59-65, TT 85-86',
                loc='left', color=INK, fontsize=12)
    a = ax[1]
    shade(a)
    a.plot(ad.oligo, 100 * ad.wt_unmarked, 'o-', color=GREY, label='WT unmarked, labeling > 0.01')
    a.plot(ad.oligo, 100 * ad.wt_unmarked_M, 'o-', color=ORANGE, label='WT unmarked, M (>= 0.6)')
    a.set_ylabel('% unmarked')
    a.legend(frameon=False, fontsize=10.5, loc='lower right')
    a = ax[2]
    shade(a)
    tf = an[an.read_len.between(*TF)]
    a.hist(tf.read_start, bins=np.arange(40.5, 96.5), color=BLUE, alpha=.8, label='start of 20-40 bp calls covering 73')
    a.hist(tf.read_end, bins=np.arange(40.5, 96.5), color=ORANGE, alpha=.8, label='end')
    a.set_ylabel('WT molecules')
    a.set_xlabel('oligo position')
    a.legend(frameon=False, fontsize=10.5)
    for a in ax:
        a.spines[['top', 'right']].set_visible(False)
        a.set_xlim(49.5, 95.5)
    fig.savefig(out / 'figures/fp_anatomy.png', dpi=100, bbox_inches='tight')
    plt.close(fig)
    # call length with hidden bases
    fig, ax = plt.subplots(figsize=(12, 4.2))
    bins = np.arange(0, 101, 2)
    for tag, lab_, c in (('read', 'as read', INK), ('no_arun', 'A run 59-65 hidden', ORANGE), ('no_arun_tt', 'A run and TT 85-86 hidden', BLUE)):
        L = an[f'{tag}_len'].dropna()
        ax.hist(L, bins=bins, histtype='step', lw=2, color=c, label=f'{lab_}: {100 * len(L) / len(an):.1f}% covered, '
                f'{100 * L.between(*TF).sum() / len(an):.1f}% by 20-40 bp')
    ax.axvspan(19.5, 40.5, color=GREY, alpha=.12, lw=0)
    ax.set_xlabel('length of the call covering oligo 73 (bp)')
    ax.set_ylabel('WT molecules')
    ax.legend(frameon=False, fontsize=10.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_title(f'{len(an):,} WT molecules (labeling >= 0.2) re-called with bases hidden', loc='left', color=INK, fontsize=12)
    fig.savefig(out / 'figures/fp_lengths.png', dpi=100, bbox_inches='tight')
    plt.close(fig)
    # labeling
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.4), sharey=False, gridspec_kw=dict(wspace=0.25))
    for a, pr in zip(ax, (PROBE, comp)):
        g = lab[lab.probe == pr]
        x = np.arange(len(g))
        a.plot(x, 100 * g.observed, 'o-', color=INK, label='WT observed')
        a.plot(x, 100 * g.simulated_naked, 'o--', color=ORANGE, label='naked DNA simulated from this sequence')
        a.set_xticks(x, g.labeling)
        a.set_xlabel('molecule labeling')
        a.set_ylabel('20-40 bp coverage (%)')
        a.set_title(f'oligo {pr}' + (' (60-86 footprint)' if pr == PROBE else ' (most-covered 20-40 bp position elsewhere)'),
                    loc='left', color=INK, fontsize=12)
        a.spines[['top', 'right']].set_visible(False)
    ax[0].legend(frameon=False, fontsize=10.5)
    fig.savefig(out / 'figures/fp_labeling.png', dpi=100, bbox_inches='tight')
    plt.close(fig)
    # variants
    off = {'A': -.27, 'C': -.09, 'G': .09, 'T': .27}
    col = {'drops its base': ORANGE, 'gains an A/T': BLUE, 'neither': INK}
    panels = [('TF_A7', '20-40 bp coverage\nover 66-84, WT - var (pts)', 'TF_A7_hits'),
              ('PIC_A7', '41-80 bp coverage\nover 66-84, WT - var (pts)', 'PIC_A7_hits'),
              ('unmarked_core_A7', 'unmarked at 70-78\nWT - var (pts)', 'unmarked_core_A7_hits'),
              ('act', 'activity LDLR\n(log2 var/WT)', None)]
    fig, ax = plt.subplots(len(panels), 1, figsize=(15, 11), sharex=True, gridspec_kw=dict(hspace=0.3))
    for a, (k, lab_, hk) in zip(ax, panels):
        shade(a)
        for r in var.itertuples():
            v = getattr(r, k)
            if not np.isfinite(v):
                continue
            x = r.oligo + off[r.snv[-1]]
            filled = hk is None or getattr(r, hk) > 0
            a.plot([x, x], [0, v], color=GREY, lw=1)
            a.scatter([x], [v], s=34, color=col[r.caller_class] if filled else 'none', edgecolor=col[r.caller_class], lw=1.3, zorder=3)
        a.axhline(0, color=GREY, lw=.8)
        a.set_ylabel(lab_, fontsize=11)
        a.spines[['top', 'right']].set_visible(False)
    h = [plt.Line2D([], [], marker='o', ls='', color=c, label=n) for n, c in col.items()]
    h.append(plt.Line2D([], [], marker='o', ls='', markerfacecolor='none', markeredgecolor=GREY, label='no significant feature there (A7-mask test)'))
    ax[0].legend(handles=h, frameon=False, fontsize=10.5, ncol=4, loc='upper left', bbox_to_anchor=(0, 1.35))
    ax[-1].set_xlabel('variant position (oligo) · alleles A, C, G, T left to right · orange band ZNF652-like, blue Repeat 1, grey A run and TT')
    fig.savefig(out / 'figures/fp_variants.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()

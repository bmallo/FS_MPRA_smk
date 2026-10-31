"""Exploratory follow-up (after Phase A): is a WT footprint in a given window real? (HepG2 run 3)

The fp60_86.py analyses for any window, with size classes chosen per footprint:
  adenines.tsv    each A/T in the window: naked relative marking (A4 lookup), run-3 p_a/p_p, evidence per unmarked
                  adenine w = ln((1-p_p)/(1-p_a)), WT unmarked fraction (all labeled, M)
  anatomy.tsv.gz  WT calls (labeling >= 0.2, sample) covering the probe: start, end, length, LLR, number of informative
                  unmarked adenines (w >= 0.5) and the share of the LLR they carry; as read and with each --hide set
                  unobserved
  labeling.tsv    coverage at each probe by labeling bin and call size (10-19, 20-40, 10-40 bp): WT observed vs naked
                  DNA simulated from this sequence (fp60_86.simulate's model)
  variants.tsv    every SNV in the window: WT - variant coverage over --span by size (A7- and +/-5 bp-mask tests),
                  unmarked change over the span's informative adenines, hits, zone, caller class, activity, and the
                  A7-predicted caller artifact of dropping its base
  summary.json, figures/

Usage: python tools/analysis_v2/fp_region.py --root results/analysis_v2_20261028 --window 100:150 --span 114:138 \
           --probe 123 --probe 132 --compare 73 --hide runs=119:127 --hide arun=124:127 \
           --zone left=110:118 --zone runs=119:127 --zone SRE=128:138 --output <new dir>
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

_spec = importlib.util.spec_from_file_location('fp', Path(__file__).with_name('fp60_86.py'))
fp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fp)
a7, a7f, a6x = fp.a7, fp.a7f, fp.a6x
INK, BODY, BLUE, ORANGE, GREY = a7f.INK, a7f.BODY, a7f.BLUE, a7f.ORANGE, a7f.GREY
TEAL = fp.TEAL
ref, olig, calls, covering = fp.ref, fp.olig, fp.calls, fp.covering
SIZES = {'10-19': (10, 19), '20-40': (20, 40), '10-40': (10, 40)}
BINS = fp.BINS
INFO_W = 0.5
N_WT, N_SIM = 40000, 20000
SEED = 20261030
MOTIFS = [('ZNF652-like', 52, 62, ORANGE), ('SP4', 82, 102, ORANGE), ('Repeat 1 (Sp1)', 89, 104, BLUE),
          ('Repeat 2 (SRE-1)', 124, 139, BLUE), ('SREBF1/2', 127, 136, ORANGE), ('Repeat 3 (Sp1)', 140, 155, BLUE)]


def span(s):
    a, b = s.split(':')
    return int(a), int(b)


def named(s):
    k, v = s.split('=')
    return k, span(v)


def size_of(length):
    return [k for k, (lo, hi) in SIZES.items() if lo <= length <= hi]


def anatomy(Xs, dist, hit, miss, pos, probe, hides):
    """Per WT molecule: the call covering the probe, as read and with each hide set unobserved."""
    info = miss >= INFO_W
    hide_cols = {k: np.flatnonzero((pos >= ref(a)) & (pos <= ref(b))) for k, (a, b) in hides}
    rows = []
    for x, d in zip(Xs, dist):
        out = dict(labeling=d)
        for tag in ['read'] + list(hide_cols):
            y = x
            if tag != 'read':
                y = x.copy()
                y[hide_cols[tag]] = -1
            c = covering(calls(y, hit, miss, pos), ref(probe))
            if c:
                s, e, l, n = c[0]
                cols = np.flatnonzero((pos >= s) & (pos < e))
                un = cols[y[cols] == 0]
                out.update({f'{tag}_start': olig(s), f'{tag}_end': olig(e - 1), f'{tag}_len': e - s, f'{tag}_llr': l,
                            f'{tag}_info_unmarked': int(info[un].sum()),
                            f'{tag}_info_share': float(miss[un][info[un]].sum() / max(miss[un].sum(), 1e-9))})
        rows.append(out)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--window', required=True, type=span)
    ap.add_argument('--span', required=True, type=span, help='footprint positions averaged for variant effects')
    ap.add_argument('--probe', type=int, action='append', required=True)
    ap.add_argument('--compare', type=int, help='a probe outside the window for the labeling comparison')
    ap.add_argument('--hide', type=named, action='append', default=[], help='name=lo:hi set unobserved in the anatomy')
    ap.add_argument('--zone', type=named, action='append', default=[], help='name=lo:hi for grouping variants')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'figures').mkdir()
    rng = np.random.default_rng(SEED)
    W0, W1 = a.window
    S0, S1 = a.span
    ext = root / 'a9_statistics/extract/HepG2_run3'
    cfg = json.load(open(root / 'a6_proxies/extract/HepG2_run3/extract.json'))
    seq_ref = a6x.a5x.ref_seq(cfg['reference'])
    pos, X, distal, _, _ = a6x.load(ext)
    cls = np.load(ext / 'meta.npz')['cls']
    em = 'results/llr_vs_fiberhmm_20261026/results/emissions.tsv.gz'
    hit, miss = a7.emission_steps(em, pos)
    t = pd.read_csv(em, sep='\t')
    t = t[t.fold == 'A']
    pa, pp = dict(zip(t.ref_pos - 1, t.p_a)), dict(zip(t.ref_pos - 1, t.p_p))
    lookup = pd.read_csv(root / 'a4_informative/visibility/naked_lookup.tsv', sep='\t', index_col=0).p_rel.to_dict()
    rel = fp.naked_rates(seq_ref, pos, lookup)
    wt = np.flatnonzero(cls == 'WT')
    summ = dict(window=a.window, span=a.span, probes=a.probe, hides=a.hide, zones=a.zone)

    # 1. adenines
    rows = []
    for c in np.flatnonzero((pos >= ref(W0)) & (pos <= ref(W1))):
        x = X[wt, c]
        k1, km = distal[wt] > 0.01, distal[wt] >= 0.6
        obs = x >= 0
        rows.append(dict(oligo=olig(int(pos[c])), base=seq_ref[pos[c]], naked_rel=round(rel[c], 3), p_a=round(pa[pos[c]], 3),
                         p_p=round(pp[pos[c]], 4), w_unmarked=round(miss[c], 2), w_marked=round(hit[c], 2),
                         wt_unmarked=round(((x == 0) & k1).sum() / max((obs & k1).sum(), 1), 4),
                         wt_unmarked_M=round(((x == 0) & km).sum() / max((obs & km).sum(), 1), 4)))
    ad = pd.DataFrame(rows)
    ad.to_csv(out / 'adenines.tsv', sep='\t', index=False)
    sp = ad[ad.oligo.between(S0, S1)]
    summ['adenines'] = dict(at_in_span=int(len(sp)), informative_in_span=int((sp.w_unmarked >= INFO_W).sum()),
                            evidence_in_span=round(float(sp.w_unmarked.sum()), 2),
                            evidence_uninformative=round(float(sp[sp.w_unmarked < INFO_W].w_unmarked.sum()), 2))
    print('adenines', summ['adenines'], flush=True)

    # 2. anatomy at the first probe
    k2 = wt[distal[wt] >= 0.2]
    samp = rng.choice(k2, min(N_WT, k2.size), replace=False)
    probe = a.probe[0]
    an = anatomy(X[samp], distal[samp], hit, miss, pos, probe, a.hide)
    an.to_csv(out / 'anatomy.tsv.gz', sep='\t', index=False, compression='gzip')
    an_s = {}
    for tag in ['read'] + [k for k, _ in a.hide]:
        L = an[f'{tag}_len'].dropna()
        d = dict(covered=round(len(L) / len(an), 4), median_len=float(L.median()) if len(L) else None)
        for k, (lo, hi) in list(SIZES.items()) + [('41-80', (41, 80)), ('81+', (81, 10 ** 6))]:
            d[k] = round(float(L.between(lo, hi).sum() / len(an)), 4)
        sm = an[an[f'{tag}_len'].between(10, 40)]
        if len(sm):
            d.update(start_mode=int(sm[f'{tag}_start'].mode()[0]), end_mode=int(sm[f'{tag}_end'].mode()[0]),
                     llr_median=round(float(sm[f'{tag}_llr'].median()), 2),
                     info_unmarked_median=float(sm[f'{tag}_info_unmarked'].median()),
                     info_unmarked_le2=round(float((sm[f'{tag}_info_unmarked'] <= 2).mean()), 3),
                     info_share_median=round(float(sm[f'{tag}_info_share'].median()), 3))
        an_s[tag] = d
    summ['anatomy'] = dict(probe=probe, molecules=len(an), **an_s)
    print('anatomy', json.dumps(summ['anatomy']), flush=True)

    # 3. labeling: observed vs simulated naked, per size class, at each probe (+ comparison)
    elcols = np.flatnonzero((pos >= fp.EL0) & (pos < fp.EL0 + 211) & ~((pos >= ref(W0)) & (pos <= ref(W1))))
    wcols = np.flatnonzero((pos >= fp.EL0 - 150) & (pos < fp.EL0 + 211 + 150))
    probes = a.probe + ([a.compare] if a.compare else [])
    lrows = []
    for lo, hi in fp.LAB_BINS:
        m = k2[(distal[k2] >= lo) & (distal[k2] < hi)]
        if m.size < 200:
            continue
        mm = rng.choice(m, min(m.size, 6000), replace=False)
        obs = {(pr, k): 0 for pr in probes for k in SIZES}
        for x in X[mm]:
            cl = calls(x, hit, miss, pos)
            for pr in probes:
                c = covering(cl, ref(pr))
                if c:
                    for k in size_of(c[0][1] - c[0][0]):
                        obs[(pr, k)] += 1
        xe = X[mm][:, elcols]
        target = float((xe == 1).sum() / max((xe >= 0).sum(), 1))
        sim = simulate(rel[wcols], hit[wcols], miss[wcols], pos[wcols], target, N_SIM // len(fp.LAB_BINS), rng, probes)
        for pr in probes:
            for k in SIZES:
                lrows.append(dict(labeling=f'{lo:.1f}-{min(hi, 1):.1f}', probe=pr, size=k, wt_molecules=int(m.size),
                                  element_marking=round(target, 3), observed=round(obs[(pr, k)] / mm.size, 4),
                                  simulated_naked=round(sim[(pr, k)], 4)))
        print('labeling', lo, hi, flush=True)
    lab = pd.DataFrame(lrows)
    lab.to_csv(out / 'labeling.tsv', sep='\t', index=False)
    summ['labeling'] = lab.to_dict('records')

    # 4. variants in the window
    e7 = pd.read_parquet(root / 'a10_effects/HepG2_run3/effects.parquet')
    e5 = pd.read_parquet(root / 'a10_effects/HepG2_run3_pm5/effects.parquet')
    fn = pd.read_csv(root / 'a11_motifs_function/function_snvs.tsv', sep='\t').set_index('snv')
    v7 = pd.read_csv(root / 'a7_snv_edits/HepG2_run3/variants.tsv', sep='\t').set_index('source_variant')
    art = np.load(root / 'a7_snv_edits/HepG2_run3/artifacts.npz')
    info_ad = ad[ad.oligo.between(S0, S1) & (ad.w_unmarked >= INFO_W)].oligo.tolist()
    snv = e7[(e7.endpoint == 'E1') & e7.snv_oligo.between(W0 + 5, W1 - 5)].drop_duplicates('snv')[['snv', 'snv_oligo', 'caller_class']]
    sl = slice(S0 - fp.OLIGO0, S1 - fp.OLIGO0 + 1)
    zone_of = lambda o: next((k for k, (lo, hi) in a.zone if lo <= o <= hi), 'flank')
    vrows = []
    for v, o, cc in snv.itertuples(index=False):
        row = dict(snv=v, oligo=o, caller_class=cc, zone=zone_of(o), act=fn.act.get(v), act2=fn.act2.get(v))
        for tag, e in (('A7', e7), ('pm5', e5)):
            g = e[(e.snv == v) & (e.endpoint == 'E2') & e.oligo.between(S0, S1)]
            for b in BINS:
                h = g[(g['bin'] == b) & g.testable]
                row[f'{b}_{tag}'] = -100 * h.D.mean() if len(h) else np.nan
                row[f'{b}_{tag}_hits'] = int(h.hit.sum())
            g = e[(e.snv == v) & (e.endpoint == 'E1') & e.oligo.isin(info_ad) & e.testable]
            row[f'unmarked_{tag}'] = -100 * g.D.mean() if len(g) else np.nan
            row[f'unmarked_{tag}_hits'] = int(g.hit.sum())
            row[f'unmarked_{tag}_n'] = int(len(g))
        if v in v7.index and int(v7.loc[v, 'edit_drop']) >= 0:
            ed = int(v7.loc[v, 'edit_drop'])
            for bi, b in enumerate(BINS):
                if b in ('sub_TF', 'TF'):
                    for ci, cn in enumerate(art['classes']):
                        row[f'predicted_drop_{b}_{cn}'] = float(-100 * art['delta'][ed, ci, bi, sl].mean())
        vrows.append(row)
    var = pd.DataFrame(vrows).sort_values(['oligo', 'snv'])
    var['small_A7'] = var.sub_TF_A7 + var.TF_A7
    var.to_csv(out / 'variants.tsv', sep='\t', index=False)
    vs = {}
    for zone, g in var.groupby('zone'):
        vs[zone] = dict(snvs=len(g), sub_TF_A7=round(g.sub_TF_A7.mean(), 2), TF_A7=round(g.TF_A7.mean(), 2),
                        small_A7_sd=round(g.small_A7.std(), 2), unmarked_A7=round(g.unmarked_A7.mean(), 2),
                        sub_TF_hits=int(g.sub_TF_A7_hits.sum()), TF_hits=int(g.TF_A7_hits.sum()),
                        act=round(g.act.mean(), 2))
    multi = var.groupby('oligo').filter(lambda g: len(g) >= 2)
    icc = float(multi.groupby('oligo').small_A7.mean().var() / multi.small_A7.var()) if len(multi) > 3 else None
    ok = var.dropna(subset=['small_A7', 'act'])
    fe = ok.assign(x=ok.small_A7 - ok.groupby('oligo').small_A7.transform('mean'), y=ok.act - ok.groupby('oligo').act.transform('mean'))
    cg, at = var[var.caller_class == 'neither'], var[var.caller_class != 'neither']
    summ['variants'] = dict(zones=vs, allele_share_of_variance=None if icc is None else round(icc, 2),
                            neither=dict(n=len(cg), small_A7=round(cg.small_A7.mean(), 2), hits=int((cg.sub_TF_A7_hits + cg.TF_A7_hits).sum())),
                            at_changing=dict(n=len(at), small_A7=round(at.small_A7.mean(), 2)),
                            activity=dict(spearman=round(float(stats.spearmanr(ok.small_A7, ok.act).correlation), 3), n=len(ok),
                                          within_position=round(float(stats.spearmanr(fe.x, fe.y).correlation), 3)),
                            size_shift=dict(corr_subTF_TF=round(float(var[['sub_TF_A7', 'TF_A7']].corr().iloc[0, 1]), 3)))
    print('variants', json.dumps(summ['variants'], default=float), flush=True)
    json.dump(summ, open(out / 'summary.json', 'w'), indent=1, default=float)
    figures(out, ad, an, lab, var, a, probe)


def simulate(rel, hit, miss, pos, target, n, rng, probes):
    """fp60_86.simulate generalised to several size classes: coverage at each probe by size."""
    lo, hi = 0.0, 5.0
    for _ in range(50):
        q = (lo + hi) / 2
        lo, hi = (q, hi) if np.minimum(1, q * rel).mean() < target else (lo, q)
    p = np.minimum(1, q * rel)
    cov = {(pr, k): 0 for pr in probes for k in SIZES}
    for _ in range(n):
        x = (rng.random(p.size) < p).astype(np.int8)
        cl = calls(x, hit, miss, pos)
        for pr in probes:
            c = covering(cl, ref(pr))
            if c:
                for k in size_of(c[0][1] - c[0][0]):
                    cov[(pr, k)] += 1
    return {k: v / n for k, v in cov.items()}


def figures(out, ad, an, lab, var, a, probe):
    a7f.a4f.style()
    W0, W1 = a.window

    def shade(ax_):
        for _, lo, hi, c in MOTIFS:
            if hi >= W0 and lo <= W1:
                ax_.axvspan(lo - .5, hi + .5, color=c, alpha=.08, lw=0)
        for _, (lo, hi) in a.hide:
            ax_.axvspan(lo - .5, hi + .5, color=GREY, alpha=.14, lw=0)

    fig, ax = plt.subplots(3, 1, figsize=(15, 9.5), sharex=True, gridspec_kw=dict(hspace=0.35))
    shade(ax[0])
    ax[0].bar(ad.oligo, ad.w_unmarked, color=[INK if b == 'A' else BLUE for b in ad.base], width=0.7)
    ax[0].axhline(INFO_W, color=GREY, lw=1, ls=':')
    ax[0].set_ylabel('evidence per\nunmarked A/T')
    ax[0].set_title(f'Oligo {W0}-{W1} · bars: ln((1-p_p)/(1-p_a)) (A dark, T blue); dotted: informative cut 0.5; grey: hidden sets; '
                    f'blue/orange bands: motifs', loc='left', color=INK, fontsize=12)
    shade(ax[1])
    ax[1].plot(ad.oligo, 100 * ad.wt_unmarked, 'o-', color=GREY, label='WT unmarked, labeling > 0.01')
    ax[1].plot(ad.oligo, 100 * ad.wt_unmarked_M, 'o-', color=ORANGE, label='WT unmarked, M (>= 0.6)')
    ax[1].set_ylabel('% unmarked')
    ax[1].legend(frameon=False, fontsize=10.5, loc='lower right')
    shade(ax[2])
    sm = an[an.read_len.between(10, 40)]
    ax[2].hist(sm.read_start, bins=np.arange(W0 - 10.5, W1 + 10.5), color=BLUE, alpha=.8, label=f'start of 10-40 bp calls covering {probe}')
    ax[2].hist(sm.read_end, bins=np.arange(W0 - 10.5, W1 + 10.5), color=ORANGE, alpha=.8, label='end')
    ax[2].set_ylabel('WT molecules')
    ax[2].set_xlabel('oligo position')
    ax[2].legend(frameon=False, fontsize=10.5)
    for x in ax:
        x.spines[['top', 'right']].set_visible(False)
        x.set_xlim(W0 - .5, W1 + .5)
    fig.savefig(out / 'figures/fp_anatomy.png', dpi=100, bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 4.2))
    cols = [INK, ORANGE, BLUE, TEAL, a7f.PURPLE, '#9A4A12']
    for (tag, lab_), c in zip([('read', 'as read')] + [(k, f'{k} {lo}-{hi} hidden') for k, (lo, hi) in a.hide], cols):
        L = an[f'{tag}_len'].dropna()
        ax.hist(L, bins=np.arange(0, 101, 2), histtype='step', lw=2, color=c,
                label=f'{lab_}: {100 * len(L) / len(an):.1f}% covered, {100 * L.between(10, 40).sum() / len(an):.1f}% by 10-40 bp')
    ax.axvspan(9.5, 40.5, color=GREY, alpha=.12, lw=0)
    ax.set_xlabel(f'length of the call covering oligo {probe} (bp)')
    ax.set_ylabel('WT molecules')
    ax.legend(frameon=False, fontsize=10.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_title(f'{len(an):,} WT molecules (labeling >= 0.2) re-called with bases hidden', loc='left', color=INK, fontsize=12)
    fig.savefig(out / 'figures/fp_lengths.png', dpi=100, bbox_inches='tight')
    plt.close(fig)

    probes = list(dict.fromkeys(lab.probe))
    fig, ax = plt.subplots(len(SIZES), len(probes), figsize=(5.2 * len(probes), 3.6 * len(SIZES)), squeeze=False,
                           gridspec_kw=dict(hspace=0.5, wspace=0.28))
    for j, pr in enumerate(probes):
        for i, k in enumerate(SIZES):
            g = lab[(lab.probe == pr) & (lab['size'] == k)]
            x = np.arange(len(g))
            a_ = ax[i, j]
            a_.plot(x, 100 * g.observed, 'o-', color=INK, label='WT observed')
            a_.plot(x, 100 * g.simulated_naked, 'o--', color=ORANGE, label='naked simulated')
            a_.set_xticks(x, g.labeling, fontsize=9)
            a_.set_title(f'oligo {pr} · {k} bp calls' + (' (comparison)' if pr == a.compare else ''), loc='left', color=INK, fontsize=11.5)
            a_.set_ylabel('coverage (%)')
            a_.spines[['top', 'right']].set_visible(False)
    ax[-1, 0].set_xlabel('molecule labeling')
    ax[0, 0].legend(frameon=False, fontsize=10)
    fig.savefig(out / 'figures/fp_labeling.png', dpi=100, bbox_inches='tight')
    plt.close(fig)

    off = {'A': -.27, 'C': -.09, 'G': .09, 'T': .27}
    col = {'drops its base': ORANGE, 'gains an A/T': BLUE, 'neither': INK}
    S0, S1 = a.span
    panels = [('sub_TF_A7', f'10-19 bp coverage\nover {S0}-{S1}, WT - var', 'sub_TF_A7_hits'),
              ('TF_A7', f'20-40 bp coverage\nover {S0}-{S1}, WT - var', 'TF_A7_hits'),
              ('unmarked_A7', 'unmarked, informative\nadenines, WT - var', 'unmarked_A7_hits'),
              ('act', 'activity LDLR\n(log2 var/WT)', None)]
    fig, ax = plt.subplots(len(panels), 1, figsize=(15, 11), sharex=True, gridspec_kw=dict(hspace=0.3))
    for a_, (k, lab_, hk) in zip(ax, panels):
        shade(a_)
        for r in var.itertuples():
            v = getattr(r, k)
            if not np.isfinite(v):
                continue
            x = r.oligo + off[r.snv[-1]]
            filled = hk is None or getattr(r, hk) > 0
            a_.plot([x, x], [0, v], color=GREY, lw=1)
            a_.scatter([x], [v], s=34, color=col[r.caller_class] if filled else 'none', edgecolor=col[r.caller_class], lw=1.3, zorder=3)
        a_.axhline(0, color=GREY, lw=.8)
        a_.set_ylabel(lab_, fontsize=11)
        a_.spines[['top', 'right']].set_visible(False)
        a_.set_xlim(W0 - .5, W1 + .5)
    h = [plt.Line2D([], [], marker='o', ls='', color=c, label=n) for n, c in col.items()]
    h.append(plt.Line2D([], [], marker='o', ls='', markerfacecolor='none', markeredgecolor=GREY, label='no significant feature (A7-mask test)'))
    ax[0].legend(handles=h, frameon=False, fontsize=10.5, ncol=4, loc='upper left', bbox_to_anchor=(0, 1.35))
    ax[-1].set_xlabel('variant position (oligo) · alleles A, C, G, T left to right · points = percentage points')
    fig.savefig(out / 'figures/fp_variants.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()

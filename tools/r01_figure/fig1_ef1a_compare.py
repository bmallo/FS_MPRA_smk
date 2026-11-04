"""Figure 1 (EF1a): compare the EF1a promoter in pEF-GFP (replicate 3) with the endogenous EEF1A1 promoter in HEK293T
genomic Fiber-seq (fig1_rules.json: ef1a). Exploratory, descriptive.

Window = the plasmid's EF1a insert on hg38 (chr6:73,520,048-73,521,236). Per dataset (fig1_profiles.py --gene EF1a
molecules), per base: share of molecules covered by a call of each Figure 1 bin (20-30, 31-50, 51-80 bp; nucleosome
>= 81 bp; molecules covering = span of observed A/T), FIRE share (fig1_fire_profiles.py --gene EF1a) and the m6A
fraction over observed A/T in a 51-bp window. Per molecule spanning the window: share of it covered by nucleosome calls.
Summaries: mean shares over the whole window, the upstream + exon 1 part (TSS side) and intron 1; Spearman of the m6A
and nucleosome profiles between datasets (positions are not independent: descriptive only).
Usage: python tools/r01_figure/fig1_ef1a_compare.py --work results/r01_figure_20261102/fig1_ef1a
"""
import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import spearmanr

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

X0, X1 = 73520048, 73521236
TSS, EXON2_END, INTRON = 73521032, 73520056, (73520057, 73520998)   # NM_001402.6 (- strand)
PARTS = {'window': (X0, X1), 'upstream + exon 1': (73520999, X1), 'intron 1': INTRON}
BINS = {'fp_20_30': (20, 30), 'fp_31_50': (31, 50), 'fp_51_80': (51, 80), 'nuc': (81, 10 ** 6)}
COL = {'nuc': '#C7C8CA', 'FIRE': '#BA7C69', 'fp_51_80': '#E77C34', 'fp_31_50': '#338B8C', 'fp_20_30': '#755B9D'}
DATASETS = {'endogenous': ('Endogenous HEK293T', 'genome'), 'pEFGFP': ('pEF-GFP plasmid (rep 3)', 'plasmid')}
SIGMA, M6A_WIN, FP_SCALE, FIRE_ALPHA = 3, 51, 2, 0.65


def smooth(raw):
    out = np.full(raw.size, np.nan)
    ok = ~np.isnan(raw)
    edges = np.flatnonzero(np.diff(np.r_[0, ok.astype(int), 0]))
    for a, b in zip(edges[::2], edges[1::2]):
        out[a:b] = gaussian_filter1d(raw[a:b], SIGMA, mode='nearest')
    return out


def tracks(mol):
    L = X1 - X0 + 1
    cov, atn, mn = np.zeros(L), np.zeros(L), np.zeros(L)
    cnt = {c: np.zeros(L) for c in BINS}
    per_mol = []
    for r in mol.itertuples():
        at = np.array(json.loads(r.at), int)
        if not at.size:
            continue
        a, b = max(at.min(), X0), min(at.max(), X1)
        if b < a:
            continue
        cov[a - X0:b - X0 + 1] += 1
        at = at[(at >= a) & (at <= b)] - X0
        atn[at] += 1
        m = np.array(json.loads(r.m6a), int)
        m = m[(m >= a) & (m <= b)] - X0
        mn[m] += 1
        hit = {c: np.zeros(L, bool) for c in BINS}
        for s, e, ln in json.loads(r.fp):
            for c, (l0, l1) in BINS.items():
                if l0 <= ln <= l1:
                    s_, e_ = max(s, a), min(e, b)
                    if e_ >= s_:
                        hit[c][s_ - X0:e_ - X0 + 1] = True
        for c in BINS:
            cnt[c] += hit[c]
        if a == X0 and b == X1:                                 # spans the whole window
            per_mol.append(hit['nuc'].mean())
    out = {'pos': np.arange(X0, X1 + 1), 'molecules': cov}
    with np.errstate(invalid='ignore', divide='ignore'):
        k = np.ones(M6A_WIN)
        sa, sm = np.convolve(atn, k, 'same'), np.convolve(mn, k, 'same')
        out['m6a'] = np.where((cov > 0) & (sa > 0), sm / np.maximum(sa, 1), np.nan)
        for c in BINS:
            out[f'{c}_share'] = np.where(cov > 0, cnt[c] / np.maximum(cov, 1), np.nan)
            out[c] = smooth(out[f'{c}_share'])
    return pd.DataFrame(out), np.array(per_mol)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--work', required=True)
    a = ap.parse_args()
    w = Path(a.work)
    out = w / 'compare'
    out.mkdir(exist_ok=True)
    fire = pd.read_csv(w / 'fire/fire_profiles.tsv.gz', sep='\t')
    T, PM, N = {}, {}, {}
    for d, (_, ftrack) in DATASETS.items():
        mol = pd.read_parquet(w / f'profiles/{d}.molecules.parquet')
        t, pm = tracks(mol)
        fr = fire[fire.track == ftrack].set_index('pos')
        t['FIRE_share'] = fr.FIRE_share.reindex(t.pos).to_numpy()
        t['FIRE'] = smooth(t.FIRE_share.to_numpy())
        T[d], PM[d], N[d] = t, pm, len(mol)
        t.to_csv(out / f'{d}.tracks.tsv.gz', sep='\t', index=False, float_format='%.4g')
    rows = []
    for part, (lo, hi) in PARTS.items():
        for d, t in T.items():
            s = t[(t.pos >= lo) & (t.pos <= hi)]
            rows.append(dict(part=part, start=lo, end=hi, dataset=d, molecules_max=int(s.molecules.max()),
                             **{f'mean_{c}': round(s[f'{c}_share'].mean(), 3) for c in BINS},
                             mean_FIRE=round(s.FIRE_share.mean(), 3), mean_m6a=round(s.m6a.mean(), 3)))
    st = pd.DataFrame(rows)
    cors = []
    for part, (lo, hi) in PARTS.items():
        for col in ('m6a', 'nuc_share', 'FIRE_share'):
            x = T['endogenous'].set_index('pos').loc[lo:hi, col].to_numpy()
            y = T['pEFGFP'].set_index('pos').loc[lo:hi, col].to_numpy()
            ok = ~np.isnan(x) & ~np.isnan(y)
            cors.append(dict(part=part, profile=col, positions=int(ok.sum()), spearman=round(spearmanr(x[ok], y[ok])[0], 2)))
    co = pd.DataFrame(cors)
    pm = pd.DataFrame([dict(dataset=d, molecules_spanning=len(v), nuc_cover_median=round(np.median(v), 3),
                            share_under_25pct=round((v < 0.25).mean(), 3), share_over_75pct=round((v > 0.75).mean(), 3))
                       for d, v in PM.items()])
    st.to_csv(out / 'summary.tsv', sep='\t', index=False)
    co.to_csv(out / 'spearman.tsv', sep='\t', index=False)
    pm.to_csv(out / 'per_molecule_nucleosome.tsv', sep='\t', index=False)
    print(st.to_string(index=False)), print(co.to_string(index=False)), print(pm.to_string(index=False))

    # figure: gene model + motif ticks, a track per dataset, m6A overlay, per-molecule nucleosome cover
    ref = json.load(open(w / 'raw/ucsc_ncbiRefSeq_eef1a1.json'))['ncbiRefSeq']
    g = next(r for r in ref if r['name'] == 'NM_001402.6')
    mot = pd.read_csv(w / 'motifs/motifs.tsv', sep='\t')
    plt.rcParams['axes.formatter.useoffset'] = False
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(5, 4, height_ratios=[0.55, 1, 1, 0.9, 0.05], hspace=0.45, wspace=0.6)
    axg = fig.add_subplot(gs[0, :3])
    axg.set_xlim(X0, X1); axg.set_ylim(0, 1); axg.axis('off')
    axg.plot([X0, X1], [0.35, 0.35], color='#1f6fb4', lw=0.8)
    for s, e in zip(map(int, g['exonStarts'].strip(',').split(',')), map(int, g['exonEnds'].strip(',').split(','))):
        if e > X0 and s < X1:
            axg.add_patch(Rectangle((max(s, X0), 0.2), min(e, X1) - max(s, X0), 0.3, color='#1f6fb4'))
    axg.plot([TSS, TSS, TSS - 50], [0.35, 0.8, 0.8], color='#1f6fb4', lw=1)
    axg.annotate('', xy=(TSS - 90, 0.8), xytext=(TSS - 45, 0.8), arrowprops=dict(arrowstyle='-|>', color='#1f6fb4', lw=1))
    axg.text(X0, 0.85, 'EEF1A1 exon 2', fontsize=5.5, color='#1f6fb4')
    axg.text(TSS + 15, 0.85, 'TSS / exon 1', fontsize=5.5, color='#1f6fb4')
    axg.text((INTRON[0] + INTRON[1]) / 2, 0.55, 'intron 1', fontsize=5.5, ha='center', color='#1f6fb4')
    for r in mot.itertuples():
        axg.plot([r.start, r.end], [0.05, 0.05], color='#B8860B', lw=2.5, solid_capstyle='butt')
    axg.set_title(f'EF1α promoter as in pEF-GFP: chr6:{X0:,}-{X1:,} (gold: {len(mot)} motif matches, P < 1e-5)',
                  fontsize=6.5)
    for i, (d, (lab, _)) in enumerate(DATASETS.items()):
        ax = fig.add_subplot(gs[1 + i, :3], sharex=axg)
        t = T[d]
        x = t.pos.to_numpy()
        for c in ['nuc', 'FIRE'] + sorted([c for c in BINS if c != 'nuc'], key=lambda c: -np.nanmax(t[c])):
            y = np.minimum(t[c].to_numpy() * (FP_SCALE if c.startswith('fp_') else 1), 1)
            ax.fill_between(x, 0, np.nan_to_num(y), color=COL[c], lw=0, alpha=FIRE_ALPHA if c == 'FIRE' else 1,
                            label={'nuc': 'nucleosome ≥81', 'FIRE': 'FIRE'}.get(c, c[3:].replace('_', '–') + f' bp ×{FP_SCALE}'))
        ax.plot(x, t.m6a, color='#0072B2', lw=0.8, label='m6A fraction (51 bp)')
        ax.axvline(TSS, color='#1f6fb4', lw=0.6, ls=':')
        ax.set_ylim(0, 1); ax.tick_params(labelsize=5.5)
        ax.set_ylabel(f'{lab}\nn = {int(t.molecules.max()):,}', fontsize=6)
        if i == 0:
            ax.legend(fontsize=5, ncol=3, frameon=False, loc='upper left', bbox_to_anchor=(0, 1.32))
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:,.0f}'))
    axo = fig.add_subplot(gs[3, :3], sharex=axg)
    for d, c in (('endogenous', '#2E2E2E'), ('pEFGFP', '#D55E00')):
        axo.plot(T[d].pos, T[d].m6a, color=c, lw=0.9, label=DATASETS[d][0])
    axo.axvline(TSS, color='#1f6fb4', lw=0.6, ls=':')
    r_ = co[(co.part == 'window') & (co.profile == 'm6a')].spearman.iloc[0]
    axo.set_ylabel('m6A fraction\n(51-bp window)', fontsize=6); axo.tick_params(labelsize=5.5)
    axo.legend(fontsize=5.5, frameon=False, title=f'Spearman {r_:.2f} (descriptive)', title_fontsize=5.5)
    axo.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:,.0f}'))
    axh = fig.add_subplot(gs[1:3, 3])
    bins = np.linspace(0, 1, 21)
    for d, c in (('endogenous', '#2E2E2E'), ('pEFGFP', '#D55E00')):
        axh.hist(PM[d], bins, weights=np.ones(len(PM[d])) / len(PM[d]), histtype='step', color=c, lw=1.1,
                 label=f'{DATASETS[d][0]} ({len(PM[d]):,})')
    axh.set_xlabel('share of window covered\nby nucleosome calls', fontsize=6)
    axh.set_ylabel('share of molecules', fontsize=6); axh.tick_params(labelsize=5.5)
    axh.legend(fontsize=4.8, frameon=False)
    fig.text(0.5, 0.03, 'Exploratory comparison, not a figure draft. Plasmid: ~3% random subsample of replicate 3 '
             'spanning molecules. Footprint shares ×2.', ha='center', fontsize=5.5)
    fig.savefig(out / 'compare.png', dpi=220, bbox_inches='tight')


if __name__ == '__main__':
    main()

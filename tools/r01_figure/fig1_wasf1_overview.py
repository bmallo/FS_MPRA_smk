"""Figure 1 (WASF1): first-look overview to choose the 300-bp construct and the layout (fig1_rules.json: wasf1).

Exploratory, descriptive. From fig1_profiles.py --gene WASF1 molecules: per base, share of molecules with a footprint of
each fp_bins_v2 bin covering it (molecules covering = span of observed A/T; 300-bp constructs limited to their insert)
and the m6A fraction over observed A/T in a 51-bp window (accessibility), for the endogenous HEK293T reads, the 4.5-kb construct and the
two 300-bp constructs. Gaussian sigma 3 bp inside covered runs (as fig1.py).
Per 300-bp region and dataset: molecules spanning it, share with >= 1 footprint 15-89 bp whose midpoint lies in it, mean
nucleosome (>= 90 bp) coverage, mean m6A fraction; Spearman of smoothed m6A between datasets (positions are not
independent: descriptive only).
Usage: python tools/r01_figure/fig1_wasf1_overview.py --work results/r01_figure_20261102/fig1_wasf1
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

X0, X1 = 110176895, 110181442
SIGMA = 3
M6A_WIN = 51                      # bp window for the m6A fraction (m6A count / observed A/T count)
BINS = {'fp_15_30': (15, 30), 'fp_31_60': (31, 60), 'fp_61_89': (61, 89), 'nuc': (90, 10 ** 6)}
BIN_C = {'fp_15_30': '#E69F00', 'fp_31_60': '#D55E00', 'fp_61_89': '#CC79A7', 'nuc': '#9a9893'}
REGIONS = {'WASF11': ('Peak 1 (WASF11)', 110178621, 110178921), 'WASF12': ('Peak 4 (WASF12)', 110179911, 110180211)}
DATASETS = [('endogenous', 'Endogenous HEK293T'), ('WASF1_full', '4.5-kb construct')]
SITES = {'GC box (SP1/SP4/KLF)': (110180054, 110180097), 'ZNF444/ZNF281/VEZF1 site': (110179973, 110179991),
         'ZNF143 site': (110180105, 110180127)}


def smooth(raw):
    out = np.full(raw.size, np.nan)
    ok = ~np.isnan(raw)
    edges = np.flatnonzero(np.diff(np.r_[0, ok.astype(int), 0]))
    for a, b in zip(edges[::2], edges[1::2]):
        out[a:b] = gaussian_filter1d(raw[a:b], SIGMA, mode='nearest')
    return out


def tracks(mol, lo=X0, hi=X1):
    L = X1 - X0 + 1
    cov, atn, mn = np.zeros(L), np.zeros(L), np.zeros(L)
    cnt = {c: np.zeros(L) for c in BINS}
    for r in mol.itertuples():
        at = np.array(json.loads(r.at), int)
        if not at.size:
            continue
        a, b = max(at.min(), lo, X0), min(at.max(), hi, X1)
        if b < a:
            continue
        cov[a - X0:b - X0 + 1] += 1
        at = at[(at >= a) & (at <= b)] - X0
        atn[at] += 1
        m = np.array(json.loads(r.m6a), int)
        m = m[(m >= a) & (m <= b)] - X0
        mn[m] += 1
        for c, (l0, l1) in BINS.items():
            hit = np.zeros(L, bool)
            for s, e, ln in json.loads(r.fp):
                if l0 <= ln <= l1:
                    s_, e_ = max(s, a), min(e, b)
                    if e_ >= s_:
                        hit[s_ - X0:e_ - X0 + 1] = True
            cnt[c] += hit
    out = {'pos': np.arange(X0, X1 + 1), 'molecules': cov}
    with np.errstate(invalid='ignore', divide='ignore'):
        k = np.ones(M6A_WIN)
        sa, sm = np.convolve(atn, k, 'same'), np.convolve(mn, k, 'same')
        out['m6a'] = np.where((cov > 0) & (sa > 0), sm / np.maximum(sa, 1), np.nan)
        for c in BINS:
            out[c] = smooth(np.where(cov > 0, cnt[c] / np.maximum(cov, 1), np.nan))
    return pd.DataFrame(out)


def region_stats(mol, lo, hi):
    n = sub = 0
    for r in mol.itertuples():
        at = np.array(json.loads(r.at), int)
        if not at.size or at.min() > lo + 20 or at.max() < hi - 20:
            continue
        n += 1
        sub += any(15 <= ln <= 89 and lo <= (s + e) / 2 <= hi for s, e, ln in json.loads(r.fp))
    return n, sub


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--work', required=True)
    a = ap.parse_args()
    w = Path(a.work)
    out = w / 'overview'
    out.mkdir(exist_ok=True)
    mol = {d: pd.read_parquet(w / f'profiles/{d}.molecules.parquet') for d in ('endogenous', 'WASF1_full', 'WASF11', 'WASF12')}
    T = {d: tracks(m) for d, m in mol.items() if d in ('endogenous', 'WASF1_full')}
    for c, (_, lo, hi) in REGIONS.items():
        T[c] = tracks(mol[c], lo, hi)
    for d, t in T.items():
        t.to_csv(out / f'{d}.tracks.tsv.gz', sep='\t', index=False, float_format='%.4g')

    rows, cors = [], []
    for c, (name, lo, hi) in REGIONS.items():
        sel = lambda t: t[(t.pos >= lo) & (t.pos <= hi)]  # noqa: E731
        for d in ('endogenous', 'WASF1_full', c):
            n, sub = region_stats(mol[d], lo, hi)
            s = sel(T[d])
            rows.append(dict(region=name, dataset=d, molecules_spanning=n, with_15_89bp_footprint=sub,
                             share_15_89=round(sub / n, 3) if n else np.nan, mean_nuc_cover=round(s.nuc.mean(), 3),
                             mean_m6a=round(s.m6a.mean(), 3)))
        for d1, d2 in (('endogenous', c), ('endogenous', 'WASF1_full'), ('WASF1_full', c)):
            x, y = sel(T[d1]).m6a.to_numpy(), sel(T[d2]).m6a.to_numpy()
            ok = ~np.isnan(x) & ~np.isnan(y)
            cors.append(dict(region=name, pair=f'{d1} vs {d2}', positions=int(ok.sum()),
                             spearman_m6a=round(spearmanr(x[ok], y[ok])[0], 2)))
    # motif sites in peak 4 (best match per site, motifs.tsv): molecules spanning the site with a 15-89 bp call
    # overlapping it by >= 10 bp
    site_rows = []
    for site, (lo, hi) in SITES.items():
        for d in ('endogenous', 'WASF1_full', 'WASF12'):
            n = k = 0
            for r in mol[d].itertuples():
                at = np.array(json.loads(r.at), int)
                if not at.size or at.min() > lo or at.max() < hi:
                    continue
                n += 1
                k += any(15 <= ln <= 89 and min(e, hi) - max(s_, lo) + 1 >= 10 for s_, e, ln in json.loads(r.fp))
            site_rows.append(dict(site=site, start=lo, end=hi, dataset=d, molecules=n, footprinted=k,
                                  share=round(k / n, 3) if n else np.nan))
    pd.DataFrame(site_rows).to_csv(out / 'site_footprints.tsv', sep='\t', index=False)
    print(pd.DataFrame(site_rows).to_string(index=False))
    st, co = pd.DataFrame(rows), pd.DataFrame(cors)
    st.to_csv(out / 'region_stats.tsv', sep='\t', index=False)
    co.to_csv(out / 'm6a_concordance.tsv', sep='\t', index=False)
    print(st.to_string(index=False))
    print(co.to_string(index=False))

    # gene model and motifs
    ref = json.load(open(w / 'raw/ucsc_ncbiRefSeq_wasf1.json'))['ncbiRefSeq']
    genes = [next(r for r in ref if r['name'] == 'NM_003931.3'), next(r for r in ref if r['name2'] == 'CDC40')]
    mot = pd.read_csv(w / 'motifs/motifs.tsv', sep='\t')

    plt.rcParams['axes.formatter.useoffset'] = False
    kb = matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:,.0f}')
    fig = plt.figure(figsize=(7.2, 8.4))
    gs = fig.add_gridspec(9, 2, height_ratios=[0.7, 1, 1, 0.25, 0.7, 1, 1, 1, 0.3], hspace=0.35, wspace=0.25)
    # whole 4.5-kb window
    axg = fig.add_subplot(gs[0, :])
    axg.set_xlim(X0, X1); axg.set_ylim(0, 1); axg.axis('off')
    axg.plot([X0, X1], [0.35, 0.35], color='#666', lw=0.6)
    for g, y in zip(genes, (0.35, 0.35)):
        for s, e in zip(map(int, g['exonStarts'].split(',')[:-1]), map(int, g['exonEnds'].split(',')[:-1])):
            if e > X0 and s < X1:
                axg.add_patch(Rectangle((max(s, X0), y - 0.12), min(e, X1) - max(s, X0), 0.24, color='#333'))
        tss = g['txEnd'] if g['strand'] == '-' else g['txStart']
        d = -1 if g['strand'] == '-' else 1
        axg.annotate('', xy=(tss + d * 250, 0.75), xytext=(tss, 0.75), arrowprops=dict(arrowstyle='-|>', lw=0.8))
        axg.plot([tss, tss], [0.35, 0.75], color='k', lw=0.8)
        axg.text(tss + d * 270, 0.75, g['name2'], ha='right' if d < 0 else 'left', va='center', fontsize=6.5, style='italic')
    for c, (name, lo, hi) in REGIONS.items():
        axg.add_patch(Rectangle((lo, 0.02), hi - lo, 0.13, color='#56B4E9'))
        axg.text((lo + hi) / 2, -0.1, name, ha='center', va='top', fontsize=6)
    axg.set_title('WASF1 / CDC40 bidirectional promoter, chr6:110,176,895-110,181,442 (4.5-kb construct insert)', fontsize=7)
    for i, (d, lab) in enumerate(DATASETS):
        ax = fig.add_subplot(gs[1 + i, :], sharex=axg)
        t = T[d]
        for c in ('nuc', 'fp_61_89', 'fp_31_60', 'fp_15_30'):
            k = 1 if c == 'nuc' else 5
            ax.fill_between(t.pos, k * t[c], color=BIN_C[c], lw=0, alpha=0.8,
                            label=f'{c.replace("fp_", "").replace("_", "-")} bp' + (' x5' if k > 1 else '') if i == 0 else None)
        ax2 = ax.twinx()
        ax2.plot(t.pos, t.m6a, color='#0072B2', lw=0.6)
        ax2.set_ylim(0, 0.8); ax2.tick_params(labelsize=5.5, colors='#0072B2')
        ax.set_ylim(0, 1); ax.tick_params(labelsize=5.5); ax.xaxis.set_major_formatter(kb)
        ax.set_ylabel(f'{lab}\n(n = {len(mol[d])})', fontsize=6)
        for _, lo, hi in REGIONS.values():
            ax.axvspan(lo, hi, color='#56B4E9', alpha=0.12, lw=0)
        if i == 0:
            ax.legend(fontsize=5, ncol=4, frameon=False, loc='upper left')
            ax2.set_ylabel('m6A fraction', fontsize=5.5, color='#0072B2')
    ax.tick_params(labelbottom=True)
    # zooms on the two 300-bp regions
    for j, (c, (name, lo, hi)) in enumerate(REGIONS.items()):
        axm = fig.add_subplot(gs[4, j])
        axm.set_xlim(lo, hi); axm.set_ylim(0, 1); axm.axis('off')
        s = mot[(mot.start >= lo) & (mot.end <= hi)].sort_values('p')
        best = s.groupby('tf').first().sort_values('p').head(6).reset_index()
        for k, r in best.iterrows():
            y = 0.85 - k * 0.14
            axm.add_patch(Rectangle((r.start, y - 0.05), r.end - r.start, 0.1, color='#6A4FC9'))
            axm.text(r.end + 3, y, f'{r.tf} ({r.p:.0e})', fontsize=4.8, va='center')
        axm.set_title(f'{name}: top motif matches (P < 1e-5; {len(s)} matches, {s.tf.nunique()} TFs)', fontsize=6)
        for i, (d, lab) in enumerate([('endogenous', 'Endogenous'), ('WASF1_full', '4.5-kb construct'),
                                      (c, '300-bp construct')]):
            ax = fig.add_subplot(gs[5 + i, j])
            t = T[d][(T[d].pos >= lo) & (T[d].pos <= hi)]
            for b in ('nuc', 'fp_61_89', 'fp_31_60', 'fp_15_30'):
                ax.fill_between(t.pos, t[b], color=BIN_C[b], lw=0, alpha=0.8)
            ax2 = ax.twinx()
            ax2.plot(t.pos, t.m6a, color='#0072B2', lw=0.7)
            ax2.set_ylim(0, 0.8); ax2.tick_params(labelsize=5, colors='#0072B2')
            ax.set_xlim(lo, hi); ax.set_ylim(0, 1); ax.tick_params(labelsize=5); ax.xaxis.set_major_formatter(kb)
            r = st[(st.region == name) & (st.dataset == d)].iloc[0]
            ax.set_title(f'{lab}: {r.molecules_spanning} molecules, {r.share_15_89:.0%} with a 15-89 bp footprint',
                         fontsize=5.5, loc='left')
            if i < 2:
                ax.tick_params(labelbottom=False)
    fig.text(0.5, 0.06, 'Shares unscaled in the zooms (x5 for sub-nucleosome bins in the 4.5-kb tracks); blue line = m6A '
             'fraction (accessibility). Exploratory first look, not a figure draft.', ha='center', fontsize=5.5)
    fig.savefig(out / 'overview.png', dpi=220, bbox_inches='tight')


if __name__ == '__main__':
    main()

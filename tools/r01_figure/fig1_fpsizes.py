"""Figure 1: which sub-nucleosome footprint sizes are present in a 300-bp region (DPAGT1 FIRE peak 1, or WASF1 FIRE peak 4
with --gene WASF1); exploratory scan to set data-driven size bins (fig1_rules.json: fp_bins_v2, wasf1.fp_size_scan).

Per dataset (full-length construct, 300-bp construct, endogenous): for each kept molecule, the share of the region
PEAK1 covered by nucleosome calls (>= 90 bp) and its sub-nucleosome footprints (< 90 bp) whose midpoint lies in
the region. "Accessible" molecule = nucleosome calls cover < ACC_MAX of the region. Writes footprints.tsv (one row per
sub-nucleosome footprint: dataset, read, accessible, start, end, length, midpoint), molecules.tsv (per molecule:
nucleosome share, footprint count), length histograms (all and accessible molecules) and position-by-length plots.
Usage: python tools/r01_figure/fig1_fpsizes.py --profiles <fig1/profiles> --motifs <fig1/motifs.tsv> --output <dir>
           [--gene DPAGT1|WASF1]
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

NUC_MIN, ACC_MAX = 90, 0.5
GENES = {   # region, datasets, motif sites (None = CTCF matches from --motifs)
    'DPAGT1': ((119100271, 119100570), ['DPAGT1_full', 'DPAGT11', 'endogenous'], None),
    'WASF1': ((110179911, 110180211), ['WASF1_full', 'WASF12', 'endogenous'],
              [('GC box', 110180054, 110180097), ('ZNF444 site', 110179973, 110179991),
               ('ZNF143 site', 110180105, 110180127)]),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--profiles', required=True)
    ap.add_argument('--motifs', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--gene', default='DPAGT1', choices=sorted(GENES))
    a = ap.parse_args()
    PEAK1, DATASETS, SITES = GENES[a.gene]
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    L = PEAK1[1] - PEAK1[0] + 1
    fps, mols = [], []
    for ds in DATASETS:
        M = pd.read_parquet(Path(a.profiles) / f'{ds}.molecules.parquet')
        for r in M.itertuples():
            f = json.loads(r.fp)
            nuc = np.zeros(L, bool)
            covered = False
            for s, e, ln in f:
                if e >= PEAK1[0] and s <= PEAK1[1]:
                    covered = True
                if ln >= NUC_MIN:
                    nuc[max(s, PEAK1[0]) - PEAK1[0]:max(min(e, PEAK1[1]) - PEAK1[0] + 1, 0)] = True
            at = np.array(json.loads(r.at), int)
            if not ((at >= PEAK1[0]).any() and (at <= PEAK1[1]).any() and at.min() <= PEAK1[0] + 20
                    and at.max() >= PEAK1[1] - 20):
                continue                                    # molecule must span the region
            share = nuc.mean()
            sub = [(s, e, ln) for s, e, ln in f if ln < NUC_MIN and PEAK1[0] <= (s + e) / 2 <= PEAK1[1]]
            mols.append(dict(dataset=ds, read=r.read, nuc_share=share, accessible=share < ACC_MAX, n_sub=len(sub)))
            fps += [dict(dataset=ds, read=r.read, accessible=share < ACC_MAX, start=s, end=e, length=ln,
                         mid=(s + e) / 2) for s, e, ln in sub]
    F, Mo = pd.DataFrame(fps), pd.DataFrame(mols)
    F.to_csv(out / 'footprints.tsv', sep='\t', index=False)
    Mo.to_csv(out / 'molecules.tsv', sep='\t', index=False, float_format='%.4g')
    s = Mo.groupby('dataset').agg(molecules=('read', 'size'), accessible=('accessible', 'sum'),
                                  accessible_with_fp=('n_sub', lambda x: int(((x > 0) & Mo.loc[x.index, 'accessible']).sum())))
    s['footprints_all'] = F.groupby('dataset').size()
    s['footprints_accessible'] = F[F.accessible].groupby('dataset').size()
    s.to_csv(out / 'summary.tsv', sep='\t')
    hist = F[F.accessible].groupby(['dataset', 'length']).size().unstack(0, fill_value=0)
    hist.to_csv(out / 'length_hist_accessible.tsv', sep='\t')      # 1-bp length histogram, accessible molecules
    print(s.to_string())
    for ds in DATASETS:
        x = F[(F.dataset == ds) & F.accessible].length
        if len(x):
            print(ds, 'accessible footprint length quantiles (5,25,50,75,95):',
                  np.percentile(x, [5, 25, 50, 75, 95]).round(1).tolist())
            h = np.bincount(x, minlength=NUC_MIN)[:NUC_MIN]
            print('  counts by 5-bp:', {f'{i}-{i + 4}': int(h[i:i + 5].sum()) for i in range(0, NUC_MIN, 5) if h[i:i + 5].sum()})
    mot = pd.read_csv(a.motifs, sep='\t')
    if SITES is None:
        mot = mot[(mot.start >= PEAK1[0]) & (mot.end <= PEAK1[1]) & mot.tf.str.contains('CTCF')]
        mot = mot.drop_duplicates(['start', 'end']).assign(name='CTCF motif', tag='ctcf')
    else:
        mot = pd.DataFrame(SITES, columns=['name', 'start', 'end']).assign(tag='site')
    # per-site view (accessible molecules): footprints >= 10 bp by 10-bp midpoint window, and those over the CTCF motif
    acc = F[F.accessible & (F.length >= 10)].copy()
    acc['win'] = ((acc.mid - PEAK1[0]) // 10 * 10 + PEAK1[0]).astype(int)
    site = acc.groupby(['dataset', 'win']).length.describe(percentiles=[0.25, 0.5, 0.75])[['count', '25%', '50%', '75%']]
    site.to_csv(out / 'sites.tsv', sep='\t', float_format='%.4g')
    for m in mot.itertuples():
        ov = acc[(np.minimum(acc.end, m.end) - np.maximum(acc.start, m.start) + 1) >= 10]
        print(f'{m.name} {m.start}-{m.end}: footprints >= 10 bp overlapping it by >= 10 bp, accessible molecules')
        for ds, g in ov.groupby('dataset'):
            n_acc = int(Mo[(Mo.dataset == ds) & Mo.accessible].shape[0])
            print(f'  {ds}: {len(g)} footprints on {g.read.nunique()} of {n_acc} accessible molecules; length '
                  f'quantiles (10,25,50,75,90) {np.percentile(g.length, [10, 25, 50, 75, 90]).round(0).tolist()}')
        ov.to_csv(out / f'{m.tag}_{m.start}_footprints.tsv', sep='\t', index=False)
    fig, axs = plt.subplots(2, len(DATASETS), figsize=(11, 6), gridspec_kw=dict(height_ratios=[1, 1.4]))
    for j, ds in enumerate(DATASETS):
        sub = F[F.dataset == ds]
        bins = np.arange(0, NUC_MIN + 1, 2)
        axs[0, j].hist(sub.length, bins, color='#BBBBBB', label='all molecules')
        axs[0, j].hist(sub[sub.accessible].length, bins, color='#2A78D6', label='accessible molecules')
        axs[0, j].set(title=f'{ds}', xlabel='footprint length (bp)', ylabel='footprints (midpoint in 300-bp region)')
        axs[0, j].legend(fontsize=7, frameon=False)
        acc = sub[sub.accessible]
        axs[1, j].scatter(acc.mid, acc.length, s=3, alpha=0.3, color='#2A78D6', lw=0)
        for m in mot.itertuples():
            axs[1, j].axvspan(m.start, m.end, color='#E77C34', alpha=0.25, lw=0)
        axs[1, j].set(xlim=PEAK1, ylim=(0, NUC_MIN), xlabel='footprint midpoint (hg38)', ylabel='length (bp)')
        axs[1, j].ticklabel_format(useOffset=False, style='plain', axis='x')
        axs[1, j].tick_params(axis='x', labelsize=6)
    axs[1, 0].text(0.01, 0.97, 'orange: ' + ', '.join(mot.name.unique()), transform=axs[1, 0].transAxes, va='top', fontsize=7)
    fig.tight_layout()
    fig.savefig(out / 'fp_sizes.png', dpi=150)


if __name__ == '__main__':
    main()

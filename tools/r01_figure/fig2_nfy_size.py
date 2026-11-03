"""Figure 2: what size is the footprint over SLC39A4 c.-169 (NF-Y)? Per-read length of the FiberHMM footprint covering
the variant base (reads.tsv.gz site_fp, 0 = none), WT vs c.-169A>G, kept reads, per replicate and pooled.

The WT-specific footprint = WT excess over variant at each length (% of reads). Reports the excess by length, its
mode, the central range holding 80% of the pooled excess, and per-replicate excess within candidate bins. Writes
size_by_allele.tsv, excess_summary.tsv and nfy_size.png.
Usage: python tools/r01_figure/fig2_nfy_size.py --project <SLC39A4_variant_slide> --output <dir>
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

MAXL = 120
CANDIDATES = [(1, 9), (10, 19), (1, 60), (61, 120), (10, 60), (20, 60), (15, 30), (31, 60), (25, 50), (30, 50)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--project', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    r = pd.read_csv(Path(a.project) / 'analysis/reads.tsv.gz', sep='\t')
    r = r[r.kept == 1]
    rows = []
    for (rep, al), g in list(r.groupby(['rep', 'allele'])) + [(('pooled', al), g) for al, g in r.groupby('allele')]:
        h = np.bincount(np.clip(g.site_fp, 0, MAXL), minlength=MAXL + 1)
        rows += [dict(rep=rep, allele=al, length=L, reads=int(h[L]), pct=100 * h[L] / len(g)) for L in range(1, MAXL + 1)]
    T = pd.DataFrame(rows)
    T.to_csv(out / 'size_by_allele.tsv', sep='\t', index=False, float_format='%.5g')
    W = T.pivot_table(index=['rep', 'length'], columns='allele', values='pct').reset_index()
    W['excess'] = W.WT - W.Var
    pe = W[W.rep == 'pooled'].set_index('length').excess
    pos = pe.clip(lower=0)
    cum = pos.cumsum() / pos.sum()
    lo, hi = int(cum.index[cum >= 0.1][0]), int(cum.index[cum >= 0.9][0])
    print(f'pooled WT excess: mode {int(pe.idxmax())} bp; central 80% of positive excess {lo}-{hi} bp; '
          f'median {int(cum.index[cum >= 0.5][0])} bp')
    summ = []
    for c0, c1 in CANDIDATES:
        for rep, g in W.groupby('rep'):
            s = g[g.length.between(c0, c1)]
            summ.append(dict(bin=f'{c0}-{c1}', rep=rep, wt_pct=s.WT.sum(), var_pct=s.Var.sum(),
                             excess_pct=s.excess.sum(), fold=s.WT.sum() / s.Var.sum()))
    S = pd.DataFrame(summ)
    S.to_csv(out / 'excess_summary.tsv', sep='\t', index=False, float_format='%.4g')
    print(S[S.rep == 'pooled'].to_string(index=False))
    fig, ax = plt.subplots(2, 1, figsize=(6, 4.5), sharex=True)
    p = T[T.rep == 'pooled']
    for al, c in (('WT', '#2E3192'), ('Var', '#E3191F')):
        q = p[p.allele == al]
        ax[0].bar(q.length, q.pct, width=1, color=c, alpha=0.5, label=al)
    ax[0].set_ylabel('% of reads, footprint\nof this length over c.-169')
    ax[0].legend(frameon=False)
    for rep, c in zip(['rep_1', 'rep_2', 'rep_3', 'pooled'], ['#CCCCCC', '#999999', '#555555', '#338B8C']):
        g = W[W.rep == rep]
        ax[1].plot(g.length, g.excess, color=c, lw=1.5 if rep == 'pooled' else 0.8, label=rep)
    ax[1].axhline(0, color='k', lw=0.5)
    ax[1].axvspan(lo, hi, color='#338B8C', alpha=0.12, lw=0)
    ax[1].set_ylabel('WT − variant\n(% of reads)')
    ax[1].set_xlabel('Footprint length (bp)')
    ax[1].legend(frameon=False, ncol=4, fontsize=7)
    ax[1].set_xlim(0, MAXL)
    fig.tight_layout()
    fig.savefig(out / 'nfy_size.png', dpi=150)


if __name__ == '__main__':
    main()

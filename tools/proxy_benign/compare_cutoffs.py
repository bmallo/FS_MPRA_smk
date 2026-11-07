"""Compare AlphaGenome proxy benign calls in the 350 bp window with the top-20% cutoff taken from the 1 kb element (v2)
or from the 350 bp window itself (rules.json v2_alphagenome_350bp.cutoff_350_comparison).

Reads <root>/v2_350bp/snv_350bp.tsv (proxy_350.py). Writes <root>/v2_350bp/snv_350bp_cutoffs.tsv,
cutoff_comparison.tsv and counts_cut350.tsv (same columns as counts.tsv, for table_figure.py --v2 --cut350).
Usage: python tools/proxy_benign/compare_cutoffs.py --root results/proxy_benign_20261106
"""
import argparse
from pathlib import Path

import pandas as pd

TOP = 0.8


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    a = ap.parse_args()
    d = Path(a.root) / 'v2_350bp'
    s = pd.read_csv(d / 'snv_350bp.tsv', sep='\t', low_memory=False)
    s['clinvar_benchmark'] = s.clinvar_benchmark.fillna('')
    out, rows, counts = [], [], []
    for g, x in s.groupby('window', sort=False):
        x = x.copy()
        for c in ('ag_rna', 'ag_dnase'):
            x[f'{c}_cut1kb'] = x[c].abs().where(x[f'{c}_large']).min()     # smallest |score| flagged in the 1 kb rule
            x[f'{c}_cut350'] = x[c].abs().quantile(TOP)
            x[f'{c}_large350'] = x[c].abs() >= x[f'{c}_cut350']
        x['pb_1kb'] = x.ag_no_change
        x['pb_350'] = ~x.ag_rna_large350 & ~x.ag_dnase_large350
        out.append(x)
        P = x[x.clinvar_benchmark == '']
        rows.append(dict(element=g, pb_1kb_cutoff=int(P.pb_1kb.sum()), pb_350_cutoff=int(P.pb_350.sum()),
                         both=int((P.pb_1kb & P.pb_350).sum()), only_1kb=int((P.pb_1kb & ~P.pb_350).sum()),
                         only_350=int((~P.pb_1kb & P.pb_350).sum()),
                         rna_cut_1kb=x.ag_rna_cut1kb.iloc[0], rna_cut_350=x.ag_rna_cut350.iloc[0],
                         dnase_cut_1kb=x.ag_dnase_cut1kb.iloc[0], dnase_cut_350=x.ag_dnase_cut350.iloc[0]))
        counts.append(dict(element=g, possible_snvs=len(x), proxy_benign=int(P.pb_350.sum()),
                           proxy_benign_pct=round(100 * P.pb_350.sum() / len(x), 1)))
    pd.concat(out).to_csv(d / 'snv_350bp_cutoffs.tsv', sep='\t', index=False, float_format='%.6g')
    r = pd.DataFrame(rows)
    r.to_csv(d / 'cutoff_comparison.tsv', sep='\t', index=False, float_format='%.4g')
    pd.DataFrame(counts).to_csv(d / 'counts_cut350.tsv', sep='\t', index=False)
    print(r.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    print('pooled: 1kb', r.pb_1kb_cutoff.sum(), '350', r.pb_350_cutoff.sum(), 'both', r.both.sum(),
          'only 1kb', r.only_1kb.sum(), 'only 350', r.only_350.sum())


if __name__ == '__main__':
    main()

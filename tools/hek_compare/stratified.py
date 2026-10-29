"""Element vs whole-plasmid protection by distal-m6A stratum, HEK293T vs HepG2 (C-class WT, LLR min_llr 5).

Separates a labeling-level effect (both cells move together within a stratum) from a cell
effect (gap persists within strata), and an element-specific effect from a plasmid-wide one.
Usage: python tools/hek_compare/stratified.py --out results/hek293t_vs_hepg2_20261028
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

CELLS = {'HepG2': (Path('results/llr_vs_fiberhmm_20261026'), 4647, (1234, 1445),
                   Path('results/rep3_site_methylation_20261021/extract/molecules.tsv.gz')),
         'HEK293T': (Path('results/hek293t_llr_20261028'), 4718, (3236, 3447),
                     Path('results/hek293t_llr_20261028/data/molecules.tsv.gz'))}
BINS = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', required=True)
    out = Path(ap.parse_args().out)
    rows = []
    for cell, (W, L, (a, b), molp) in CELLS.items():
        ids = [l.strip() for l in open(W / 'config/wt_chrom.read_ids.txt') if l.strip()]
        df = pd.read_parquet(W / 'results/calls/wt_chrom.llr_global_foldA_pp0_dd25.parquet', columns=['min_llr', 'read_name', 'start', 'end'])
        df = df[np.isclose(df.min_llr, 5.0)]
        s0, e = df.start.values - 1, df.end.values
        df = df.assign(whole=e - s0, elem=np.clip(np.minimum(e, b) - np.maximum(s0, a), 0, None))
        g = df.groupby('read_name')[['whole', 'elem']].sum().reindex(ids).fillna(0)
        mol = pd.read_csv(molp, sep='\t', usecols=['read_id', 'distal_m6a_fraction']).set_index('read_id')
        g['distal'] = mol.distal_m6a_fraction.reindex(g.index).values
        g['whole_frac'], g['elem_frac'] = g.whole / L, g.elem / (b - a)
        g['bin'] = pd.cut(g.distal, BINS)
        for k, d in g.groupby('bin', observed=True):
            rows.append(dict(cell=cell, distal_bin=str(k), reads=len(d), element_protected=d.elem_frac.mean(),
                             plasmid_protected=d.whole_frac.mean()))
        rows.append(dict(cell=cell, distal_bin='all C', reads=len(g), element_protected=g.elem_frac.mean(),
                         plasmid_protected=g.whole_frac.mean()))
    t = pd.DataFrame(rows)
    t.to_csv(out / 'stratified_protection.tsv', sep='\t', index=False)
    print(t.round(3).to_string(index=False))


if __name__ == '__main__':
    main()

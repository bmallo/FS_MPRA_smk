"""Figure 1: how many molecules carry a CTCF-sized footprint (fp_bins_v2 31-60 bp) over the DPAGT1 CTCF motif, per dataset.

A molecule counts if it spans the motif (observed A/T on both sides within 20 bp) and has a 31-60 bp call overlapping
the motif (chr11:119,100,436-455) by >= MIN_OVERLAP bp. Also reports calls of any sub-nucleosome size over the motif
and the length/position of the 31-60 bp calls. Writes ctcf_reads.tsv (one row per qualifying call) and summary.tsv.
Usage: python tools/r01_figure/fig1_ctcf_reads.py --profiles <fig1/profiles> --output <dir>
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

CTCF = (119100436, 119100455)
MIN_OVERLAP, SLACK = 10, 20
DATASETS = ['DPAGT1_full', 'DPAGT11', 'endogenous']


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--profiles', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    rows, summ = [], []
    for ds in DATASETS:
        M = pd.read_parquet(Path(a.profiles) / f'{ds}.molecules.parquet')
        span = ctcf = anysub = 0
        for r in M.itertuples():
            at = np.array(json.loads(r.at), int)
            if not at.size or at.min() > CTCF[0] - SLACK or at.max() < CTCF[1] + SLACK:
                continue
            span += 1
            over = [(s, e, ln) for s, e, ln in json.loads(r.fp)
                    if min(e, CTCF[1]) - max(s, CTCF[0]) + 1 >= MIN_OVERLAP]
            anysub += any(15 <= ln < 90 for *_, ln in over)
            hit = [(s, e, ln) for s, e, ln in over if 31 <= ln <= 60]
            ctcf += bool(hit)
            rows += [dict(dataset=ds, read=r.read, start=s, end=e, length=ln) for s, e, ln in hit]
        summ.append(dict(dataset=ds, molecules=len(M), spanning_motif=span, with_15_89_over_motif=anysub,
                         with_31_60_over_motif=ctcf, share=ctcf / span if span else np.nan))
    S = pd.DataFrame(summ)
    S.to_csv(out / 'summary.tsv', sep='\t', index=False, float_format='%.4g')
    pd.DataFrame(rows).to_csv(out / 'ctcf_reads.tsv', sep='\t', index=False)
    print(S.to_string(index=False))
    R = pd.DataFrame(rows)
    for ds, g in R.groupby('dataset'):
        print(ds, 'lengths', np.percentile(g.length, [10, 50, 90]).tolist(), 'start median', g.start.median(),
              'end median', g.end.median())


if __name__ == '__main__':
    main()

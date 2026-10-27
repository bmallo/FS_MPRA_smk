"""Pick WT read sets for the plasmid-LLR vs FiberHMM comparison.

Same-reaction proxy controls, classified on distal m6A only (as in the
site-methylation stage), from analysis-eligible WT molecules:
  wt_unmeth : U class (distal <= 0.01)  -> stands in for naked, no Hia5   (p_p)
  wt_meth   : M class (distal >= 0.60)  -> stands in for naked + Hia5    (p_a)
  wt_chrom  : C class (in between)      -> the chromatin sample to call

Usage: python tools/llr_compare/select_reads.py --workdir results/llr_vs_fiberhmm_20261026
"""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

MOLECULES = Path('results/rep3_site_methylation_20261021/extract/molecules.tsv.gz')
U_MAX, M_MIN = 0.01, 0.60
SEED = 20261026


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workdir', required=True)
    ap.add_argument('--n-control', type=int, default=93030)
    ap.add_argument('--n-chrom', type=int, default=60000)
    a = ap.parse_args()
    wd = Path(a.workdir)
    (wd / 'config').mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(MOLECULES, sep='\t', usecols=['read_id', 'observed_variant_id', 'footprint_analysis_eligible', 'distal_m6a_fraction'])
    variants = json.load(open('results/rep3_site_methylation_20261021/extract/variants.json'))
    wt_id = next(v['variant_id'] for v in variants if v['class'] == 'WT')
    wt = df[(df.observed_variant_id == wt_id) & (df.footprint_analysis_eligible == 1)]
    rng = np.random.default_rng(SEED)
    sets = {'wt_unmeth': wt[wt.distal_m6a_fraction <= U_MAX], 'wt_meth': wt[wt.distal_m6a_fraction >= M_MIN],
            'wt_chrom': wt[(wt.distal_m6a_fraction > U_MAX) & (wt.distal_m6a_fraction < M_MIN)]}
    n = {'wt_unmeth': a.n_control, 'wt_meth': a.n_control, 'wt_chrom': a.n_chrom}
    summary = {}
    for name, d in sets.items():
        take = d.read_id.values if len(d) <= n[name] else rng.choice(d.read_id.values, n[name], replace=False)
        with open(wd / 'config' / f'{name}.read_ids.txt', 'w') as h:
            h.write('\n'.join(sorted(take)) + '\n')
        summary[name] = dict(available=int(len(d)), selected=int(len(take)))
    with open(wd / 'config' / 'samples.tsv', 'w') as h:
        h.write('sample\trole\nwt_meth\tcontrol\nwt_unmeth\tcontrol\nwt_chrom\tchromatin\n')
    json.dump(dict(source=str(MOLECULES), wt_variant_id=wt_id, u_max=U_MAX, m_min=M_MIN, seed=SEED, sets=summary),
              open(wd / 'config' / 'selection.json', 'w'), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()

"""Are SNV molecules comparable to WT molecules on read-level characteristics?

For every SNV cohort (HepG2, eligible + evaluated molecules) and each read covariate, the
standardized mean difference SMD = (mean_SNV - mean_WT) / SD_WT. Technical covariates (passes,
read quality, mapq, circular) should not differ; labeling (distal m6A) and class mix are handled
by the U/M/C-strata standardization but are reported too. Each SMD is then correlated with the
SNV's chromatin z (msp, nuc, local A_norm) to see whether imbalance tracks the effect calls.

Usage: python tools/variant_effects/covariate_balance.py --dir results/variant_effects_20261028 --output <dir>
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

COVS = ['passes', 'read_quality', 'mapq', 'circular', 'distal_m6a', 'whole_m6a', 'distal_nucleosome']


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    V, out = Path(a.dir), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    ext = V / 'hepg2/extract_llr'
    m = np.load(ext / 'meta.npz')
    variants = json.load(open(ext / 'variants.json'))
    wt_i = next(i for i, v in enumerate(variants) if v['class'] == 'WT')
    keep = m['eligible'].astype(bool) & m['evaluated'].astype(bool)
    vi = m['variant_index'].astype(int)
    df = pd.DataFrame({c: m[c].astype(float) for c in COVS})
    df['vi'], df['cluster'] = vi, m['cluster']
    df = df[keep]
    wt = df[df.vi == wt_i]
    mu, sd = wt[COVS].mean(), wt[COVS].std()
    g = df[df.vi != wt_i].groupby('vi')
    smd = ((g[COVS].mean() - mu) / sd)
    smd['molecules'] = g.size()
    smd['clusters'] = g.cluster.nunique()
    smd['source_variant'] = [variants[i].get('source_variant') for i in smd.index]
    smd = smd[smd.molecules >= 300]
    cls = pd.read_csv(V / 'classes/variants.tsv', sep='\t').set_index('source_variant')
    smd = smd.join(cls[['oligo_pos1', 'hepg2_msp_z', 'hepg2_nuc_z', 'hepg2_tf_local_z', 'chromatin_hit', 'functional_hit']], on='source_variant')
    smd = smd[smd.oligo_pos1.notna()]
    smd.to_csv(out / 'smd_per_snv.tsv', sep='\t', index=False)
    rows = []
    for c in COVS:
        r = dict(covariate=c, wt_mean=mu[c], wt_sd=sd[c], smd_median=smd[c].median(), smd_p05=smd[c].quantile(.05),
                 smd_p95=smd[c].quantile(.95), snvs_abs_smd_gt_0_1=int((smd[c].abs() > .1).sum()), snvs=len(smd))
        for z in ('hepg2_msp_z', 'hepg2_nuc_z', 'hepg2_tf_local_z'):
            ok = smd[c].notna() & smd[z].notna()
            r[f'rho_{z}'] = stats.spearmanr(smd.loc[ok, c], smd.loc[ok, z])[0]
        h = smd.chromatin_hit == True
        r['smd_median_chromatin_hits'] = smd.loc[h, c].median()
        r['smd_median_non_hits'] = smd.loc[~h, c].median()
        rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(out / 'balance_summary.tsv', sep='\t', index=False)
    pd.set_option('display.width', 250)
    print(t.round(3).to_string(index=False))
    print('clusters per SNV: median', smd.clusters.median(), 'p10', smd.clusters.quantile(.1))


if __name__ == '__main__':
    main()

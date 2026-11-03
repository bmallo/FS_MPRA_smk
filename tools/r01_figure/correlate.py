"""Correlations between evidence types for LDLR promoter SNVs (pairs fixed in correlation_rules.json).

Writes pairs.tsv (Spearman rho, bootstrap CI, P, BH q, n per pair x subset; pairs whose columns are missing are
skipped and listed) and clinical.tsv (AUROC P/LP vs B/LB per score, plus each control's value) to --output.
Usage: python tools/r01_figure/correlate.py --evidence <evidence.tsv> --output <dir>
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).parent
RULES = json.load(open(HERE / 'correlation_rules.json'))
SEED, NBOOT = 20261102, 2000
PLP = ['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic']
BLB = ['Benign', 'Likely benign']
# sign so that larger = more pathogenic-like (loss of expression, gain of protection)
CLIN_SCORES = {'act': -1, 'nuc_D_pts': 1, 'local_inacc_pts': 1, 'C1': 1, 'ag_rna_seq_hepg2': -1,
               'ag_rna_seq_liver': -1, 'ag_cage_hepg2': -1, 'ag_dnase_hepg2': -1, 'pai_score': -1}


def spearman_boot(x, y, rng):
    r = stats.spearmanr(x, y)
    n = len(x)
    bs = []
    for _ in range(NBOOT):
        i = rng.integers(0, n, n)
        bs.append(stats.spearmanr(x[i], y[i]).correlation)
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    return r.correlation, lo, hi, r.pvalue


def auroc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    if not len(pos) or not len(neg):
        return np.nan
    return float(((pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum())
                 / (len(pos) * len(neg)))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--evidence', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    t = pd.read_csv(a.evidence, sep='\t')
    rng = np.random.default_rng(SEED)
    subsets = {'all': t.index == t.index, 'sp1_site_143_151': t.oligo.between(143, 151),
               'outside_138_157': ~t.oligo.between(138, 157)}
    rows, skipped = [], []
    for group, pairs in RULES['pairs'].items():
        for x, y, label in pairs:
            if x not in t or y not in t:
                skipped.append(label)
                continue
            for sname, m in subsets.items():
                s = t[m]
                if 'hek' in x or 'hek' in y:
                    s = s[s.hek_clones >= 3]
                s = s[[x, y]].dropna()
                if len(s) < 8:
                    continue
                rho, lo, hi, p = spearman_boot(s[x].values, s[y].values, rng)
                rows.append(dict(group=group, pair=label, x=x, y=y, subset=sname, n=len(s), rho=rho, lo=lo, hi=hi, p=p))
    r = pd.DataFrame(rows)
    if len(r):
        ok = r.p.notna()
        r.loc[ok, 'q'] = stats.false_discovery_control(r.loc[ok, 'p'])
    r.to_csv(out / 'pairs.tsv', sep='\t', index=False, float_format='%.4g')

    crow = []
    for col, sgn in CLIN_SCORES.items():
        if col not in t:
            continue
        p = t[t.clinvar.isin(PLP)][['snv', col]].dropna()
        b = t[t.clinvar.isin(BLB)][['snv', col]].dropna()
        crow.append(dict(score=col, n_plp=len(p), n_blb=len(b), auroc=auroc(sgn * p[col], sgn * b[col]),
                         plp=';'.join(f'{v}:{x:.3g}' for v, x in zip(p.snv, p[col])),
                         blb=';'.join(f'{v}:{x:.3g}' for v, x in zip(b.snv, b[col]))))
    pd.DataFrame(crow).to_csv(out / 'clinical.tsv', sep='\t', index=False, float_format='%.4g')
    json.dump(dict(skipped=skipped), open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

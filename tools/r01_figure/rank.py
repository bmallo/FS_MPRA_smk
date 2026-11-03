"""R01 figure, step 2 (shortlist): score SNV x footprint pairs for the activator and repressor panels using the rules
fixed in tools/r01_figure/selection_rules.json before this was run. No new tests: uses stored fp_features and A10 results.

Writes ranked.tsv (all pool pairs with every criterion and the score) and summary.json to --output.
Usage: python tools/r01_figure/rank.py --root results/analysis_v2_20261028 --step1 <step1 dir> --output <dir>
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

TF_FOOTPRINTS = ['F04', 'F07', 'F10', 'F13', 'F15', 'F18', 'F21', 'F22']


def caller_free(root, pairs, sub):
    """Mean A10 E1 D (inaccessible fraction, variant - WT) over testable adenines inside each footprint span."""
    snvs = pairs.snv.unique().tolist()
    e = pd.read_parquet(root / f'a10_effects/{sub}/effects.parquet', columns=['snv', 'endpoint', 'oligo', 'D', 'p'],
                        filters=[('endpoint', '==', 'E1'), ('snv', 'in', snvs)])
    e = e[e.p.notna()]
    out = []
    for r in pairs.itertuples():
        g = e[(e.snv == r.snv) & e.oligo.between(r.span_lo, r.span_hi)]
        out.append((100 * g.D.mean() if len(g) else np.nan, len(g)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--step1', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    rules = json.load(open(Path(__file__).with_name('selection_rules.json')))

    p = pd.read_csv(Path(a.step1) / 'snv_footprint.tsv', sep='\t')
    n0 = len(p)
    p = p[p.footprint.isin(TF_FOOTPRINTS)]
    n1 = len(p)
    p = p[(p.act_p < 0.05) & (np.sign(p.act) == np.sign(p.act2)) & (p.act != 0)].copy()
    n2 = len(p)
    p['panel'] = np.where(np.sign(p.D_A7_confident) == np.sign(p.act), 'activator', 'repressor')
    p['D_pts'] = 100 * p.D_A7_confident
    p['D_pm5_pts'] = 100 * p.D_pm5_confident

    p['pt_activity'] = np.where(p.act.abs() >= 1, 2, np.where(p.act.abs() >= 0.5, 1, 0))
    p['pt_footprint'] = np.where(p.D_pts.abs() >= 2, 2, np.where(p.D_pts.abs() >= 1, 1, 0))
    p['pt_confident'] = (p.confident_share >= 0.8).astype(int)
    same = np.sign(p.D_pm5_confident) == np.sign(p.D_A7_confident)
    pm5sig = p.sig_pm5_confident.fillna(False).astype(bool)
    p['pt_mask'] = np.where(pm5sig & same, 2, np.where(same & (p.p_pm5_confident < 0.1), 1, 0))
    p['pt_distance'] = (p.distance_to_span > 10).astype(int)
    p['pt_keeps_base'] = (~p.drops_its_base.astype(bool)).astype(int)
    key = p.groupby(['panel', 'footprint', 'pos']).snv.transform('count') - 1
    p['agreeing_alleles'] = key
    p['pt_alleles'] = key.clip(upper=2)
    pts = [c for c in p.columns if c.startswith('pt_')]
    p['score'] = p[pts].sum(axis=1)

    for sub, tag in (('HepG2_run3', 'A7'), ('HepG2_run3_pm5', 'pm5')):
        cf = caller_free(root, p, sub)
        p[f'cf_inacc_{tag}_pts'] = [x[0] for x in cf]
        p[f'cf_adenines_{tag}'] = [x[1] for x in cf]

    p = p.sort_values(['panel', 'score', 'p_A7_confident'], ascending=[True, False, True])
    p['rank'] = p.groupby('panel').cumcount() + 1
    cols = ['panel', 'rank', 'score', 'snv', 'pos', 'footprint', 'span_lo', 'span_hi', 'distance_to_span',
            'drops_its_base', 'D_pts', 'p_A7_confident', 'D_pm5_pts', 'p_pm5_confident', 'sig_pm5_confident',
            'act', 'act_p', 'act2', 'agreeing_alleles', 'cf_inacc_A7_pts', 'cf_adenines_A7', 'cf_inacc_pm5_pts',
            'cf_adenines_pm5', 'C1', 'hek_molecules', 'hek_clones'] + pts
    p[cols].to_csv(out / 'ranked.tsv', sep='\t', index=False, float_format='%.4g')
    json.dump(dict(rules=rules, attrition=dict(sig_pairs=n0, tf_sized=n1, activity_gate=n2),
                   by_panel=p.panel.value_counts().to_dict()), open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

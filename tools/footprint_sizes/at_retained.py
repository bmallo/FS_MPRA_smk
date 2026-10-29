"""Post hoc: footprint occupancy tests with A/T-changing SNVs removed where they can touch a call.

Variant molecules are scored with WT emissions, and the substituted base no longer matches the
reference, so the caller drops it. An SNV that removes (or adds) an A/T where a counted call can
lie changes the evidence the caller sees, independently of any protein (inside footprint spans:
residual -0.26 pp for A/T loss, -0.08 retained, +0.05 gain, with near-naked molecules moving the
same way). This rerun keeps the estimates in hepg2_effects unchanged (same deltas, SEs, slopes,
calibration gate) and only changes which SNV x feature pairs form the primary family:

  call reach of a feature = [mid_lo - len_hi/2, mid_hi + len_hi/2] (element index; every call
      counted in the box lies inside it)
  primary   = A/T-retaining SNVs anywhere + A/T-changing SNVs outside the call reach
  secondary = A/T-changing SNVs inside the call reach (reported, BH separately)

Both families are restricted to the 12 features whose occupancy-residual gate passed. The pooled
within-span test is repeated on A/T-retaining SNVs only. Post hoc: chosen after the A/T pattern was
seen in the per-SNV effects (2026-10-28).

Usage: python tools/footprint_sizes/at_retained.py --dir results/footprint_sizes_20261028 \
           --classes results/variant_effects_20261028/classes
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from fsmpra.site_effects import bh


def pooled(e):
    d, se = e.occupancy_resid_delta.to_numpy(), e.occupancy_resid_se.to_numpy()
    pse = np.sqrt((se ** 2).sum()) / len(d)
    z = d.mean() / pse
    return dict(k=len(d), mean_pp=100 * d.mean(), z=z, p=2 * stats.norm.sf(abs(z)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--classes', required=True)
    ap.add_argument('--alpha', type=float, default=0.05)
    a = ap.parse_args()
    R = Path(a.dir)
    out = R / 'at_retained'
    out.mkdir(exist_ok=False)
    feats = pd.read_csv(R / 'hepg2_map/features.tsv', sep='\t').set_index('feature')
    cal = pd.read_csv(R / 'hepg2_effects/calibration_summary.tsv', sep='\t')
    ok = set(cal[(cal.endpoint == 'occupancy_resid') & (cal.gate_pass == True)].feature) - {'all'}
    cls = pd.read_csv(Path(a.classes) / 'variants.tsv', sep='\t').set_index('source_variant')
    e = pd.read_csv(R / 'hepg2_effects/effects.tsv.gz', sep='\t')
    e = e[e.feature.isin(ok) & e.occupancy_resid_p.notna()].copy()
    f = feats.loc[e.feature]
    p0 = e.oligo_pos1.to_numpy() - 21
    reach_lo = f.mid_lo.to_numpy() - f.len_hi.to_numpy() / 2
    reach_hi = f.mid_hi.to_numpy() + f.len_hi.to_numpy() / 2
    e['in_reach'] = (p0 >= reach_lo) & (p0 <= reach_hi)
    e['at_changed'] = e.at_change != 'retained_count'
    e['family'] = np.where(e.at_changed & e.in_reach, 'secondary_at_changed_in_reach', 'primary')
    e['q_all_gated'] = bh(e.occupancy_resid_p.to_numpy())  # the original 7,560-test family, for comparison
    e['q_family'] = np.nan
    for fam, idx in e.groupby('family').groups.items():
        e.loc[idx, 'q_family'] = bh(e.loc[idx, 'occupancy_resid_p'].to_numpy())
    e = e.join(cls[['activity', 'functional_hit']], on='source_variant')
    cols = ['source_variant', 'oligo_pos1', 'at_change', 'feature', 'snv_in_footprint', 'in_reach', 'family', 'molecules',
            'occupancy_C_delta', 'occupancy_M_delta', 'slope_b', 'occupancy_resid_delta', 'occupancy_resid_z', 'q_all_gated', 'q_family',
            'activity', 'functional_hit']
    e[cols].to_csv(out / 'tests.tsv.gz', sep='\t', index=False)

    old = e[e.q_all_gated < a.alpha]
    new = e[(e.family == 'primary') & (e.q_family < a.alpha)]
    sec = e[(e.family != 'primary') & (e.q_family < a.alpha)]
    key = lambda d: set(zip(d.source_variant, d.feature))
    status = []
    for r in e[(e.q_all_gated < a.alpha) | (e.q_family < a.alpha)].itertuples():
        k = (r.source_variant, r.feature)
        status.append(dict(source_variant=r.source_variant, oligo_pos1=r.oligo_pos1, feature=r.feature, at_change=r.at_change,
                           in_reach=r.in_reach, family=r.family, delta_pp=100 * r.occupancy_resid_delta, z=r.occupancy_resid_z,
                           q_original=r.q_all_gated, q_rerun=r.q_family, original_hit=k in key(old),
                           rerun=('primary hit' if k in key(new) else 'secondary hit' if k in key(sec) else
                                  'not significant (primary)' if r.family == 'primary' else 'not significant (secondary)'),
                           activity=r.activity, functional_hit=r.functional_hit))
    st = pd.DataFrame(status).sort_values(['feature', 'oligo_pos1', 'source_variant'])
    st.to_csv(out / 'hit_comparison.tsv', sep='\t', index=False)

    prow = []
    for ft in sorted(ok):
        x = e[(e.feature == ft) & e.snv_in_footprint]
        for lab, sub in (('all inside', x), ('A/T-retaining inside', x[~x.at_changed]), ('A/T-changing inside', x[x.at_changed])):
            if len(sub):
                prow.append(dict(feature=ft, centre_oligo=int(feats.loc[ft, 'mid_peak_oligo']), length=int(feats.loc[ft, 'len_peak']), set=lab, **pooled(sub)))
    pt = pd.DataFrame(prow)
    for lab in pt.set.unique():
        k = pt.set == lab
        pt.loc[k, 'q'] = bh(pt.loc[k, 'p'].to_numpy())
    pt.to_csv(out / 'pooled_inside.tsv', sep='\t', index=False)

    summ = dict(tests_original=len(e), hits_original=len(old), tests_primary=int((e.family == 'primary').sum()),
                hits_primary=len(new), tests_secondary=int((e.family != 'primary').sum()), hits_secondary=len(sec),
                original_hits_kept_primary=len(key(old) & key(new)), original_hits_moved_secondary=int((st.original_hit & (st.family != 'primary')).sum()),
                original_hits_lost_primary=int((st.original_hit & (st.family == 'primary') & (st.rerun != 'primary hit')).sum()),
                new_primary_hits=len(key(new) - key(old)), primary_hits_by_feature=new.feature.value_counts().to_dict(),
                rule=__doc__.split('\n\n')[1].strip().split('\n')[0], post_hoc=True)
    json.dump(summ, open(out / 'summary.json', 'w'), indent=1)
    pd.set_option('display.width', 250)
    print(json.dumps(summ, indent=1))
    print(st.round(4).to_string(index=False))
    print(pt.round(4).to_string(index=False))


if __name__ == '__main__':
    main()

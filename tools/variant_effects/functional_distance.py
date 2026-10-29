"""How many functional SNVs change chromatin anywhere, and where the change sits relative to the SNV.

One family per SNV x test, BH across all of them, over every HepG2 chromatin test whose WT
pseudo-variant gate passes:
  site m6A (A_norm): local (5-10 bp from the SNV) and flank (11-30 bp) windows, and the fixed
      regions (WT_P1-4, P1_core, Repeat1-3, adjacent_143_150), each excluding +-4 bp of the SNV
  LLR footprint coverage (observability-corrected 'net'): TF- and nucleosome-sized, fp_local / fp_flank
  element scalars: MSP accessibility and 85-250 bp protection (C_std)
  size-resolved footprints: occupancy residual of the 12 gated features
Distance = bp from the SNV to the tested window/region/footprint span (0 = the SNV lies in it).
Bins: overlap (0, or the local windows), near (1-30 bp, or the flank windows), distal (> 30 bp),
element-wide (MSP / 85-250 bp scalars, no position).

Usage: python tools/variant_effects/functional_distance.py --dir results/variant_effects_20261028 \
           --sizes results/footprint_sizes_20261028 --config config/ldlr_footprints.json --output <dir>
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fsmpra.site_effects import bh

NEAR = 30


def gated(path):
    c = pd.read_csv(path, sep='\t')
    return set(c.loc[c.gate_pass.astype(str) == 'True', 'endpoint'])


def dist(p0, lo, hi):
    """bp from p0 to the closed interval [lo, hi]."""
    return np.where(p0 < lo, lo - p0, np.where(p0 > hi, p0 - hi, 0))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ('--dir', '--sizes', '--config', '--output'):
        ap.add_argument(k, required=True)
    ap.add_argument('--alpha', type=float, default=0.05)
    a = ap.parse_args()
    V, S, out = Path(a.dir), Path(a.sizes), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    regions = json.load(open(a.config))['analysis']['regions0']
    cls = pd.read_csv(V / 'classes/variants.tsv', sep='\t')
    p0 = cls.set_index('source_variant').oligo_pos1 - 21

    parts = []
    site = pd.read_csv(V / 'hepg2/effects/effects.tsv.gz', sep='\t')
    g = gated(V / 'hepg2/effects/calibration_summary.tsv')
    s = site[(site.measure == 'A_norm') & site.window.isin(['local', 'flank', *regions]) & site.endpoint.isin(g)].copy()
    s['family'] = 'site_m6A'
    s['distance'] = np.where(s.window == 'local', 0, np.where(s.window == 'flank', 11, -1))
    for name, (lo, hi) in regions.items():
        k = s.window == name
        s.loc[k, 'distance'] = dist(s.loc[k, 'oligo_pos1'].to_numpy() - 21, lo, hi - 1)
    s['bin_hint'] = np.where(s.window == 'local', 'overlap', np.where(s.window == 'flank', 'near', ''))
    parts.append(s)
    e = site[site.endpoint.isin(['msp_C_std', 'nuc_85_250_C_std']) & site.endpoint.isin(g)].copy()
    e['family'], e['distance'], e['bin_hint'] = 'element_scalar', np.nan, 'element-wide'
    parts.append(e)
    fp = pd.read_csv(V / 'hepg2/footprints/effects.tsv.gz', sep='\t')
    gf = gated(V / 'hepg2/footprints/calibration_summary.tsv')
    f = fp[(fp.measure == 'net') & fp.endpoint.isin(gf)].copy()
    f['family'] = 'llr_window'
    f['distance'] = np.where(f.window == 'fp_local', 0, 11)
    f['bin_hint'] = np.where(f.window == 'fp_local', 'overlap', 'near')
    parts.append(f)
    feats = pd.read_csv(S / 'hepg2_map/features.tsv', sep='\t').set_index('feature')
    fe = pd.read_csv(S / 'hepg2_effects/effects.tsv.gz', sep='\t')
    gs = pd.read_csv(S / 'hepg2_effects/calibration_summary.tsv', sep='\t')
    ok = set(gs[(gs.endpoint == 'occupancy_resid') & (gs.gate_pass == True)].feature) - {'all'}
    fe = fe[fe.feature.isin(ok) & fe.occupancy_resid_p.notna()].copy()
    span_lo = feats.mid_peak - feats.len_peak / 2
    span_hi = feats.mid_peak + feats.len_peak / 2
    fe['distance'] = np.ceil(dist(fe.oligo_pos1.to_numpy() - 21, span_lo[fe.feature].to_numpy(), span_hi[fe.feature].to_numpy()))
    fe = fe.rename(columns={'occupancy_resid_delta': 'delta', 'occupancy_resid_z': 'z', 'occupancy_resid_p': 'p'})
    fe['endpoint'], fe['family'], fe['bin_hint'] = 'size:' + fe.feature, 'size_feature', ''
    parts.append(fe)

    cols = ['source_variant', 'oligo_pos1', 'family', 'endpoint', 'distance', 'bin_hint', 'delta', 'z', 'p']
    t = pd.concat([x[cols] for x in parts], ignore_index=True)
    t = t[t.p.notna() & t.source_variant.isin(cls.source_variant)]
    t['q_all'] = bh(t.p.to_numpy())
    t['bin'] = np.where(t.bin_hint != '', t.bin_hint, np.where(t.distance == 0, 'overlap', np.where(t.distance <= NEAR, 'near', 'distal')))
    t = t.merge(cls[['source_variant', 'activity', 'functional_hit']], on='source_variant')
    t.to_csv(out / 'tests.tsv.gz', sep='\t', index=False)

    hit = t[t.q_all < a.alpha]
    bins = ['overlap', 'near', 'distal', 'element-wide']
    per = pd.DataFrame(index=cls.source_variant)
    for b in bins:
        per[b] = per.index.isin(hit[hit.bin == b].source_variant)
    per['any'] = per[bins].any(axis=1)
    per['nearest'] = np.select([per[b] for b in bins], bins, 'none')
    per['tests'] = t.groupby('source_variant').size().reindex(per.index).fillna(0).astype(int)
    per['hits'] = hit.groupby('source_variant').size().reindex(per.index).fillna(0).astype(int)
    per = per.join(cls.set_index('source_variant')[['oligo_pos1', 'activity', 'functional_hit', 'hepg2_molecules']])
    per.to_csv(out / 'per_snv.tsv', sep='\t')

    rows = []
    for lab, sel in (('functional', per.functional_hit == True), ('non-functional', per.functional_hit == False),
                     ('functional, activity down', (per.functional_hit == True) & (per.activity < 0)),
                     ('functional, activity up', (per.functional_hit == True) & (per.activity > 0)),
                     ('functional, |log2| >= 1', (per.functional_hit == True) & (per.activity.abs() >= 1))):
        p = per[sel]
        r = dict(group=lab, snvs=len(p), any=int(p['any'].sum()), any_pct=100 * p['any'].mean())
        for b in bins:
            r[f'{b}'] = int(p[b].sum())
        for b in bins + ['none']:
            r[f'nearest_{b}'] = int((p.nearest == b).sum())
        rows.append(r)
    summ = pd.DataFrame(rows)
    summ.to_csv(out / 'summary.tsv', sep='\t', index=False)
    fam = t.groupby(['family', 'bin']).agg(tests=('p', 'size'), hits=('q_all', lambda q: int((q < a.alpha).sum())))
    fam.to_csv(out / 'tests_by_family.tsv', sep='\t')
    hb = hit.merge(per[['nearest']], left_on='source_variant', right_index=True)
    dd = hit[hit.distance.notna()].assign(func=lambda d: d.functional_hit == True)
    dist_tab = dd.groupby([pd.cut(dd.distance, [-1, 0, 10, 30, 60, 120, 250]), 'func'], observed=False).size().unstack(fill_value=0)
    dist_tab.to_csv(out / 'hit_distance_histogram.tsv', sep='\t')
    json.dump(dict(tests=len(t), hits=len(hit), alpha=a.alpha, near_bp=NEAR, families=t.family.value_counts().to_dict()),
              open(out / 'summary.json', 'w'), indent=1)
    pd.set_option('display.width', 250)
    print(summ.round(1).to_string(index=False))
    print(fam.to_string())
    print(dist_tab.to_string())


if __name__ == '__main__':
    main()

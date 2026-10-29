"""Join the size-resolved footprint features with replication, calibration, variant effects, activity and motifs.

Reads results/footprint_sizes_20261028/{hepg2_map, hepg2_effects, hepg2_map_subset_dd25, hepg2_map_subset_dd10,
hek293t_map} and results/variant_effects_20261028/classes. Per-SNV q and the pooled within-footprint q are
recomputed over the features whose occupancy-residual calibration gate passes (the gate rule: z/p/q only
where it passes). Writes <dir>/summary/{feature_summary.tsv, snv_hits.tsv, summary.json}.

Usage: python tools/footprint_sizes/summarize.py --dir results/footprint_sizes_20261028 \
           --classes results/variant_effects_20261028/classes
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from fsmpra.footprint_sizes import box_mass
from fsmpra.site_effects import bh

LO = 10


def match(f, other, dmid=3, dlen=4):
    """Nearest peak in another feature table within +-dmid bp midpoint and +-dlen bp length."""
    c = other[(abs(other.mid_peak_oligo - f.mid_peak_oligo) <= dmid) & (abs(other.len_peak - f.len_peak) <= dlen)]
    if c.empty:
        return ''
    c = c.assign(d=abs(c.mid_peak_oligo - f.mid_peak_oligo) + abs(c.len_peak - f.len_peak)).sort_values('d')
    return f'{int(c.mid_peak_oligo.iloc[0])}/{int(c.len_peak.iloc[0])}'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--classes', required=True)
    a = ap.parse_args()
    R, K = Path(a.dir), Path(a.classes)
    out = R / 'summary'
    out.mkdir(exist_ok=True)
    feats = pd.read_csv(R / 'hepg2_map/features.tsv', sep='\t')
    maps = {k: np.load(R / k / 'map.npz') for k in ('hepg2_map', 'hepg2_map_subset_dd25', 'hepg2_map_subset_dd10', 'hek293t_map')}
    other = {k: pd.read_csv(R / k / 'features.tsv', sep='\t') for k in ('hepg2_map_subset_dd10', 'hek293t_map')}
    cal = pd.read_csv(R / 'hepg2_effects/calibration_summary.tsv', sep='\t')
    gate = cal[cal.endpoint == 'occupancy_resid'].set_index('feature')
    eff = pd.read_csv(R / 'hepg2_effects/effects.tsv.gz', sep='\t')
    cls = pd.read_csv(K / 'variants.tsv', sep='\t').set_index('source_variant')
    motifs = pd.read_csv(K / 'motif_sites_ranked.tsv', sep='\t')
    motifs[['m_lo', 'm_hi']] = motifs.site_oligo.str.split('-', expand=True).astype(int)

    ok = set(gate[gate.gate_pass == True].index) - {'all'}
    u = eff[eff.feature.isin(ok) & eff.occupancy_resid_p.notna()].copy()
    u['q_gated'] = bh(u.occupancy_resid_p.values)
    u = u.join(cls[['activity', 'functional_hit', 'hepg2_msp_delta', 'hepg2_msp_q', 'hepg2_nuc_delta', 'class']], on='source_variant')

    rows = []
    for _, f in feats.iterrows():
        d = f.to_dict()
        span = (f.mid_peak_oligo - f.len_peak / 2, f.mid_peak_oligo + f.len_peak / 2)
        e = eff[eff.feature == f.feature].join(cls[['activity']], on='source_variant')
        inside = e[e.snv_in_footprint & e.occupancy_resid_delta.notna()]
        r = dict(feature=f.feature, centre_oligo=f.mid_peak_oligo, window_oligo=f'{f.mid_lo_oligo}-{f.mid_hi_oligo}',
                 length_peak=f.len_peak, length_range=f'{f.len_lo}-{f.len_hi}', span_oligo=f'{span[0]:.0f}-{span[1]:.0f}',
                 flag=f.flag, c_over_m=f.c_over_m,
                 occ_discovery=box_mass(maps['hepg2_map']['discovery_C'], d, LO), occ_reference=box_mass(maps['hepg2_map']['reference_C'], d, LO),
                 occ_subset_dd25=box_mass(maps['hepg2_map_subset_dd25']['all_C'], d, LO),
                 occ_subset_dd10=box_mass(maps['hepg2_map_subset_dd10']['all_C'], d, LO), dd10_peak=match(f, other['hepg2_map_subset_dd10']),
                 occ_hek293t=box_mass(maps['hek293t_map']['all_C'], d, LO), hek293t_peak=match(f, other['hek293t_map']),
                 gate_pass=f.feature in ok, gate_sd_z=gate.loc[f.feature, 'sd_z'], gate_rate_high=gate.loc[f.feature, 'rate_ci_high'],
                 slope_b=e.slope_b.iloc[0], snvs_inside=len(inside))
        if len(inside):
            dd, se = inside.occupancy_resid_delta.to_numpy(), inside.occupancy_resid_se.to_numpy()
            pse = np.sqrt((se ** 2).sum()) / len(dd)
            r.update(pooled_inside_delta=dd.mean(), pooled_inside_se=pse, pooled_inside_z=dd.mean() / pse,
                     pooled_inside_p=2 * stats.norm.sf(abs(dd.mean() / pse)), inside_fraction_negative=(dd < 0).mean())
            fin = inside.activity.notna()
            if fin.sum() >= 10:
                rho, p = stats.spearmanr(inside.occupancy_resid_delta[fin], inside.activity[fin])
                r.update(rho_occupancy_activity_inside=rho, p_rho_inside=p)
        h = u[u.feature == f.feature]
        r.update(snv_hits=int((h.q_gated < .05).sum()), snv_hits_loss=int(((h.q_gated < .05) & (h.occupancy_resid_delta < 0)).sum()),
                 snv_hits_inside=int(((h.q_gated < .05) & h.snv_in_footprint).sum()), powered_snvs=int((e.expected_events >= 10).sum()))
        ov = motifs[(motifs.m_hi >= span[0]) & (motifs.m_lo <= span[1])].sort_values('rho_activity', ascending=False)
        r['top_motifs'] = '; '.join(f'{m.name} {m.site_oligo} (rho {m.rho_activity:.2f})' for m in ov.drop_duplicates('name').head(3).itertuples())
        rows.append(r)
    fs = pd.DataFrame(rows)
    g = fs.gate_pass & fs.pooled_inside_p.notna()
    fs.loc[g, 'pooled_inside_q_gated'] = bh(fs.loc[g, 'pooled_inside_p'].values)
    fs.to_csv(out / 'feature_summary.tsv', sep='\t', index=False)
    hits = u[u.q_gated < .05].sort_values(['feature', 'oligo_pos1'])
    cols = ['source_variant', 'oligo_pos1', 'feature', 'snv_in_footprint', 'distance_to_feature', 'at_change', 'molecules', 'expected_events',
            'occupancy_C_delta', 'occupancy_M_delta', 'slope_b', 'occupancy_resid_delta', 'occupancy_resid_z', 'q_gated',
            'activity', 'functional_hit', 'hepg2_msp_delta', 'hepg2_msp_q', 'hepg2_nuc_delta', 'class']
    hits[cols].to_csv(out / 'snv_hits.tsv', sep='\t', index=False)
    hek_r = np.corrcoef(fs.occ_discovery, fs.occ_hek293t)[0, 1]
    summ = dict(features=len(fs), gate_pass=int(fs.gate_pass.sum()), tests_gated=len(u), snv_hits=len(hits),
                snv_hits_by_feature=hits.feature.value_counts().to_dict(),
                pooled_inside_q_lt_05=fs.loc[fs.pooled_inside_q_gated < .05, 'feature'].tolist(),
                halves_r=float(np.corrcoef(fs.occ_discovery, fs.occ_reference)[0, 1]), hek293t_r=float(hek_r),
                hek293t_spearman=float(stats.spearmanr(fs.occ_discovery, fs.occ_hek293t)[0]),
                hek293t_matched=int((fs.hek293t_peak != '').sum()), dd10_matched=int((fs.dd10_peak != '').sum()),
                hits_functional=int(hits.functional_hit.fillna(False).astype(bool).sum()))
    json.dump(summ, open(out / 'summary.json', 'w'), indent=1)
    pd.set_option('display.width', 260)
    print(fs[['feature', 'centre_oligo', 'length_peak', 'span_oligo', 'flag', 'occ_discovery', 'occ_hek293t', 'hek293t_peak', 'dd10_peak', 'gate_pass',
              'pooled_inside_delta', 'pooled_inside_z', 'pooled_inside_q_gated', 'rho_occupancy_activity_inside', 'snv_hits', 'top_motifs']].round(4).to_string(index=False))
    print(json.dumps(summ, indent=1))


if __name__ == '__main__':
    main()

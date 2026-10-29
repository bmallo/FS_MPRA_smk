"""Collect every per-position and per-SNV layer into one JSON for the element track browser.

Inputs (all existing outputs): results/variant_effects_20261028/{classes, hepg2/effects, hepg2/footprints,
hek293t/footprints, functional_distance} and results/footprint_sizes_20261028/{hepg2_map, hek293t_map,
hepg2_effects, summary, shuffle_null (optional)}. Size maps are smoothed as on the deck slide (sigma 1.5 bp
midpoint, 2 bp length) and stored as 0-255 on a square-root scale with their maxima.

Usage: python tools/footprint_sizes/browser_data.py --output results/footprint_sizes_20261028/browser/data.json
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage

from fsmpra.site_effects import bh

V = Path('results/variant_effects_20261028')
S = Path('results/footprint_sizes_20261028')
REGIONS = [('P1 core', 54, 62), ('Repeat 1 (Sp1)', 89, 104), ('Repeat 2 (SRE-1)', 124, 139), ('Repeat 3 (Sp1)', 140, 155)]
MOTIFS = [(52, 62, 'ZNF652-like'), (82, 102, 'SP4'), (127, 136, 'SREBF1/2'), (140, 156, 'SP1/SP3/KLF'), (182, 192, 'CREB1/CREM')]


def r(x, d):
    return None if x is None or not np.isfinite(x) else round(float(x), d)


def qmap(M):
    M = ndimage.gaussian_filter(np.clip(M, 0, None), sigma=(1.5, 2.0))
    m = float(M.max())
    return dict(max=m, q=np.rint(255 * np.sqrt(M / m)).astype(int).T.ravel().tolist())  # length-major rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    variants = json.load(open(V / 'hepg2/extract_llr/variants.json'))
    seq = next(v for v in variants if v['class'] == 'WT')['sequence']
    wt = pd.read_csv(V / 'hepg2/footprints/wt_profile.tsv', sep='\t')
    hw = pd.read_csv(V / 'hek293t/footprints/wt_profile.tsv', sep='\t').set_index('site_oligo_pos1').reindex(wt.site_oligo_pos1)
    hs = pd.read_csv(V / 'hek293t/effects/effects.tsv.gz', sep='\t')
    hgate = pd.read_csv(V / 'hek293t/effects/calibration_summary.tsv', sep='\t').set_index('endpoint').gate_pass.astype(str).eq('True')
    cls = pd.read_csv(V / 'classes/variants.tsv', sep='\t')
    site = pd.read_csv(V / 'hepg2/effects/effects.tsv.gz', sep='\t')
    loc = site[site.endpoint == 'local:A_norm'].set_index('source_variant')
    fl = site[site.endpoint == 'flank:A_norm'].set_index('source_variant')
    fd = pd.read_csv(V / 'functional_distance/per_snv.tsv', sep='\t').set_index('source_variant')

    feats = pd.read_csv(S / 'summary/feature_summary.tsv', sep='\t')
    fmap = pd.read_csv(S / 'hepg2_map/features.tsv', sep='\t').set_index('feature')
    eff = pd.read_csv(S / 'hepg2_effects/effects.tsv.gz', sep='\t')
    ok = set(feats.loc[feats.gate_pass, 'feature'])
    g = eff[eff.feature.isin(ok) & eff.occupancy_resid_p.notna()].copy()
    at = pd.read_csv(S / 'at_retained/tests.tsv.gz', sep='\t')[['source_variant', 'feature', 'family', 'q_family']]
    at = at.rename(columns={'q_family': 'q'}).assign(sec=lambda d: d.family != 'primary')
    eff = eff.merge(at[['source_variant', 'feature', 'q', 'sec']], on=['source_variant', 'feature'], how='left')
    shuf = S / 'shuffle_local/shuffle_local.tsv'
    shuf = pd.read_csv(shuf, sep='\t').set_index('feature') if shuf.exists() else None

    # Shared per-SNV fields are sequence and the published MPRA activity only; chromatin values live per cell.
    snvs, hep_vals = [], {}
    for v in cls.itertuples():
        s = v.source_variant
        snvs.append(dict(id=s, pos=int(v.oligo_pos1), ref=v.ref, alt=v.alt, at=v.at_change if isinstance(v.at_change, str) else '',
                         act=r(v.activity, 3), func=bool(v.functional_hit)))
        hep_vals[s] = dict(msp=r(100 * v.hepg2_msp_delta, 2), msp_q=r(v.hepg2_msp_q, 4), nuc=r(100 * v.hepg2_nuc_delta, 2), nuc_q=r(v.hepg2_nuc_q, 4),
                           loc_z=r(loc.z.get(s, np.nan), 2), loc_q=r(loc.q_bh_within_endpoint.get(s, np.nan), 4),
                           near=fd.nearest.get(s, 'none'), mol=int(v.hepg2_molecules) if np.isfinite(v.hepg2_molecules) else None)
    hek_vals = {}
    for s, g_ in hs.groupby('source_variant'):
        e_ = g_.set_index('endpoint')
        get = lambda ep, c: e_[c].get(ep, np.nan)
        hek_vals[s] = dict(msp=r(100 * get('msp_C_std', 'delta'), 2), msp_q=r(get('msp_C_std', 'q_bh_primary_family'), 4),
                           nuc=r(100 * get('nuc_85_250_C_std', 'delta'), 2), nuc_q=r(get('nuc_85_250_C_std', 'q_bh_primary_family'), 4),
                           loc_z=r(get('local:A_norm', 'z'), 2), loc_q=r(get('local:A_norm', 'q_bh_within_endpoint'), 4),
                           mol=int(g_.molecules.iloc[0]), clones=int(g_.clones.iloc[0]))
    idx = {d['id']: i for i, d in enumerate(snvs)}
    fout = []
    for f in feats.itertuples():
        m = fmap.loc[f.feature]
        e = eff[eff.feature == f.feature]
        z = [None] * len(snvs)
        dl = [None] * len(snvs)
        q = [None] * len(snvs)
        sec = [0] * len(snvs)
        for row in e.itertuples():
            i = idx.get(row.source_variant)
            if i is None:
                continue
            z[i], dl[i], q[i] = r(row.occupancy_resid_z, 2), r(100 * row.occupancy_resid_delta, 3), r(row.q, 4)
            sec[i] = int(row.sec == True)
        d = dict(id=f.feature, mid=int(f.centre_oligo), len=int(f.length_peak), mid_lo=int(m.mid_lo_oligo), mid_hi=int(m.mid_hi_oligo),
                 len_lo=int(m.len_lo), len_hi=int(m.len_hi), flag=f.flag, c_over_m=r(f.c_over_m, 2), occ=r(100 * f.occ_discovery, 3),
                 gate=bool(f.gate_pass), pooled_z=r(f.pooled_inside_z, 2), pooled_q=r(f.pooled_inside_q_gated, 4),
                 hits=int(sum(1 for i in range(len(snvs)) if q[i] is not None and q[i] < .05 and not sec[i])), sec=sec, motifs=f.top_motifs if isinstance(f.top_motifs, str) else '', z=z, d=dl, q=q)
        if shuf is not None and f.feature in shuf.index:
            sr = shuf.loc[f.feature]
            d.update(shuf_c=r(sr.C_fold, 2), shuf_m=r(sr.M_fold, 2), shuf_obs_m=r(100 * sr.observed_M_free, 3), shuf_null_m=r(100 * sr.null_M_free, 3),
                     shuf_obs_c=r(100 * sr.observed_C_std_free, 3), shuf_null_c=r(100 * sr.null_C_std_free, 3))
        fout.append(d)

    hep = np.load(S / 'hepg2_map/map.npz')
    hek = np.load(S / 'hek293t_map/map.npz')
    hf = pd.read_csv(S / 'hek293t_map/features.tsv', sep='\t')
    hek_feats = [dict(id='K' + f.feature[1:], mid=int(f.mid_peak_oligo), len=int(f.len_peak), mid_lo=int(f.mid_lo_oligo), mid_hi=int(f.mid_hi_oligo),
                      len_lo=int(f.len_lo), len_hi=int(f.len_hi), flag=f.flag, c_over_m=r(f.c_over_m, 2), occ=r(100 * f.occupancy_C_std, 3),
                      occ_m=r(100 * f.occupancy_M, 3), gate=False) for f in hf.itertuples()]
    wtd = lambda p: dict(msp=[r(100 * x, 2) for x in p.wt_msp_C_std], msp_m=[r(100 * x, 2) for x in p.wt_msp_M],
                         nuc=[r(100 * x, 2) for x in p.wt_nuc_C_std], nuc_m=[r(100 * x, 2) for x in p.wt_nuc_M],
                         tf=[r(100 * x, 3) for x in p.wt_tf_C_std], tf_m=[r(100 * x, 3) for x in p.wt_tf_M])
    cells = {
        'HepG2': dict(wt=wtd(wt), vals=hep_vals, features=fout, effects=True, gates=dict(msp=True, nuc=True, loc_z=True),
                      maps=dict(C=qmap(hep['discovery_C']), M=qmap(hep['discovery_M'])),
                      source='replicate 3 · WT discovery half for the map'),
        'HEK293T': dict(wt=wtd(hw), vals=hek_vals, features=hek_feats, effects=False,
                        gates=dict(msp=bool(hgate.get('msp_C_std', False)), nuc=bool(hgate.get('nuc_85_250_C_std', False)),
                                   loc_z=bool(hgate.get('local:A_norm', False))),
                        maps=dict(C=qmap(hek['all_C']), M=qmap(hek['all_M'])),
                        source='all WT molecules for the map · per-SNV values are single clones vs the WT clone null'),
    }
    out = dict(oligo0=21, seq=seq, regions=REGIONS, motifs=MOTIFS, snvs=snvs, cells=cells, lengths=[int(x) for x in hep['lengths']])
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(a.output, 'w'), separators=(',', ':'))
    print(Path(a.output).stat().st_size, 'bytes;', len(snvs), 'SNVs;', len(fout), 'features; shuffle', shuf is not None)


if __name__ == '__main__':
    main()

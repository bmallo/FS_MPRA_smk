"""Per-variant footprint effects from plasmid-LLR coverage: per-base profiles and calibrated local windows.

Reads a site extraction whose coverage bits come from plasmid-LLR (fsmpra.llr_coverage):
bit 1 TF-sized call (10-60 bp), 2 MSP, 4 nucleosome-sized call (85-250 bp). Each bit is treated
like a per-site "call" at every element position, so site_effects.measures standardizes it the
same way as m6A: per class U, M, C_raw, C_std (over 10 WT labeling strata).

Outputs per SNV:
- per-base profiles (point estimates): delta vs WT of C_std occupancy for TF, nucleosome and MSP
  at all 211 positions, and of M (the near-naked proxy);
- calibrated windows around the substitution, excluding +-mask_radius bp:
  fp_local (5-20 bp) and fp_flank (21-40 bp), for TF and nucleosome tracks; endpoints C_std, M and
  net = C_std - M. The substitution changes which adenines can be observed nearby (the
  observability kernel), and LLR calls inherit that; M molecules carry the same sequence but
  little protein, so subtracting their change removes the part that is measurement.

Uncertainty and calibration follow site_effects: barcode-cluster Poisson bootstrap per SNV, WT
pseudo-variants drawn as whole clusters at random positions, gate z SD 0.85-1.20 and upper exact
bound of the |z| > 1.96 rate <= 0.10. BH over the prespecified footprint family.

Usage: python -m fsmpra.footprint_effects --config config/ldlr_footprints.json --extract <dir> --output <dir>
"""
import argparse
from pathlib import Path

import numpy as np
from scipy import stats

from .io import atomic_output, read_json, source_info, write_json, write_table
from .site_effects import (COHORTS, GENOMIC_OFFSET, assign_groups, bh, bootstrap, contrast, exact_ci, grouped_sum,
                           homopolymer_context, load_activity, measures)

TRACKS = {'tf': 1, 'msp': 2, 'nuc': 4}


def footprint_masks(p0, windows, mask_radius):
    """All element positions (footprints are defined at every base) in each distance window from p0."""
    d = np.abs(np.arange(211) - p0)
    names = list(windows)
    rows = [(d >= lo) & (d <= hi) & (d > mask_radius) for lo, hi in windows.values()]
    return names, np.array(rows, dtype=float)


def endpoints(site):
    """C_std, M and net = C_std - M from a measures() site dict."""
    return {'C_std': site['C_std'], 'M': site['M'], 'net': site['C_std'] - site['M']}


def run(config, extract, output):
    cfg = read_json(config)
    spec, fp = cfg['analysis'], cfg['footprints']
    extract = Path(extract)
    audit = read_json(extract / 'audit.json')
    if 'llr_coverage' not in audit:
        raise ValueError('Expected an extraction with plasmid-LLR coverage (fsmpra.llr_coverage)')
    variants = read_json(extract / 'variants.json')
    wt_i = next(i for i, v in enumerate(variants) if v['class'] == 'WT')
    seq = variants[wt_i]['sequence']
    meta = dict(np.load(extract / 'meta.npz'))
    cohort = meta[COHORTS[spec['cohort']]].astype(bool) & meta.get('evaluated', np.ones(len(meta['variant_index']), bool)).astype(bool)
    vi_all = meta['variant_index'].astype(int)
    groups_all, weights, edges = assign_groups(meta['distal_m6a'].astype(float), (vi_all == wt_i) & cohort,
                                               spec['unmethylated_max_distal_m6a'], spec['methylated_min_distal_m6a'],
                                               spec['labeling_strata'])
    n_groups = 2 + len(weights)
    keep = np.flatnonzero(cohort & (groups_all >= 0))
    cov_all = np.load(extract / 'coverage.npy', mmap_mode='r')
    bits = {t: np.zeros((len(keep), 211), dtype=bool) for t in TRACKS}
    for a in range(0, len(keep), 200000):
        cov = cov_all[keep[a:a + 200000]]
        for t, bit in TRACKS.items():
            bits[t][a:a + len(cov)] = (cov & bit) > 0
    ones = np.ones((len(keep), 211), dtype=bool)
    noscal = np.zeros((len(keep), 3), dtype=np.float32)
    vi, groups, clusters = vi_all[keep], groups_all[keep].astype(int), meta['cluster'][keep]
    keys = vi * n_groups + groups
    nv = len(variants)
    cnt = np.bincount(keys, minlength=nv * n_groups).reshape(nv, n_groups).astype(float)
    vs = grouped_sum(keys, ones, nv * n_groups).reshape(nv, n_groups, 211)
    cs = {t: grouped_sum(keys, bits[t], nv * n_groups).reshape(nv, n_groups, 211) for t in TRACKS}
    ss = np.zeros((nv, n_groups, 3))
    kw = dict(min_class=spec['min_class_molecules'], min_stratum=spec['min_stratum_molecules'])
    windows = {k: tuple(v) for k, v in fp['windows'].items()}
    rng = np.random.default_rng(fp['seed'])
    by_variant = {i: np.flatnonzero(vi == i) for i in range(nv)}
    act = load_activity({k: Path(v) for k, v in cfg['published_activity'].items()})
    print(f'Classified {len(keep):,} molecules', flush=True)

    # Per-base profiles for every variant (point estimates).
    eye = np.eye(211)
    profiles = {}
    for t in TRACKS:
        site, _ = measures(cs[t], vs, ss, cnt, eye, weights, **kw)
        profiles[f'{t}_C_std'] = site['C_std']
        profiles[f'{t}_M'] = site['M']

    effects, per_variant = [], {}
    for i, v in enumerate(variants):
        if v['class'] != 'SNV':
            continue
        idx = by_variant[i]
        if len(idx) < fp['min_molecules']:
            continue
        p0 = v['oligo_pos1'] - 21
        names, masks = footprint_masks(p0, windows, spec['mask_radius'])
        base = dict(source_variant=v['source_variant'], variant_id=v['variant_id'], oligo_pos1=v['oligo_pos1'],
                    genomic_pos1=GENOMIC_OFFSET + v['oligo_pos1'], ref=v['ref'], alt=v['alt'],
                    at_change=('retained_count' if (v['ref'] in 'AT') == (v['alt'] in 'AT') else 'loss' if v['ref'] in 'AT' else 'gain'),
                    context=homopolymer_context(seq, p0), molecules=len(idx))
        for src in act:
            base[f'activity_{src}'] = act[src].get((GENOMIC_OFFSET + v['oligo_pos1'], v['ref'], v['alt']))
        per_variant[i] = base
        for t in fp['tracks']:
            pt, _ = measures(cs[t][i], vs[i], ss[i], cnt[i], masks, weights, **kw)
            rf, _ = measures(cs[t][wt_i], vs[wt_i], ss[wt_i], cnt[wt_i], masks, weights, **kw)
            bt, _ = bootstrap(bits[t][idx], ones[idx], noscal[idx], groups[idx], clusters[idx], masks, weights,
                              n_groups, spec['bootstrap'], rng, **kw)
            pe, re_, be = endpoints(pt), endpoints(rf), endpoints(bt)
            for m in pe:
                delta, se, z = contrast(pe[m], re_[m], be[m], len(idx) / len(by_variant[wt_i]))
                for w, name in enumerate(names):
                    effects.append(dict(base, endpoint=f'{t}:{name}:{m}', track=t, window=name, measure=m,
                                        variant_value=pe[m][w], wt_value=re_[m][w], delta=delta[w], se=se[w], z=z[w]))
    for e in effects:
        e['p'] = float(2 * stats.norm.sf(abs(e['z']))) if np.isfinite(e['z']) else np.nan
    eps = sorted({e['endpoint'] for e in effects})
    for ep in eps:
        rows = [e for e in effects if e['endpoint'] == ep and np.isfinite(e['p'])]
        for e, q in zip(rows, bh([e['p'] for e in rows])):
            e['q_bh_within_endpoint'] = float(q)
    fam = [e for e in effects if e['endpoint'] in fp['primary_endpoints'] and np.isfinite(e['p'])]
    for e, q in zip(fam, bh([e['p'] for e in fam])):
        e['q_bh_footprint_family'] = float(q)
    print(f'Variant contrasts: {len(per_variant)} SNVs', flush=True)

    # WT pseudo-variants: whole barcode clusters matched to the median SNV cohort, random positions.
    wt_idx = by_variant[wt_i]
    target = int(np.median([len(by_variant[i]) for i in per_variant]))
    wt_cl, wt_inv = np.unique(clusters[wt_idx], return_inverse=True)
    members = np.split(np.argsort(wt_inv, kind='stable'), np.cumsum(np.bincount(wt_inv))[:-1])
    positions = sorted({v['oligo_pos1'] - 21 for v in variants if v['class'] == 'SNV'})
    calibration = []
    for draw in range(fp['aa_draws']):
        chosen, total = [], 0
        for c in rng.permutation(len(wt_cl)):
            chosen.append(members[c])
            total += len(members[c])
            if total >= target:
                break
        sel = wt_idx[np.concatenate(chosen)]
        p0 = int(rng.choice(positions))
        names, masks = footprint_masks(p0, windows, spec['mask_radius'])
        k = groups[sel]
        pcnt = np.bincount(k, minlength=n_groups).astype(float)
        pvs = grouped_sum(k, ones[sel], n_groups)
        for t in fp['tracks']:
            pcs = grouped_sum(k, bits[t][sel], n_groups)
            pt, _ = measures(pcs, pvs, ss[wt_i] * 0, pcnt, masks, weights, **kw)
            rf, _ = measures(cs[t][wt_i] - pcs, vs[wt_i] - pvs, ss[wt_i], cnt[wt_i] - pcnt, masks, weights, **kw)
            bt, _ = bootstrap(bits[t][sel], ones[sel], noscal[sel], k, clusters[sel], masks, weights, n_groups,
                              spec['bootstrap'], rng, **kw)
            pe, re_, be = endpoints(pt), endpoints(rf), endpoints(bt)
            for m in pe:
                delta, se, z = contrast(pe[m], re_[m], be[m], len(sel) / (len(wt_idx) - len(sel)))
                for w, name in enumerate(names):
                    calibration.append(dict(draw=draw, pseudo_oligo_pos1=p0 + 21, molecules=len(sel),
                                            endpoint=f'{t}:{name}:{m}', delta=delta[w], z=z[w]))
    summary_rows = []
    for ep in eps:
        z = np.array([c['z'] for c in calibration if c['endpoint'] == ep], float)
        z = z[np.isfinite(z)]
        if len(z) > 1:
            lo, hi = exact_ci(int((np.abs(z) > 1.96).sum()), len(z))
            summary_rows.append(dict(endpoint=ep, draws=len(z), mean_z=z.mean(), sd_z=z.std(ddof=1),
                                     rate_abs_z_gt_1_96=float((np.abs(z) > 1.96).mean()), rate_ci_low=lo, rate_ci_high=hi,
                                     gate_pass=bool(.85 <= z.std(ddof=1) <= 1.2 and hi <= .10)))
    print(f'Calibration: {fp["aa_draws"]} draws', flush=True)

    def clean(rows):
        return [{k: (None if isinstance(x, (float, np.floating)) and not np.isfinite(x) else
                     float(x) if isinstance(x, (float, np.floating)) else x) for k, x in r.items()} for r in rows]
    profile_rows = []
    for i, b in per_variant.items():
        for j in range(211):
            profile_rows.append(dict(source_variant=b['source_variant'], oligo_pos1=b['oligo_pos1'], site_oligo_pos1=j + 21,
                                     **{f'delta_{k}': float(arr[i, j] - arr[wt_i, j]) for k, arr in profiles.items()}))
    wt_profile = [dict(site_oligo_pos1=j + 21, **{f'wt_{k}': float(arr[wt_i, j]) for k, arr in profiles.items()}) for j in range(211)]
    with atomic_output(output) as out:
        write_table(out / 'effects.tsv.gz', clean(effects), list(dict.fromkeys(list(effects[0]) + ['p', 'q_bh_within_endpoint', 'q_bh_footprint_family'])))
        write_table(out / 'calibration_draws.tsv.gz', clean(calibration), list(calibration[0]))
        write_table(out / 'calibration_summary.tsv', clean(summary_rows), list(summary_rows[0]))
        write_table(out / 'profiles.tsv.gz', clean(profile_rows), list(profile_rows[0]))
        write_table(out / 'wt_profile.tsv', clean(wt_profile), list(wt_profile[0]))
        np.savez_compressed(out / 'profiles.npz', variants=np.array([variants[i]['source_variant'] for i in range(nv)]),
                            wt_index=wt_i, **profiles)
        write_json(out / 'summary.json', dict(
            schema_version='fsmpra.footprint_effects.v1', sources=[source_info(Path(config), True), source_info(extract / 'audit.json', True)],
            molecules=int(len(keep)), variants_tested=len(per_variant), footprints=fp, mask_radius=spec['mask_radius'],
            labeling_strata_edges=edges.tolist(), aa_target_molecules=target,
            interpretation='Footprint calls are plasmid-LLR (min_llr 5) assay annotations, not identified proteins. net = C_std - M '
                           'removes changes also seen on near-naked molecules (observability). Use z/p/q only for endpoints '
                           'whose calibration gate passes.'))
    print('Complete', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', required=True)
    ap.add_argument('--extract', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    run(a.config, a.extract, a.output)


if __name__ == '__main__':
    main()

"""Per-variant effects when each variant is (essentially) one plasmid clone: a WT clone-level null.

In the HEK293T library the top barcode holds a median 99% of a variant's molecules, so a
barcode-cluster bootstrap has one unit and would ignore clone-to-clone variation. Here a
variant's effect is compared with how much single WT clones deviate from the WT pool:

  d_v      = measure(variant molecules) - measure(all WT)          (same windows/standardization
                                                                     as site_effects.measures)
  s2_v     = within-clone sampling variance of d_v, molecule-level Poisson bootstrap
  tau2     = between-clone variance: mean over deep WT clones c of d_c^2 - s2_c, floored at 0,
             with d_c = clone c minus the WT pool without c, at the variant's own window masks
  z_v      = d_v / sqrt(tau2 + s2_v);  p from Student t with (WT clones - 1) df

Calibration gate (same rule as site_effects): each deep WT clone is a pseudo-variant at random
substituted positions, scored against tau2 re-estimated without it. Usable endpoints need
z SD within 0.85-1.20 and the upper exact 95% bound of the |z| > 1.96 rate <= 0.10, with the
bound computed on the number of clones (the independent unit), not the number of draws.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy import stats

from .io import atomic_output, read_json, source_info, write_json, write_table
from .processing import check_sources
from .site_effects import (COHORTS, GENOMIC_OFFSET, SCALARS, SITE_MEASURES, assign_groups, bh, exact_ci, grouped_sum,
                           homopolymer_context, load_activity, measures, window_masks)


def sums(idx, keys_local, called, valid, scal, n_groups, w=None):
    """Per-group sums for one molecule set; w (B, n) gives bootstrap-weighted sums with a leading B axis."""
    k = keys_local[idx]
    if w is None:
        return (grouped_sum(k, called[idx], n_groups), grouped_sum(k, valid[idx], n_groups),
                grouped_sum(k, scal[idx], n_groups), np.bincount(k, minlength=n_groups).astype(float))
    b = len(w)
    cs, vs = np.zeros((b, n_groups, 211)), np.zeros((b, n_groups, 211))
    ss, cnt = np.zeros((b, n_groups, len(SCALARS))), np.zeros((b, n_groups))
    for g in range(n_groups):
        sel = k == g
        if sel.any():
            wg = w[:, sel]
            cs[:, g], vs[:, g] = wg @ called[idx][sel], wg @ valid[idx][sel]
            ss[:, g], cnt[:, g] = wg @ scal[idx][sel], wg.sum(1)
    return cs, vs, ss, cnt


def run(config, extract, output):
    cfg = read_json(config)
    spec = cfg['analysis']
    cl = cfg['clone_null']
    u_max, m_min = spec['unmethylated_max_distal_m6a'], spec['methylated_min_distal_m6a']
    extract = Path(extract)
    audit = read_json(extract / 'audit.json')
    if audit.get('schema_version') != 'fsmpra.site_extract.v1':
        raise ValueError('Expected a site extraction directory')
    act_paths = {k: Path(v) for k, v in cfg['published_activity'].items()}
    sources = [source_info(Path(config), True), source_info(extract / 'audit.json', True), source_info(extract / 'meta.npz', True)]
    sources += [source_info(p, True) for p in act_paths.values()]
    sources += [source_info(Path(__file__).with_name(n), True) for n in ('site_effects_clone.py', 'site_effects.py', 'io.py')]
    variants = read_json(extract / 'variants.json')
    wt_i = next(i for i, v in enumerate(variants) if v['class'] == 'WT')
    seq = variants[wt_i]['sequence']
    meta = dict(np.load(extract / 'meta.npz'))
    cohort = meta[COHORTS[spec['cohort']]].astype(bool) & meta['evaluated'].astype(bool)
    vi_all = meta['variant_index'].astype(int)
    groups_all, weights, edges = assign_groups(meta['distal_m6a'].astype(float), (vi_all == wt_i) & cohort, u_max, m_min,
                                               spec['labeling_strata'])
    n_groups = 2 + len(weights)
    keep = np.flatnonzero(cohort & (groups_all >= 0))
    arrays = {n: np.load(extract / f'{n}.npy', mmap_mode='r') for n in ('scores', 'known', 'qualities', 'coverage')}
    v_ = arrays['known'][keep] & (arrays['qualities'][keep] >= spec['baseq_min'])
    called, valid = v_ & (arrays['scores'][keep] >= spec['ml_min']), v_
    cov = arrays['coverage'][keep]
    scal = np.stack([((cov & bit) > 0).mean(1) for bit in (2, 4, 1)], 1).astype(np.float32)
    vi, groups, clusters = vi_all[keep], groups_all[keep].astype(int), meta['cluster'][keep]
    kw = dict(min_class=spec['min_class_molecules'], min_stratum=spec['min_stratum_molecules'])
    windows = {k: tuple(v) for k, v in spec['windows'].items()}
    regions = {k: tuple(v) for k, v in spec['regions0'].items()}
    rng = np.random.default_rng(spec['seed'])
    B = cl['within_bootstrap']
    act = load_activity(act_paths)
    print(f'Cohort {cohort.sum():,}; classified {len(keep):,}', flush=True)

    wt_idx = np.flatnonzero(vi == wt_i)
    wt_sum = sums(wt_idx, groups, called, valid, scal, n_groups)
    cid, inv = np.unique(clusters[wt_idx], return_inverse=True)
    size = np.bincount(inv)
    deep = np.flatnonzero(size >= cl['min_wt_clone_molecules'])
    clone_idx = [wt_idx[inv == c] for c in deep]
    K = len(clone_idx)
    if K < 20:
        raise ValueError(f'Only {K} deep WT clones')
    # Per-clone point sums, leave-one-out WT pool sums and molecule bootstrap sums.
    c_pt = [sums(ix, groups, called, valid, scal, n_groups) for ix in clone_idx]
    c_ref = [tuple(a - b for a, b in zip(wt_sum, cs)) for cs in c_pt]
    c_bt = [sums(ix, groups, called, valid, scal, n_groups, rng.poisson(1., (B, len(ix))).astype(float)) for ix in clone_idx]
    print(f'{K} WT clones with >= {cl["min_wt_clone_molecules"]} molecules ({sum(map(len, clone_idx)):,} molecules)', flush=True)

    def clone_null(masks):
        """d (K, ...), s2 (K, ...) per endpoint for the WT clones at these masks."""
        d, s2 = {}, {}
        for k in range(K):
            pt_site, pt_scal = measures(*c_pt[k], masks, weights, **kw)
            rf_site, rf_scal = measures(*c_ref[k], masks, weights, **kw)
            bt_site, bt_scal = measures(*c_bt[k], masks, weights, **kw)
            for m in SITE_MEASURES:
                d.setdefault(m, []).append(pt_site[m] - rf_site[m])
                s2.setdefault(m, []).append(np.nanvar(bt_site[m], axis=0, ddof=1))
            for m in pt_scal:
                d.setdefault(m, []).append(np.asarray(pt_scal[m] - rf_scal[m]))
                s2.setdefault(m, []).append(np.asarray(np.nanvar(bt_scal[m], ddof=1)))
        return {m: np.array(x) for m, x in d.items()}, {m: np.array(x) for m, x in s2.items()}

    def tau2(d, s2, drop=None):
        sel = np.ones(len(d), bool) if drop is None else np.arange(len(d)) != drop
        with np.errstate(invalid='ignore'):
            return np.maximum(np.nanmean(d[sel] ** 2 - s2[sel], axis=0), 0.)

    df_t = K - 1
    effects, per_variant = [], {}
    tau_rows = []
    for i, v in enumerate(variants):
        if v['class'] != 'SNV':
            continue
        idx = np.flatnonzero(vi == i)
        if len(idx) < cl['min_variant_molecules']:
            continue
        p0 = v['oligo_pos1'] - 21
        names, masks = window_masks(seq, p0, windows, regions, spec['mask_radius'])
        pt_site, pt_scal = measures(*sums(idx, groups, called, valid, scal, n_groups), masks, weights, **kw)
        rf_site, rf_scal = measures(*wt_sum, masks, weights, **kw)
        bt_site, bt_scal = measures(*sums(idx, groups, called, valid, scal, n_groups,
                                          rng.poisson(1., (B, len(idx))).astype(float)), masks, weights, **kw)
        d_null, s2_null = clone_null(masks)
        vc = np.bincount(np.unique(clusters[idx], return_inverse=True)[1])
        base = dict(source_variant=v['source_variant'], variant_id=v['variant_id'], oligo_pos1=v['oligo_pos1'],
                    genomic_pos1=GENOMIC_OFFSET + v['oligo_pos1'], ref=v['ref'], alt=v['alt'],
                    at_change=('retained_count' if (v['ref'] in 'AT') == (v['alt'] in 'AT') else 'loss' if v['ref'] in 'AT' else 'gain'),
                    context=homopolymer_context(seq, p0), molecules=len(idx), clones=int((vc >= 2).sum()),
                    top_clone_fraction=float(vc.max() / len(idx)), molecules_U=int((groups[idx] == 0).sum()),
                    molecules_M=int((groups[idx] == 1).sum()), molecules_C=int((groups[idx] >= 2).sum()))
        for src in act:
            base[f'activity_{src}'] = act[src].get((GENOMIC_OFFSET + v['oligo_pos1'], v['ref'], v['alt']))
        per_variant[i] = base

        def add(m, point, ref, boot, w=None, name=None):
            t2 = tau2(d_null[m], s2_null[m])
            delta = point - ref
            with np.errstate(invalid='ignore'):
                s2 = np.nanvar(boot, axis=0, ddof=1)
                se = np.sqrt(t2 + s2)
            if w is not None:
                delta, se, t2, s2, point, ref = delta[w], se[w], t2[w], s2[w], point[w], ref[w]
            z = delta / se if se > 1e-9 else np.nan
            effects.append(dict(base, endpoint=name or m, window=name.split(':')[0] if name else 'element', measure=m,
                                sites=int(masks[w].sum()) if w is not None else 211, variant_value=float(point),
                                wt_value=float(ref), delta=float(delta), se=float(se), tau=float(np.sqrt(t2)),
                                se_within=float(np.sqrt(s2)), z=float(z),
                                p=float(2 * stats.t.sf(abs(z), df_t)) if np.isfinite(z) else np.nan))
        for m in SITE_MEASURES:
            for w, name in enumerate(names):
                add(m, pt_site[m], rf_site[m], bt_site[m], w, f'{name}:{m}')
        for m in pt_scal:
            add(m, np.asarray(pt_scal[m]), np.asarray(rf_scal[m]), bt_scal[m])
    print(f'Variant contrasts: {len(per_variant)} SNVs with >= {cl["min_variant_molecules"]} molecules', flush=True)
    endpoints = sorted({e['endpoint'] for e in effects})
    for ep in endpoints:
        rows = [e for e in effects if e['endpoint'] == ep and np.isfinite(e['p'])]
        if rows:
            for e, q in zip(rows, bh([e['p'] for e in rows])):
                e['q_bh_within_endpoint'] = float(q)
    primary = [e for e in effects if e['endpoint'] in spec['primary_endpoints'] and np.isfinite(e['p'])]
    if primary:
        for e, q in zip(primary, bh([e['p'] for e in primary])):
            e['q_bh_primary_family'] = float(q)

    # Calibration: every deep WT clone as a pseudo-variant at random SNV positions, tau2 without it.
    positions = sorted({v['oligo_pos1'] - 21 for v in variants if v['class'] == 'SNV'})
    calibration = []
    for draw in range(cl['calibration_positions']):
        p0 = int(rng.choice(positions))
        names, masks = window_masks(seq, p0, windows, regions, spec['mask_radius'])
        d_null, s2_null = clone_null(masks)
        for k in range(K):
            for m, dm in d_null.items():
                t2 = tau2(dm, s2_null[m], drop=k)
                with np.errstate(invalid='ignore', divide='ignore'):
                    sd = np.sqrt(t2 + s2_null[m][k])
                    z = np.where(sd > 1e-9, dm[k] / np.where(sd > 1e-9, sd, 1), np.nan)
                if m in SITE_MEASURES:
                    for w, name in enumerate(names):
                        calibration.append(dict(draw=draw, clone=int(cid[deep[k]]), clone_molecules=len(clone_idx[k]),
                                                pseudo_oligo_pos1=p0 + 21, endpoint=f'{name}:{m}', delta=float(dm[k][w]), z=float(z[w])))
                elif draw == 0:
                    calibration.append(dict(draw=draw, clone=int(cid[deep[k]]), clone_molecules=len(clone_idx[k]),
                                            pseudo_oligo_pos1=None, endpoint=m, delta=float(dm[k]), z=float(z)))
    summary_rows = []
    for ep in endpoints:
        z = np.array([c['z'] for c in calibration if c['endpoint'] == ep], float)
        z = z[np.isfinite(z)]
        if len(z) > 1:
            rate = float((np.abs(z) > 1.96).mean())
            lo, hi = exact_ci(int(round(rate * K)), K)
            summary_rows.append(dict(endpoint=ep, draws=len(z), clones=K, mean_z=z.mean(), sd_z=z.std(ddof=1),
                                     rate_abs_z_gt_1_96=rate, rate_ci_low=lo, rate_ci_high=hi,
                                     gate_pass=bool(.85 <= z.std(ddof=1) <= 1.2 and hi <= .10)))
    print(f'Calibration: {K} clones x {cl["calibration_positions"]} positions', flush=True)

    def clean(rows):
        return [{k: (None if isinstance(x, (float, np.floating)) and not np.isfinite(x) else
                     float(x) if isinstance(x, (float, np.floating)) else x) for k, x in r.items()} for r in rows]
    fields = list(dict.fromkeys(list(effects[0]) + ['q_bh_within_endpoint', 'q_bh_primary_family']))
    check_sources(sources)
    with atomic_output(output) as out:
        write_table(out / 'effects.tsv.gz', clean(effects), fields)
        write_table(out / 'calibration_draws.tsv.gz', clean(calibration), list(calibration[0]))
        write_table(out / 'calibration_summary.tsv', clean(summary_rows), list(summary_rows[0]))
        write_json(out / 'summary.json', dict(
            schema_version='fsmpra.site_effects_clone.v1', sources=sources, cohort=spec['cohort'],
            thresholds=dict(unmethylated_max_distal_m6a=u_max, methylated_min_distal_m6a=m_min,
                            labeling_strata_edges=edges.tolist(), labeling_strata_wt_weights=weights.tolist()),
            wt_clones=K, wt_clone_molecules=[len(x) for x in clone_idx], variants_tested=len(per_variant),
            t_df=df_t, clone_null=cl, primary_endpoints=spec['primary_endpoints'], windows=windows, regions0=regions,
            mask_radius=spec['mask_radius'],
            interpretation='Each variant is compared with the spread of single WT clones (between-clone tau plus within-clone '
                           'sampling). Variant and clone identity are confounded; z/p/q are usable only for endpoints whose '
                           'calibration_summary gate passes.'))
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

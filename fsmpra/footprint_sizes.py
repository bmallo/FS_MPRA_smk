"""Footprint sizes along the element: a WT position x length map, a frozen feature list, per-feature variant tests.

The 10-60 bp "TF" bit lumps every sub-nucleosomal call together. Here each plasmid-LLR call keeps
its midpoint (oligo order) and length, so footprints of different sizes at different places can
be told apart.

map: WT molecules only (variant-blind). Calls per chromatin molecule are binned by midpoint (1 bp)
  x length (1 bp, 10-90), standardized over the 10 WT labeling strata (C_std) and, separately,
  for near-naked (M) molecules. Peaks of the smoothed C_std density become features (a midpoint
  window plus a length range); see config 'map' for the exact rule. Discovery uses half 0 of the
  WT barcode clusters; the other half is kept for the variant tests.

effects: for every SNV with enough molecules and every frozen feature,
  occupancy = fraction of molecules with a call inside the feature box, and
  size = mean length of 10-90 bp calls centred in the feature window (molecules with one).
  Both use site_effects.measures (C_std, M) against the reference half of WT, a barcode-cluster
  Poisson bootstrap, and WT pseudo-variants from the reference half for the calibration gate.
  The occupancy endpoint is the observability residual dC_std - b_f * dM with a per-feature slope.

Usage:
  python -m fsmpra.footprint_sizes map --config config/ldlr_footprint_sizes.json --extract <dir> \\
      --calls <a.parquet> [...] --element 1234 1445 [--reverse] [--all-wt] [--reads ids.txt] --output <dir>
  python -m fsmpra.footprint_sizes effects --config config/ldlr_footprint_sizes.json --extract <dir> \\
      --calls <a.parquet> [...] --element 1234 1445 --features <map>/features.tsv --output <dir>
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy import ndimage, stats

from .io import atomic_output, open_text, read_json, source_info, write_json, write_table
from .site_effects import (COHORTS, GENOMIC_OFFSET, SCALARS, assign_groups, bh, contrast, exact_ci, grouped_sum,
                           load_activity, measures, ratio)

WIDTH = 211
OLIGO0 = 21
SIZE_MIN_CALLED = 10


def load_calls(extract, calls, element, reverse, min_llr, length_range, reads=None):
    """(row, midpoint index 0..210 in oligo order, length) for calls of the given length centred in the element."""
    import pandas as pd
    with open_text(Path(extract) / 'molecules.tsv.gz') as h:
        index = {r['read_id']: int(r['row']) for r in csv.DictReader(h, delimiter='\t')}
    if reads is not None:
        index = {k: v for k, v in index.items() if k in reads}
    lo, hi = length_range
    parts = []
    for p in calls:
        df = pd.read_parquet(p, columns=['min_llr', 'read_name', 'start', 'end'])
        df = df[np.isclose(df.min_llr, min_llr)]
        s0, e = df.start.to_numpy() - 1, df.end.to_numpy()
        j = midpoint_index(s0, e, element, reverse)
        length = e - s0
        ok = (length >= lo) & (length <= hi) & (j >= 0) & (j < WIDTH)
        rows = df.read_name[ok].map(index).to_numpy()
        has = ~pd.isna(rows)
        parts.append((rows[has].astype(np.int64), j[ok][has], length[ok][has].astype(np.int64)))
    rows, mid, length = (np.concatenate(x) for x in zip(*parts))
    return rows, mid, length, np.array(sorted(index.values()), dtype=np.int64)


def midpoint_index(s0, e, element, reverse=False):
    """Element index of a call midpoint ([s0, e) 0-based reference); reverse maps a minus-strand element to oligo order."""
    a, b = element
    mid = (np.asarray(s0) + np.asarray(e) - 1) / 2
    return np.rint((b - 1) - mid if reverse else mid - a).astype(np.int64)


def classify(extract, spec):
    variants = read_json(extract / 'variants.json')
    wt_i = next(i for i, v in enumerate(variants) if v['class'] == 'WT')
    meta = dict(np.load(extract / 'meta.npz'))
    cohort = meta[COHORTS[spec['cohort']]].astype(bool) & meta['evaluated'].astype(bool)
    vi = meta['variant_index'].astype(int)
    groups, weights, edges = assign_groups(meta['distal_m6a'].astype(float), (vi == wt_i) & cohort,
                                           spec['unmethylated_max_distal_m6a'], spec['methylated_min_distal_m6a'],
                                           spec['labeling_strata'])
    groups = np.where(cohort, groups, -1).astype(int)
    return variants, wt_i, vi, groups, weights, edges, meta['cluster'].astype(np.int64)


def wt_halves(clusters, is_wt, seed):
    """0/1 per WT molecule by a seeded permutation of WT clusters; -1 elsewhere."""
    cl = np.unique(clusters[is_wt])
    h = np.zeros(len(cl), dtype=np.int8)
    h[np.random.default_rng(seed).permutation(len(cl))[len(cl) // 2:]] = 1
    half = np.full(len(clusters), -1, dtype=np.int8)
    half[is_wt] = h[np.searchsorted(cl, clusters[is_wt])]
    return half


def density(rows, mid, length, groups, select, n_groups, length_range):
    """Calls per molecule by group x midpoint x length, for molecules in select; and molecules per group."""
    lo, hi = length_range
    L = hi - lo + 1
    k = select[rows]
    idx = (groups[rows[k]].astype(np.int64) * WIDTH + mid[k]) * L + (length[k] - lo)
    counts = np.bincount(idx, minlength=n_groups * WIDTH * L).reshape(n_groups, WIDTH, L).astype(float)
    mol = np.bincount(groups[select], minlength=n_groups).astype(float)
    return counts / np.maximum(mol, 1)[:, None, None], mol


def standardize(dens, weights):
    return np.einsum('s,sml->ml', weights, dens[2:])


def find_features(C, M, mp, length_range):
    """Peaks of the smoothed C_std density above a local median background -> boxes (0-based midpoint index, lengths in bp)."""
    lo = length_range[0]
    sig, nb = mp['smoothing_sigma_bp'], mp['peak_neighbourhood_bp']
    sm = ndimage.gaussian_filter(C, (sig['midpoint'], sig['length']), mode='nearest')
    mx = ndimage.maximum_filter(sm, size=(2 * nb['midpoint'] + 1, 2 * nb['length'] + 1), mode='nearest')
    cand = np.argwhere((sm == mx) & (sm > 0))

    def extent(profile, i, lim, half):
        a = i
        while a > 0 and i - (a - 1) <= lim[1] and profile[a - 1] >= half:
            a -= 1
        b = i
        while b < len(profile) - 1 and (b + 1) - i <= lim[1] and profile[b + 1] >= half:
            b += 1
        return min(a, max(0, i - lim[0])), max(b, min(len(profile) - 1, i + lim[0]))

    bn = mp['background_neighbourhood_bp']
    bg = ndimage.median_filter(sm, size=(2 * bn['midpoint'] + 1, 2 * bn['length'] + 1), mode='nearest')
    boxes = []
    for m, l in cand:
        if sm[m, l] <= bg[m, l]:
            continue
        half = bg[m, l] + (sm[m, l] - bg[m, l]) / 2
        m0, m1 = extent(sm[:, l], m, mp['box_midpoint_halfwidth_bp'], half)
        l0, l1 = extent(sm[m, :], l, mp['box_length_halfwidth_bp'], half)
        box = (slice(m0, m1 + 1), slice(l0, l1 + 1))
        mc, mm = float((C[box] - bg[box]).sum()), float(M[box].sum())
        boxes.append(dict(mid_peak=int(m), len_peak=int(l + lo), mid_lo=int(m0), mid_hi=int(m1), len_lo=int(l0 + lo),
                          len_hi=int(l1 + lo), peak_density=float(sm[m, l]), background_density=float(bg[m, l]),
                          occupancy_C_std=float(C[box].sum()), occupancy_above_background=mc, occupancy_M=mm))
    boxes.sort(key=lambda b: -b['occupancy_above_background'])
    kept = []
    for b in boxes:
        if b['occupancy_above_background'] < mp['min_box_occupancy_C_std']:
            continue
        area = (b['mid_hi'] - b['mid_lo'] + 1) * (b['len_hi'] - b['len_lo'] + 1)
        if any(overlap(b, k) > mp['max_box_overlap'] * area for k in kept):
            continue
        kept.append(b)
    kept.sort(key=lambda b: (b['mid_peak'], b['len_peak']))
    for i, b in enumerate(kept):
        b['feature'] = f'F{i + 1:02d}'
        b['c_over_m'] = b['occupancy_C_std'] / b['occupancy_M'] if b['occupancy_M'] > 0 else np.inf
        b['flag'] = 'chromatin' if b['c_over_m'] >= mp['near_naked_ratio_flag'] else 'shared_with_near_naked'
    return kept, sm, bg


def overlap(a, b):
    dm = min(a['mid_hi'], b['mid_hi']) - max(a['mid_lo'], b['mid_lo']) + 1
    dl = min(a['len_hi'], b['len_hi']) - max(a['len_lo'], b['len_lo']) + 1
    return max(dm, 0) * max(dl, 0)


def box_mass(D, f, lo):
    return float(D[f['mid_lo']:f['mid_hi'] + 1, f['len_lo'] - lo:f['len_hi'] - lo + 1].sum())


def run_map(config, extract, calls, element, reverse, output, all_wt=False, reads_file=None):
    cfg = read_json(config)
    spec = read_json(Path(cfg['base_config']))['analysis']
    lr, mp = tuple(cfg['calls']['length_range']), cfg['map']
    extract = Path(extract)
    variants, wt_i, vi, groups, weights, edges, clusters = classify(extract, spec)
    n_groups = 2 + len(weights)
    reads = set(Path(reads_file).read_text().split()) if reads_file else None
    rows, mid, length, rows_in = load_calls(extract, calls, element, reverse, cfg['calls']['min_llr'], lr, reads)
    is_wt = (vi == wt_i) & (groups >= 0)
    half = wt_halves(clusters, is_wt, cfg['split']['seed'])
    restrict = np.zeros(len(vi), bool)
    restrict[rows_in] = True
    sets = {'all': is_wt & restrict} if all_wt else {'discovery': (half == 0) & restrict, 'reference': (half == 1) & restrict}
    maps, mols = {}, {}
    for name, sel in sets.items():
        dens, mol = density(rows, mid, length, groups, sel, n_groups, lr)
        maps[name] = dict(C=standardize(dens, weights), M=dens[1], U=dens[0])
        mols[name] = mol
    disc = 'all' if all_wt else 'discovery'
    feats, smooth, background = find_features(maps[disc]['C'], maps[disc]['M'], mp, lr)
    for f in feats:
        for name in maps:
            f[f'occupancy_C_std_{name}'] = box_mass(maps[name]['C'], f, lr[0])
            f[f'occupancy_M_{name}'] = box_mass(maps[name]['M'], f, lr[0])
    rows_out = [dict(feature=f['feature'], mid_peak_oligo=f['mid_peak'] + OLIGO0, mid_lo_oligo=f['mid_lo'] + OLIGO0,
                     mid_hi_oligo=f['mid_hi'] + OLIGO0, **{k: v for k, v in f.items() if k not in ('feature',)}) for f in feats]
    with atomic_output(output) as out:
        write_table(out / 'features.tsv', rows_out, list(rows_out[0]))
        np.savez_compressed(out / 'map.npz', lengths=np.arange(lr[0], lr[1] + 1), oligo=np.arange(WIDTH) + OLIGO0,
                            smoothed_discovery_C=smooth, background_discovery_C=background, **{f'{n}_{k}': v for n, d in maps.items() for k, v in d.items()})
        write_json(out / 'summary.json', dict(
            schema_version='fsmpra.footprint_sizes_map.v1',
            sources=[source_info(Path(config), True), source_info(extract / 'audit.json', True)] + [source_info(Path(p)) for p in calls],
            element=list(element), reversed=reverse, calls_in_range=int(len(rows)), molecules_restricted=int(restrict.sum()),
            sets={n: dict(molecules=int(sets[n].sum()), per_group=mols[n].tolist(), clusters=int(len(np.unique(clusters[sets[n]]))))
                  for n in sets},
            features=len(feats), labeling_strata_edges=edges.tolist(), map_rule=mp,
            interpretation='Features are frozen from the discovery set only; occupancy_C_std is the fraction of chromatin '
                           'molecules with a call in the box (one call per molecule per box assumed).'))
    print(f'{len(feats)} features from {int(sets[disc].sum()):,} WT molecules', flush=True)


def read_features(path):
    with open(path) as h:
        rows = list(csv.DictReader(h, delimiter='\t'))
    return [dict(feature=r['feature'], mid_lo=int(r['mid_lo']), mid_hi=int(r['mid_hi']), mid_peak=int(r['mid_peak']),
                 len_lo=int(r['len_lo']), len_hi=int(r['len_hi']), len_peak=int(r['len_peak']), flag=r['flag']) for r in rows]


def feature_values(rows, mid, length, n, feats):
    """Per molecule x feature: occupancy (bool), size sum and count of window calls."""
    F = len(feats)
    occ = np.zeros((n, F), dtype=bool)
    size_sum, size_n = np.zeros((n, F)), np.zeros((n, F))
    for f, ft in enumerate(feats):
        win = (mid >= ft['mid_lo']) & (mid <= ft['mid_hi'])
        box = win & (length >= ft['len_lo']) & (length <= ft['len_hi'])
        occ[rows[box], f] = True
        size_sum[:, f] = np.bincount(rows[win], weights=length[win], minlength=n)
        size_n[:, f] = np.bincount(rows[win], minlength=n)
    has = size_n > 0
    return {'occupancy': (occ.astype(float), np.ones((n, F))), 'size': (np.where(has, size_sum / np.maximum(size_n, 1), 0.), has.astype(float))}


def bootstrap_f(called, valid, groups, clusters, weights, n_groups, b, rng, **kw):
    """Barcode-cluster Poisson bootstrap for an (n, F) endpoint; returns measures() site dict of (B, F) arrays."""
    _, inv = np.unique(clusters, return_inverse=True)
    w = rng.poisson(1., (b, inv.max() + 1)).astype(np.float64)[:, inv]
    F = called.shape[1]
    cs, vs, cnt = np.zeros((b, n_groups, F)), np.zeros((b, n_groups, F)), np.zeros((b, n_groups))
    for g in range(n_groups):
        sel = groups == g
        if sel.any():
            wg = w[:, sel]
            cs[:, g], vs[:, g], cnt[:, g] = wg @ called[sel], wg @ valid[sel], wg.sum(1)
    site, _ = measures(cs, vs, np.zeros((b, n_groups, len(SCALARS))), cnt, np.eye(F), weights, **kw)
    return site


def point(cs, vs, cnt, weights, n_groups, **kw):
    F = cs.shape[-1]
    site, _ = measures(cs, vs, np.zeros(cnt.shape + (len(SCALARS),)), cnt, np.eye(F), weights, **kw)
    return site


def slope(dM, dC):
    ok = np.isfinite(dM) & np.isfinite(dC)
    if ok.sum() < 10 or np.ptp(dM[ok]) == 0:
        return 0.
    return float(np.clip(np.polyfit(dM[ok], dC[ok], 1)[0], 0., 1.))


def gate(z):
    z = z[np.isfinite(z)]
    if len(z) < 2:
        return dict(draws=len(z))
    lo, hi = exact_ci(int((np.abs(z) > 1.96).sum()), len(z))
    return dict(draws=len(z), mean_z=float(z.mean()), sd_z=float(z.std(ddof=1)), rate_abs_z_gt_1_96=float((np.abs(z) > 1.96).mean()),
                rate_ci_low=lo, rate_ci_high=hi, gate_pass=bool(.85 <= z.std(ddof=1) <= 1.2 and hi <= .10))


def pooled_tests(effects, feats):
    """Per feature: mean occupancy residual over SNVs inside vs outside the footprint span; BH across features per set."""
    out = []
    for inside in (True, False):
        rows = []
        for ft in feats:
            e = [r for r in effects if r['feature'] == ft['feature'] and r['snv_in_footprint'] == inside
                 and np.isfinite(r['occupancy_resid_delta']) and np.isfinite(r['occupancy_resid_se'])]
            if not e:
                continue
            d = np.array([r['occupancy_resid_delta'] for r in e])
            se = np.sqrt(np.sum(np.array([r['occupancy_resid_se'] for r in e]) ** 2)) / len(e)
            z = d.mean() / se if se > 0 else np.nan
            rows.append(dict(feature=ft['feature'], snvs='inside' if inside else 'outside', k=len(e), mean_delta=float(d.mean()),
                             se=float(se), z=float(z), p=float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan,
                             fraction_negative=float((d < 0).mean())))
        fin = [r for r in rows if np.isfinite(r['p'])]
        for r, q in zip(fin, bh([r['p'] for r in fin])):
            r['q_bh'] = float(q)
        out += rows
    return out


def run_effects(config, extract, calls, element, reverse, features, output):
    cfg = read_json(config)
    base = read_json(Path(cfg['base_config']))
    spec, ef, lr = base['analysis'], cfg['effects'], tuple(cfg['calls']['length_range'])
    extract = Path(extract)
    variants, wt_i, vi, groups, weights, edges, clusters = classify(extract, spec)
    n_groups, nv = 2 + len(weights), len(variants)
    feats = read_features(features)
    F = len(feats)
    rows, mid, length, _ = load_calls(extract, calls, element, reverse, cfg['calls']['min_llr'], lr)
    half = wt_halves(clusters, (vi == wt_i) & (groups >= 0), cfg['split']['seed'])
    keep = np.flatnonzero((groups >= 0) & ((vi != wt_i) | (half == 1)))
    pos = np.full(len(vi), -1, dtype=np.int64)
    pos[keep] = np.arange(len(keep))
    ok = pos[rows] >= 0
    vals = feature_values(pos[rows[ok]], mid[ok], length[ok], len(keep), feats)
    vi_k, g_k, cl_k = vi[keep], groups[keep], clusters[keep]
    keys = vi_k * n_groups + g_k
    cnt = np.bincount(keys, minlength=nv * n_groups).reshape(nv, n_groups).astype(float)
    kw = dict(min_class=spec['min_class_molecules'], min_stratum=spec['min_stratum_molecules'])
    pts = {}
    for e, (c, v) in vals.items():
        cs = grouped_sum(keys, c, nv * n_groups).reshape(nv, n_groups, F)
        vs = grouped_sum(keys, v, nv * n_groups).reshape(nv, n_groups, F)
        pts[e] = dict(sums=(cs, vs), site=point(cs, vs, cnt, weights, n_groups, **kw))
    by_variant = {i: np.flatnonzero(vi_k == i) for i in range(nv)}
    print(f'{len(keep):,} molecules (WT reference half {len(by_variant[wt_i]):,}), {F} features', flush=True)

    # Per-SNV x feature size maps (C_std calls per molecule by midpoint x length), for descriptive figures.
    L = lr[1] - lr[0] + 1
    rk = pos[rows[ok]]
    idx = (keys[rk].astype(np.int64) * WIDTH + mid[ok]) * L + (length[ok] - lr[0])
    counts = np.bincount(idx, minlength=nv * n_groups * WIDTH * L).reshape(nv, n_groups, WIDTH, L).astype(np.float32)
    counts /= np.maximum(cnt, 1)[:, :, None, None].astype(np.float32)
    size_maps_C = np.einsum('s,vsml->vml', weights.astype(np.float32), counts[:, 2:])
    size_maps_M = counts[:, 1].copy()
    del counts

    rng = np.random.default_rng(ef['seed'])
    act = load_activity({k: Path(v) for k, v in base['published_activity'].items()})
    seq = variants[wt_i]['sequence']
    snvs = [i for i, v in enumerate(variants) if v['class'] == 'SNV' and len(by_variant[i]) >= ef['min_molecules']]
    n_ref = len(by_variant[wt_i])
    boots = {}
    for i in snvs:
        idx_i = by_variant[i]
        boots[i] = {e: bootstrap_f(vals[e][0][idx_i], vals[e][1][idx_i], g_k[idx_i], cl_k[idx_i], weights, n_groups,
                                   ef['bootstrap'], rng, **kw) for e in vals}
    print(f'Bootstrapped {len(snvs)} SNVs', flush=True)
    occ = pts['occupancy']['site']
    dC = np.array([occ['C_std'][i] - occ['C_std'][wt_i] for i in snvs])
    dM = np.array([occ['M'][i] - occ['M'][wt_i] for i in snvs])
    b = np.array([slope(dM[:, f], dC[:, f]) for f in range(F)])

    def endpoint_set(site_v, site_r, boot, ratio_, n_chrom, n_size):
        """Occupancy SE is floored at the binomial SE of WT occupancy at this cohort's chromatin depth (a resample with
        no chromatin events otherwise gives SE ~ 0). Size uses pooled chromatin (C_raw), only with >= 10 molecules with a call."""
        out = {}
        occ_v, occ_r, occ_b = site_v['occupancy'], site_r['occupancy'], boot['occupancy']
        p_ref = np.clip(occ_r['C_std'], 0, 1)
        floor = np.sqrt(p_ref * (1 - p_ref) / max(n_chrom, 1)) * np.sqrt(1 + ratio_)
        for name, v_, r_, b_ in (('occupancy_resid', occ_v['C_std'] - b * occ_v['M'], occ_r['C_std'] - b * occ_r['M'],
                                  occ_b['C_std'] - b * occ_b['M']),
                                 ('occupancy_C', occ_v['C_std'], occ_r['C_std'], occ_b['C_std'])):
            d, se, _ = contrast(v_, r_, b_, ratio_)
            se = np.fmax(se, floor)
            out[name] = (d, se, np.where(se > 0, d / np.where(se > 0, se, 1), np.nan))
        out['occupancy_M'] = contrast(occ_v['M'], occ_r['M'], occ_b['M'], ratio_)
        d, se, z = contrast(site_v['size']['C_raw'], site_r['size']['C_raw'], boot['size']['C_raw'], ratio_)
        few = n_size < SIZE_MIN_CALLED
        out['size'] = tuple(np.where(few, np.nan, x) for x in (d, se, z))
        return out

    sv = lambda i: {e: {m: pts[e]['site'][m][i] for m in ('C_std', 'M', 'C_raw')} for e in pts}
    n_chrom = cnt[:, 2:].sum(1)
    n_size = pts['size']['sums'][1][:, 2:].sum(1)
    effects = []
    for i in snvs:
        v = variants[i]
        p0 = v['oligo_pos1'] - OLIGO0
        res = endpoint_set(sv(i), sv(wt_i), boots[i], len(by_variant[i]) / n_ref, n_chrom[i], n_size[i])
        common = dict(source_variant=v['source_variant'], oligo_pos1=v['oligo_pos1'], ref=v['ref'], alt=v['alt'],
                      at_change=('retained_count' if (v['ref'] in 'AT') == (v['alt'] in 'AT') else 'loss' if v['ref'] in 'AT' else 'gain'),
                      molecules=len(by_variant[i]),
                      **{f'activity_{s}': act[s].get((GENOMIC_OFFSET + v['oligo_pos1'], v['ref'], v['alt'])) for s in act})
        for f, ft in enumerate(feats):
            span = (ft['mid_peak'] - ft['len_peak'] / 2, ft['mid_peak'] + ft['len_peak'] / 2)
            row = dict(common, feature=ft['feature'], feature_flag=ft['flag'], distance_to_feature=p0 - ft['mid_peak'],
                       snv_in_footprint=bool(span[0] <= p0 <= span[1]), slope_b=b[f],
                       wt_occupancy_C_std=pts['occupancy']['site']['C_std'][wt_i][f],
                       wt_size_C_raw=pts['size']['site']['C_raw'][wt_i][f], size_called_molecules=n_size[i][f],
                       expected_events=pts['occupancy']['site']['C_std'][wt_i][f] * cnt[i, 2:].sum())
            for name, (d, se, z) in res.items():
                row.update({f'{name}_delta': d[f], f'{name}_se': se[f], f'{name}_z': z[f],
                            f'{name}_p': float(2 * stats.norm.sf(abs(z[f]))) if np.isfinite(z[f]) else np.nan})
            effects.append(row)
    for name, fam in (('occupancy_resid', 'q_primary'), ('size', 'q_secondary')):
        fin = [e for e in effects if np.isfinite(e[f'{name}_p'])]
        for e, q in zip(fin, bh([e[f'{name}_p'] for e in fin])):
            e[fam] = float(q)
    print(f'{len(effects):,} SNV x feature contrasts', flush=True)
    pooled = pooled_tests(effects, feats)

    # WT pseudo-variants from the reference half, median SNV cohort size.
    wt_idx = by_variant[wt_i]
    target = int(np.median([len(by_variant[i]) for i in snvs]))
    wt_cl, wt_inv = np.unique(cl_k[wt_idx], return_inverse=True)
    members = np.split(np.argsort(wt_inv, kind='stable'), np.cumsum(np.bincount(wt_inv))[:-1])
    calibration = []
    for draw in range(ef['aa_draws']):
        chosen, total = [], 0
        for c in rng.permutation(len(wt_cl)):
            chosen.append(members[c])
            total += len(members[c])
            if total >= target:
                break
        sel = wt_idx[np.concatenate(chosen)]
        k = g_k[sel]
        pc = np.bincount(k, minlength=n_groups).astype(float)
        site_p, site_r, boot = {}, {}, {}
        for e, (c, v) in vals.items():
            ps, pv = grouped_sum(k, c[sel], n_groups), grouped_sum(k, v[sel], n_groups)
            cs, vs = pts[e]['sums']
            site_p[e] = point(ps, pv, pc, weights, n_groups, **kw)
            site_r[e] = point(cs[wt_i] - ps, vs[wt_i] - pv, cnt[wt_i] - pc, weights, n_groups, **kw)
            boot[e] = bootstrap_f(c[sel], v[sel], k, cl_k[sel], weights, n_groups, ef['bootstrap'], rng, **kw)
        res = endpoint_set(site_p, site_r, boot, len(sel) / (n_ref - len(sel)), pc[2:].sum(),
                           grouped_sum(k, vals['size'][1][sel], n_groups)[2:].sum(0))
        for f, ft in enumerate(feats):
            calibration.append(dict(draw=draw, molecules=len(sel), feature=ft['feature'],
                                    **{f'{name}_z': z[f] for name, (_, _, z) in res.items()},
                                    **{f'{name}_delta': d[f] for name, (d, _, _) in res.items()}))
    summary_rows = []
    for name in ('occupancy_resid', 'occupancy_C', 'occupancy_M', 'size'):
        allz = np.array([c[f'{name}_z'] for c in calibration], float)
        summary_rows.append(dict(endpoint=name, feature='all', **gate(allz)))
        for ft in feats:
            summary_rows.append(dict(endpoint=name, feature=ft['feature'],
                                     **gate(np.array([c[f'{name}_z'] for c in calibration if c['feature'] == ft['feature']], float))))
    print(f'Calibration: {ef["aa_draws"]} draws', flush=True)

    def clean(rs):
        return [{k: (None if isinstance(x, (float, np.floating)) and not np.isfinite(x) else
                     float(x) if isinstance(x, (float, np.floating)) else
                     int(x) if isinstance(x, (np.integer,)) else bool(x) if isinstance(x, np.bool_) else x) for k, x in r.items()} for r in rs]
    cols = list(dict.fromkeys(list(effects[0]) + ['q_primary', 'q_secondary']))
    keys_cal = list(dict.fromkeys(k for r in summary_rows for k in r))
    with atomic_output(output) as out:
        write_table(out / 'effects.tsv.gz', clean(effects), cols)
        write_table(out / 'feature_tests.tsv', clean(pooled), list(pooled[0]))
        write_table(out / 'calibration_draws.tsv.gz', clean(calibration), list(calibration[0]))
        write_table(out / 'calibration_summary.tsv', clean(summary_rows), keys_cal)
        np.savez_compressed(out / 'size_maps.npz', variants=np.array([v['source_variant'] for v in variants]), wt_index=wt_i,
                            lengths=np.arange(lr[0], lr[1] + 1), C_std=size_maps_C, M=size_maps_M,
                            molecules=np.array([len(by_variant[i]) for i in range(nv)]))
        write_json(out / 'summary.json', dict(
            schema_version='fsmpra.footprint_sizes_effects.v1',
            sources=[source_info(Path(config), True), source_info(extract / 'audit.json', True), source_info(Path(features), True)]
                    + [source_info(Path(p)) for p in calls],
            molecules=int(len(keep)), wt_reference_molecules=int(n_ref), snvs_tested=len(snvs), features=F,
            slopes={ft['feature']: float(b[f]) for f, ft in enumerate(feats)}, aa_target_molecules=target, effects_rule=ef,
            interpretation='Exploratory size-resolved footprints (config interpretation). occupancy_resid is the primary '
                           'endpoint (BH over all SNV x feature, q_primary); size is secondary and not corrected for '
                           'observability (q_secondary). Use z/p/q only where the calibration gate passes.'))
    print('Complete', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('command', choices=['map', 'effects'])
    ap.add_argument('--config', required=True)
    ap.add_argument('--extract', required=True)
    ap.add_argument('--calls', nargs='+', required=True)
    ap.add_argument('--element', nargs=2, type=int, required=True)
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--all-wt', action='store_true', help='map: use every WT molecule (no discovery/reference split)')
    ap.add_argument('--reads', help='map: restrict to these read ids (one per line)')
    ap.add_argument('--features')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    if a.command == 'map':
        run_map(a.config, a.extract, a.calls, tuple(a.element), a.reverse, a.output, a.all_wt, a.reads)
    else:
        run_effects(a.config, a.extract, a.calls, tuple(a.element), a.reverse, a.features, a.output)


if __name__ == '__main__':
    main()

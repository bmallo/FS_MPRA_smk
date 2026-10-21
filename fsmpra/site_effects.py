"""Proxy-control-normalized per-site m6A effects for every designed SNV.

Molecules are classed from distal (non-element) retained m6A only, so a variant
cannot move its own molecules between classes:

- U, distal m6A <= u_max: plasmids that apparently never met Hia5. Element calls
  estimate sequence-dependent false-positive calling.
- M, distal m6A >= m_min: heavily labeled, presumably unchromatinized plasmids.
  Element calls estimate sequence-dependent labeling/calling near saturation.
- C, in between: the chromatin signal, directly standardized over WT strata of
  distal labeling so labeling composition does not masquerade as an effect.

Normalized accessibility is A = (C - U) / (M - U) per window. U and M are
proxies from the same reaction, not separately prepared controls; M molecules
may include open nuclear plasmids. Uncertainty is a barcode-cluster Poisson
bootstrap of each SNV with WT fixed (WT has ~250x more molecules), calibrated
by WT pseudo-variants drawn as whole barcode clusters.
"""
import argparse
import csv
import warnings
from pathlib import Path
import numpy as np
from scipy import stats
from .io import atomic_output, open_text, read_json, source_info, write_json, write_table
from .library_stats import circular_shift
from .processing import check_sources

GENOMIC_OFFSET = 11089264  # genomic_pos1 = offset + oligo_pos1 (GRCh38, forward strand)
SITE_MEASURES = ['U', 'M', 'C_raw', 'C_std', 'ALL', 'A_norm']
SCALARS = ['msp', 'nuc_85_250', 'tf_small']
COHORTS = {'eligible': 'eligible', 'eligible_with_barcode': 'eligible_with_barcode', 'sequence_qualified': None}


def assign_groups(distal, is_wt, u_max, m_min, strata):
    """Group 0 = U, 1 = M, 2.. = C strata from WT C-molecule quantiles; -1 = unclassifiable."""
    finite = np.isfinite(distal)
    u, m = finite & (distal <= u_max), finite & (distal >= m_min)
    c = finite & ~u & ~m
    if not (c & is_wt).any():
        raise ValueError('No WT molecules in the chromatin class')
    edges = np.unique(np.quantile(distal[c & is_wt], np.linspace(0, 1, strata + 1)[1:-1]))
    groups = np.full(len(distal), -1, dtype=np.int16)
    groups[u], groups[m] = 0, 1
    groups[c] = 2 + np.searchsorted(edges, distal[c], side='right')
    n_strata = len(edges) + 1
    wt_c = np.bincount(groups[c & is_wt] - 2, minlength=n_strata).astype(float)
    return groups, wt_c / wt_c.sum(), edges


def grouped_sum(keys, values, n_keys):
    """Sum rows of values by integer key without np.add.at."""
    out = np.zeros((n_keys,) + values.shape[1:], dtype=np.float64)
    if not len(keys):
        return out
    order = np.argsort(keys, kind='stable')
    k = keys[order]
    starts = np.r_[0, np.flatnonzero(k[1:] != k[:-1]) + 1]
    out[k[starts]] = np.add.reduceat(values[order].astype(np.float64), starts, axis=0)
    return out


def ratio(num, den):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)


def measures(cs, vs, ss, cnt, masks, weights, min_class=30, min_stratum=10):
    """Window and scalar measures from per-group sums.

    cs, vs: (..., G, 211) called / valid sums; ss: (..., G, 3) scalar sums;
    cnt: (..., G) molecule weights; masks: (W, 211) float site indicators.
    Returns {name: (..., W)} for site measures and {name: (...)} for scalars."""
    wc, wv = cs @ masks.T, vs @ masks.T                  # (..., G, W)
    enough_u = (cnt[..., 0] >= min_class)[..., None]
    enough_m = (cnt[..., 1] >= min_class)[..., None]
    strata_ok = (cnt[..., 2:] >= min_stratum).all(-1)[..., None]
    out = {'U': np.where(enough_u, ratio(wc[..., 0, :], wv[..., 0, :]), np.nan),
           'M': np.where(enough_m, ratio(wc[..., 1, :], wv[..., 1, :]), np.nan),
           'C_raw': ratio(wc[..., 2:, :].sum(-2), wv[..., 2:, :].sum(-2)),
           'C_std': np.where(strata_ok, np.einsum('...sw,s->...w', ratio(wc[..., 2:, :], wv[..., 2:, :]), weights), np.nan),
           'ALL': ratio(wc.sum(-2), wv.sum(-2))}
    out['A_norm'] = ratio(out['C_std'] - out['U'], out['M'] - out['U'])
    scal = {}
    for k, name in enumerate(SCALARS):
        scal[f'{name}_C_std'] = np.where(strata_ok[..., 0], (ratio(ss[..., 2:, k], cnt[..., 2:]) * weights).sum(-1), np.nan)
        scal[f'{name}_ALL'] = ratio(ss[..., :, k].sum(-1), cnt.sum(-1))
    total = cnt.sum(-1)
    scal['fraction_U'], scal['fraction_M'] = ratio(cnt[..., 0], total), ratio(cnt[..., 1], total)
    return out, scal


def homopolymer_context(seq, p0, min_run=3):
    """'homopolymer' if the substituted base is in, or adjacent to, a WT A or T run of >= min_run."""
    i = 0
    while i < len(seq):
        j = i
        while j + 1 < len(seq) and seq[j + 1] == seq[i]:
            j += 1
        if seq[i] in 'AT' and j - i + 1 >= min_run and i - 1 <= p0 <= j + 1:
            return 'homopolymer'
        i = j + 1
    return 'other'


def window_masks(seq, p0, windows, regions, mask_radius):
    """Unchanged A/T sites only; the substituted base is always excluded."""
    j = np.arange(211)
    d = np.abs(j - p0)
    unchanged = np.array([b in 'AT' for b in seq]) & (j != p0)
    names, rows = [], []
    for name, (lo, hi) in windows.items():
        names.append(name)
        rows.append(unchanged & (d >= lo) & (d <= hi))
    names.append('element')
    rows.append(unchanged & (d > mask_radius))
    for name, (a, b) in regions.items():
        names.append(name)
        rows.append(unchanged & (d > mask_radius) & (j >= a) & (j < b))
    return names, np.array(rows, dtype=float)


def bootstrap(called, valid, scal, groups, clusters, masks, weights, n_groups, b, rng, **kw):
    """Barcode-cluster Poisson bootstrap of one cohort; returns (B, ...) measures."""
    _, inv = np.unique(clusters, return_inverse=True)
    w = rng.poisson(1., (b, inv.max() + 1)).astype(np.float64)[:, inv]
    cs, vs = np.zeros((b, n_groups, 211)), np.zeros((b, n_groups, 211))
    ss, cnt = np.zeros((b, n_groups, len(SCALARS))), np.zeros((b, n_groups))
    for g in range(n_groups):
        sel = groups == g
        if sel.any():
            wg = w[:, sel]
            cs[:, g], vs[:, g] = wg @ called[sel], wg @ valid[sel]
            ss[:, g], cnt[:, g] = wg @ scal[sel], wg.sum(1)
    return measures(cs, vs, ss, cnt, masks, weights, **kw)


def bh(p):
    """Benjamini-Hochberg q-values (SciPy-version independent)."""
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.empty_like(p)
    q[order] = np.minimum(1., np.minimum.accumulate(ranked[::-1])[::-1])
    return q


def exact_ci(k, n):
    lo = stats.beta.ppf(.025, k, n - k + 1) if k > 0 else 0.
    hi = stats.beta.ppf(.975, k + 1, n - k) if k < n else 1.
    return float(lo), float(hi)


def contrast(point, ref, boot, ref_ratio=0.):
    """Delta, SE and z. The reference is not resampled; ref_ratio = n_variant / n_reference
    adds its variance assuming similar per-molecule variance (about 0.3% for real WT)."""
    delta = point - ref
    with np.errstate(invalid='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        finite = np.isfinite(boot).mean(0)
        se = np.where(finite >= .9, np.nanstd(boot, axis=0, ddof=1), np.nan) if len(boot) > 1 else np.full(np.shape(point), np.nan)
        se = se * np.sqrt(1 + ref_ratio)
        z = np.where(se > 0, delta / np.where(se > 0, se, 1), np.nan)
    return delta, se, z


def load_activity(paths):
    act = {}
    for name, path in paths.items():
        delim = ',' if str(path).endswith('.csv') else '\t'
        with open(path) as h:
            for r in csv.DictReader(h, delimiter=delim):
                if r['Alt'] in ('A', 'C', 'G', 'T'):
                    act.setdefault(name, {})[(int(r['Position']), r['Ref'], r['Alt'])] = float(r['Value'])
    return act


def run(config, extract, output, u_max=None, m_min=None, draws=None, cohort=None):
    cfg = read_json(config)
    spec = cfg['analysis']
    u_max = spec['unmethylated_max_distal_m6a'] if u_max is None else u_max
    m_min = spec['methylated_min_distal_m6a'] if m_min is None else m_min
    draws = spec['aa_draws'] if draws is None else draws
    spec = dict(spec, cohort=cohort or spec['cohort'])
    if not 0 <= u_max < m_min <= 1:
        raise ValueError('Require 0 <= u_max < m_min <= 1')
    extract = Path(extract)
    audit = read_json(extract / 'audit.json')
    if audit.get('schema_version') != 'fsmpra.site_extract.v1':
        raise ValueError('Expected a site extraction directory')
    check_sources([s for s in audit['sources'] if not s['path'].endswith(('.bam', '.bai'))])
    act_paths = {k: Path(v) for k, v in cfg['published_activity'].items()}
    sources = [source_info(Path(config), True), source_info(extract / 'audit.json', True), source_info(extract / 'meta.npz', True)]
    sources += [source_info(p, True) for p in act_paths.values()]
    sources += [source_info(Path(__file__).with_name(n), True) for n in ('site_effects.py', 'library_stats.py', 'io.py')]
    variants = read_json(extract / 'variants.json')
    wt_i = next(i for i, v in enumerate(variants) if v['class'] == 'WT')
    seq = variants[wt_i]['sequence']
    meta = dict(np.load(extract / 'meta.npz'))
    cohort_key = COHORTS[spec['cohort']]
    cohort = meta[cohort_key].astype(bool) if cohort_key else np.ones(len(meta['variant_index']), bool)
    # Unevaluated FiberHMM tracks are missing, not zero coverage: never analysed.
    cohort &= meta.get('evaluated', np.ones(len(cohort), bool)).astype(bool)
    vi_all = meta['variant_index'].astype(int)
    groups_all, weights, edges = assign_groups(meta['distal_m6a'].astype(float), (vi_all == wt_i) & cohort, u_max, m_min,
                                               spec['labeling_strata'])
    n_groups = 2 + len(weights)
    keep = np.flatnonzero(cohort & (groups_all >= 0))
    print(f'Cohort {cohort.sum():,} molecules; classified {len(keep):,}', flush=True)

    # Dense per-molecule outcomes for the analysed cohort (about 0.5 GB per million molecules).
    arrays = {n: np.load(extract / f'{n}.npy', mmap_mode='r') for n in ('scores', 'known', 'qualities', 'coverage')}
    called = np.zeros((len(keep), 211), dtype=bool)
    valid = np.zeros((len(keep), 211), dtype=bool)
    scal = np.zeros((len(keep), len(SCALARS)), dtype=np.float32)
    for a in range(0, len(keep), 200000):
        rows = keep[a:a + 200000]
        v = arrays['known'][rows] & (arrays['qualities'][rows] >= spec['baseq_min'])
        valid[a:a + len(rows)] = v
        called[a:a + len(rows)] = v & (arrays['scores'][rows] >= spec['ml_min'])
        cov = arrays['coverage'][rows]
        scal[a:a + len(rows)] = np.stack([((cov & bit) > 0).mean(1) for bit in (2, 4, 1)], 1)
    vi, groups, clusters = vi_all[keep], groups_all[keep].astype(int), meta['cluster'][keep]
    keys = vi * n_groups + groups
    nv = len(variants)
    cs = grouped_sum(keys, called, nv * n_groups).reshape(nv, n_groups, 211)
    vs = grouped_sum(keys, valid, nv * n_groups).reshape(nv, n_groups, 211)
    ss = grouped_sum(keys, scal, nv * n_groups).reshape(nv, n_groups, len(SCALARS))
    cnt = np.bincount(keys, minlength=nv * n_groups).reshape(nv, n_groups).astype(float)
    kw = dict(min_class=spec['min_class_molecules'], min_stratum=spec['min_stratum_molecules'])
    windows = {k: tuple(v) for k, v in spec['windows'].items()}
    regions = {k: tuple(v) for k, v in spec['regions0'].items()}
    rng = np.random.default_rng(spec['seed'])
    by_variant = {i: np.flatnonzero(vi == i) for i in range(nv)}
    act = load_activity(act_paths)

    # Per-SNV contrasts against the full WT reference.
    effects, per_variant = [], {}
    for i, v in enumerate(variants):
        if v['class'] != 'SNV':
            continue
        p0 = v['oligo_pos1'] - 21
        if seq[p0] != v['ref'] or v['sequence'][p0] != v['alt']:
            raise ValueError('Catalog allele disagreement for ' + v['source_variant'])
        names, masks = window_masks(seq, p0, windows, regions, spec['mask_radius'])
        ref_site, ref_scal = measures(cs[wt_i], vs[wt_i], ss[wt_i], cnt[wt_i], masks, weights, **kw)
        pt_site, pt_scal = measures(cs[i], vs[i], ss[i], cnt[i], masks, weights, **kw)
        idx = by_variant[i]
        if len(idx) >= 2 * spec['min_class_molecules']:
            bt_site, bt_scal = bootstrap(called[idx], valid[idx], scal[idx], groups[idx], clusters[idx], masks, weights,
                                         n_groups, spec['bootstrap'], rng, **kw)
        else:
            bt_site = {k: np.full((1,) + np.shape(x), np.nan) for k, x in pt_site.items()}
            bt_scal = {k: np.full((1,), np.nan) for k in pt_scal}
        base = dict(source_variant=v['source_variant'], variant_id=v['variant_id'], oligo_pos1=v['oligo_pos1'],
                    genomic_pos1=GENOMIC_OFFSET + v['oligo_pos1'], ref=v['ref'], alt=v['alt'],
                    at_change=('retained_count' if (v['ref'] in 'AT') == (v['alt'] in 'AT') else 'loss' if v['ref'] in 'AT' else 'gain'),
                    context=homopolymer_context(seq, p0), molecules=int(cnt[i].sum()), molecules_U=int(cnt[i, 0]),
                    molecules_M=int(cnt[i, 1]), molecules_C=int(cnt[i, 2:].sum()))
        for src in act:
            base[f'activity_{src}'] = act[src].get((GENOMIC_OFFSET + v['oligo_pos1'], v['ref'], v['alt']))
        per_variant[i] = base
        for m in SITE_MEASURES:
            delta, se, z = contrast(pt_site[m], ref_site[m], bt_site[m], len(idx) / len(by_variant[wt_i]))
            for w, name in enumerate(names):
                effects.append(dict(base, endpoint=f'{name}:{m}', window=name, measure=m, sites=int(masks[w].sum()),
                                    variant_value=pt_site[m][w], wt_value=ref_site[m][w], delta=delta[w], se=se[w], z=z[w]))
        for m in pt_scal:
            delta, se, z = contrast(pt_scal[m], ref_scal[m], bt_scal[m], len(idx) / len(by_variant[wt_i]))
            effects.append(dict(base, endpoint=m, window='element', measure=m, sites=211, variant_value=float(pt_scal[m]),
                                wt_value=float(ref_scal[m]), delta=float(delta), se=float(se), z=float(z)))
    for e in effects:
        e['p'] = float(2 * stats.norm.sf(abs(e['z']))) if np.isfinite(e['z']) else np.nan
    endpoints = sorted({e['endpoint'] for e in effects})
    for ep in endpoints:
        rows = [e for e in effects if e['endpoint'] == ep and np.isfinite(e['p'])]
        if rows:
            q = bh([e['p'] for e in rows])
            for e, value in zip(rows, q):
                e['q_bh_within_endpoint'] = float(value)
    # One prespecified primary family across endpoints; within-endpoint q-values alone would let
    # ~5% of the many endpoint families show a false discovery under a global null.
    primary = [e for e in effects if e['endpoint'] in spec['primary_endpoints'] and np.isfinite(e['p'])]
    if primary:
        for e, value in zip(primary, bh([e['p'] for e in primary])):
            e['q_bh_primary_family'] = float(value)
    print(f'Variant contrasts complete: {len(per_variant)} SNVs', flush=True)

    # WT pseudo-variants: whole barcode clusters, matched to the median SNV cohort size.
    wt_idx = by_variant[wt_i]
    target = int(np.median([len(by_variant[i]) for i in per_variant]))
    wt_clusters, wt_inv = np.unique(clusters[wt_idx], return_inverse=True)
    members = np.split(np.argsort(wt_inv, kind='stable'), np.cumsum(np.bincount(wt_inv))[:-1])
    positions = sorted({v['oligo_pos1'] - 21 for v in variants if v['class'] == 'SNV'})
    calibration = []
    for draw in range(draws):
        chosen, total = [], 0
        for c in rng.permutation(len(wt_clusters)):
            chosen.append(members[c])
            total += len(members[c])
            if total >= target:
                break
        sel = wt_idx[np.concatenate(chosen)]
        p0 = int(rng.choice(positions))
        names, masks = window_masks(seq, p0, windows, regions, spec['mask_radius'])
        k = groups[sel]
        pcs, pvs = grouped_sum(k, called[sel], n_groups), grouped_sum(k, valid[sel], n_groups)
        pss, pcnt = grouped_sum(k, scal[sel], n_groups), np.bincount(k, minlength=n_groups).astype(float)
        ref_site, ref_scal = measures(cs[wt_i] - pcs, vs[wt_i] - pvs, ss[wt_i] - pss, cnt[wt_i] - pcnt, masks, weights, **kw)
        pt_site, pt_scal = measures(pcs, pvs, pss, pcnt, masks, weights, **kw)
        bt_site, bt_scal = bootstrap(called[sel], valid[sel], scal[sel], k, clusters[sel], masks, weights,
                                     n_groups, spec['bootstrap'], rng, **kw)
        for m in SITE_MEASURES:
            delta, se, z = contrast(pt_site[m], ref_site[m], bt_site[m], len(sel) / (len(wt_idx) - len(sel)))
            for w, name in enumerate(names):
                calibration.append(dict(draw=draw, pseudo_oligo_pos1=p0 + 21, molecules=len(sel), endpoint=f'{name}:{m}', delta=delta[w], se=se[w], z=z[w]))
        for m in pt_scal:
            delta, se, z = contrast(pt_scal[m], ref_scal[m], bt_scal[m], len(sel) / (len(wt_idx) - len(sel)))
            calibration.append(dict(draw=draw, pseudo_oligo_pos1=p0 + 21, molecules=len(sel), endpoint=m,
                                    delta=float(delta), se=float(se), z=float(z)))
    calibration_summary = []
    for ep in endpoints:
        z = np.array([c['z'] for c in calibration if c['endpoint'] == ep], dtype=float)
        z = z[np.isfinite(z)]
        if len(z):
            lo, hi = exact_ci(int((np.abs(z) > 1.96).sum()), len(z))
            calibration_summary.append(dict(endpoint=ep, draws=len(z), mean_z=z.mean(), sd_z=z.std(ddof=1) if len(z) > 1 else np.nan,
                                            rate_abs_z_gt_1_96=float((np.abs(z) > 1.96).mean()), rate_ci_low=lo, rate_ci_high=hi,
                                            rate_lower_tail=float((z < -1.96).mean()), rate_upper_tail=float((z > 1.96).mean()),
                                            gate_pass=bool(len(z) > 1 and .85 <= z.std(ddof=1) <= 1.2 and hi <= .10)))
    print(f'WT pseudo-variant calibration complete: {draws} draws', flush=True)

    # Library kernel: per-site change by offset from the substitution, per class.
    eye = np.eye(211)
    site_all, _ = measures(cs, vs, ss, cnt, eye, weights, **kw)        # {m: (V, 211)}
    kernel_rows, kmax = [], spec['kernel_max_offset']
    for m in ('U', 'M', 'C_std', 'ALL', 'A_norm'):
        deltas = {}
        for i, base in per_variant.items():
            p0 = base['oligo_pos1'] - 21
            for j in range(max(0, p0 - kmax), min(211, p0 + kmax + 1)):
                if j != p0 and seq[j] in 'AT':
                    x = site_all[m][i, j] - site_all[m][wt_i, j]
                    if np.isfinite(x):
                        for cls in (base['at_change'], 'all'):
                            for ctx in (base['context'], 'all'):
                                deltas.setdefault((cls, ctx, j - p0), []).append(x)
        for (cls, ctx, d), xs in sorted(deltas.items()):
            xs = np.array(xs)
            kernel_rows.append(dict(measure=m, at_change=cls, context=ctx, offset=d, pairs=len(xs), mean_delta=xs.mean(),
                                    se=xs.std(ddof=1) / np.sqrt(len(xs)) if len(xs) > 1 else np.nan, median_delta=np.median(xs)))

    # Library association with published activity and observability regression.
    association = []
    by_ep = {}
    for e in effects:
        by_ep.setdefault(e['endpoint'], []).append(e)
    # One SNV set for every endpoint: tiny cohorts would otherwise enter only the unstandardized endpoints.
    common_set = {b['source_variant'] for b in per_variant.values() if b['molecules'] >= spec.get('min_association_molecules', 300)}
    for ep, rows in by_ep.items():
        rows = [r for r in rows if r['source_variant'] in common_set]
        for src in act:
            x = np.array([r['delta'] for r in rows], float)
            y = np.array([np.nan if r.get(f'activity_{src}') is None else r[f'activity_{src}'] for r in rows], float)
            pos = np.array([r['oligo_pos1'] for r in rows])
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() >= 50 and len(np.unique(pos[ok])) > 40 and np.std(x[ok]) > 0:
                association.append(dict(endpoint=ep, activity_source=src, variants=int(ok.sum()),
                                        variant_pearson=float(np.corrcoef(x[ok], y[ok])[0, 1]), **circular_shift(pos[ok], x[ok], y[ok])))
    observability = []
    for w in list(windows) + ['element'] + list(regions):
        xm = {e['source_variant']: e['delta'] for e in by_ep.get(f'{w}:M', [])}
        xc = {e['source_variant']: e['delta'] for e in by_ep.get(f'{w}:C_std', [])}
        common = [k for k in xm if k in xc and np.isfinite(xm[k]) and np.isfinite(xc[k])]
        if len(common) > 10 and np.std([xm[k] for k in common]) > 0:
            a, b = np.array([xm[k] for k in common]), np.array([xc[k] for k in common])
            fit = stats.linregress(a, b)
            observability.append(dict(window=w, variants=len(common), pearson_M_vs_C=fit.rvalue, slope_C_on_M=fit.slope,
                                      sd_delta_M=a.std(ddof=1), sd_delta_C_std=b.std(ddof=1)))

    # Distal methylation distribution (for choosing/confirming proxy thresholds).
    bins = np.linspace(0, 1, 101)
    dm = meta['distal_m6a'].astype(float)
    wt_h = np.histogram(dm[cohort & (vi_all == wt_i) & np.isfinite(dm)], bins)[0]
    snv_h = np.histogram(dm[cohort & (vi_all != wt_i) & np.isfinite(dm)], bins)[0]
    hist = [dict(bin_low=float(bins[k]), bin_high=float(bins[k + 1]), wt_molecules=int(wt_h[k]), snv_molecules=int(snv_h[k]))
            for k in range(100)]
    class_rows = [dict(source_variant=b['source_variant'], molecules=b['molecules'], molecules_U=b['molecules_U'],
                       molecules_M=b['molecules_M'], molecules_C=b['molecules_C']) for b in per_variant.values()]
    n_mol = np.array([r['molecules'] for r in class_rows], float)
    dispersion = {}
    for name, col, k in (('U', 'molecules_U', 0), ('M', 'molecules_M', 1)):
        p_wt = cnt[wt_i, k] / cnt[wt_i].sum()
        f = np.array([r[col] for r in class_rows]) / np.where(n_mol > 0, n_mol, 1)
        z = (f - p_wt) / np.sqrt(p_wt * (1 - p_wt) / np.where(n_mol > 0, n_mol, 1))
        dispersion[f'fraction_{name}_pearson_dispersion'] = float(np.mean(z[n_mol > 0] ** 2))
    dispersion['interpretation'] = ('Mean squared standardized deviation of each SNV class fraction from WT (Pearson chi-square/df). '
                                    '~1 means binomial sampling only; well above 1 indicates barcode clustering or genotype-dependent '
                                    'class membership. Robust to SNVs with very few molecules, unlike a variance ratio.')

    def clean(rows):
        return [{k: (None if isinstance(v, (float, np.floating)) and not np.isfinite(v) else
                     float(v) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()} for r in rows]
    fields = list(effects[0]) + ['p', 'q_bh_within_endpoint', 'q_bh_primary_family']
    fields = list(dict.fromkeys(fields))
    check_sources(sources)
    with atomic_output(output) as out:
        write_table(out / 'effects.tsv.gz', clean(effects), fields)
        write_table(out / 'calibration_draws.tsv.gz', clean(calibration), list(calibration[0]))
        # Sparse smoke runs can leave these empty; an absent file means no estimable rows.
        if calibration_summary:
            write_table(out / 'calibration_summary.tsv', clean(calibration_summary), list(calibration_summary[0]))
        if kernel_rows:
            write_table(out / 'kernel.tsv', clean(kernel_rows), list(kernel_rows[0]))
        if association:
            write_table(out / 'activity_association.tsv', clean(association), list(association[0]))
        if observability:
            write_table(out / 'observability_regression.tsv', clean(observability), list(observability[0]))
        write_table(out / 'distal_m6a_histogram.tsv', clean(hist), list(hist[0]))
        write_table(out / 'class_counts.tsv', class_rows, list(class_rows[0]))
        write_json(out / 'summary.json', dict(
            schema_version='fsmpra.site_effects.v1', sources=sources, cohort=spec['cohort'],
            thresholds=dict(unmethylated_max_distal_m6a=u_max, methylated_min_distal_m6a=m_min,
                            labeling_strata_edges=edges.tolist(), labeling_strata_wt_weights=weights.tolist()),
            molecules=dict(cohort=int(cohort.sum()), classified=int(len(keep)),
                           wt=dict(U=int(cnt[wt_i, 0]), M=int(cnt[wt_i, 1]), C=int(cnt[wt_i, 2:].sum())),
                           snv=dict(U=int(cnt[:, 0].sum() - cnt[wt_i, 0]), M=int(cnt[:, 1].sum() - cnt[wt_i, 1]),
                                    C=int(cnt[:, 2:].sum() - cnt[wt_i, 2:].sum()))),
            class_dispersion=dispersion, bootstrap=spec['bootstrap'], aa_draws=draws, aa_target_molecules=target,
            primary_endpoints=spec['primary_endpoints'], windows=windows, regions0=regions, mask_radius=spec['mask_radius'], baseq_min=spec['baseq_min'], ml_min=spec['ml_min'],
            interpretation='U/M are same-reaction proxies, not prepared controls. z/p/q are only as good as the WT '
                           'pseudo-variant calibration in calibration_summary.tsv; check it per endpoint before use.'))
        if spec.get('figures', True) and kernel_rows:
            figures(out, kernel_rows, calibration, effects, hist, u_max, m_min)
    print('Complete', flush=True)


def figures(out, kernel_rows, calibration, effects, hist, u_max, m_min):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = {'U': '#7a7a7a', 'M': '#c0504d', 'C_std': '#2f6db3', 'A_norm': '#3c8d4f'}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
    for ax, cls in zip(axes, ('gain', 'loss', 'retained_count')):
        for m in ('U', 'M', 'C_std'):
            rows = [r for r in kernel_rows if r['measure'] == m and r['at_change'] == cls and r['context'] == 'all']
            d = np.array([r['offset'] for r in rows])
            y = 100 * np.array([r['mean_delta'] for r in rows])
            e = 100 * np.array([r['se'] for r in rows], dtype=float)
            ax.plot(d, y, color=colors[m], lw=1.6, label=m)
            ax.fill_between(d, y - 2 * e, y + 2 * e, color=colors[m], alpha=.15, lw=0)
        ax.axhline(0, color='k', lw=.6)
        ax.axvline(0, color='k', lw=.6, ls=':')
        ax.set(title=f'A/T {cls.replace("_", " ")}', xlabel='Offset of unchanged A/T site from substitution (bp)')
        ax.grid(alpha=.15)
    axes[0].set_ylabel('Mean change in retained m6A vs WT (pp)')
    axes[0].legend(title='Molecule class')
    fig.suptitle('Library kernel: how far substitutions change per-site calls in each proxy class (±2 SE)')
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(out / f'kernel.{ext}', dpi=150)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.array([h['bin_low'] for h in hist]) + .005
    for key, color in (('wt_molecules', '#273f4b'), ('snv_molecules', '#cc7b32')):
        y = np.array([h[key] for h in hist], float)
        ax.plot(x, y / y.sum(), color=color, label=key.replace('_molecules', ''))
    for t in (u_max, m_min):
        ax.axvline(t, color='k', ls='--', lw=1)
    ax.set(xlabel='Distal retained m6A fraction', ylabel='Fraction of molecules', yscale='log',
           title='Proxy-control thresholds on distal methylation')
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / 'distal_m6a_histogram.png', dpi=150)
    plt.close(fig)
    eps = ['local:C_std', 'local:M', 'local:A_norm', 'near:M', 'msp_C_std', 'nuc_85_250_C_std']
    fig, axes = plt.subplots(1, len(eps), figsize=(3 * len(eps), 3.2), sharey=True)
    grid = np.linspace(-5, 5, 41)
    for ax, ep in zip(axes, eps):
        z = np.array([c['z'] for c in calibration if c['endpoint'] == ep], float)
        z = z[np.isfinite(z)]
        if len(z):
            ax.hist(z, bins=grid, density=True, color='#8aa5c2')
        ax.plot(grid, stats.norm.pdf(grid), color='k', lw=1)
        ax.set(title=f'{ep}\nSD {z.std(ddof=1):.2f}' if len(z) > 1 else ep, xlabel='WT pseudo-variant z')
    fig.tight_layout()
    fig.savefig(out / 'calibration_z.png', dpi=150)
    plt.close(fig)
    xm = {e['source_variant']: e['delta'] for e in effects if e['endpoint'] == 'local:M'}
    xc = {e['source_variant']: e['delta'] for e in effects if e['endpoint'] == 'local:C_std'}
    keys = [k for k in xm if k in xc and np.isfinite(xm[k]) and np.isfinite(xc[k])]
    if keys:
        fig, ax = plt.subplots(figsize=(5.5, 5))
        ax.scatter([100 * xm[k] for k in keys], [100 * xc[k] for k in keys], s=6, color='#9aa7b1')
        for k in keys:
            if k.startswith('125'):
                ax.scatter(100 * xm[k], 100 * xc[k], s=30, color='#c0504d')
                ax.annotate(k, (100 * xm[k], 100 * xc[k]), fontsize=8)
        ax.axhline(0, color='k', lw=.5)
        ax.axvline(0, color='k', lw=.5)
        ax.set(xlabel='Local Δ m6A in M (heavily labeled proxy), pp', ylabel='Local Δ m6A in C (standardized), pp',
               title='Observability (x) vs in-cell change (y), sites 3-10 bp from substitution')
        fig.tight_layout()
        fig.savefig(out / 'observability_scatter.png', dpi=150)
        plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--config', default='config/ldlr_site_methylation.json')
    p.add_argument('--extract', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--u-max', type=float, help='Override unmethylated distal m6A ceiling (sensitivity runs)')
    p.add_argument('--m-min', type=float, help='Override methylated distal m6A floor (sensitivity runs)')
    p.add_argument('--aa-draws', type=int)
    p.add_argument('--cohort', choices=sorted(COHORTS), help='Override the analysed cohort (default from config)')
    a = p.parse_args()
    run(a.config, a.extract, a.output, a.u_max, a.m_min, a.aa_draws, a.cohort)


if __name__ == '__main__':
    main()

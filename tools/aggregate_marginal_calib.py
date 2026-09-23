#!/usr/bin/env python3
"""
aggregate_marginal_calib.py — verdict for the marginal-null calibration.

Pools the per-seed JSONs written by marginal_permutation_test.py (one per
permutation seed) and produces the calibration verdict.

Each permutation seed is one independent Monte-Carlo replicate, so the
across-seed mean and standard error of a metric are honest even though
the ~900 per-variant p-values WITHIN one seed are correlated (shared
reads across overlapping variant pools). This is why the verdict uses
the seed-level spread, not a within-seed KS test.

Metrics
-------
* frac(p<=0.05): expected 0.05 under H0. Reported overall and per
  variant-read-count quartile. A rise with read count is the
  self-contamination signature (a variant's reads are in the all-reads
  null pool; contamination ~ N_variant / N_pool).
* realized FDR: under the global null every BH discovery is false, so
  mean n_significant_bh should be ~0. This is the metric that decides
  whether marginal-mode calls are trustworthy at FDR 10%.

Verdict
-------
  PASS            — overall + all quartiles within ~2 SE of 0.05,
                    mean BH discoveries ~0.
  PASS_WITH_CAVEAT— low/mid quartiles calibrated and BH discoveries ~0,
                    but the top read-count quartile is elevated (mild
                    self-contamination). Marginal calls are usable at
                    FDR 10%; the leave-one-variant-out null (plan
                    Phase 2) removes the residual bias.
  FAIL            — overall exceedance materially > 0.05 or BH discovers
                    real false positives under the global null.
"""

import argparse
import glob
import json
import os

import numpy as np


def mean_se(values):
    a = np.asarray(values, dtype=float)
    n = len(a)
    if n == 0:
        return float('nan'), float('nan')
    se = a.std(ddof=1) / np.sqrt(n) if n > 1 else float('nan')
    return float(a.mean()), float(se)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--calib-dir', required=True,
                   help='Directory containing seed_*/marginal_calib.*.json')
    p.add_argument('--output', default=None,
                   help='Aggregate JSON path (default: <calib-dir>/'
                        'marginal_calib.aggregate.json)')
    p.add_argument('--alpha', type=float, default=0.05)
    args = p.parse_args()

    paths = sorted(glob.glob(os.path.join(
        args.calib_dir, 'seed_*', 'marginal_calib.*.json')))
    if not paths:
        raise SystemExit(f"No per-seed JSONs under {args.calib_dir}")

    seeds = []
    for path in paths:
        with open(path) as fh:
            seeds.append(json.load(fh))

    n_seeds = len(seeds)
    fracs = [s['frac_p_le_alpha'] for s in seeds]
    n_sigs = [s['n_significant_bh'] for s in seeds]
    n_variants = int(np.median([s['n_variants'] for s in seeds]))

    frac_mean, frac_se = mean_se(fracs)
    nsig_mean, nsig_se = mean_se(n_sigs)

    # Per-quartile: align quartiles across seeds by index (each seed
    # splits its own variants into 4 equal read-count bins).
    n_q = max(len(s.get('by_n_quartile', [])) for s in seeds)
    quartiles = []
    for qi in range(n_q):
        qf = [s['by_n_quartile'][qi]['frac_p_le_alpha']
              for s in seeds if len(s.get('by_n_quartile', [])) > qi]
        rngs = [s['by_n_quartile'][qi]['n_reads_range']
                for s in seeds if len(s.get('by_n_quartile', [])) > qi]
        m, se = mean_se(qf)
        quartiles.append({
            'quartile': qi + 1,
            'n_reads_range_typical': rngs[len(rngs) // 2] if rngs else None,
            'frac_p_le_alpha_mean': m,
            'frac_p_le_alpha_se': se,
        })

    # Verdict. ~2 SE band around alpha; BH discoveries must be ~0.
    def elevated(mean, se):
        if not np.isfinite(se) or se == 0:
            return mean > args.alpha + 0.01
        return (mean - args.alpha) > 2 * se

    overall_elevated = elevated(frac_mean, frac_se)
    top_elevated = (len(quartiles) > 0
                    and elevated(quartiles[-1]['frac_p_le_alpha_mean'],
                                 quartiles[-1]['frac_p_le_alpha_se']))
    lowmid_ok = all(
        not elevated(q['frac_p_le_alpha_mean'], q['frac_p_le_alpha_se'])
        for q in quartiles[:-1]) if len(quartiles) > 1 else True
    fdr_ok = nsig_mean <= max(1.0, 0.01 * n_variants)

    if not fdr_ok or (overall_elevated and not lowmid_ok):
        verdict = 'FAIL'
    elif top_elevated and fdr_ok:
        verdict = 'PASS_WITH_CAVEAT'
    elif not overall_elevated and fdr_ok:
        verdict = 'PASS'
    else:
        verdict = 'PASS_WITH_CAVEAT'

    agg = {
        'n_seeds': n_seeds,
        'n_variants_typical': n_variants,
        'alpha': args.alpha,
        'frac_p_le_alpha_mean': frac_mean,
        'frac_p_le_alpha_se': frac_se,
        'n_significant_bh_mean': nsig_mean,
        'n_significant_bh_se': nsig_se,
        'n_significant_bh_per_seed': n_sigs,
        'by_n_quartile': quartiles,
        'verdict': verdict,
        'seed_json_files': paths,
    }

    out = args.output or os.path.join(args.calib_dir,
                                      'marginal_calib.aggregate.json')
    with open(out, 'w') as fh:
        json.dump(agg, fh, indent=2)

    print(f"=== Marginal-null calibration — {n_seeds} seeds, "
          f"~{n_variants} variants/seed ===")
    print(f"  overall frac(p<=0.05) : {frac_mean:.4f} +/- {frac_se:.4f} "
          f"(expected {args.alpha})")
    print(f"  BH discoveries/seed   : {nsig_mean:.2f} +/- {nsig_se:.2f} "
          f"(expected ~0 under global null)")
    print("  frac(p<=0.05) by variant read-count quartile:")
    for q in quartiles:
        print(f"    Q{q['quartile']} N~{q['n_reads_range_typical']}: "
              f"{q['frac_p_le_alpha_mean']:.4f} +/- "
              f"{q['frac_p_le_alpha_se']:.4f}")
    print(f"\n  VERDICT: {verdict}")
    print(f"  JSON: {out}")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
wt_vs_wt_calibration.py — negative-control calibration harness.

Splits WT into a reference pool (used to build the null + WT mean) and a
disjoint pseudo-variant pool. Draws K pseudo-"variants" from the pseudo
pool (each a sample of N WT reads) and runs them through the real Stage 3
testing path. Because pseudo-variants are WT, under a correctly calibrated
test:
  * best_cluster_p should be ~Uniform(0,1)
  * P(best_cluster_p <= alpha) should be ~= alpha
  * cross-variant BH should call ~0 variants significant

This is the negative control for Issues 1-3 AND the deterministic
regression oracle for the de-duplication refactor (fixed --seed
--threads 1 -> identical best_cluster_p list before/after).

Read-only w.r.t. the pipeline; writes a small JSON snapshot for diffing.
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..',
                                'workflow', 'scripts'))
from common import (  # noqa: E402
    build_analysis_bins, parse_region, get_bam_ref_info, resolve_target_chrom,
    parse_bam, group_variants, compute_ground_truth_nc,
    run_null_calibration, run_variant_testing_parallel,
    benjamini_hochberg, setup_logging,
)


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bam', required=True)
    p.add_argument('--target-region', required=True)
    p.add_argument('--promoter-region', default=None)
    p.add_argument('--n-pseudo', type=int, default=300,
                   help='# pseudo-variants drawn from the WT pseudo pool')
    p.add_argument('--pseudo-n', type=int, default=None,
                   help='Reads per pseudo-variant (default: median '
                        'testable real-variant N)')
    p.add_argument('--n-null-iterations', type=int, default=2000)
    p.add_argument('--ref-fraction', type=float, default=0.5,
                   help='Fraction of WT used as the null/reference pool')
    p.add_argument('--threads', type=int, default=1,
                   help='1 = deterministic serial path (regression oracle)')
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--stratify', action='store_true', default=False,
                   help='Reuse nulls across similar pseudo-variants '
                        '(default: off = independent per-variant null)')
    p.add_argument('--null-strata-n-tol', type=float, default=0.10)
    p.add_argument('--null-strata-nc-dist', type=float, default=0.30)
    p.add_argument('--out', default='results/phase0/wt_vs_wt_calib.json')
    args = p.parse_args()
    setup_logging(quiet=True)

    analysis_bins = build_analysis_bins()
    bin_labels = [b[0] for b in analysis_bins]

    target_chrom, t_start, t_end = parse_region(args.target_region)
    ref_info = get_bam_ref_info(args.bam)
    ref_name, ref_length = resolve_target_chrom(ref_info, target_chrom)
    if args.promoter_region:
        ps = args.promoter_region.split('-')
        prom_s, prom_e = int(ps[0]), int(ps[1]) + 1
    elif t_start is not None:
        prom_s, prom_e = t_start, t_end
    else:
        prom_s, prom_e = 0, ref_length
    analysis_region = (prom_s, prom_e)

    rd, ref_length, ref_name, _ = parse_bam(
        args.bam, analysis_bins, target_chrom=target_chrom,
        analysis_region=analysis_region)

    # Target N from the real testable-variant distribution
    real_groups = group_variants(rd, include_multi=False, min_reads=50)
    real_Ns = sorted(len(v) for v in real_groups.values())
    pseudo_n = args.pseudo_n or (int(np.median(real_Ns)) if real_Ns else 130)

    # Split WT into disjoint reference and pseudo pools (seeded)
    rng = np.random.default_rng(args.seed)
    wt_all = rd.wt_indices.copy()
    rng.shuffle(wt_all)
    n_ref = int(len(wt_all) * args.ref_fraction)
    wt_ref = np.sort(wt_all[:n_ref])
    wt_pseudo = wt_all[n_ref:]

    print(f"\n=== WT-vs-WT calibration ===")
    print(f"WT total={len(wt_all):,}  ref={len(wt_ref):,}  "
          f"pseudo_pool={len(wt_pseudo):,}")
    print(f"pseudo-variants: K={args.n_pseudo}  N={pseudo_n}  "
          f"B={args.n_null_iterations}  threads={args.threads}  "
          f"seed={args.seed}")

    # K pseudo-variants: each N draws from the pseudo pool (with
    # replacement across variants, without within) — each is an
    # independent genuine null sample vs the reference-WT null.
    pseudo_groups = {}
    for k in range(args.n_pseudo):
        sub = rng.choice(wt_pseudo, size=pseudo_n, replace=False)
        pseudo_groups[f"pseudo_{k:04d}"] = np.sort(sub)

    # Per-variant nulls (Issue 3 fix): each pseudo-variant builds its
    # own null from the reference WT pool at its exact N, NC-matched to
    # itself, Option-A reference. Built inside run_variant_testing_parallel.
    results = run_variant_testing_parallel(
        rd, wt_ref, pseudo_groups, analysis_region,
        prom_s, prom_e, random_seed=args.seed,
        n_null_iterations=args.n_null_iterations,
        n_workers=args.threads,
        stratify=args.stratify,
        n_tol=args.null_strata_n_tol,
        nc_dist=args.null_strata_nc_dist)

    best_p = np.array(sorted(r['best_cluster_p'] for r in results),
                      dtype=np.float64)
    n = len(best_p)
    if n == 0:
        print("ERROR: no pseudo-variants returned a result")
        sys.exit(1)

    # P-value uniformity: KS distance to U(0,1) (here KS-vs-uniform IS
    # the right tool — we are testing p-value calibration, not NC)
    ecdf = np.arange(1, n + 1) / n
    ks = float(np.max(np.abs(ecdf - best_p)))
    qs = benjamini_hochberg(best_p)

    print(f"\n=== Calibration result (n={n} pseudo-variants tested) ===")
    print(f"best_cluster_p: min={best_p.min():.4f} "
          f"median={np.median(best_p):.4f} mean={best_p.mean():.4f}")
    print(f"KS distance to Uniform(0,1) = {ks:.4f} "
          f"(smaller = better calibrated)")
    print("Type-I error vs nominal (calibrated => observed ~= nominal):")
    for a in (0.01, 0.05, 0.10, 0.20):
        obs = float(np.mean(best_p <= a))
        flag = "  <-- inflated" if obs > a * 1.5 else ""
        print(f"  alpha={a:<4}  observed P(p<=alpha)={obs:.4f}{flag}")
    for q in (0.05, 0.10, 0.20):
        print(f"  cross-variant BH q<{q}: {int(np.sum(qs < q))} "
              f"(expect ~0 under the global null)")

    snap = {
        'n': n, 'pseudo_n': pseudo_n, 'B': args.n_null_iterations,
        'seed': args.seed, 'threads': args.threads, 'ks_uniform': ks,
        'best_cluster_p_sorted': [round(x, 8) for x in best_p.tolist()],
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as fh:
        json.dump(snap, fh, indent=1)
    print(f"\nsnapshot -> {args.out} "
          f"(diff before/after de-dup; must be identical)")


if __name__ == '__main__':
    main()

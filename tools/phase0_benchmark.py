#!/usr/bin/env python3
"""
phase0_benchmark.py — measured inputs for the Stage 3 redesign sizing.

Reports, on a real tagged BAM:
  1. Testable-variant count and reads-per-variant (N) distribution.
  2. Per-variant nucleosome-count (NC) distribution spread vs WT
     (Wasserstein-1), plus pairwise spread among variants — sets a
     defensible default for --null-strata-nc-dist (null stratification).
  3. Real per-null-iteration wall time (setup vs per-iteration separated
     by a two-point fit), single-core and parallel, extrapolated to
     B = 2k / 10k for per-variant vs stratified nulls.

This replaces the timing estimates in docs/stage3_redesign_plan.md.
Read-only; does not modify the pipeline.
"""

import argparse
import os
import sys
import time

import numpy as np
from scipy.stats import wasserstein_distance

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..',
                                'workflow', 'scripts'))
from common import (  # noqa: E402
    build_analysis_bins, parse_region, get_bam_ref_info, resolve_target_chrom,
    parse_bam, group_variants, compute_ground_truth_nc,
    run_null_calibration, setup_logging,
)


def _variant_nc_samples(rd, indices, min_nuc, max_nuc):
    idx = rd.get_indices_filtered(indices, min_nuc=min_nuc, max_nuc=max_nuc)
    nc = rd.nuc_counts[idx]
    return nc[nc >= 0]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bam', required=True, help='Tagged BAM (Stage 2 output)')
    p.add_argument('--target-region', required=True,
                   help='CHROM:START-END, e.g. LDLR:3184-3501')
    p.add_argument('--promoter-region', default=None,
                   help='START-END within target (default: = target)')
    p.add_argument('--min-reads', type=int, default=50)
    p.add_argument('--time-iters-small', type=int, default=30)
    p.add_argument('--time-iters-large', type=int, default=300)
    p.add_argument('--threads', type=int, default=8,
                   help='Workers for the parallel timing point')
    p.add_argument('--pairwise-sample', type=int, default=200,
                   help='# variants sampled for pairwise NC spread')
    p.add_argument('--seed', type=int, default=42)
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

    print(f"\n=== Parse ===")
    print(f"ref={ref_name} ({ref_length} bp)  "
          f"analysis_region={analysis_region} "
          f"({prom_e - prom_s} bp)")
    t0 = time.time()
    rd, ref_length, ref_name, pstats = parse_bam(
        args.bam, analysis_bins, target_chrom=target_chrom,
        analysis_region=analysis_region)
    print(f"parsed {pstats['parsed']:,} reads in {time.time() - t0:.1f}s")

    wt_idx = rd.wt_indices
    gtn = compute_ground_truth_nc(rd, wt_idx)
    variant_groups = group_variants(rd, include_multi=False,
                                    min_reads=args.min_reads)

    # ---- 1. Coverage distribution ----
    Ns = np.array(sorted(len(v) for v in variant_groups.values()))
    print(f"\n=== 1. Testable variants (min_reads={args.min_reads}) ===")
    print(f"n_variants={len(Ns)}  WT_reads={len(wt_idx):,}")
    if len(Ns):
        for q in (0, 5, 25, 50, 75, 95, 100):
            print(f"  N p{q:>3} = {int(np.percentile(Ns, q)):>6,}")
    median_N = int(np.median(Ns)) if len(Ns) else 200

    # ---- 2. NC-distribution spread ----
    wt_nc = _variant_nc_samples(rd, wt_idx, None, None)
    var_ids = list(variant_groups.keys())
    w_vs_wt = []
    nc_samples = {}
    for vid in var_ids:
        s = _variant_nc_samples(rd, variant_groups[vid], None, None)
        if len(s) >= 2:
            nc_samples[vid] = s
            w_vs_wt.append(wasserstein_distance(s, wt_nc))
    w_vs_wt = np.array(w_vs_wt)
    print(f"\n=== 2. NC spread: Wasserstein-1 variant vs WT "
          f"(WT mean NC={wt_nc.mean():.2f}) ===")
    if len(w_vs_wt):
        for q in (50, 75, 90, 95, 99, 100):
            print(f"  p{q:>3} = {np.percentile(w_vs_wt, q):.3f} nucleosomes")

    sample_ids = var_ids[:args.pairwise_sample]
    sample_ids = [v for v in sample_ids if v in nc_samples]
    pw = []
    for i in range(len(sample_ids)):
        for j in range(i + 1, len(sample_ids)):
            pw.append(wasserstein_distance(nc_samples[sample_ids[i]],
                                           nc_samples[sample_ids[j]]))
    pw = np.array(pw)
    print(f"  pairwise among {len(sample_ids)} variants "
          f"({len(pw):,} pairs):")
    if len(pw):
        for q in (50, 75, 90, 95, 99):
            print(f"    p{q:>3} = {np.percentile(pw, q):.3f}")
        print(f"  -> a --null-strata-nc-dist near the pairwise p50 "
              f"({np.percentile(pw, 50):.2f}) groups ~half of all "
              f"variant pairs into shared strata")

    # ---- 3. Per-iteration timing (two-point fit) ----
    print(f"\n=== 3. Null timing at N={median_N} (median variant) ===")

    def _time_null(n_iter, workers):
        t = time.time()
        run_null_calibration(
            rd, wt_idx, gtn, coverage_level=median_N,
            n_iterations=n_iter, random_seed=args.seed,
            analysis_region=analysis_region, n_workers=workers)
        return time.time() - t

    s, L = args.time_iters_small, args.time_iters_large
    t_s = _time_null(s, 1)
    t_L = _time_null(L, 1)
    per_iter = max((t_L - t_s) / (L - s), 1e-9)
    setup = max(t_s - per_iter * s, 0.0)
    print(f"  single-core: setup~{setup:.2f}s  "
          f"per-iteration~{per_iter * 1000:.2f} ms "
          f"(from {s} iter={t_s:.1f}s, {L} iter={t_L:.1f}s)")
    t_par = _time_null(L, args.threads)
    print(f"  {args.threads} workers: {L} iter in {t_par:.1f}s "
          f"(speedup {t_L / max(t_par, 1e-9):.1f}x)")

    n_var = max(len(Ns), 1)
    print(f"\n=== Extrapolation (per-iteration {per_iter*1000:.2f} ms, "
          f"single-core CPU-time) ===")
    for B in (2000, 10000):
        per_null = setup + per_iter * B
        print(f"  B={B:>5}: 1 null ~{per_null:.1f}s | "
              f"per-variant x{n_var} ~{per_null * n_var / 3600:.2f} CPU-h | "
              f"50 strata ~{per_null * 50 / 60:.1f} CPU-min")
    print()


if __name__ == '__main__':
    main()

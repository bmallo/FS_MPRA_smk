#!/usr/bin/env python3
"""
marginal_power_benchmark.py — power / sensitivity benchmark for
--null-model marginal (validation test 3c).

Question
--------
Calibration (marginal_permutation_test.py) showed marginal mode does not
produce false positives. This asks the complementary question: when a
variant DOES change chromatin, does marginal mode detect it, and how big
must the effect be?

Design
------
1. Start from a PERMUTED background: variant labels are shuffled across
   reads (the H0 state proven null by the calibration test). On this
   background every variant is, by construction, a true null — so any
   detection after we inject an effect is attributable to the injection,
   and un-injected variants are genuine specificity controls.

2. Inject known effects. For a grid of variants spanning a range of read
   counts N, knock out a footprint: for a fraction `f` (penetrance) of
   the variant's reads, zero the chosen coverage bin (default TF) inside
   a randomly placed window. Windows are scattered so the library-mean
   baseline is essentially unperturbed.

3. The injected effect is measured EMPIRICALLY per variant as the change
   in that variant's mean occupancy over its window (original minus
   modified) — exact ground truth, robust to any read-sharing crosstalk.

4. Run the real marginal testing path and compare recovered effect and
   detection (FDR q<0.10) against the injected truth, and against the
   pipeline's own predicted minimum detectable effect (MDE).

Outputs a per-variant TSV and a JSON summary (power curve, recovered-vs-
injected fit, specificity, MDE concordance).
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'workflow', 'scripts'))
from common import (  # noqa: E402
    build_analysis_bins, parse_region, get_bam_ref_info, resolve_target_chrom,
    parse_bam, group_variants, run_variant_testing_parallel,
    benjamini_hochberg, setup_logging,
)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from marginal_permutation_test import shuffle_bam  # noqa: E402

# Fraction of a variant's reads whose footprint is knocked out. The
# injected occupancy delta is ~ penetrance x bin-baseline-occupancy, so
# this range x the nuc bin (~0.46 occupancy) spans ~0.02-0.30 — the
# realistic effect-size range for this assay.
PENETRANCES = [0.05, 0.10, 0.20, 0.35, 0.50, 0.65]


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bam', required=True, help='Tagged BAM (Stage 2 output)')
    p.add_argument('--target-region', required=True)
    p.add_argument('--promoter-region', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--n-spikes', type=int, default=150,
                   help='Number of variants to inject an effect into')
    p.add_argument('--window-width', type=int, default=40,
                   help='Knockout window width in bp (default: 40)')
    p.add_argument('--bin', default='nuc',
                   help='Coverage bin to knock out (default: nuc — the '
                        'only bin with occupancy high enough to inject '
                        'effects across the realistic 0.02-0.30 range)')
    p.add_argument('--n-null-iterations', type=int, default=2000)
    p.add_argument('--min-reads', type=int, default=50)
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--fdr', type=float, default=0.10)
    p.add_argument('--skip-shuffle', action='store_true',
                   help='Treat --bam as already permuted (H0 background)')
    args = p.parse_args()
    setup_logging(quiet=True)
    os.makedirs(args.output_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # --- H0 background: permute variant labels across reads ------------
    if args.skip_shuffle:
        bg_bam = args.bam
    else:
        bg_bam = os.path.join(args.output_dir, 'background.bam')
        n = shuffle_bam(args.bam, bg_bam, args.seed)
        print(f"Permuted background: {n:,} reads relabelled")

    # --- Parse -------------------------------------------------------
    analysis_bins = build_analysis_bins()
    target_chrom, t_start, t_end = parse_region(args.target_region)
    ref_info = get_bam_ref_info(bg_bam)
    ref_name, ref_length = resolve_target_chrom(ref_info, target_chrom)
    ps = args.promoter_region.split('-')
    prom_s, prom_e = int(ps[0]), int(ps[1]) + 1
    analysis_region = (prom_s, prom_e)
    analysis_length = prom_e - prom_s

    rd, ref_length, ref_name, _ = parse_bam(
        bg_bam, analysis_bins, target_chrom=target_chrom,
        analysis_region=analysis_region)

    if args.bin not in rd.coverage_matrices:
        raise SystemExit(f"bin {args.bin} not in {list(rd.coverage_matrices)}")

    # --- Marginal-mode variant pools ---------------------------------
    variant_groups = group_variants(rd, include_multi=True,
                                    min_reads=args.min_reads,
                                    multi_mode='all')
    # Tiebreak by variant ID for fully deterministic spike selection
    # independent of dict insertion order.
    vids = sorted(variant_groups,
                  key=lambda v: (len(variant_groups[v]), v))
    if not vids:
        raise SystemExit("no testable variants")

    # --- Inject known effects ----------------------------------------
    # Spike variants are spread evenly across the read-count range.
    pick = np.unique(np.linspace(0, len(vids) - 1,
                                 min(args.n_spikes, len(vids))).astype(int))
    spike_vids = [vids[i] for i in pick]
    orig_bin = rd.coverage_matrices[args.bin].copy()
    spec = {}  # vid -> (w0, w1, penetrance)
    for i, vid in enumerate(spike_vids):
        f = PENETRANCES[i % len(PENETRANCES)]
        idx = variant_groups[vid]
        n_zero = int(round(f * len(idx)))
        if n_zero == 0:
            continue
        zero_reads = rng.choice(idx, size=n_zero, replace=False)
        w0 = int(rng.integers(0, analysis_length - args.window_width))
        w1 = w0 + args.window_width
        rd.coverage_matrices[args.bin][np.ix_(zero_reads,
                                              np.arange(w0, w1))] = 0
        spec[vid] = (w0, w1, f)
    print(f"Injected effects into {len(spec)} variants "
          f"({len(vids)} total; {len(vids) - len(spec)} null controls)")

    # --- Run the real marginal testing path --------------------------
    null_idx = np.arange(rd.n_reads)
    results = run_variant_testing_parallel(
        rd, null_idx, variant_groups, analysis_region, prom_s, prom_e,
        random_seed=args.seed, n_null_iterations=args.n_null_iterations,
        n_workers=args.threads, stratify=False)

    qmap = dict(zip([r['variant_id'] for r in results],
                    benjamini_hochberg(
                        np.array([r['best_cluster_p'] for r in results]))))

    # --- Per-variant injected vs recovered ---------------------------
    rows = []
    mod_bin = rd.coverage_matrices[args.bin]
    for r in results:
        vid = r['variant_id']
        idx = variant_groups[vid]
        spiked = vid in spec
        if spiked:
            w0, w1, f = spec[vid]
        else:
            w0, w1, f = 0, analysis_length, 0.0
        injected = float(mod_bin[idx][:, w0:w1].mean()
                          - orig_bin[idx][:, w0:w1].mean())
        delta = r[args.bin]['delta_obs']
        recovered = float(np.mean(delta[w0:w1]))
        q = float(qmap[vid])
        rows.append({
            'variant_id': vid, 'spiked': spiked, 'n_reads': len(idx),
            'penetrance': f, 'injected_effect': injected,
            'recovered_effect': recovered,
            'best_cluster_p': float(r['best_cluster_p']), 'q': q,
            'detected': bool(q < args.fdr),
            'mde_median': float(r.get('mde_median', float('nan'))),
        })

    # --- Summarise ---------------------------------------------------
    tsv = os.path.join(args.output_dir, 'power_per_variant.tsv')
    cols = ['variant_id', 'spiked', 'n_reads', 'penetrance',
            'injected_effect', 'recovered_effect', 'best_cluster_p', 'q',
            'detected', 'mde_median']
    with open(tsv, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for row in rows:
            fh.write('\t'.join(str(row[c]) for c in cols) + '\n')

    spiked = [r for r in rows if r['spiked']]
    nulls = [r for r in rows if not r['spiked']]
    inj = np.array([abs(r['injected_effect']) for r in spiked])
    rec = np.array([r['recovered_effect'] for r in spiked])
    det = np.array([r['detected'] for r in spiked])

    # Power curve: detection rate by injected-effect magnitude bin.
    edges = [0.0, 0.02, 0.05, 0.10, 0.15, 0.22, 1.0]
    power_curve = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (inj >= lo) & (inj < hi)
        if sel.sum() == 0:
            continue
        power_curve.append({
            'injected_abs_range': [lo, hi], 'n': int(sel.sum()),
            'detection_rate': float(det[sel].mean())})

    # Recovered vs injected (signed) linear fit.
    inj_s = np.array([r['injected_effect'] for r in spiked])
    if len(inj_s) > 2 and inj_s.std() > 0:
        slope, intercept = np.polyfit(inj_s, rec, 1)
        corr = float(np.corrcoef(inj_s, rec)[0, 1])
    else:
        slope = intercept = corr = float('nan')

    # MDE concordance: detection rate above vs below the predicted MDE.
    mde = np.array([r['mde_median'] for r in spiked])
    above = inj > mde
    mde_check = {
        'detection_rate_injected_gt_mde':
            float(det[above].mean()) if above.sum() else None,
        'n_gt_mde': int(above.sum()),
        'detection_rate_injected_le_mde':
            float(det[~above].mean()) if (~above).sum() else None,
        'n_le_mde': int((~above).sum()),
    }

    specificity = {
        'n_null_controls': len(nulls),
        'false_positive_rate': float(np.mean([r['detected']
                                              for r in nulls])),
        'expected': args.fdr,
    }

    summary = {
        'n_variants': len(rows), 'n_spiked': len(spiked),
        'bin': args.bin, 'window_width': args.window_width,
        'n_null_iterations': args.n_null_iterations, 'seed': args.seed,
        'recovered_vs_injected': {'slope': float(slope),
                                  'intercept': float(intercept),
                                  'pearson_r': corr},
        'power_curve': power_curve,
        'mde_concordance': mde_check,
        'specificity': specificity,
    }
    out_json = os.path.join(args.output_dir, 'power_summary.json')
    with open(out_json, 'w') as fh:
        json.dump(summary, fh, indent=2)

    print(f"\n=== Power benchmark ({len(spiked)} spiked, "
          f"{len(nulls)} null controls) ===")
    print(f"recovered vs injected: slope={slope:.3f} r={corr:.3f} "
          f"(ideal slope~1, r~1)")
    print("power curve (detection rate by |injected effect|):")
    for pc in power_curve:
        lo, hi = pc['injected_abs_range']
        print(f"  |Δ| {lo:.2f}-{hi:.2f}: {pc['detection_rate']:.2f} "
              f"(n={pc['n']})")
    print(f"MDE concordance: detect rate  |inj|>MDE = "
          f"{mde_check['detection_rate_injected_gt_mde']} "
          f"(n={mde_check['n_gt_mde']})  |inj|<=MDE = "
          f"{mde_check['detection_rate_injected_le_mde']} "
          f"(n={mde_check['n_le_mde']})")
    print(f"specificity: false-positive rate among {len(nulls)} null "
          f"controls = {specificity['false_positive_rate']:.4f} "
          f"(expected ~{args.fdr})")
    print(f"\nTSV:  {tsv}\nJSON: {out_json}")


if __name__ == '__main__':
    main()

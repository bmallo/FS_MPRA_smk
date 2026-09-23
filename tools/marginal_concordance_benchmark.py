#!/usr/bin/env python3
"""
marginal_concordance_benchmark.py — controlled concordance benchmark for
--null-model marginal (validation test 3e, simulation-only).

Question
--------
Does marginal mode recover the same answer the standard WT-null gives
when the WT-null is valid (abundant WT) — and does it keep working when
the WT pool is too small for the WT-null, the regime that motivated
marginal mode in the first place?

Design
------
A single simulated dataset, three testing conditions:

  1. WT-null, ABUNDANT WT  — null sampled from a large pseudo-WT pool.
                             The standard method in its valid regime;
                             serves as the gold-standard reference.
  2. WT-null, SCARCE WT    — null sampled from a small pseudo-WT pool
                             (~the real saturating-library condition).
  3. MARGINAL              — null sampled from all reads (WT-abundance
                             irrelevant).

Only the null pool varies; the variant pools (marginal all-pooling) and
the injected truth are identical across conditions, so any difference is
attributable to the null pool alone.

Background and injection follow marginal_power_benchmark.py: variant
labels are permuted (H0), then known footprint knockouts are injected
into a grid of variants. The pseudo-WT pools are drawn from reads that
received no injected knockout, so they are a clean WT proxy; SCARCE is a
nested subsample of ABUNDANT.

Expected: condition 1 ~= condition 3 (concordance — both recover the
spikes); condition 2 degrades (power loss and/or specificity break);
condition 3 is unaffected by WT scarcity.
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

# See marginal_power_benchmark.py: penetrance x nuc-bin occupancy (~0.46)
# spans the realistic ~0.02-0.30 injected-effect range.
PENETRANCES = [0.05, 0.10, 0.20, 0.35, 0.50, 0.65]


def test_condition(rd, pool_idx, variant_groups, analysis_region,
                   prom_s, prom_e, seed, n_iter, threads):
    """Run the testing path with a given null pool; return {vid: result}."""
    results = run_variant_testing_parallel(
        rd, pool_idx, variant_groups, analysis_region, prom_s, prom_e,
        random_seed=seed, n_null_iterations=n_iter, n_workers=threads,
        stratify=False)
    qs = benjamini_hochberg(np.array([r['best_cluster_p'] for r in results]))
    return {r['variant_id']: {'result': r, 'q': float(q)}
            for r, q in zip(results, qs)}


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bam', required=True)
    p.add_argument('--target-region', required=True)
    p.add_argument('--promoter-region', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--n-spikes', type=int, default=150)
    p.add_argument('--window-width', type=int, default=40)
    p.add_argument('--bin', default='nuc')
    p.add_argument('--wt-abundant', type=int, default=8000,
                   help='Pseudo-WT pool size, ABUNDANT condition')
    p.add_argument('--wt-scarce', type=int, default=500,
                   help='Pseudo-WT pool size, SCARCE condition')
    p.add_argument('--n-null-iterations', type=int, default=2000)
    p.add_argument('--min-reads', type=int, default=50)
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--fdr', type=float, default=0.10)
    p.add_argument('--skip-shuffle', action='store_true')
    args = p.parse_args()
    setup_logging(quiet=True)
    os.makedirs(args.output_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # --- H0 background -----------------------------------------------
    if args.skip_shuffle:
        bg_bam = args.bam
    else:
        bg_bam = os.path.join(args.output_dir, 'background.bam')
        shuffle_bam(args.bam, bg_bam, args.seed)

    # --- Parse -------------------------------------------------------
    analysis_bins = build_analysis_bins()
    target_chrom, _, _ = parse_region(args.target_region)
    ref_info = get_bam_ref_info(bg_bam)
    ref_name, ref_length = resolve_target_chrom(ref_info, target_chrom)
    ps = args.promoter_region.split('-')
    prom_s, prom_e = int(ps[0]), int(ps[1]) + 1
    analysis_region = (prom_s, prom_e)
    analysis_length = prom_e - prom_s

    rd, ref_length, ref_name, _ = parse_bam(
        bg_bam, analysis_bins, target_chrom=target_chrom,
        analysis_region=analysis_region)
    variant_groups = group_variants(rd, include_multi=True,
                                    min_reads=args.min_reads,
                                    multi_mode='all')
    # Tiebreak by variant ID for fully deterministic spike selection
    # independent of dict insertion order.
    vids = sorted(variant_groups,
                  key=lambda v: (len(variant_groups[v]), v))

    # --- Inject known effects ----------------------------------------
    pick = np.unique(np.linspace(0, len(vids) - 1,
                                 min(args.n_spikes, len(vids))).astype(int))
    spike_vids = [vids[i] for i in pick]
    orig_bin = rd.coverage_matrices[args.bin].copy()
    spec, zeroed = {}, set()
    for i, vid in enumerate(spike_vids):
        f = PENETRANCES[i % len(PENETRANCES)]
        idx = variant_groups[vid]
        n_zero = int(round(f * len(idx)))
        if n_zero == 0:
            continue
        zr = rng.choice(idx, size=n_zero, replace=False)
        w0 = int(rng.integers(0, analysis_length - args.window_width))
        w1 = w0 + args.window_width
        rd.coverage_matrices[args.bin][np.ix_(zr, np.arange(w0, w1))] = 0
        spec[vid] = (w0, w1, f)
        zeroed.update(zr.tolist())

    # --- Pseudo-WT pools (clean of injected reads; SCARCE nested) ----
    clean = np.array([i for i in range(rd.n_reads) if i not in zeroed])
    rng.shuffle(clean)
    wt_abundant = np.sort(clean[:args.wt_abundant])
    wt_scarce = np.sort(wt_abundant[:args.wt_scarce])
    all_reads = np.arange(rd.n_reads)
    print(f"Injected {len(spec)} variants; pseudo-WT abundant="
          f"{len(wt_abundant)} scarce={len(wt_scarce)}; "
          f"marginal pool={len(all_reads)}")

    # --- Three conditions --------------------------------------------
    conditions = {
        'wt_abundant': wt_abundant,
        'wt_scarce': wt_scarce,
        'marginal': all_reads,
    }
    cond_res = {}
    for name, pool in conditions.items():
        print(f"  running condition: {name} (pool={len(pool):,})...")
        cond_res[name] = test_condition(
            rd, pool, variant_groups, analysis_region, prom_s, prom_e,
            args.seed, args.n_null_iterations, args.threads)

    # --- Per-variant injected truth ----------------------------------
    mod_bin = rd.coverage_matrices[args.bin]
    injected = {}
    for vid in variant_groups:
        idx = variant_groups[vid]
        w0, w1, _ = spec.get(vid, (0, analysis_length, 0.0))
        injected[vid] = float(mod_bin[idx][:, w0:w1].mean()
                               - orig_bin[idx][:, w0:w1].mean())

    # --- Summarise per condition -------------------------------------
    def metrics(cmap):
        spiked_det, null_det, recov = [], [], {}
        for vid, d in cmap.items():
            detected = d['q'] < args.fdr
            if vid in spec:
                spiked_det.append(detected)
                w0, w1, _ = spec[vid]
                recov[vid] = float(np.mean(
                    d['result'][args.bin]['delta_obs'][w0:w1]))
            else:
                null_det.append(detected)
        return {
            'power': float(np.mean(spiked_det)) if spiked_det else None,
            'n_spiked': len(spiked_det),
            'false_positive_rate':
                float(np.mean(null_det)) if null_det else None,
            'n_null': len(null_det),
            '_recovered': recov,
        }

    per_cond = {name: metrics(cmap) for name, cmap in cond_res.items()}

    # --- Concordance: abundant-WT (gold standard) vs marginal --------
    gold = per_cond['wt_abundant']['_recovered']
    marg = per_cond['marginal']['_recovered']
    common = sorted(set(gold) & set(marg))
    if len(common) > 2:
        g = np.array([gold[v] for v in common])
        m = np.array([marg[v] for v in common])
        concord_r = float(np.corrcoef(g, m)[0, 1])
    else:
        concord_r = float('nan')

    def detected_set(name):
        return {v for v, d in cond_res[name].items()
                if d['q'] < args.fdr and v in spec}
    d_gold = detected_set('wt_abundant')
    d_marg = detected_set('marginal')
    jac = (len(d_gold & d_marg) / len(d_gold | d_marg)
           if (d_gold | d_marg) else float('nan'))

    summary = {
        'seed': args.seed, 'bin': args.bin,
        'n_null_iterations': args.n_null_iterations,
        'wt_abundant_n': len(wt_abundant), 'wt_scarce_n': len(wt_scarce),
        'conditions': {name: {k: v for k, v in m.items()
                              if not k.startswith('_')}
                       for name, m in per_cond.items()},
        'concordance_wt_abundant_vs_marginal': {
            'recovered_effect_pearson_r': concord_r,
            'detected_set_jaccard': jac,
            'n_detected_wt_abundant': len(d_gold),
            'n_detected_marginal': len(d_marg),
        },
    }
    out_json = os.path.join(args.output_dir, 'concordance_summary.json')
    with open(out_json, 'w') as fh:
        json.dump(summary, fh, indent=2)

    print("\n=== Concordance benchmark ===")
    print(f"{'condition':<14} {'power':>8} {'FP-rate':>9} {'n_spiked':>9}")
    for name in ('wt_abundant', 'wt_scarce', 'marginal'):
        m = per_cond[name]
        pw = f"{m['power']:.3f}" if m['power'] is not None else 'NA'
        fp = (f"{m['false_positive_rate']:.4f}"
              if m['false_positive_rate'] is not None else 'NA')
        print(f"{name:<14} {pw:>8} {fp:>9} {m['n_spiked']:>9}")
    print(f"\nconcordance (WT-abundant vs marginal): "
          f"recovered-effect r={concord_r:.3f}, "
          f"detected-set Jaccard={jac:.3f}")
    print(f"JSON: {out_json}")


if __name__ == '__main__':
    main()

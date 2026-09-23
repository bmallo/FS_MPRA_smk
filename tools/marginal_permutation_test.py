#!/usr/bin/env python3
"""
marginal_permutation_test.py — calibration harness for --null-model marginal.

This is the negative control for the marginal null model (the analog of
tools/wt_vs_wt_calibration.py for the standard WT null).

Idea
----
Marginal mode tests each variant's reads against a null sampled from ALL
reads. To check that null is calibrated we destroy every real variant
effect while preserving the data's structure, then confirm the test
reports no signal.

We do this by PERMUTING the per-read variant labels: every primary read
keeps its own footprint data (MA/ns/nl tags, coverage) but is reassigned
the PV tag of a randomly chosen other read. This preserves
  * the variant-frequency spectrum (which combos exist, how often),
  * the per-read variant-count distribution,
  * the footprint data itself (real, not simulated),
and destroys any real variant<->footprint association. Under H0 every
"variant" is now a random set of reads, so a correctly calibrated test
must yield:
  * best_cluster_p ~ Uniform(0,1)
  * P(best_cluster_p <= alpha) ~= alpha
  * cross-variant BH calls ~0 variants significant

The shuffled BAM is run through the real 03_analyze_library.py with
--null-model marginal, so this also validates the flag wiring end to end.

Writes a small JSON snapshot for diffing / record-keeping.
"""

import argparse
import json
import os
import subprocess
import sys
import csv

import numpy as np
import pysam
from scipy import stats


def vc_from_pv(pv):
    """Variant count consistent with a PV tag value."""
    if pv == 'WT' or not pv:
        return 0
    try:
        return len(json.loads(pv))
    except (ValueError, TypeError):
        return pv.count(';') + 1


def shuffle_bam(in_bam, out_bam, seed):
    """Write out_bam: every primary read keeps its alignment + footprint
    tags but gets a PV/VC pair drawn (by permutation) from another read.
    Returns the number of reads whose labels were permuted."""
    # Pass 1: collect PV tags of primary reads, in file order.
    pvs = []
    with pysam.AlignmentFile(in_bam, 'rb') as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            if read.has_tag('PV'):
                pvs.append(read.get_tag('PV'))

    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(pvs))
    shuffled = [pvs[i] for i in perm]

    # Pass 2: reassign PV (and a consistent VC) to primary reads in the
    # same file order; everything else is copied verbatim.
    with pysam.AlignmentFile(in_bam, 'rb') as bam:
        with pysam.AlignmentFile(out_bam, 'wb', template=bam) as out:
            c = 0
            for read in bam:
                if (not read.is_secondary and not read.is_supplementary
                        and read.has_tag('PV')):
                    pv = shuffled[c]
                    c += 1
                    read.set_tag('PV', pv, 'Z')
                    read.set_tag('VC', vc_from_pv(pv), 'i')
                out.write(read)
    pysam.index(out_bam)
    return len(pvs)


def analyze_summary(tsv_path, alpha=0.05, fdr=0.10, n_quartiles=4):
    """Read a Stage 3 summary TSV; collect per-seed calibration metrics.

    NOTE: the KS test is reported as a DIAGNOSTIC only. It assumes the
    per-variant p-values are independent, which marginal all-pooling
    violates — every read sits in ~2-3 variant pools, so neighbouring
    variants' tests share reads and their p-values are correlated. A KS
    test on correlated values over-rejects. The authoritative metrics
    are the realized FDR and the 0.05-exceedance fraction, AGGREGATED
    across many permutation seeds (see aggregate_marginal_calib.py),
    where each seed is one independent Monte-Carlo replicate.

    The per-N-quartile exceedance is a self-contamination probe: a
    variant's reads are themselves in the all-reads null pool, so any
    residual anti-conservatism should grow with variant read count.
    """
    ids, nreads, pvals, qvals = [], [], [], []
    with open(tsv_path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            try:
                p = float(row['best_cluster_p'])
                q = float(row['variant_fdr_q'])
                nr = int(float(row['n_reads']))
            except (KeyError, ValueError):
                continue
            ids.append(row['variant_id'])
            nreads.append(nr)
            pvals.append(p)
            qvals.append(q)
    pvals = np.asarray(pvals, dtype=float)
    qvals = np.asarray(qvals, dtype=float)
    nreads = np.asarray(nreads, dtype=int)
    n = len(pvals)
    if n == 0:
        return {'n_variants': 0, 'error': 'no variants in summary TSV'}

    ks = stats.kstest(pvals, 'uniform')
    n_sig = int(np.sum(qvals < fdr))

    # Per-N-quartile 0.05-exceedance (self-contamination probe).
    order = np.argsort(nreads)
    by_n = []
    for qi in range(n_quartiles):
        sel = order[qi * n // n_quartiles:(qi + 1) * n // n_quartiles]
        if len(sel) == 0:
            continue
        by_n.append({
            'n_reads_range': [int(nreads[sel].min()),
                              int(nreads[sel].max())],
            'n_variants': int(len(sel)),
            'frac_p_le_alpha': float(np.mean(pvals[sel] <= alpha)),
        })

    return {
        'n_variants': n,
        'alpha': alpha,
        'fdr_threshold': fdr,
        'frac_p_le_alpha': float(np.mean(pvals <= alpha)),
        'n_significant_bh': n_sig,
        'p_min': float(np.min(pvals)),
        'p_median': float(np.median(pvals)),
        'by_n_quartile': by_n,
        'ks_statistic_diagnostic': float(ks.statistic),
        'ks_pvalue_diagnostic': float(ks.pvalue),
    }


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bam', required=True, help='Tagged BAM (Stage 2 output)')
    p.add_argument('--target-region', required=True)
    p.add_argument('--promoter-region', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--seed', type=int, default=42,
                   help='Permutation seed (default: 42)')
    p.add_argument('--n-null-iterations', type=int, default=2000)
    p.add_argument('--min-reads', type=int, default=50)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--analyze-script', default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..', 'workflow', 'scripts', '03_analyze_library.py'))
    p.add_argument('--python', default=sys.executable)
    p.add_argument('--skip-shuffle', action='store_true',
                   help='Reuse an existing shuffled.bam in --output-dir')
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    shuffled = os.path.join(args.output_dir, 'shuffled.bam')

    if not args.skip_shuffle:
        print(f"[1/3] Permuting variant labels (seed={args.seed})...")
        n = shuffle_bam(args.bam, shuffled, args.seed)
        print(f"      {n:,} primary reads relabelled -> {shuffled}")
    else:
        print(f"[1/3] Reusing existing {shuffled}")

    print("[2/3] Running 03_analyze_library.py --null-model marginal...")
    cmd = [
        args.python, args.analyze_script,
        '--bam', shuffled,
        '-o', args.output_dir, '-n', 'permtest',
        '--target-region', args.target_region,
        '--promoter-region', args.promoter_region,
        '--null-model', 'marginal',
        '--min-reads', str(args.min_reads),
        '--n-null-iterations', str(args.n_null_iterations),
        '--threads', str(args.threads),
        '-q',
    ]
    r = subprocess.run(cmd)
    if r.returncode != 0:
        sys.exit(f"03_analyze_library.py failed (exit {r.returncode})")

    print("[3/3] Assessing p-value calibration...")
    summary_tsv = os.path.join(args.output_dir, 'permtest_summary.tsv')
    result = analyze_summary(summary_tsv)
    result['seed'] = args.seed
    result['bam'] = os.path.abspath(args.bam)
    result['target_region'] = args.target_region
    result['n_null_iterations'] = args.n_null_iterations

    out_json = os.path.join(args.output_dir,
                            f'marginal_calib.seed{args.seed}.json')
    with open(out_json, 'w') as fh:
        json.dump(result, fh, indent=2)

    print(f"\n=== Marginal-null calibration: seed {args.seed} "
          f"(one Monte-Carlo replicate) ===")
    for k in ('n_variants', 'frac_p_le_alpha', 'n_significant_bh',
              'p_median', 'ks_statistic_diagnostic'):
        print(f"  {k:28s} {result.get(k)}")
    print("  by_n_quartile (frac p<=0.05 vs variant read count):")
    for q in result.get('by_n_quartile', []):
        rng = q['n_reads_range']
        print(f"    N={rng[0]:>4d}-{rng[1]:<5d} "
              f"({q['n_variants']:>3d} vars)  "
              f"frac={q['frac_p_le_alpha']:.4f}")
    print(f"\nJSON: {out_json}")
    print("This is ONE permutation. Run several seeds and aggregate with "
          "aggregate_marginal_calib.py for the calibration verdict — "
          "a single seed's KS test is invalid (correlated p-values).")


if __name__ == '__main__':
    main()

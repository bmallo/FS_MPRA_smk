#!/usr/bin/env python3
"""Aggregate the WT-vs-WT calibration sweep by arm.

Per arm: pool best_cluster_p across its cells for the MARGINAL check
(KS vs Uniform, P(p<=alpha), median). FDR is judged PER CELL — each
cell's disjoint pseudo-variants are internally independent, so BH
within a cell is the valid all-null FDR test — then aggregated as
total false rejections / total pseudo-variants across the arm.

Usage: aggregate_sweep.py [sweep_dir]
  default sweep_dir = results/phase0_hmm_full/sweep
"""
import glob
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "workflow", "scripts"))
from common import benjamini_hochberg  # noqa: E402

sweep = (sys.argv[1] if len(sys.argv) > 1
         else os.path.join(_HERE, "..", "..", "results",
                           "phase0_hmm_full", "sweep"))

cells = {}
for fp in sorted(glob.glob(os.path.join(sweep, "cell_*_*.json"))):
    arm = os.path.basename(fp).split("_", 2)[2].rsplit("_seed", 1)[0]
    d = json.load(open(fp))
    cells.setdefault(arm, []).append(d.get("best_cluster_p_sorted", []))

if not cells:
    sys.exit(f"no cell_*.json under {sweep}")

print(f"{'arm':<12}{'cells':>6}{'pooled_n':>9}{'median':>8}"
      f"{'KS':>7}{'KScrit':>8}{'P<=.01':>8}{'P<=.05':>8}"
      f"{'FP q<.05':>10}{'FP q<.10':>10}")
for arm in sorted(cells):
    arrs = [np.asarray(a, float) for a in cells[arm] if len(a)]
    if not arrs:
        continue
    pooled = np.sort(np.concatenate(arrs))
    n = len(pooled)
    ks = float(np.max(np.abs(np.arange(1, n + 1) / n - pooled)))
    fp05 = fp10 = 0
    for a in arrs:
        q = benjamini_hochberg(np.sort(a))
        fp05 += int(np.sum(q < 0.05))
        fp10 += int(np.sum(q < 0.10))
    print(f"{arm:<12}{len(arrs):>6}{n:>9}{np.median(pooled):>8.3f}"
          f"{ks:>7.3f}{1.36/np.sqrt(n):>8.3f}"
          f"{float(np.mean(pooled<=0.01)):>8.3f}"
          f"{float(np.mean(pooled<=0.05)):>8.3f}"
          f"{f'{fp05}/{n}':>10}{f'{fp10}/{n}':>10}")
print("\nPASS = per-cell BH false positives ~0 (FP q<.05/.10), "
      "P(p<=a) ~ a, KS within ~KScrit — across ALL arms.")

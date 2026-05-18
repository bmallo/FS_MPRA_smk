#!/usr/bin/env python3
"""Generate the calibration-sweep cell manifest (deterministic).

Each cell = one independent disjoint WT-vs-WT run:
  A_prod     production config (median N, stratify ON) — headline verdict
  B_stratoff stratify OFF at production N — validates the stratification
             approximation at scale
  C_lowN/C_highN  calibration across the coverage range (p10 / p90)

Writes cells.tsv next to this script. Re-run to regenerate.
"""
import collections
import os

ROWS = []
_seed = 100


def add(arm, N, stratify, k):
    global _seed
    for _ in range(k):
        ROWS.append((arm, N, stratify, _seed))
        _seed += 1


# (arm, pseudo-variant N, stratify, n_seeds)
add("A_prod",     1248, "on",  20)
add("B_stratoff", 1248, "off",  8)
add("C_lowN",      200, "on",   4)
add("C_highN",    2900, "on",  10)

out = os.path.join(os.path.dirname(__file__), "cells.tsv")
with open(out, "w") as f:
    f.write("idx\tarm\tN\tstratify\tseed\n")
    for i, (a, N, s, sd) in enumerate(ROWS, 1):
        f.write(f"{i}\t{a}\t{N}\t{s}\t{sd}\n")
print(f"{len(ROWS)} cells -> {out}")
print(dict(collections.Counter(r[0] for r in ROWS)))

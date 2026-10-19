# Expanded replicate-3 WT calibration

This diagnostic uses the fixed contract in
`config/ldlr_heldout_calibration.json`. It excludes the 10,000 WT barcodes in
the earlier footprint-quality audit and retains all 948,189 eligible molecules
from the remaining 258,283 WT barcodes. These barcodes are disjoint from that
recent audit sample. They are **not globally untouched**: the earlier full-data
calibration drew from the entire WT pool, and all barcodes came from the same
biological preparation. This stage does not release variant significance.

Run on a Slurm compute node, using one CPU and the allocated job memory. The
workflow reads the existing joined footprint TSV and aggregate arrays; it does
not reopen the annotation BAM. It writes three atomic directories under
`results/rep3_heldout_calibration_20261018`:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
export MPLCONFIGDIR=/tmp/fsmpra-matplotlib
/mmfs1/home/bmallo/miniconda3/envs/snakemake/bin/snakemake \
  -s workflow/heldout_calibration.smk --cores 1 --resources mem_mb=4096
```

The equivalent direct commands, in dependency order, are:

```bash
/mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.heldout_outcomes \
  --config config/ldlr_heldout_calibration.json \
  --output results/rep3_heldout_calibration_20261018/outcomes
/mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.heldout_calibration \
  --config config/ldlr_heldout_calibration.json \
  --outcomes results/rep3_heldout_calibration_20261018/outcomes \
  --output results/rep3_heldout_calibration_20261018/calibration
/mmfs1/home/bmallo/miniconda3/bin/python tools/plot_heldout_calibration.py \
  --calibration results/rep3_heldout_calibration_20261018/calibration \
  --output results/rep3_heldout_calibration_20261018/plots
```

`outcomes/endpoint_matrix.npz` has sorted barcode IDs, depths, 22 fixed endpoint
IDs and coordinates, float64 `presence` and `coverage` matrices, and
`valid_molecules`. Presence is calculated per molecule before averaging within
barcode. An evaluated empty call contributes zero. Every held-out molecule
must align across all 211 element bases. The builder checks every barcode's
depth against the aggregate metadata and every coverage endpoint against the
existing per-base count arrays, with maximum allowed discrepancy `1e-12`.
`outcomes/audit.json` records conservation totals and frozen source identities.

`calibration/summary.tsv` contains one row per endpoint, scenario, and method.
The fixed WT A/A grid has 21 TF presence endpoints, 21 TF coverage endpoints,
and one whole-element MSP coverage endpoint, at pseudo-variant sizes 50, 100,
400 and 1,000. Each size has 1,000 observed draws compared with its entire WT
complement. Welch and studentized permutation results use the same observed
draws; the permutation reference bank has 9,999 draws per size and is shared
across endpoints and observed contrasts. `calibration/iterations.tsv.gz` holds
the individual contrasts. The minimum possible permutation p-value is 0.0001.
The summary gives rejection rates at 0.05 and 0.01 with exact binomial 95%
intervals. For WT A/A, those intervals condition on both this fixed WT pool and
the shared reference bank. They do not include biological-preparation or
reference-bank Monte Carlo uncertainty.

The same calibration output includes 32 synthetic equal-mean stress scenarios
and 18 alternative scenarios, each with 200 independent iterations and 999
label permutations per contrast. The synthetic reference arm has 10,000
barcodes. The model crosses sparse and common means, unequal dispersion, and
independent or propensity-associated depths. It is a bounded approximation to
the full WT reference size. Welch and Hoeffding 95% intervals are evaluated
separately for truth coverage and width; neither is a permutation interval.

`plots/wt_coverage_rates.png/.svg` and `wt_presence_rates.png/.svg` display the
WT A/A 0.05 rejection percentages for both methods across all sizes and fixed
endpoints. Their dashed colorbar mark is the nominal 5% target, not an
automatic pass/fail threshold. `synthetic_null_rates.png/.svg` shows rates and
exact binomial intervals for all 32 null scenarios. The compact
`effect_and_interval_summary.tsv` lists 0.05 and 0.01 power, directional
rejections, and Welch/Hoeffding interval coverage and width for the 18
alternative scenarios. These figures and tables report diagnostics, including
failures; they do not compute variant p/q-values.

Each stage checks source identities before and after computation and publishes
only a completed directory. Existing outputs are never overwritten. The
calibration accepts `--pilot-iterations` and `--pilot-permutations` only for a
separate output directory. Pilot results carry `partial: true`; plotting rejects
them unless explicitly called with `--allow-pilot`, and pilot figures are
labeled. The workflow runs the frozen complete settings only. For the fixed
design and interpretation details, see
`docs/rep3_heldout_calibration_design.md`.

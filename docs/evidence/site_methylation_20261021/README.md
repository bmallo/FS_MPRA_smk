# Site-methylation stage, replicate 3 (run 2026-10-20)

Prespecified run of `fsmpra.site_extract` / `fsmpra.site_effects` (see `docs/site_methylation_usage.md`). Numbers only; interpretation rules are applied mechanically below.

## Provenance

- Git commit: `fc4731a891a3c56b9377db47c6fb720aabda54ec` (branch `marginal-null-model`); no local changes at run time.
- Source directories (Klone, repo root `/mmfs1/gscratch/stergachislab/bmallo/large_home/code/mine/FS_MPRA_smk`):
  - `results/rep3_site_methylation_20261021/extract` → `extract_audit.json`
  - `results/rep3_site_methylation_20261021/effects` → `effects/` (primary)
  - `results/rep3_site_methylation_20261021/sensitivity/u0.02_m0.7` → `u0.02_m0.7/`
  - `results/rep3_site_methylation_20261021/sensitivity/u0.005_m0.5` → `u0.005_m0.5/`
  - `results/rep3_site_methylation_20261021/effects_supported_barcodes` → `effects_supported_barcodes/`
- Slurm: `--account=stergachislab --partition=cpu-g2`, 1 CPU, 16G, BLAS/OpenMP threads = 1.
  - Smoke/tests: job 42090480 (n3480), 7:11 elapsed, 132 tests OK (1 skipped). Smoke output: `results/site_smoke_42090480/`.
  - Full: job 42090780 (n3480), 59:17 elapsed, 0 requeues, sacct MaxRSS 8,619,636 KB. Log: `results/site_methylation_jobs/full_42090780.log`.

| Stage | Wall clock | Max RSS (KB) |
|---|---|---|
| site_extract | 53:38.66 | 2,559,184 |
| site_effects (primary, 200 AA draws) | 2:24.43 | 8,539,884 |
| sensitivity u0.02_m0.7 (100 AA draws) | 0:51.36 | 8,539,880 |
| sensitivity u0.005_m0.5 (100 AA draws) | 0:50.98 | 8,540,080 |
| effects_supported_barcodes (100 AA draws) | 0:54.94 | 7,359,440 |

Extraction audit: `complete` True; `molecules` 2,659,690; `barcodes` 568,386; `bam_primary_rows_scanned` 8,591,169. Nonzero `retained_positive_ml_histogram` bins: 240: 52,457; 241: 278,701; 242: 231,341; 243: 324,507; 244: 288,115; 245: 412,692; 246: 365,544; 247: 411,645; 248: 462,843; 249: 692,560; 250: 664,702; 251: 817,600; 252: 1,070,592; 253: 1,513,042; 254: 2,573,229; 255: 82,008,006.

Snakemake dry run (`workflow/site_methylation.smk -n`): "Nothing to be done"; 5 jobs flagged with missing provenance metadata (site_effects, site_effects_sensitivity, site_effects_supported_barcodes, site_extract), as expected for direct runs.

## Thresholds

Primary (config defaults, unchanged): U distal m6A ≤ 0.01, M distal m6A ≥ 0.60. Sensitivities: 0.02/0.70 and 0.005/0.50. Labeling strata edges are in section 1.

Conventions: deltas and SEs in percentage points (×100). Region membership of an SNV uses the code's mask rule `a ≤ oligo_pos1 − 21 < b` on `regions0`. Region SE = sqrt(Σ se²)/n. "≈ 0" predictions are judged as |z| ≤ 2. `p_circular_shift` = 0.005 is the floor for 200 rotations.

## Rule outcomes (mechanical)

- Rule 1 (calibration gate): 56/73 non-U endpoints pass. Primary endpoints: flank:A_norm pass, msp_C_std pass, nuc_85_250_C_std pass, **local:A_norm fails** (rate_ci_high 0.121 > 0.10; sd_z 1.053).
- Rule 2 (observability reach): M kernel has |mean| > 2 SE at |offset| > 2 for offsets −14, +6, +10 (at_change = all, ±15 bp examined).
- Rule 3 (position 125): observability clause met (WT_P3:M same sign and ≥ half of WT_P3:C_std for 3/3 substitutions). Chromatin clause not met (WT_P3:A_norm z ≤ 3 for 125A>C and 125A>G under u0.005_m0.5: z = 2.82, 2.56).
- Rule 4 (regions): P1_core:A_norm > 0 agree; P1_core:M no change disagree (+2.82 pp, z 15.9); Repeat1 and Repeat3 A_norm < 0, msp_C_std < 0, nuc_85_250_C_std > 0 all agree; SRE-1 msp_C_std ≈ 0 disagree (+0.24 pp, z 2.41).
- Rule 5 (association): r(msp_C_std) > r(msp_ALL) > 0 and r(nuc_85_250_C_std) < 0 hold for both activity sources, variant- and position-level, in both the primary and supported-barcode cohorts.

## Figures (descriptive)

- `effects/distal_m6a_histogram.png`: WT and SNV distributions overlap. The U cutoff (0.01) sits on the left-edge spike at ~0; the M cutoff (0.60) sits on the rising shoulder below a mode near 0.70, after which the distribution falls steeply to ~0.8.
- `effects/kernel.png`: in A/T-gain and A/T-loss panels, M and C_std both show a large dip/peak (M ~±30 pp, C_std ~±16 pp) concentrated within about ±5 bp; U is flat; the retained-count panel is flat apart from wider bands near 0.
- `effects/calibration_z.png`: WT pseudo-variant z histograms for six endpoints, roughly matching the standard normal; SDs 0.96–1.07.
- `effects/observability_scatter.png`: local Δ in M (x) vs Δ in C_std (y) shows a positive linear trend; the three 125 substitutions sit in the upper-right quadrant (M ≈ +4 to +6 pp, C ≈ +4 to +5 pp).

## Step-4 tables
## 1. Molecules and classes (primary)

```json
{
 "molecules": {
  "cohort": 2594368,
  "classified": 2594368,
  "wt": {
   "U": 93030,
   "M": 266758,
   "C": 884103
  },
  "snv": {
   "U": 101223,
   "M": 289915,
   "C": 959339
  }
 },
 "thresholds": {
  "unmethylated_max_distal_m6a": 0.01,
  "methylated_min_distal_m6a": 0.6,
  "labeling_strata_edges": [
   0.04260780289769173,
   0.07854209095239639,
   0.12012320011854172,
   0.18326488137245178,
   0.2643737196922302,
   0.34240245819091797,
   0.4111909568309784,
   0.4758726954460144,
   0.5385010242462158
  ],
  "labeling_strata_wt_weights": [
   0.09951216091337774,
   0.09939905192042103,
   0.10021230557977973,
   0.10030053059428597,
   0.1001206872954848,
   0.0999317952772471,
   0.09992953309738797,
   0.10057425435724118,
   0.09937190576211143,
   0.10064777520266303
  ]
 },
 "class_dispersion": {
  "fraction_U_variance_ratio": 1.6210649833237256,
  "fraction_M_variance_ratio": 0.597508226258629,
  "interpretation": "Observed between-SNV variance of class fractions over the molecule-level binomial expectation; values well above 1 reflect barcode clustering and/or genotype-dependent class membership."
 }
}
```

WT classified molecules 1,243,891: U fraction 0.0748, M fraction 0.2145, C fraction 0.7108

- u0.02_m0.7: WT U 121,914 (0.0980), M 79,521 (0.0639), C 1,042,456
- u0.005_m0.5: WT U 75,302 (0.0605), M 409,699 (0.3294), C 758,890

## 2. Calibration gate (primary; non-U endpoints)

Passing: 56 / 73

Primary endpoints:

| endpoint | draws | sd_z | rate_abs_z_gt_1_96 | rate_ci_high | rate_lower_tail | rate_upper_tail | pass |
|---|---|---|---|---|---|---|---|
| flank:A_norm | 200 | 0.973 | 0.055 | 0.096 | 0.025 | 0.030 | yes |
| local:A_norm | 200 | 1.053 | 0.075 | 0.121 | 0.035 | 0.040 | NO |
| msp_C_std | 200 | 0.962 | 0.050 | 0.090 | 0.025 | 0.025 | yes |
| nuc_85_250_C_std | 200 | 0.997 | 0.040 | 0.077 | 0.010 | 0.030 | yes |

Failing non-U endpoints:

| endpoint | draws | sd_z | rate_abs_z_gt_1_96 | rate_ci_high | rate_lower_tail | rate_upper_tail |
|---|---|---|---|---|---|---|
| P1_core:M | 200 | 0.997 | 0.065 | 0.109 | 0.030 | 0.035 |
| Repeat1_Sp1:M | 200 | 1.051 | 0.075 | 0.121 | 0.035 | 0.040 |
| Repeat2_SRE1:A_norm | 200 | 0.968 | 0.060 | 0.102 | 0.035 | 0.025 |
| Repeat2_SRE1:C_std | 200 | 0.986 | 0.070 | 0.115 | 0.040 | 0.030 |
| WT_P2:M | 200 | 1.060 | 0.065 | 0.109 | 0.030 | 0.035 |
| WT_P3:ALL | 200 | 1.031 | 0.070 | 0.115 | 0.035 | 0.035 |
| WT_P3:M | 200 | 1.017 | 0.060 | 0.102 | 0.030 | 0.030 |
| WT_P4:A_norm | 200 | 1.014 | 0.060 | 0.102 | 0.035 | 0.025 |
| WT_P4:M | 200 | 1.071 | 0.060 | 0.102 | 0.010 | 0.050 |
| adjacent_143_150:M | 199 | 1.137 | 0.070 | 0.115 | 0.010 | 0.060 |
| element:M | 200 | 1.076 | 0.080 | 0.127 | 0.035 | 0.045 |
| flank:M | 200 | 1.004 | 0.060 | 0.102 | 0.015 | 0.045 |
| local:ALL | 200 | 1.006 | 0.060 | 0.102 | 0.025 | 0.035 |
| local:A_norm | 200 | 1.053 | 0.075 | 0.121 | 0.035 | 0.040 |
| local:M | 200 | 1.069 | 0.080 | 0.127 | 0.025 | 0.055 |
| near:A_norm | 184 | 1.023 | 0.082 | 0.131 | 0.033 | 0.049 |
| tf_small_ALL | 200 | 1.011 | 0.060 | 0.102 | 0.030 | 0.030 |

## 3. Observability kernel (primary, at_change = all; pp)

| measure | offset | pairs | mean_delta_pp | se_pp | |mean|>2SE & |off|>2 |
|---|---|---|---|---|---|
| C_std | -15 | 281 | -0.15 | 0.11 |  |
| C_std | -14 | 285 | -0.22 | 0.10 | * |
| C_std | -13 | 284 | -0.14 | 0.12 |  |
| C_std | -12 | 284 | -0.25 | 0.12 | * |
| C_std | -11 | 285 | -0.15 | 0.13 |  |
| C_std | -10 | 286 | -0.08 | 0.17 |  |
| C_std | -9 | 290 | 0.11 | 0.17 |  |
| C_std | -8 | 293 | 0.06 | 0.34 |  |
| C_std | -7 | 292 | -0.27 | 0.23 |  |
| C_std | -6 | 293 | 0.13 | 0.20 |  |
| C_std | -5 | 296 | -0.03 | 0.24 |  |
| C_std | -4 | 295 | -0.26 | 0.25 |  |
| C_std | -3 | 295 | 0.24 | 0.32 |  |
| C_std | -2 | 299 | -0.58 | 0.48 |  |
| C_std | -1 | 300 | -0.55 | 1.07 |  |
| C_std | 1 | 301 | -1.93 | 0.97 |  |
| C_std | 2 | 299 | 0.38 | 0.50 |  |
| C_std | 3 | 297 | 0.34 | 0.33 |  |
| C_std | 4 | 291 | -0.28 | 0.28 |  |
| C_std | 5 | 293 | -0.04 | 0.21 |  |
| C_std | 6 | 292 | -0.33 | 0.20 |  |
| C_std | 7 | 289 | 0.08 | 0.24 |  |
| C_std | 8 | 285 | 0.16 | 0.25 |  |
| C_std | 9 | 287 | 0.28 | 0.18 |  |
| C_std | 10 | 286 | -0.44 | 0.14 | * |
| C_std | 11 | 285 | -0.15 | 0.14 |  |
| C_std | 12 | 284 | -0.32 | 0.11 | * |
| C_std | 13 | 280 | -0.11 | 0.11 |  |
| C_std | 14 | 281 | -0.16 | 0.10 |  |
| C_std | 15 | 281 | -0.12 | 0.12 |  |
| M | -15 | 281 | -0.26 | 0.20 |  |
| M | -14 | 285 | -0.44 | 0.18 | * |
| M | -13 | 284 | -0.05 | 0.20 |  |
| M | -12 | 284 | -0.30 | 0.22 |  |
| M | -11 | 285 | -0.11 | 0.26 |  |
| M | -10 | 286 | -0.18 | 0.32 |  |
| M | -9 | 290 | 0.09 | 0.35 |  |
| M | -8 | 293 | -0.61 | 0.69 |  |
| M | -7 | 292 | -0.34 | 0.41 |  |
| M | -6 | 293 | 0.36 | 0.41 |  |
| M | -5 | 296 | -0.05 | 0.43 |  |
| M | -4 | 295 | -0.56 | 0.49 |  |
| M | -3 | 295 | -0.19 | 0.66 |  |
| M | -2 | 299 | -1.26 | 1.03 |  |
| M | -1 | 300 | -1.88 | 2.37 |  |
| M | 1 | 301 | -5.51 | 2.16 |  |
| M | 2 | 299 | -0.80 | 1.04 |  |
| M | 3 | 297 | 0.76 | 0.63 |  |
| M | 4 | 291 | -0.16 | 0.48 |  |
| M | 5 | 293 | 0.12 | 0.44 |  |
| M | 6 | 292 | -0.91 | 0.40 | * |
| M | 7 | 289 | 0.04 | 0.52 |  |
| M | 8 | 285 | 0.53 | 0.47 |  |
| M | 9 | 287 | 0.55 | 0.38 |  |
| M | 10 | 286 | -0.61 | 0.26 | * |
| M | 11 | 285 | -0.28 | 0.24 |  |
| M | 12 | 284 | -0.22 | 0.16 |  |
| M | 13 | 280 | -0.07 | 0.15 |  |
| M | 14 | 281 | -0.14 | 0.16 |  |
| M | 15 | 281 | -0.20 | 0.17 |  |

Flagged offsets: M ['-14', '6', '10']; C_std ['-14', '-12', '10', '12']

## 4. Position 125

### primary u0.01_m0.6

| variant | endpoint | delta_pp | se_pp | z | molecules_M |
|---|---|---|---|---|---|
| 125A>C | WT_P3:M | 4.37 | 0.64 | 6.82 | 400 |
| 125A>C | WT_P3:C_std | 4.07 | 0.72 | 5.62 | 400 |
| 125A>C | WT_P3:A_norm | 3.82 | 1.14 | 3.34 | 400 |
| 125A>C | local:M | 4.39 | 0.65 | 6.77 | 400 |
| 125A>C | local:C_std | 3.89 | 0.71 | 5.48 | 400 |
| 125A>C | local:A_norm | 3.50 | 1.15 | 3.05 | 400 |
| 125A>G | WT_P3:M | 6.12 | 0.71 | 8.58 | 451 |
| 125A>G | WT_P3:C_std | 4.38 | 0.59 | 7.36 | 451 |
| 125A>G | WT_P3:A_norm | 3.11 | 0.96 | 3.22 | 451 |
| 125A>G | local:M | 6.20 | 0.70 | 8.82 | 451 |
| 125A>G | local:C_std | 3.95 | 0.58 | 6.86 | 451 |
| 125A>G | local:A_norm | 2.36 | 0.94 | 2.52 | 451 |
| 125A>T | WT_P3:M | 6.17 | 0.64 | 9.58 | 430 |
| 125A>T | WT_P3:C_std | 5.46 | 0.69 | 7.89 | 430 |
| 125A>T | WT_P3:A_norm | 4.52 | 1.06 | 4.27 | 430 |
| 125A>T | local:M | 5.83 | 0.68 | 8.58 | 430 |
| 125A>T | local:C_std | 5.35 | 0.66 | 8.08 | 430 |
| 125A>T | local:A_norm | 4.56 | 1.07 | 4.25 | 430 |

### u0.02_m0.7

| variant | endpoint | delta_pp | se_pp | z | molecules_M |
|---|---|---|---|---|---|
| 125A>C | WT_P3:M | 2.92 | 1.07 | 2.74 | 122 |
| 125A>C | WT_P3:C_std | 4.35 | 0.64 | 6.80 | 122 |
| 125A>C | WT_P3:A_norm | 4.39 | 1.13 | 3.89 | 122 |
| 125A>C | local:M | 2.77 | 1.12 | 2.48 | 122 |
| 125A>C | local:C_std | 4.19 | 0.63 | 6.69 | 122 |
| 125A>C | local:A_norm | 4.24 | 1.17 | 3.64 | 122 |
| 125A>G | WT_P3:M | 4.51 | 1.24 | 3.63 | 129 |
| 125A>G | WT_P3:C_std | 4.93 | 0.54 | 9.19 | 129 |
| 125A>G | WT_P3:A_norm | 4.16 | 1.19 | 3.50 | 129 |
| 125A>G | local:M | 4.73 | 1.20 | 3.95 | 129 |
| 125A>G | local:C_std | 4.58 | 0.52 | 8.78 | 129 |
| 125A>G | local:A_norm | 3.54 | 1.15 | 3.07 | 129 |
| 125A>T | WT_P3:M | 4.58 | 0.98 | 4.66 | 127 |
| 125A>T | WT_P3:C_std | 5.81 | 0.62 | 9.39 | 127 |
| 125A>T | WT_P3:A_norm | 5.13 | 1.14 | 4.50 | 127 |
| 125A>T | local:M | 4.44 | 1.03 | 4.31 | 127 |
| 125A>T | local:C_std | 5.65 | 0.60 | 9.48 | 127 |
| 125A>T | local:A_norm | 5.05 | 1.19 | 4.26 | 127 |

### u0.005_m0.5

| variant | endpoint | delta_pp | se_pp | z | molecules_M |
|---|---|---|---|---|---|
| 125A>C | WT_P3:M | 5.34 | 0.76 | 7.04 | 570 |
| 125A>C | WT_P3:C_std | 3.69 | 0.74 | 5.00 | 570 |
| 125A>C | WT_P3:A_norm | 3.41 | 1.21 | 2.82 | 570 |
| 125A>C | local:M | 5.47 | 0.75 | 7.27 | 570 |
| 125A>C | local:C_std | 3.41 | 0.71 | 4.79 | 570 |
| 125A>C | local:A_norm | 2.82 | 1.16 | 2.42 | 570 |
| 125A>G | WT_P3:M | 6.39 | 0.74 | 8.64 | 696 |
| 125A>G | WT_P3:C_std | 3.68 | 0.61 | 6.00 | 696 |
| 125A>G | WT_P3:A_norm | 2.69 | 1.05 | 2.56 | 696 |
| 125A>G | local:M | 6.41 | 0.73 | 8.82 | 696 |
| 125A>G | local:C_std | 3.20 | 0.59 | 5.39 | 696 |
| 125A>G | local:A_norm | 1.77 | 1.02 | 1.73 | 696 |
| 125A>T | WT_P3:M | 6.22 | 0.72 | 8.69 | 648 |
| 125A>T | WT_P3:C_std | 5.04 | 0.73 | 6.94 | 648 |
| 125A>T | WT_P3:A_norm | 4.77 | 1.14 | 4.18 | 648 |
| 125A>T | local:M | 6.03 | 0.71 | 8.49 | 648 |
| 125A>T | local:C_std | 4.89 | 0.70 | 7.00 | 648 |
| 125A>T | local:A_norm | 4.56 | 1.15 | 3.98 | 648 |

- 125A>C: WT_P3:M 4.37 pp vs WT_P3:C_std 4.07 pp -> same direction & >= half: True
- 125A>G: WT_P3:M 6.12 pp vs WT_P3:C_std 4.38 pp -> same direction & >= half: True
- 125A>T: WT_P3:M 6.17 pp vs WT_P3:C_std 5.46 pp -> same direction & >= half: True
- Rule 3 observability clause (>= 2 of 3, primary thresholds): True (3/3; 0 missing)
- Rule 3 chromatin clause (WT_P3:A_norm z > 3 for all 3 x 3 thresholds): False; primary u0.01_m0.6/125A>C z=3.34, primary u0.01_m0.6/125A>G z=3.22, primary u0.01_m0.6/125A>T z=4.27, u0.02_m0.7/125A>C z=3.89, u0.02_m0.7/125A>G z=3.50, u0.02_m0.7/125A>T z=4.50, u0.005_m0.5/125A>C z=2.82, u0.005_m0.5/125A>G z=2.56, u0.005_m0.5/125A>T z=4.18

## 5. Prespecified regions (primary; mean per-SNV delta, pp)

| region | endpoint | SNVs (rows/with value) | mean_delta_pp | se_pp | z | prediction | agree/disagree |
|---|---|---|---|---|---|---|---|
| P1_core | P1_core:A_norm | 27/27 | 2.96 | 0.28 | 10.69 | > 0 | agree |
| P1_core | P1_core:M | 27/27 | 2.82 | 0.18 | 15.92 | ~0 (no M change) | disagree (|z|>2) |
| P1_core | msp_C_std | 27/27 | -0.39 | 0.14 | -2.79 |  |  |
| P1_core | nuc_85_250_C_std | 27/27 | 0.25 | 0.16 | 1.55 |  |  |
| Repeat1_Sp1 | Repeat1_Sp1:A_norm | 48/48 | -1.42 | 0.19 | -7.62 | < 0 | agree |
| Repeat1_Sp1 | Repeat1_Sp1:M | 48/48 | -0.64 | 0.15 | -4.18 |  |  |
| Repeat1_Sp1 | msp_C_std | 48/48 | -1.00 | 0.11 | -9.45 | < 0 | agree |
| Repeat1_Sp1 | nuc_85_250_C_std | 48/48 | 0.85 | 0.12 | 6.89 | > 0 | agree |
| Repeat2_SRE1 | Repeat2_SRE1:A_norm | 48/48 | 0.92 | 0.17 | 5.27 |  |  |
| Repeat2_SRE1 | Repeat2_SRE1:M | 48/48 | 1.85 | 0.12 | 15.77 |  |  |
| Repeat2_SRE1 | msp_C_std | 48/48 | 0.24 | 0.10 | 2.41 | ~0 | disagree (|z|>2) |
| Repeat2_SRE1 | nuc_85_250_C_std | 48/48 | -0.18 | 0.12 | -1.54 |  |  |
| Repeat3_Sp1 | Repeat3_Sp1:A_norm | 48/48 | -2.24 | 0.20 | -11.28 | < 0 | agree |
| Repeat3_Sp1 | Repeat3_Sp1:M | 48/48 | 0.47 | 0.14 | 3.44 |  |  |
| Repeat3_Sp1 | msp_C_std | 48/48 | -2.02 | 0.10 | -20.55 | < 0 | agree |
| Repeat3_Sp1 | nuc_85_250_C_std | 48/48 | 1.92 | 0.12 | 16.33 | > 0 | agree |

"~0" judged mechanically as |z| <= 2.

## 6. Primary-family discoveries (primary)

| endpoint | tested | q < 0.05 |
|---|---|---|
| local:A_norm | 627 | 106 |
| flank:A_norm | 627 | 78 |
| msp_C_std | 630 | 80 |
| nuc_85_250_C_std | 630 | 39 |

| source_variant | endpoint | delta_pp | z | q_bh_primary_family | activity_LDLR | activity_LDLR_2 |
|---|---|---|---|---|---|---|
| 157A>G | msp_C_std | 7.71 | 11.77 | 1.36e-28 | 0.43 | 0.5 |
| 157A>G | local:A_norm | 8.73 | 9.33 | 1.34e-17 | 0.43 | 0.5 |
| 147T>C | msp_C_std | -4.73 | -8.02 | 9.16e-13 | -2.36 | -2.75 |
| 145C>T | flank:A_norm | -6.21 | -7.88 | 2.02e-12 | -2.87 | -3.55 |
| 157A>G | flank:A_norm | 7.09 | 7.69 | 6.14e-12 | 0.43 | 0.5 |
| 157A>G | nuc_85_250_C_std | -5.31 | -7.71 | 6.14e-12 | 0.43 | 0.5 |
| 191G>T | local:A_norm | 5.13 | 7.37 | 6.3e-11 | 0.11 | 0.1 |
| 145C>T | msp_C_std | -5.11 | -7.31 | 8.21e-11 | -2.87 | -3.55 |
| 146C>T | msp_C_std | -4.28 | -7.28 | 9.28e-11 | -2.74 | -3.42 |
| 148C>A | msp_C_std | -4.19 | -7.17 | 1.88e-10 | -2.38 | -3.48 |
| 146C>G | msp_C_std | -4.00 | -7.02 | 5.25e-10 | -2.7 | -3.31 |
| 80G>C | msp_C_std | 3.77 | 6.43 | 2.61e-08 | 0.03 | 0.17 |
| 98C>A | msp_C_std | -3.29 | -6.34 | 4.03e-08 | -2.43 | -2.38 |
| 149C>G | msp_C_std | -2.93 | -6.35 | 4.03e-08 | -2.57 | -3.06 |
| 146C>T | nuc_85_250_C_std | 4.15 | 6.30 | 4.86e-08 | -2.74 | -3.42 |
| 135A>C | local:A_norm | -6.31 | -6.27 | 5.57e-08 | -2.37 | -2.56 |
| 143C>T | msp_C_std | -4.16 | -6.08 | 1.66e-07 | -2.28 | -3.09 |
| 148C>T | flank:A_norm | -4.69 | -6.09 | 1.66e-07 | -2.59 | -3.52 |
| 56G>T | local:A_norm | 6.20 | 5.98 | 2.95e-07 | -1.18 | -1.32 |
| 147T>C | flank:A_norm | -4.39 | -5.87 | 5.52e-07 | -2.36 | -2.75 |

## 7. Library association

### primary

| endpoint | source | variants | variant_pearson | positions | position_pearson | null_abs_max | p_circular_shift |
|---|---|---|---|---|---|---|---|
| local:A_norm | LDLR | 628 | 0.244 | 211 | 0.258 | 0.247 | 0.005 |
| local:A_norm | LDLR_2 | 628 | 0.254 | 211 | 0.283 | 0.236 | 0.005 |
| flank:A_norm | LDLR | 628 | 0.354 | 211 | 0.374 | 0.307 | 0.005 |
| flank:A_norm | LDLR_2 | 628 | 0.370 | 211 | 0.409 | 0.299 | 0.005 |
| msp_C_std | LDLR | 630 | 0.500 | 211 | 0.555 | 0.308 | 0.005 |
| msp_C_std | LDLR_2 | 630 | 0.520 | 211 | 0.585 | 0.311 | 0.005 |
| nuc_85_250_C_std | LDLR | 630 | -0.477 | 211 | -0.568 | 0.295 | 0.005 |
| nuc_85_250_C_std | LDLR_2 | 630 | -0.500 | 211 | -0.603 | 0.301 | 0.005 |
| msp_ALL | LDLR | 633 | 0.151 | 211 | 0.209 | 0.269 | 0.036 |
| msp_ALL | LDLR_2 | 633 | 0.161 | 211 | 0.225 | 0.278 | 0.031 |
| nuc_85_250_ALL | LDLR | 633 | -0.355 | 211 | -0.475 | 0.233 | 0.005 |
| nuc_85_250_ALL | LDLR_2 | 633 | -0.371 | 211 | -0.499 | 0.259 | 0.005 |
| local:C_std | LDLR | 630 | 0.070 | 211 | 0.048 | 0.326 | 0.756 |
| local:C_std | LDLR_2 | 630 | 0.084 | 211 | 0.070 | 0.315 | 0.585 |
| local:M | LDLR | 630 | -0.139 | 211 | -0.177 | 0.294 | 0.088 |
| local:M | LDLR_2 | 630 | -0.128 | 211 | -0.164 | 0.279 | 0.124 |
| element:A_norm | LDLR | 628 | 0.486 | 211 | 0.546 | 0.307 | 0.005 |
| element:A_norm | LDLR_2 | 628 | 0.501 | 211 | 0.574 | 0.305 | 0.005 |

- primary / LDLR / variant_pearson: r(msp_C_std)=0.500 > r(msp_ALL)=0.151 > 0: True; r(nuc_85_250_C_std)=-0.477 < 0: True
- primary / LDLR / position_pearson: r(msp_C_std)=0.555 > r(msp_ALL)=0.209 > 0: True; r(nuc_85_250_C_std)=-0.568 < 0: True
- primary / LDLR_2 / variant_pearson: r(msp_C_std)=0.520 > r(msp_ALL)=0.161 > 0: True; r(nuc_85_250_C_std)=-0.500 < 0: True
- primary / LDLR_2 / position_pearson: r(msp_C_std)=0.585 > r(msp_ALL)=0.225 > 0: True; r(nuc_85_250_C_std)=-0.603 < 0: True

### effects_supported_barcodes

| endpoint | source | variants | variant_pearson | positions | position_pearson | null_abs_max | p_circular_shift |
|---|---|---|---|---|---|---|---|
| local:A_norm | LDLR | 627 | 0.245 | 211 | 0.262 | 0.237 | 0.005 |
| local:A_norm | LDLR_2 | 627 | 0.253 | 211 | 0.286 | 0.233 | 0.005 |
| flank:A_norm | LDLR | 627 | 0.353 | 211 | 0.372 | 0.291 | 0.005 |
| flank:A_norm | LDLR_2 | 627 | 0.367 | 211 | 0.405 | 0.283 | 0.005 |
| msp_C_std | LDLR | 630 | 0.492 | 211 | 0.549 | 0.284 | 0.005 |
| msp_C_std | LDLR_2 | 630 | 0.509 | 211 | 0.578 | 0.289 | 0.005 |
| nuc_85_250_C_std | LDLR | 630 | -0.465 | 211 | -0.562 | 0.305 | 0.005 |
| nuc_85_250_C_std | LDLR_2 | 630 | -0.487 | 211 | -0.596 | 0.315 | 0.005 |
| msp_ALL | LDLR | 632 | 0.298 | 211 | 0.378 | 0.198 | 0.005 |
| msp_ALL | LDLR_2 | 632 | 0.314 | 211 | 0.409 | 0.204 | 0.005 |
| nuc_85_250_ALL | LDLR | 632 | -0.409 | 211 | -0.522 | 0.255 | 0.005 |
| nuc_85_250_ALL | LDLR_2 | 632 | -0.429 | 211 | -0.549 | 0.262 | 0.005 |
| local:C_std | LDLR | 630 | 0.076 | 211 | 0.053 | 0.316 | 0.731 |
| local:C_std | LDLR_2 | 630 | 0.088 | 211 | 0.074 | 0.306 | 0.570 |
| local:M | LDLR | 630 | -0.139 | 211 | -0.178 | 0.297 | 0.088 |
| local:M | LDLR_2 | 630 | -0.126 | 211 | -0.165 | 0.279 | 0.119 |
| element:A_norm | LDLR | 627 | 0.475 | 211 | 0.536 | 0.296 | 0.005 |
| element:A_norm | LDLR_2 | 627 | 0.488 | 211 | 0.563 | 0.294 | 0.005 |

- effects_supported_barcodes / LDLR / variant_pearson: r(msp_C_std)=0.492 > r(msp_ALL)=0.298 > 0: True; r(nuc_85_250_C_std)=-0.465 < 0: True
- effects_supported_barcodes / LDLR / position_pearson: r(msp_C_std)=0.549 > r(msp_ALL)=0.378 > 0: True; r(nuc_85_250_C_std)=-0.562 < 0: True
- effects_supported_barcodes / LDLR_2 / variant_pearson: r(msp_C_std)=0.509 > r(msp_ALL)=0.314 > 0: True; r(nuc_85_250_C_std)=-0.487 < 0: True
- effects_supported_barcodes / LDLR_2 / position_pearson: r(msp_C_std)=0.578 > r(msp_ALL)=0.409 > 0: True; r(nuc_85_250_C_std)=-0.596 < 0: True

## 8. Observability regression (primary)

| window | variants | pearson_M_vs_C | slope_C_on_M | sd_delta_M | sd_delta_C_std |
|---|---|---|---|---|---|
| near | 585 | 0.9123 | 0.4575 | 0.2239 | 0.1123 |
| local | 630 | 0.7168 | 0.4281 | 0.0392 | 0.0234 |
| flank | 630 | 0.3344 | 0.4032 | 0.0101 | 0.0122 |
| element | 630 | 0.3693 | 0.4858 | 0.0065 | 0.0086 |
| WT_P1 | 630 | 0.5537 | 0.4290 | 0.0129 | 0.0100 |
| WT_P2 | 630 | 0.4473 | 0.2822 | 0.0201 | 0.0127 |
| WT_P3 | 630 | 0.5824 | 0.4757 | 0.0124 | 0.0101 |
| WT_P4 | 630 | 0.4115 | 0.3402 | 0.0172 | 0.0142 |
| adjacent_143_150 | 624 | 0.6002 | 0.6302 | 0.0189 | 0.0198 |
| P1_core | 630 | 0.7620 | 0.6856 | 0.0128 | 0.0115 |
| Repeat1_Sp1 | 630 | 0.3806 | 0.3422 | 0.0154 | 0.0139 |
| Repeat2_SRE1 | 630 | 0.5709 | 0.3812 | 0.0157 | 0.0105 |
| Repeat3_Sp1 | 630 | 0.3315 | 0.2228 | 0.0151 | 0.0102 |


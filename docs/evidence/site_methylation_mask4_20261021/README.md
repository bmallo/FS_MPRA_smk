# Site-methylation mask-4 follow-up, replicate 3 (run 2026-10-20)

Prespecified rerun of the analysis stage only (`docs/rep3_site_methylation_20261021.md`, "Follow-up run, prespecified before it is executed"; instructions `docs/klone_agent_prompt_site_methylation_mask4.md`). Numbers only; rules applied mechanically.

## Provenance

- Git commit: `f8840e8` (branch `marginal-null-model`); `git diff fc4731a -- config/ldlr_site_methylation.json fsmpra/site_extract.py fsmpra/io.py fsmpra/processing.py` empty.
- Config: `config/ldlr_site_methylation_mask4.json` (mask ±4; near 1–4, local 5–10, flank 11–30; 1,000 WT pseudo-variants; `min_association_molecules` 300).
- Extraction reused: `results/rep3_site_methylation_20261021/extract` (job 42090780; audit `complete` True, 2,659,690 molecules) → `extract_audit.json`.
- Sources → bundle subdirectories:
  - `results/rep3_site_methylation_20261021/mask4/effects` → `effects/` (primary, U ≤ 0.01, M ≥ 0.60)
  - `.../mask4/sensitivity/u0.02_m0.7` → `u0.02_m0.7/`
  - `.../mask4/sensitivity/u0.005_m0.5` → `u0.005_m0.5/`
  - `.../mask4/effects_supported_barcodes` → `effects_supported_barcodes/`
- Slurm job 42093363 (`--account=stergachislab --partition=cpu-g2`, 1 CPU, 16G, threads = 1), node n3433, 6:39 elapsed, sacct MaxRSS 8,578,848 KB. Tests: `Ran 133 tests`, `OK (skipped=1)`. Log: `results/site_methylation_jobs/mask4_42093363.log`. A duplicate submission (42093366) created by an `sbatch` socket timeout was cancelled before it started.

| Run | Wall clock | Max RSS (KB) |
|---|---|---|
| effects (primary) | 1:15.46 | 8,539,944 |
| sensitivity u0.02_m0.7 | 1:14.32 | 8,541,756 |
| sensitivity u0.005_m0.5 | 1:13.50 | 8,539,976 |
| effects_supported_barcodes | 0:50.50 | 7,357,744 |

Conventions: deltas/SEs in pp (×100). Region membership `a ≤ oligo_pos1 − 21 < b` on `regions0`; element bins in table 5 are inclusive oligo positions. Rule-3 SE uses sample variance of per-SNV deltas, SNVs with ≥ 300 molecules. `p_circular_shift` 0.005 is the floor for 200 rotations.

## Rule outcomes (mechanical)

- Rule 1 (calibration): primary 73/73 non-U endpoints pass `gate_pass`, all four primary endpoints included. u0.02_m0.7 72/73 (fails: `adjacent_143_150:M`); u0.005_m0.5 73/73; supported barcodes 73/73.
- Rule 2 (homopolymer reach): triggered. Homopolymer-context M kernel has |z| > 2 at |offset| ≥ 5 for A/T gain at offsets −14 (z −2.18), −8 (z −2.04), −5 (z −3.10). Local endpoints for homopolymer-context SNVs are "unresolved".
- Rule 3 (residual artifact): `local:A_norm` gain − loss = −0.59 pp (SE 0.22, z −2.72; gate pass) → usable (|z| < 3). Context "other" only: −0.24 pp (SE 0.26, z −0.91) → usable.
- Rule 4 (position 125): observability clause met (3/3). Chromatin clause not met (`WT_P3:A_norm` z ≤ 3 for 125A>C primary 2.76, u0.02_m0.7 3.00, u0.005_m0.5 2.37; 125A>G u0.005_m0.5 2.44). All three 125 substitutions are homopolymer context.
- Library association: r(msp_C_std) > r(msp_ALL) > 0 and r(nuc_85_250_C_std) < 0 hold in all 8 source × level × cohort combinations.

## Step-2 tables
Primary run: mask_radius 4, windows {'near': [1, 4], 'local': [5, 10], 'flank': [11, 30]}, aa_draws 1000, thresholds U<=0.01 M>=0.6

molecules: {"cohort": 2594368, "classified": 2594368, "wt": {"U": 93030, "M": 266758, "C": 884103}, "snv": {"U": 101223, "M": 289915, "C": 959339}}

class_dispersion: {"fraction_U_pearson_dispersion": 1.0232669623288997, "fraction_M_pearson_dispersion": 1.0779811249613078}

## 1. Calibration (rule 1)

### primary u0.01_m0.6: gate_pass 73 / 73 non-U endpoints

| endpoint | draws | sd_z | rate_abs_z_gt_1_96 | rate_ci_high | rate_lower_tail | rate_upper_tail | gate_pass |
|---|---|---|---|---|---|---|---|
| flank:A_norm | 1000 | 0.999 | 0.050 | 0.065 | 0.027 | 0.023 | True |
| local:A_norm | 996 | 0.971 | 0.044 | 0.059 | 0.025 | 0.019 | True |
| msp_C_std | 1000 | 0.993 | 0.050 | 0.065 | 0.023 | 0.027 | True |
| nuc_85_250_C_std | 1000 | 1.029 | 0.059 | 0.075 | 0.033 | 0.026 | True |

Failing non-U endpoints:

none

### u0.02_m0.7: gate_pass 72 / 73 non-U endpoints

| endpoint | draws | sd_z | rate_abs_z_gt_1_96 | rate_ci_high | rate_lower_tail | rate_upper_tail | gate_pass |
|---|---|---|---|---|---|---|---|
| flank:A_norm | 1000 | 1.003 | 0.050 | 0.065 | 0.031 | 0.019 | True |
| local:A_norm | 996 | 0.973 | 0.036 | 0.050 | 0.023 | 0.013 | True |
| msp_C_std | 1000 | 0.995 | 0.052 | 0.068 | 0.026 | 0.026 | True |
| nuc_85_250_C_std | 1000 | 1.034 | 0.064 | 0.081 | 0.037 | 0.027 | True |

Failing non-U endpoints:

| endpoint | draws | sd_z | rate_abs_z_gt_1_96 | rate_ci_high | rate_lower_tail | rate_upper_tail | gate_pass |
|---|---|---|---|---|---|---|---|
| adjacent_143_150:M | 922 | 1.236 | 0.104 | 0.126 | 0.001 | 0.103 | False |

### u0.005_m0.5: gate_pass 73 / 73 non-U endpoints

| endpoint | draws | sd_z | rate_abs_z_gt_1_96 | rate_ci_high | rate_lower_tail | rate_upper_tail | gate_pass |
|---|---|---|---|---|---|---|---|
| flank:A_norm | 1000 | 0.984 | 0.046 | 0.061 | 0.026 | 0.020 | True |
| local:A_norm | 996 | 0.991 | 0.044 | 0.059 | 0.021 | 0.023 | True |
| msp_C_std | 1000 | 1.005 | 0.052 | 0.068 | 0.030 | 0.022 | True |
| nuc_85_250_C_std | 1000 | 1.034 | 0.067 | 0.084 | 0.038 | 0.029 | True |

Failing non-U endpoints:

none

### effects_supported_barcodes: gate_pass 73 / 73 non-U endpoints

| endpoint | draws | sd_z | rate_abs_z_gt_1_96 | rate_ci_high | rate_lower_tail | rate_upper_tail | gate_pass |
|---|---|---|---|---|---|---|---|
| flank:A_norm | 1000 | 0.982 | 0.048 | 0.063 | 0.025 | 0.023 | True |
| local:A_norm | 995 | 1.006 | 0.048 | 0.063 | 0.020 | 0.028 | True |
| msp_C_std | 1000 | 1.002 | 0.054 | 0.070 | 0.022 | 0.032 | True |
| nuc_85_250_C_std | 1000 | 1.022 | 0.064 | 0.081 | 0.027 | 0.037 | True |

Failing non-U endpoints:

none

## 2. Homopolymer reach (rule 2; primary kernel, offsets -14..+14, pp)

### M / gain / homopolymer

| offset | pairs | mean_delta_pp | se_pp | z | |z|>2 & |off|>=5 |
|---|---|---|---|---|---|
| -14 | 12 | -1.87 | 0.86 | -2.18 | * |
| -13 | 6 | -0.21 | 0.79 | -0.26 |  |
| -12 | 8 | -2.23 | 1.32 | -1.69 |  |
| -11 | 12 | 0.97 | 0.87 | 1.12 |  |
| -10 | 4 | 0.02 | 1.51 | 0.02 |  |
| -9 | 6 | -0.73 | 1.23 | -0.59 |  |
| -8 | 6 | -1.69 | 0.83 | -2.04 | * |
| -7 | 10 | -1.08 | 0.61 | -1.79 |  |
| -6 | 10 | -1.04 | 0.97 | -1.07 |  |
| -5 | 8 | -6.07 | 1.96 | -3.10 | * |
| -4 | 6 | -5.30 | 3.85 | -1.38 |  |
| -3 | 20 | -10.46 | 2.59 | -4.04 |  |
| -2 | 22 | -16.02 | 4.03 | -3.98 |  |
| -1 | 18 | -46.06 | 7.51 | -6.13 |  |
| 1 | 18 | -41.08 | 7.59 | -5.41 |  |
| 2 | 20 | -15.88 | 4.47 | -3.56 |  |
| 3 | 16 | -10.29 | 2.59 | -3.98 |  |
| 4 | 12 | -6.35 | 2.68 | -2.37 |  |
| 5 | 10 | -1.08 | 1.99 | -0.54 |  |
| 6 | 12 | -1.16 | 0.78 | -1.50 |  |
| 7 | 12 | 0.56 | 1.48 | 0.38 |  |
| 8 | 10 | -0.27 | 0.27 | -0.99 |  |
| 9 | 8 | 4.10 | 2.51 | 1.63 |  |
| 10 | 14 | -3.89 | 2.08 | -1.87 |  |
| 11 | 12 | 1.04 | 1.10 | 0.95 |  |
| 12 | 12 | 0.35 | 0.81 | 0.44 |  |
| 13 | 12 | 0.38 | 0.94 | 0.40 |  |
| 14 | 10 | -1.29 | 1.01 | -1.28 |  |

### M / gain / other

| offset | pairs | mean_delta_pp | se_pp | z | |z|>2 & |off|>=5 |
|---|---|---|---|---|---|
| -14 | 86 | 0.12 | 0.26 | 0.47 |  |
| -13 | 88 | -0.11 | 0.29 | -0.40 |  |
| -12 | 91 | -0.25 | 0.29 | -0.88 |  |
| -11 | 94 | -0.33 | 0.42 | -0.80 |  |
| -10 | 97 | -0.20 | 0.44 | -0.45 |  |
| -9 | 96 | 0.69 | 0.60 | 1.15 |  |
| -8 | 93 | 2.17 | 0.93 | 2.34 | * |
| -7 | 101 | 0.63 | 0.61 | 1.04 |  |
| -6 | 97 | 1.71 | 0.72 | 2.36 | * |
| -5 | 88 | 1.63 | 0.71 | 2.29 | * |
| -4 | 103 | -2.64 | 0.90 | -2.92 |  |
| -3 | 81 | -0.99 | 1.10 | -0.90 |  |
| -2 | 79 | -8.81 | 2.08 | -4.24 |  |
| -1 | 91 | -28.88 | 3.25 | -8.89 |  |
| 1 | 92 | -31.92 | 3.33 | -9.59 |  |
| 2 | 81 | -8.07 | 2.30 | -3.50 |  |
| 3 | 86 | -2.52 | 1.03 | -2.44 |  |
| 4 | 95 | -2.66 | 0.86 | -3.08 |  |
| 5 | 83 | 0.94 | 1.00 | 0.95 |  |
| 6 | 95 | 0.62 | 0.56 | 1.11 |  |
| 7 | 97 | 1.11 | 0.95 | 1.17 |  |
| 8 | 85 | 0.49 | 0.74 | 0.67 |  |
| 9 | 92 | 0.97 | 0.76 | 1.28 |  |
| 10 | 87 | 0.22 | 0.51 | 0.44 |  |
| 11 | 94 | -0.24 | 0.29 | -0.81 |  |
| 12 | 88 | -0.28 | 0.24 | -1.20 |  |
| 13 | 79 | -0.13 | 0.21 | -0.64 |  |
| 14 | 86 | 0.01 | 0.29 | 0.04 |  |

### M / loss / homopolymer

| offset | pairs | mean_delta_pp | se_pp | z | |z|>2 & |off|>=5 |
|---|---|---|---|---|---|
| -14 | 40 | -1.01 | 0.63 | -1.61 |  |
| -13 | 37 | -0.39 | 0.43 | -0.90 |  |
| -12 | 38 | -0.60 | 0.70 | -0.85 |  |
| -11 | 32 | 1.36 | 0.81 | 1.69 |  |
| -10 | 31 | 0.91 | 0.84 | 1.07 |  |
| -9 | 29 | -0.16 | 0.94 | -0.17 |  |
| -8 | 30 | -3.41 | 2.29 | -1.49 |  |
| -7 | 29 | -1.33 | 0.97 | -1.37 |  |
| -6 | 36 | 0.99 | 0.89 | 1.11 |  |
| -5 | 39 | 1.24 | 0.70 | 1.78 |  |
| -4 | 39 | 4.43 | 0.99 | 4.48 |  |
| -3 | 39 | 8.74 | 1.45 | 6.04 |  |
| -2 | 50 | 7.40 | 2.31 | 3.20 |  |
| -1 | 58 | 44.00 | 4.26 | 10.33 |  |
| 1 | 56 | 28.98 | 3.73 | 7.77 |  |
| 2 | 48 | 10.78 | 2.04 | 5.27 |  |
| 3 | 38 | 5.27 | 1.25 | 4.23 |  |
| 4 | 30 | 2.45 | 0.91 | 2.69 |  |
| 5 | 36 | 0.23 | 1.10 | 0.21 |  |
| 6 | 33 | -1.88 | 1.42 | -1.32 |  |
| 7 | 33 | -3.34 | 2.49 | -1.34 |  |
| 8 | 34 | 0.55 | 1.25 | 0.44 |  |
| 9 | 37 | 0.56 | 0.78 | 0.73 |  |
| 10 | 29 | 0.23 | 0.55 | 0.41 |  |
| 11 | 30 | -1.63 | 1.19 | -1.36 |  |
| 12 | 25 | -0.57 | 0.61 | -0.93 |  |
| 13 | 29 | 0.52 | 0.46 | 1.13 |  |
| 14 | 33 | -0.74 | 0.41 | -1.78 |  |

### M / loss / other

| offset | pairs | mean_delta_pp | se_pp | z | |z|>2 & |off|>=5 |
|---|---|---|---|---|---|
| -14 | 52 | -0.34 | 0.40 | -0.84 |  |
| -13 | 58 | 0.58 | 0.55 | 1.07 |  |
| -12 | 52 | 0.16 | 0.61 | 0.26 |  |
| -11 | 52 | -0.13 | 0.64 | -0.21 |  |
| -10 | 58 | -0.41 | 0.53 | -0.78 |  |
| -9 | 62 | -0.77 | 0.66 | -1.17 |  |
| -8 | 66 | -2.91 | 1.73 | -1.68 |  |
| -7 | 54 | -0.47 | 0.74 | -0.64 |  |
| -6 | 52 | -2.40 | 0.93 | -2.57 | * |
| -5 | 62 | -1.64 | 1.11 | -1.48 |  |
| -4 | 48 | -0.51 | 1.17 | -0.44 |  |
| -3 | 56 | -2.75 | 1.69 | -1.63 |  |
| -2 | 48 | 4.67 | 1.82 | 2.57 |  |
| -1 | 32 | 6.28 | 3.87 | 1.62 |  |
| 1 | 34 | 17.60 | 4.55 | 3.87 |  |
| 2 | 50 | 5.92 | 2.36 | 2.51 |  |
| 3 | 58 | 3.77 | 1.64 | 2.30 |  |
| 4 | 56 | 4.51 | 1.18 | 3.81 |  |
| 5 | 66 | -0.64 | 0.88 | -0.74 |  |
| 6 | 54 | -4.40 | 1.33 | -3.32 | * |
| 7 | 50 | 0.04 | 0.53 | 0.07 |  |
| 8 | 60 | -0.28 | 0.79 | -0.35 |  |
| 9 | 54 | -0.62 | 0.92 | -0.67 |  |
| 10 | 60 | -1.07 | 0.44 | -2.46 | * |
| 11 | 54 | 0.04 | 0.45 | 0.09 |  |
| 12 | 64 | -0.01 | 0.42 | -0.03 |  |
| 13 | 66 | -0.40 | 0.34 | -1.19 |  |
| 14 | 58 | -0.30 | 0.35 | -0.84 |  |

### M / retained_count / homopolymer

| offset | pairs | mean_delta_pp | se_pp | z | |z|>2 & |off|>=5 |
|---|---|---|---|---|---|
| -14 | 26 | -1.34 | 0.83 | -1.61 |  |
| -13 | 22 | -1.14 | 0.83 | -1.38 |  |
| -12 | 23 | -0.82 | 0.81 | -1.01 |  |
| -11 | 22 | -0.70 | 0.76 | -0.93 |  |
| -10 | 18 | -0.34 | 0.61 | -0.55 |  |
| -9 | 18 | -0.53 | 0.76 | -0.69 |  |
| -8 | 18 | 1.11 | 4.34 | 0.26 |  |
| -7 | 20 | 1.40 | 3.04 | 0.46 |  |
| -6 | 23 | 1.71 | 1.22 | 1.41 |  |
| -5 | 24 | 0.04 | 1.24 | 0.03 |  |
| -4 | 23 | -0.80 | 1.62 | -0.49 |  |
| -3 | 30 | 2.20 | 1.68 | 1.30 |  |
| -2 | 36 | -1.40 | 2.03 | -0.69 |  |
| -1 | 39 | 11.79 | 6.19 | 1.91 |  |
| 1 | 38 | -7.43 | 5.23 | -1.42 |  |
| 2 | 34 | -2.46 | 2.09 | -1.18 |  |
| 3 | 27 | -0.16 | 2.33 | -0.07 |  |
| 4 | 22 | -2.80 | 1.77 | -1.58 |  |
| 5 | 23 | -0.83 | 1.17 | -0.71 |  |
| 6 | 23 | 0.60 | 0.94 | 0.63 |  |
| 7 | 23 | -2.81 | 2.55 | -1.10 |  |
| 8 | 23 | 0.56 | 1.26 | 0.45 |  |
| 9 | 23 | -0.11 | 1.03 | -0.11 |  |
| 10 | 22 | -1.30 | 1.08 | -1.20 |  |
| 11 | 21 | -0.97 | 1.53 | -0.63 |  |
| 12 | 19 | -0.96 | 0.73 | -1.32 |  |
| 13 | 21 | 0.04 | 0.79 | 0.05 |  |
| 14 | 22 | 0.75 | 0.57 | 1.31 |  |

### M / retained_count / other

| offset | pairs | mean_delta_pp | se_pp | z | |z|>2 & |off|>=5 |
|---|---|---|---|---|---|
| -14 | 69 | -0.30 | 0.30 | -0.99 |  |
| -13 | 73 | 0.05 | 0.39 | 0.12 |  |
| -12 | 72 | -0.15 | 0.49 | -0.30 |  |
| -11 | 73 | -0.44 | 0.60 | -0.74 |  |
| -10 | 78 | -0.38 | 0.87 | -0.43 |  |
| -9 | 79 | 0.34 | 0.81 | 0.42 |  |
| -8 | 80 | -1.18 | 1.21 | -0.97 |  |
| -7 | 78 | -1.48 | 0.88 | -1.69 |  |
| -6 | 75 | 0.00 | 0.91 | 0.00 |  |
| -5 | 75 | -0.78 | 0.96 | -0.81 |  |
| -4 | 76 | 0.13 | 0.88 | 0.15 |  |
| -3 | 69 | -0.26 | 1.25 | -0.21 |  |
| -2 | 64 | 2.01 | 2.14 | 0.94 |  |
| -1 | 62 | -5.13 | 3.12 | -1.64 |  |
| 1 | 63 | 1.22 | 2.93 | 0.42 |  |
| 2 | 66 | 0.05 | 1.43 | 0.03 |  |
| 3 | 72 | 2.68 | 1.09 | 2.45 |  |
| 4 | 76 | 0.23 | 0.67 | 0.34 |  |
| 5 | 75 | 0.27 | 0.79 | 0.34 |  |
| 6 | 75 | -0.31 | 0.63 | -0.49 |  |
| 7 | 74 | 0.95 | 0.66 | 1.44 |  |
| 8 | 73 | 1.33 | 1.30 | 1.02 |  |
| 9 | 73 | 0.70 | 0.69 | 1.02 |  |
| 10 | 74 | -0.73 | 0.41 | -1.79 |  |
| 11 | 74 | -0.06 | 0.41 | -0.14 |  |
| 12 | 76 | -0.11 | 0.28 | -0.41 |  |
| 13 | 73 | -0.05 | 0.26 | -0.18 |  |
| 14 | 72 | -0.05 | 0.27 | -0.18 |  |

### C_std / gain / homopolymer

| offset | pairs | mean_delta_pp | se_pp | z | |z|>2 & |off|>=5 |
|---|---|---|---|---|---|
| -14 | 12 | -1.90 | 0.44 | -4.34 | * |
| -13 | 6 | -1.99 | 0.87 | -2.28 | * |
| -12 | 8 | -2.63 | 0.77 | -3.41 | * |
| -11 | 12 | -0.97 | 0.42 | -2.33 | * |
| -10 | 4 | -2.96 | 1.68 | -1.76 |  |
| -9 | 6 | -0.28 | 2.11 | -0.13 |  |
| -8 | 6 | -1.26 | 0.97 | -1.30 |  |
| -7 | 10 | -3.05 | 1.06 | -2.87 | * |
| -6 | 10 | -1.01 | 0.45 | -2.24 | * |
| -5 | 8 | -2.53 | 0.95 | -2.67 | * |
| -4 | 6 | -3.15 | 3.26 | -0.97 |  |
| -3 | 20 | -5.19 | 1.04 | -4.97 |  |
| -2 | 22 | -7.16 | 1.65 | -4.34 |  |
| -1 | 18 | -19.02 | 3.03 | -6.27 |  |
| 1 | 18 | -16.46 | 2.74 | -6.00 |  |
| 2 | 20 | -7.83 | 2.00 | -3.90 |  |
| 3 | 16 | -5.79 | 1.43 | -4.06 |  |
| 4 | 12 | -5.74 | 1.23 | -4.67 |  |
| 5 | 10 | -1.66 | 0.84 | -1.97 |  |
| 6 | 12 | -1.12 | 0.75 | -1.49 |  |
| 7 | 12 | -0.82 | 1.03 | -0.79 |  |
| 8 | 10 | 0.27 | 0.77 | 0.35 |  |
| 9 | 8 | -1.04 | 1.48 | -0.70 |  |
| 10 | 14 | -2.35 | 0.82 | -2.86 | * |
| 11 | 12 | -0.27 | 0.33 | -0.82 |  |
| 12 | 12 | -1.25 | 0.50 | -2.50 | * |
| 13 | 12 | -0.99 | 0.56 | -1.78 |  |
| 14 | 10 | -0.97 | 0.45 | -2.14 | * |

### C_std / gain / other

| offset | pairs | mean_delta_pp | se_pp | z | |z|>2 & |off|>=5 |
|---|---|---|---|---|---|
| -14 | 86 | -0.22 | 0.13 | -1.65 |  |
| -13 | 88 | -0.63 | 0.16 | -3.88 | * |
| -12 | 91 | -0.26 | 0.17 | -1.49 |  |
| -11 | 94 | -0.56 | 0.21 | -2.69 | * |
| -10 | 97 | -0.54 | 0.24 | -2.22 | * |
| -9 | 96 | 0.32 | 0.28 | 1.16 |  |
| -8 | 93 | 1.02 | 0.51 | 2.00 | * |
| -7 | 101 | 0.19 | 0.30 | 0.63 |  |
| -6 | 97 | 0.79 | 0.37 | 2.16 | * |
| -5 | 88 | 0.41 | 0.44 | 0.95 |  |
| -4 | 103 | -2.03 | 0.44 | -4.65 |  |
| -3 | 81 | -1.23 | 0.58 | -2.13 |  |
| -2 | 79 | -6.38 | 0.87 | -7.31 |  |
| -1 | 91 | -15.72 | 1.24 | -12.72 |  |
| 1 | 92 | -16.07 | 1.20 | -13.35 |  |
| 2 | 81 | -4.45 | 0.86 | -5.18 |  |
| 3 | 86 | -1.99 | 0.53 | -3.78 |  |
| 4 | 95 | -2.41 | 0.47 | -5.18 |  |
| 5 | 83 | -0.03 | 0.47 | -0.05 |  |
| 6 | 95 | 0.15 | 0.28 | 0.52 |  |
| 7 | 97 | 0.48 | 0.40 | 1.18 |  |
| 8 | 85 | 0.40 | 0.41 | 0.97 |  |
| 9 | 92 | 0.17 | 0.32 | 0.53 |  |
| 10 | 87 | -0.39 | 0.26 | -1.47 |  |
| 11 | 94 | -0.79 | 0.18 | -4.36 | * |
| 12 | 88 | -0.68 | 0.18 | -3.86 | * |
| 13 | 79 | -0.38 | 0.21 | -1.82 |  |
| 14 | 86 | -0.33 | 0.19 | -1.75 |  |

### C_std / loss / homopolymer

| offset | pairs | mean_delta_pp | se_pp | z | |z|>2 & |off|>=5 |
|---|---|---|---|---|---|
| -14 | 40 | -0.06 | 0.37 | -0.16 |  |
| -13 | 37 | 0.13 | 0.37 | 0.36 |  |
| -12 | 38 | 0.24 | 0.34 | 0.71 |  |
| -11 | 32 | 0.96 | 0.41 | 2.35 | * |
| -10 | 31 | 1.27 | 0.51 | 2.51 | * |
| -9 | 29 | 0.33 | 0.49 | 0.66 |  |
| -8 | 30 | -0.43 | 1.01 | -0.42 |  |
| -7 | 29 | -0.47 | 0.76 | -0.62 |  |
| -6 | 36 | 0.25 | 0.54 | 0.46 |  |
| -5 | 39 | 0.70 | 0.41 | 1.69 |  |
| -4 | 39 | 3.04 | 0.48 | 6.38 |  |
| -3 | 39 | 4.69 | 0.66 | 7.14 |  |
| -2 | 50 | 4.71 | 0.90 | 5.24 |  |
| -1 | 58 | 21.33 | 1.72 | 12.37 |  |
| 1 | 56 | 15.06 | 1.46 | 10.30 |  |
| 2 | 48 | 6.96 | 0.91 | 7.63 |  |
| 3 | 38 | 3.60 | 0.76 | 4.74 |  |
| 4 | 30 | 2.75 | 0.96 | 2.85 |  |
| 5 | 36 | 1.54 | 0.64 | 2.40 | * |
| 6 | 33 | 0.06 | 0.93 | 0.06 |  |
| 7 | 33 | -0.41 | 1.12 | -0.37 |  |
| 8 | 34 | 0.04 | 0.62 | 0.07 |  |
| 9 | 37 | 1.07 | 0.46 | 2.31 | * |
| 10 | 29 | 0.74 | 0.38 | 1.97 |  |
| 11 | 30 | 0.33 | 0.56 | 0.58 |  |
| 12 | 25 | -0.20 | 0.28 | -0.72 |  |
| 13 | 29 | 0.49 | 0.31 | 1.61 |  |
| 14 | 33 | 0.16 | 0.30 | 0.55 |  |

### C_std / loss / other

| offset | pairs | mean_delta_pp | se_pp | z | |z|>2 & |off|>=5 |
|---|---|---|---|---|---|
| -14 | 52 | 0.11 | 0.23 | 0.47 |  |
| -13 | 58 | 0.88 | 0.29 | 3.00 | * |
| -12 | 52 | 0.17 | 0.33 | 0.50 |  |
| -11 | 52 | 0.34 | 0.26 | 1.31 |  |
| -10 | 58 | 0.48 | 0.34 | 1.41 |  |
| -9 | 62 | -0.56 | 0.36 | -1.55 |  |
| -8 | 66 | -1.04 | 0.75 | -1.39 |  |
| -7 | 54 | -0.33 | 0.49 | -0.68 |  |
| -6 | 52 | -0.78 | 0.54 | -1.43 |  |
| -5 | 62 | -0.41 | 0.63 | -0.65 |  |
| -4 | 48 | 0.60 | 0.61 | 0.98 |  |
| -3 | 56 | 0.57 | 0.77 | 0.74 |  |
| -2 | 48 | 4.78 | 0.95 | 5.06 |  |
| -1 | 32 | 9.26 | 2.00 | 4.63 |  |
| 1 | 34 | 13.13 | 2.02 | 6.51 |  |
| 2 | 50 | 6.22 | 1.22 | 5.11 |  |
| 3 | 58 | 2.36 | 0.69 | 3.40 |  |
| 4 | 56 | 3.08 | 0.54 | 5.73 |  |
| 5 | 66 | -0.57 | 0.41 | -1.39 |  |
| 6 | 54 | -1.87 | 0.53 | -3.52 | * |
| 7 | 50 | -0.33 | 0.34 | -0.97 |  |
| 8 | 60 | -0.51 | 0.44 | -1.15 |  |
| 9 | 54 | -0.37 | 0.37 | -1.01 |  |
| 10 | 60 | -0.20 | 0.23 | -0.84 |  |
| 11 | 54 | 0.70 | 0.32 | 2.15 | * |
| 12 | 64 | 0.42 | 0.26 | 1.61 |  |
| 13 | 66 | 0.21 | 0.23 | 0.93 |  |
| 14 | 58 | -0.09 | 0.22 | -0.41 |  |

### C_std / retained_count / homopolymer

| offset | pairs | mean_delta_pp | se_pp | z | |z|>2 & |off|>=5 |
|---|---|---|---|---|---|
| -14 | 26 | -0.62 | 0.41 | -1.51 |  |
| -13 | 22 | -0.62 | 0.40 | -1.54 |  |
| -12 | 23 | -0.71 | 0.44 | -1.64 |  |
| -11 | 22 | 0.20 | 0.50 | 0.40 |  |
| -10 | 18 | -0.02 | 0.54 | -0.03 |  |
| -9 | 18 | 0.87 | 0.72 | 1.22 |  |
| -8 | 18 | 2.78 | 2.32 | 1.20 |  |
| -7 | 20 | 0.94 | 1.76 | 0.54 |  |
| -6 | 23 | 0.80 | 0.64 | 1.25 |  |
| -5 | 24 | -0.61 | 0.51 | -1.19 |  |
| -4 | 23 | 0.06 | 0.74 | 0.08 |  |
| -3 | 30 | 1.48 | 0.65 | 2.26 |  |
| -2 | 36 | -0.40 | 0.77 | -0.52 |  |
| -1 | 39 | 4.51 | 2.27 | 1.98 |  |
| 1 | 38 | -2.02 | 1.78 | -1.13 |  |
| 2 | 34 | -0.80 | 1.03 | -0.78 |  |
| 3 | 27 | 0.66 | 1.05 | 0.63 |  |
| 4 | 22 | -0.37 | 0.75 | -0.49 |  |
| 5 | 23 | -0.37 | 0.56 | -0.67 |  |
| 6 | 23 | 0.51 | 0.71 | 0.72 |  |
| 7 | 23 | -0.75 | 1.27 | -0.59 |  |
| 8 | 23 | 0.32 | 1.01 | 0.31 |  |
| 9 | 23 | 0.64 | 0.67 | 0.96 |  |
| 10 | 22 | -0.27 | 0.52 | -0.52 |  |
| 11 | 21 | -0.37 | 0.63 | -0.58 |  |
| 12 | 19 | -1.00 | 0.40 | -2.48 | * |
| 13 | 21 | -0.52 | 0.40 | -1.31 |  |
| 14 | 22 | -0.10 | 0.24 | -0.42 |  |

### C_std / retained_count / other

| offset | pairs | mean_delta_pp | se_pp | z | |z|>2 & |off|>=5 |
|---|---|---|---|---|---|
| -14 | 69 | -0.13 | 0.20 | -0.65 |  |
| -13 | 73 | -0.23 | 0.25 | -0.92 |  |
| -12 | 72 | -0.38 | 0.25 | -1.50 |  |
| -11 | 73 | -0.44 | 0.26 | -1.69 |  |
| -10 | 78 | -0.34 | 0.36 | -0.95 |  |
| -9 | 79 | 0.15 | 0.35 | 0.43 |  |
| -8 | 80 | -0.49 | 0.58 | -0.85 |  |
| -7 | 78 | -0.69 | 0.38 | -1.81 |  |
| -6 | 75 | -0.19 | 0.35 | -0.54 |  |
| -5 | 75 | -0.16 | 0.50 | -0.32 |  |
| -4 | 76 | 0.03 | 0.39 | 0.08 |  |
| -3 | 69 | 0.19 | 0.64 | 0.30 |  |
| -2 | 64 | 0.58 | 0.91 | 0.63 |  |
| -1 | 62 | -1.64 | 1.36 | -1.21 |  |
| 1 | 63 | -0.30 | 1.52 | -0.20 |  |
| 2 | 66 | 0.19 | 0.70 | 0.27 |  |
| 3 | 72 | 1.02 | 0.63 | 1.62 |  |
| 4 | 76 | -0.41 | 0.35 | -1.16 |  |
| 5 | 75 | -0.04 | 0.38 | -0.12 |  |
| 6 | 75 | -0.13 | 0.30 | -0.44 |  |
| 7 | 74 | 0.48 | 0.35 | 1.35 |  |
| 8 | 73 | 0.43 | 0.64 | 0.67 |  |
| 9 | 73 | 0.51 | 0.36 | 1.40 |  |
| 10 | 74 | -0.84 | 0.27 | -3.10 | * |
| 11 | 74 | -0.05 | 0.26 | -0.19 |  |
| 12 | 76 | -0.26 | 0.22 | -1.17 |  |
| 13 | 73 | -0.09 | 0.22 | -0.40 |  |
| 14 | 72 | -0.05 | 0.22 | -0.22 |  |

Rule 2: homopolymer-context M rows with |z| > 2 at |offset| >= 5: 3 -> local endpoints for homopolymer-context SNVs UNRESOLVED; gain -14 (z=-2.18); gain -8 (z=-2.04); gain -5 (z=-3.10)

## 3. Residual artifact (rule 3; local:A_norm, SNVs with >= 300 molecules)

| SNV context | n_gain | mean_gain_pp | n_loss | mean_loss_pp | gain-loss_pp | se_pp | z | gate_pass | local:A_norm |
|---|---|---|---|---|---|---|---|---|---|
| all | 217 | -0.25 | 200 | 0.34 | -0.59 | 0.22 | -2.72 | True | usable |
| other | 191 | -0.12 | 124 | 0.12 | -0.24 | 0.26 | -0.91 | True | usable |

## 4. Position 125 (rule 4 = usage-doc rule 3)

### primary u0.01_m0.6

| variant | context | endpoint | sites | delta_pp | se_pp | z | molecules_M |
|---|---|---|---|---|---|---|---|
| 125A>C | homopolymer | WT_P3:M | 6 | 5.09 | 0.77 | 6.61 | 400 |
| 125A>C | homopolymer | WT_P3:C_std | 6 | 4.46 | 0.91 | 4.88 | 400 |
| 125A>C | homopolymer | WT_P3:A_norm | 6 | 3.19 | 1.16 | 2.76 | 400 |
| 125A>C | homopolymer | local:M | 6 | 5.13 | 0.78 | 6.56 | 400 |
| 125A>C | homopolymer | local:C_std | 6 | 4.18 | 0.89 | 4.71 | 400 |
| 125A>C | homopolymer | local:A_norm | 6 | 2.80 | 1.16 | 2.41 | 400 |
| 125A>G | homopolymer | WT_P3:M | 6 | 3.79 | 0.72 | 5.23 | 451 |
| 125A>G | homopolymer | WT_P3:C_std | 6 | 3.90 | 0.73 | 5.33 | 451 |
| 125A>G | homopolymer | WT_P3:A_norm | 6 | 3.00 | 0.99 | 3.02 | 451 |
| 125A>G | homopolymer | local:M | 6 | 3.91 | 0.71 | 5.51 | 451 |
| 125A>G | homopolymer | local:C_std | 6 | 3.26 | 0.70 | 4.66 | 451 |
| 125A>G | homopolymer | local:A_norm | 6 | 2.09 | 0.96 | 2.18 | 451 |
| 125A>T | homopolymer | WT_P3:M | 6 | 6.32 | 0.77 | 8.21 | 430 |
| 125A>T | homopolymer | WT_P3:C_std | 6 | 5.90 | 0.87 | 6.76 | 430 |
| 125A>T | homopolymer | WT_P3:A_norm | 6 | 4.07 | 1.08 | 3.76 | 430 |
| 125A>T | homopolymer | local:M | 6 | 5.81 | 0.83 | 6.97 | 430 |
| 125A>T | homopolymer | local:C_std | 6 | 5.74 | 0.82 | 6.98 | 430 |
| 125A>T | homopolymer | local:A_norm | 6 | 4.14 | 1.09 | 3.78 | 430 |

### u0.02_m0.7

| variant | context | endpoint | sites | delta_pp | se_pp | z | molecules_M |
|---|---|---|---|---|---|---|---|
| 125A>C | homopolymer | WT_P3:M | 6 | 3.82 | 1.30 | 2.93 | 122 |
| 125A>C | homopolymer | WT_P3:C_std | 6 | 4.78 | 0.81 | 5.90 | 122 |
| 125A>C | homopolymer | WT_P3:A_norm | 6 | 3.53 | 1.18 | 3.00 | 122 |
| 125A>C | homopolymer | local:M | 6 | 3.59 | 1.38 | 2.60 | 122 |
| 125A>C | homopolymer | local:C_std | 6 | 4.55 | 0.79 | 5.76 | 122 |
| 125A>C | homopolymer | local:A_norm | 6 | 3.35 | 1.22 | 2.74 | 122 |
| 125A>G | homopolymer | WT_P3:M | 6 | 1.81 | 1.26 | 1.43 | 129 |
| 125A>G | homopolymer | WT_P3:C_std | 6 | 4.14 | 0.64 | 6.43 | 129 |
| 125A>G | homopolymer | WT_P3:A_norm | 6 | 4.05 | 1.14 | 3.55 | 129 |
| 125A>G | homopolymer | local:M | 6 | 2.13 | 1.20 | 1.77 | 129 |
| 125A>G | homopolymer | local:C_std | 6 | 3.61 | 0.62 | 5.83 | 129 |
| 125A>G | homopolymer | local:A_norm | 6 | 3.27 | 1.10 | 2.97 | 129 |
| 125A>T | homopolymer | WT_P3:M | 6 | 4.95 | 1.19 | 4.15 | 127 |
| 125A>T | homopolymer | WT_P3:C_std | 6 | 6.16 | 0.78 | 7.89 | 127 |
| 125A>T | homopolymer | WT_P3:A_norm | 6 | 4.32 | 1.17 | 3.69 | 127 |
| 125A>T | homopolymer | local:M | 6 | 4.74 | 1.24 | 3.82 | 127 |
| 125A>T | homopolymer | local:C_std | 6 | 5.92 | 0.74 | 7.99 | 127 |
| 125A>T | homopolymer | local:A_norm | 6 | 4.22 | 1.20 | 3.51 | 127 |

### u0.005_m0.5

| variant | context | endpoint | sites | delta_pp | se_pp | z | molecules_M |
|---|---|---|---|---|---|---|---|
| 125A>C | homopolymer | WT_P3:M | 6 | 6.04 | 0.92 | 6.56 | 570 |
| 125A>C | homopolymer | WT_P3:C_std | 6 | 4.08 | 0.94 | 4.36 | 570 |
| 125A>C | homopolymer | WT_P3:A_norm | 6 | 2.92 | 1.23 | 2.37 | 570 |
| 125A>C | homopolymer | local:M | 6 | 6.24 | 0.91 | 6.83 | 570 |
| 125A>C | homopolymer | local:C_std | 6 | 3.65 | 0.90 | 4.08 | 570 |
| 125A>C | homopolymer | local:A_norm | 6 | 2.20 | 1.18 | 1.86 | 570 |
| 125A>G | homopolymer | WT_P3:M | 6 | 4.26 | 0.81 | 5.26 | 696 |
| 125A>G | homopolymer | WT_P3:C_std | 6 | 3.45 | 0.76 | 4.54 | 696 |
| 125A>G | homopolymer | WT_P3:A_norm | 6 | 2.66 | 1.09 | 2.44 | 696 |
| 125A>G | homopolymer | local:M | 6 | 4.28 | 0.79 | 5.40 | 696 |
| 125A>G | homopolymer | local:C_std | 6 | 2.74 | 0.73 | 3.75 | 696 |
| 125A>G | homopolymer | local:A_norm | 6 | 1.56 | 1.06 | 1.48 | 696 |
| 125A>T | homopolymer | WT_P3:M | 6 | 6.24 | 0.86 | 7.28 | 648 |
| 125A>T | homopolymer | WT_P3:C_std | 6 | 5.56 | 0.92 | 6.07 | 648 |
| 125A>T | homopolymer | WT_P3:A_norm | 6 | 4.37 | 1.18 | 3.70 | 648 |
| 125A>T | homopolymer | local:M | 6 | 5.96 | 0.85 | 7.00 | 648 |
| 125A>T | homopolymer | local:C_std | 6 | 5.33 | 0.87 | 6.14 | 648 |
| 125A>T | homopolymer | local:A_norm | 6 | 4.13 | 1.18 | 3.51 | 648 |

- 125A>C: WT_P3:M 5.09 pp vs WT_P3:C_std 4.46 pp -> same direction & >= half: True
- 125A>G: WT_P3:M 3.79 pp vs WT_P3:C_std 3.90 pp -> same direction & >= half: True
- 125A>T: WT_P3:M 6.32 pp vs WT_P3:C_std 5.90 pp -> same direction & >= half: True
- Observability clause (>= 2 of 3, primary thresholds): True (3/3; 0 missing)
- Chromatin clause (WT_P3:A_norm positive with z > 3, all 3 x 3): False; primary u0.01_m0.6/125A>C z=2.76, primary u0.01_m0.6/125A>G z=3.02, primary u0.01_m0.6/125A>T z=3.76, u0.02_m0.7/125A>C z=3.00, u0.02_m0.7/125A>G z=3.55, u0.02_m0.7/125A>T z=3.69, u0.005_m0.5/125A>C z=2.37, u0.005_m0.5/125A>G z=2.44, u0.005_m0.5/125A>T z=3.70

## 5. Primary family (primary run)

| endpoint | tested | q<0.05 | up | down | sign agrees w/ activity_LDLR | P1 core 54-62 | Repeat 1 89-104 | SRE-1 124-139 | Repeat 3 140-155 | other pos | homopolymer | other ctx |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| local:A_norm | 624 | 93 | 42 | 51 | 67/90 (0.74) | 6 | 7 | 9 | 15 | 56 | 28 | 65 |
| flank:A_norm | 627 | 75 | 24 | 51 | 57/73 (0.78) | 5 | 7 | 5 | 20 | 38 | 13 | 62 |
| msp_C_std | 630 | 80 | 18 | 62 | 69/79 (0.87) | 2 | 11 | 7 | 25 | 35 | 16 | 64 |
| nuc_85_250_C_std | 630 | 38 | 31 | 7 | 38/38 (1.00) | 0 | 5 | 1 | 21 | 11 | 7 | 31 |

Sign agreement: delta sign equals activity sign (opposite for nuc_85_250_C_std); hits with blank/zero activity excluded from the denominator.

| source_variant | context | endpoint | delta_pp | z | q_bh_primary_family | activity_LDLR | activity_LDLR_2 |
|---|---|---|---|---|---|---|---|
| 157A>G | other | msp_C_std | 7.71 | 11.77 | 1.36e-28 | 0.43 | 0.5 |
| 147T>C | other | msp_C_std | -4.73 | -8.02 | 1.37e-12 | -2.36 | -2.75 |
| 145C>T | other | flank:A_norm | -6.21 | -7.88 | 2.69e-12 | -2.87 | -3.55 |
| 157A>G | other | flank:A_norm | 7.09 | 7.69 | 7.36e-12 | 0.43 | 0.5 |
| 157A>G | other | nuc_85_250_C_std | -5.31 | -7.71 | 7.36e-12 | 0.43 | 0.5 |
| 157A>G | other | local:A_norm | 6.75 | 7.56 | 1.63e-11 | 0.43 | 0.5 |
| 191G>T | other | local:A_norm | 5.13 | 7.37 | 6.29e-11 | 0.11 | 0.1 |
| 145C>T | other | msp_C_std | -5.11 | -7.31 | 8.2e-11 | -2.87 | -3.55 |
| 146C>T | other | msp_C_std | -4.28 | -7.28 | 9.27e-11 | -2.74 | -3.42 |
| 148C>A | other | msp_C_std | -4.19 | -7.17 | 1.88e-10 | -2.38 | -3.48 |
| 146C>G | other | msp_C_std | -4.00 | -7.02 | 5.24e-10 | -2.7 | -3.31 |
| 80G>C | other | msp_C_std | 3.77 | 6.43 | 2.61e-08 | 0.03 | 0.17 |
| 98C>A | other | msp_C_std | -3.29 | -6.34 | 4.03e-08 | -2.43 | -2.38 |
| 149C>G | other | msp_C_std | -2.93 | -6.35 | 4.03e-08 | -2.57 | -3.06 |
| 146C>T | other | nuc_85_250_C_std | 4.15 | 6.30 | 4.85e-08 | -2.74 | -3.42 |
| 135A>C | other | local:A_norm | -6.31 | -6.27 | 5.56e-08 | -2.37 | -2.56 |
| 143C>T | homopolymer | msp_C_std | -4.16 | -6.08 | 1.65e-07 | -2.28 | -3.09 |
| 148C>T | other | flank:A_norm | -4.69 | -6.09 | 1.65e-07 | -2.59 | -3.52 |
| 147T>C | other | flank:A_norm | -4.39 | -5.87 | 5.8e-07 | -2.36 | -2.75 |
| 149C>G | other | flank:A_norm | -3.60 | -5.78 | 9.15e-07 | -2.57 | -3.06 |

## 6. Library association

### primary

| endpoint | source | variants | variant_pearson | positions | position_pearson | null_abs_max | p_circular_shift |
|---|---|---|---|---|---|---|---|
| local:A_norm | LDLR | 625 | 0.237 | 210 | 0.258 | 0.233 | 0.005 |
| local:A_norm | LDLR_2 | 625 | 0.239 | 210 | 0.276 | 0.225 | 0.005 |
| flank:A_norm | LDLR | 628 | 0.354 | 211 | 0.374 | 0.307 | 0.005 |
| flank:A_norm | LDLR_2 | 628 | 0.370 | 211 | 0.409 | 0.299 | 0.005 |
| msp_C_std | LDLR | 630 | 0.500 | 211 | 0.555 | 0.308 | 0.005 |
| msp_C_std | LDLR_2 | 630 | 0.520 | 211 | 0.585 | 0.311 | 0.005 |
| nuc_85_250_C_std | LDLR | 630 | -0.477 | 211 | -0.568 | 0.295 | 0.005 |
| nuc_85_250_C_std | LDLR_2 | 630 | -0.500 | 211 | -0.603 | 0.301 | 0.005 |
| msp_ALL | LDLR | 630 | 0.336 | 211 | 0.418 | 0.220 | 0.005 |
| msp_ALL | LDLR_2 | 630 | 0.355 | 211 | 0.444 | 0.219 | 0.005 |
| nuc_85_250_ALL | LDLR | 630 | -0.456 | 211 | -0.553 | 0.254 | 0.005 |
| nuc_85_250_ALL | LDLR_2 | 630 | -0.477 | 211 | -0.585 | 0.252 | 0.005 |
| element:C_std | LDLR | 630 | 0.500 | 211 | 0.542 | 0.326 | 0.005 |
| element:C_std | LDLR_2 | 630 | 0.518 | 211 | 0.574 | 0.329 | 0.005 |
| local:C_std | LDLR | 627 | 0.047 | 210 | 0.026 | 0.323 | 0.859 |
| local:C_std | LDLR_2 | 627 | 0.048 | 210 | 0.035 | 0.302 | 0.771 |

- primary / LDLR / variant_pearson: r(msp_C_std)=0.500 > r(msp_ALL)=0.336 > 0: True; r(nuc_85_250_C_std)=-0.477 < 0: True
- primary / LDLR / position_pearson: r(msp_C_std)=0.555 > r(msp_ALL)=0.418 > 0: True; r(nuc_85_250_C_std)=-0.568 < 0: True
- primary / LDLR_2 / variant_pearson: r(msp_C_std)=0.520 > r(msp_ALL)=0.355 > 0: True; r(nuc_85_250_C_std)=-0.500 < 0: True
- primary / LDLR_2 / position_pearson: r(msp_C_std)=0.585 > r(msp_ALL)=0.444 > 0: True; r(nuc_85_250_C_std)=-0.603 < 0: True

### effects_supported_barcodes

| endpoint | source | variants | variant_pearson | positions | position_pearson | null_abs_max | p_circular_shift |
|---|---|---|---|---|---|---|---|
| local:A_norm | LDLR | 624 | 0.238 | 210 | 0.261 | 0.244 | 0.005 |
| local:A_norm | LDLR_2 | 624 | 0.237 | 210 | 0.278 | 0.239 | 0.005 |
| flank:A_norm | LDLR | 627 | 0.353 | 211 | 0.372 | 0.291 | 0.005 |
| flank:A_norm | LDLR_2 | 627 | 0.367 | 211 | 0.405 | 0.283 | 0.005 |
| msp_C_std | LDLR | 629 | 0.492 | 211 | 0.548 | 0.288 | 0.005 |
| msp_C_std | LDLR_2 | 629 | 0.510 | 211 | 0.577 | 0.293 | 0.005 |
| nuc_85_250_C_std | LDLR | 629 | -0.466 | 211 | -0.563 | 0.306 | 0.005 |
| nuc_85_250_C_std | LDLR_2 | 629 | -0.488 | 211 | -0.597 | 0.315 | 0.005 |
| msp_ALL | LDLR | 629 | 0.323 | 211 | 0.405 | 0.205 | 0.005 |
| msp_ALL | LDLR_2 | 629 | 0.343 | 211 | 0.437 | 0.206 | 0.005 |
| nuc_85_250_ALL | LDLR | 629 | -0.447 | 211 | -0.547 | 0.252 | 0.005 |
| nuc_85_250_ALL | LDLR_2 | 629 | -0.465 | 211 | -0.579 | 0.258 | 0.005 |
| element:C_std | LDLR | 629 | 0.494 | 211 | 0.536 | 0.315 | 0.005 |
| element:C_std | LDLR_2 | 629 | 0.512 | 211 | 0.567 | 0.318 | 0.005 |
| local:C_std | LDLR | 626 | 0.053 | 210 | 0.031 | 0.314 | 0.807 |
| local:C_std | LDLR_2 | 626 | 0.053 | 210 | 0.040 | 0.295 | 0.729 |

- effects_supported_barcodes / LDLR / variant_pearson: r(msp_C_std)=0.492 > r(msp_ALL)=0.323 > 0: True; r(nuc_85_250_C_std)=-0.466 < 0: True
- effects_supported_barcodes / LDLR / position_pearson: r(msp_C_std)=0.548 > r(msp_ALL)=0.405 > 0: True; r(nuc_85_250_C_std)=-0.563 < 0: True
- effects_supported_barcodes / LDLR_2 / variant_pearson: r(msp_C_std)=0.510 > r(msp_ALL)=0.343 > 0: True; r(nuc_85_250_C_std)=-0.488 < 0: True
- effects_supported_barcodes / LDLR_2 / position_pearson: r(msp_C_std)=0.577 > r(msp_ALL)=0.437 > 0: True; r(nuc_85_250_C_std)=-0.597 < 0: True

## 7. Regions (primary; mean per-SNV delta, pp; SNV in region if a <= oligo_pos1-21 < b)

| region | endpoint | SNVs (rows/with value) | sites (min-max) | SNVs with 0 sites | mean_delta_pp | se_pp | z |
|---|---|---|---|---|---|---|---|
| P1_core | P1_core:A_norm | 27/16 | 0-4 | 3 | 2.84 | 1.25 | 2.28 |
| P1_core | P1_core:M | 27/24 | 0-4 | 3 | 2.32 | 0.21 | 10.89 |
| P1_core | msp_C_std | 27/27 | 211-211 | 0 | -0.39 | 0.14 | -2.79 |
| P1_core | nuc_85_250_C_std | 27/27 | 211-211 | 0 | 0.25 | 0.16 | 1.55 |
| Repeat1_Sp1 | Repeat1_Sp1:A_norm | 48/48 | 3-6 | 0 | -0.84 | 0.20 | -4.27 |
| Repeat1_Sp1 | Repeat1_Sp1:M | 48/48 | 3-6 | 0 | -0.52 | 0.16 | -3.16 |
| Repeat1_Sp1 | msp_C_std | 48/48 | 211-211 | 0 | -1.00 | 0.11 | -9.45 |
| Repeat1_Sp1 | nuc_85_250_C_std | 48/48 | 211-211 | 0 | 0.85 | 0.12 | 6.89 |
| Repeat2_SRE1 | Repeat2_SRE1:A_norm | 48/48 | 2-6 | 0 | 0.73 | 0.20 | 3.76 |
| Repeat2_SRE1 | Repeat2_SRE1:M | 48/48 | 2-6 | 0 | 1.63 | 0.13 | 12.29 |
| Repeat2_SRE1 | msp_C_std | 48/48 | 211-211 | 0 | 0.24 | 0.10 | 2.41 |
| Repeat2_SRE1 | nuc_85_250_C_std | 48/48 | 211-211 | 0 | -0.18 | 0.12 | -1.54 |
| Repeat3_Sp1 | Repeat3_Sp1:A_norm | 48/48 | 1-5 | 0 | -1.90 | 0.21 | -9.00 |
| Repeat3_Sp1 | Repeat3_Sp1:M | 48/48 | 1-5 | 0 | 1.03 | 0.15 | 6.85 |
| Repeat3_Sp1 | msp_C_std | 48/48 | 211-211 | 0 | -2.02 | 0.10 | -20.55 |
| Repeat3_Sp1 | nuc_85_250_C_std | 48/48 | 211-211 | 0 | 1.92 | 0.12 | 16.33 |

Usage-doc rule 4 predictions: P1_core:A_norm > 0 with no M change; Sp1 repeats A_norm < 0, msp_C_std < 0, nuc_85_250_C_std > 0; SRE-1 msp_C_std ~ 0 (|z| <= 2).


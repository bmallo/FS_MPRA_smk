# Klone agent instructions: run the site-methylation stage on replicate 3

Hand this file to an agent running on Klone. Fill in the four settings first.

```
REPO       = /mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/FS_MPRA_smk
SLURM      = --account=<ACCOUNT> --partition=<PARTITION>   (e.g. --account=stergachislab-ckpt --partition=ckpt-g2; ckpt is preemptible)
THRESHOLDS = <blank = config defaults: U distal m6A <= 0.01, M distal m6A >= 0.60>
PUSH       = <yes/no: may the agent commit the small summary bundle to docs/evidence/ and push branch marginal-null-model?>
```

---

You are running a prespecified analysis on Klone (UW HPC) in the FS_MPRA_smk repository at `REPO`. Your job is to execute and report, not to redesign. Read these first: `docs/site_methylation_usage.md` (what the stage does, interpretation rules, outputs) and the "Recommended plan" section of `docs/INDEPENDENT_REVIEW_20261020.md`. `CLAUDE.md` still describes a legacy pipeline; ignore its pipeline description for this task.

## Hard rules

- **Compute only inside Slurm jobs**, with one CPU and BLAS/OpenMP threads set to 1. The login node is for git, `ls`, reading small files and `sbatch`/`squeue`/`sacct` only.
- **Do not edit code, tests, audits, or source data to make something pass.** If a test, source guard (`Input changed during processing`, `differs from frozen ... source`) or run fails: stop, and report the exact command, the full traceback and the job log path. The only allowed fixes are environmental (a wrong path in a job script you wrote); report any you make.
- **Never overwrite or delete existing `results/` directories.** Outputs are atomic and refuse to overwrite. If a target directory already exists, stop and report it. You may delete only staging directories (`.extract-*`, `.effects-*`) left inside output directories that your own jobs created.
- **Do not change thresholds, windows, cohorts or endpoints after seeing any result.** If `THRESHOLDS` is filled in, set them in `config/ldlr_site_methylation.json` (`unmethylated_max_distal_m6a`, `methylated_min_distal_m6a`) before step 3 and report the change. Otherwise leave the config untouched.
- **Report numbers, not conclusions.** State outcomes against the prespecified rules mechanically, and keep any commentary clearly separate.
- Do not commit or push unless `PUSH = yes`, and then only step 6's bundle, never read IDs or `results/`.

## Step 0 — preflight (login node)

```bash
cd "$REPO"
git branch --show-current   # must be marginal-null-model; otherwise stop and report
git status --short          # report any local changes; do not discard them
git fetch origin marginal-null-model
git merge --ff-only origin/marginal-null-model   # if this is not a fast-forward, stop and report
git log --oneline -3        # must include "Add library-wide per-site m6A stage with same-reaction proxy controls"
PY=/mmfs1/home/bmallo/miniconda3/bin/python
$PY -c "import json; c=json.load(open('config/ldlr_site_methylation.json')); [print(k, c[k]) for k in ('features','feature_audit','catalog','footprint_contract','output_dir')]"
ls -l results/rep3_heterogeneity_20261017/features/molecule_heterogeneity.tsv.gz results/rep3_heterogeneity_20261017/features/audit.json \
      results/rep3_ledger_20261016/catalog/catalog.json config/rep3_existing_footprints.json
$PY -c "import json; b=json.load(open('config/rep3_existing_footprints.json'))['bam']; print(b)" | xargs -I{} ls -l {} {}.bai
ls -d results/rep3_site_methylation_20261021 2>/dev/null && echo "OUTPUT EXISTS: stop and report"
mkdir -p results/site_methylation_jobs
```

Report every missing file. Do not substitute other paths.

## Step 1 — tests and smoke run (Slurm job, ~30 min)

Write `results/site_methylation_jobs/smoke.sbatch` (fill in `SLURM` and `<REPO>`) and submit it with `sbatch` **from `$REPO`**, because log paths are relative. If `/usr/bin/time` does not exist on the compute node, drop it from the commands and report memory from `sacct ... --format=MaxRSS` instead.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=site_smoke
#SBATCH <SLURM account and partition options>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=results/site_methylation_jobs/smoke_%j.log
set -euo pipefail
cd "<REPO>"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
export MPLCONFIGDIR=/tmp/fsmpra-matplotlib-$SLURM_JOB_ID
PY=/mmfs1/home/bmallo/miniconda3/bin/python
echo "host=$(hostname) job=$SLURM_JOB_ID start=$(date)"
$PY -c "import sys,numpy,scipy,pysam,matplotlib; print(sys.version); print('numpy',numpy.__version__,'scipy',scipy.__version__,'pysam',pysam.__version__,'mpl',matplotlib.__version__)"
$PY -m unittest discover -s tests -p 'test_*.py'
OUT=results/site_smoke_$SLURM_JOB_ID
/usr/bin/time -v $PY -m fsmpra.site_extract --output $OUT/extract --limit-molecules 20000
/usr/bin/time -v $PY -m fsmpra.site_effects --extract $OUT/extract --output $OUT/effects --aa-draws 20
echo "SMOKE DONE $(date)"
```

Wait with `squeue -j <id>` / `sacct -j <id> --format=JobID,State,Elapsed,MaxRSS`, polling every few minutes, not continuously. Success means:
- all tests pass (expect about 132 tests, one optional skip);
- the extraction prints `Complete: 20,000 molecules`;
- the effects run prints `Complete`;
- the log ends with `SMOKE DONE`.

Report the test summary line, both commands' "Elapsed (wall clock)" and "Maximum resident set size", and the smoke `summary.json` `molecules` block. Smoke statistics are not results; do not interpret them. If anything fails, stop here and report.

## Step 2 — confirm thresholds (login node)

If `THRESHOLDS` is filled in, edit the config now (two numbers only) and show `git diff config/ldlr_site_methylation.json`. Otherwise confirm that the defaults are unchanged.

## Step 3 — full run (Slurm job, ~2–3 h)

Write `results/site_methylation_jobs/full.sbatch`. It is idempotent, so a requeue after preemption skips completed stages:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=site_full
#SBATCH <SLURM account and partition options>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=05:00:00
#SBATCH --requeue
#SBATCH --output=results/site_methylation_jobs/full_%j.log
set -euo pipefail
cd "<REPO>"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
export MPLCONFIGDIR=/tmp/fsmpra-matplotlib-$SLURM_JOB_ID
PY=/mmfs1/home/bmallo/miniconda3/bin/python
OUT=results/rep3_site_methylation_20261021
echo "host=$(hostname) job=$SLURM_JOB_ID restart=${SLURM_RESTART_COUNT:-0} start=$(date)"
[ -d $OUT/extract ] || /usr/bin/time -v $PY -m fsmpra.site_extract --output $OUT/extract
[ -d $OUT/effects ] || /usr/bin/time -v $PY -m fsmpra.site_effects --extract $OUT/extract --output $OUT/effects
for pair in "0.02 0.7" "0.005 0.5"; do
  set -- $pair
  [ -d $OUT/sensitivity/u$1_m$2 ] || /usr/bin/time -v $PY -m fsmpra.site_effects --extract $OUT/extract \
      --output $OUT/sensitivity/u$1_m$2 --u-max $1 --m-min $2 --aa-draws 100
done
[ -d $OUT/effects_supported_barcodes ] || /usr/bin/time -v $PY -m fsmpra.site_effects --extract $OUT/extract \
    --output $OUT/effects_supported_barcodes --cohort eligible_with_barcode --aa-draws 100
echo "FULL DONE $(date)"
```

If `THRESHOLDS` changed the primary cutoffs, keep the two sensitivity pairs as written. The extraction prints progress every 500,000 BAM records; expect about 8.7M primary records. When the job finishes, run a Snakemake dry run as a consistency check only. It will report missing execution metadata because the stages ran directly; that is expected:

```bash
/mmfs1/home/bmallo/miniconda3/envs/snakemake/bin/snakemake -s workflow/site_methylation.smk -n --cores 1
```

Report:
- each stage's wall time and max RSS;
- `extract/audit.json` fields `complete`, `molecules`, `barcodes`, `bam_primary_rows_scanned`, and the nonzero bins of `retained_positive_ml_histogram`;
- any requeues.

## Step 4 — report the prespecified readouts (any node; small files)

For `effects/` (primary) and, where noted, the three sensitivity directories, extract and tabulate the following. Use short Python reading the TSV/JSON outputs; do not plot new figures. Values in `effects.tsv.gz` and `kernel.tsv` are fractions: multiply `delta`, `se`, `mean_delta`, `variant_value` and `wt_value` by 100 to report percentage points. Blank cells are unavailable measurements, not zeros.

1. **Molecules and classes** (`summary.json`): `molecules`, `thresholds` (including the labeling strata edges), and `class_dispersion`. Also the WT fraction of U and M molecules.
2. **Calibration gate** (`calibration_summary.tsv`): the number of endpoints passing (`sd_z` in 0.85–1.20 and `rate_ci_high` ≤ 0.10) out of the total, excluding `:U` endpoints. List every failing non-U endpoint with `sd_z`, `rate_abs_z_gt_1_96`, `rate_lower_tail`, `rate_upper_tail`. Give the four primary endpoints (`local:A_norm`, `flank:A_norm`, `msp_C_std`, `nuc_85_250_C_std`) explicitly.
3. **Observability kernel** (`kernel.tsv`, `at_change == 'all'`): for measures `M` and `C_std`, the offsets from −15 to +15 with `mean_delta` (×100, pp) and `se` (pp). Flag offsets with |mean| > 2·SE and |offset| > 2.
4. **Position 125** (`effects.tsv.gz`, `source_variant` in `125A>C`, `125A>G`, `125A>T`): for endpoints `WT_P3:M`, `WT_P3:C_std`, `WT_P3:A_norm`, `local:M`, `local:C_std`, `local:A_norm`, give `delta` (pp), `se`, `z`, `molecules_M`. Do this for the primary run and both threshold sensitivities. Then state rule 3 of `docs/site_methylation_usage.md` mechanically:
   - Does `WT_P3:M` move in the same direction as `WT_P3:C_std` with at least half its magnitude, for at least two of the three substitutions?
   - Is `WT_P3:A_norm` z > 3 for all three, under all three threshold pairs?
5. **Prespecified regions**: average the per-SNV `delta` across SNVs whose position lies inside the region (from `regions0` in the config: `P1_core`, `Repeat1_Sp1`, `Repeat2_SRE1`, `Repeat3_Sp1`). Do this for endpoints `<region>:A_norm`, `<region>:M`, `msp_C_std`, `nuc_85_250_C_std`. SE is sqrt(sum se²)/n. Show the predictions from rule 4 next to the numbers, as agree/disagree.
6. **Primary-family discoveries**: count rows with `q_bh_primary_family` < 0.05 per primary endpoint. List the top 20 by q with `source_variant`, `endpoint`, `delta`, `z`, `q_bh_primary_family`, and both activity columns.
7. **Library association** (`activity_association.tsv`): for the primary endpoints plus `msp_ALL`, `nuc_85_250_ALL`, `local:C_std`, `local:M` and `element:A_norm`, give `variant_pearson`, `position_pearson`, `null_abs_max` and `p_circular_shift` for both activity sources. Do this for the primary run and `effects_supported_barcodes`. State rule 5: is r(`msp_C_std`) > r(`msp_ALL`) > 0, and is r(`nuc_85_250_C_std`) < 0?
8. **Observability regression** (`observability_regression.tsv`): all rows.

## Step 5 — look at the figures

Open and briefly describe what you see in each, without inferring biology:
- `effects/distal_m6a_histogram.png`: where the thresholds fall relative to the distribution;
- `kernel.png`;
- `calibration_z.png`;
- `observability_scatter.png`.

## Step 6 — summary bundle

Create `docs/evidence/site_methylation_20261021/` (stop if it exists). Copy into it, from each of `effects/`, `sensitivity/u0.02_m0.7/`, `sensitivity/u0.005_m0.5/` and `effects_supported_barcodes/` (in subdirectories of the same names):
- `summary.json`, `calibration_summary.tsv`, `kernel.tsv`, `activity_association.tsv`, `observability_regression.tsv`, `class_counts.tsv`, `distal_m6a_histogram.tsv`, `effects.tsv.gz`, and the PNG figures.
- `extract/audit.json` (as `extract_audit.json`).

Do **not** copy `molecules.tsv.gz`, any `.npy`/`.npz`, or `calibration_draws.tsv.gz`. Check the total size (expect < 40 MB) and report it. Write a `README.md` in the bundle listing:
- the source directories, the job IDs, and the git commit used;
- the thresholds;
- the tables from step 4, verbatim.

If `PUSH = yes`, run `git add docs/evidence/site_methylation_20261021`. Commit only that path, with a message naming the run, and run `git push origin marginal-null-model`. If `PUSH = no`, tell the user the bundle path.

## Final message to the user

Keep it short:
1. Did everything run? Give job IDs, times and memory.
2. The step-4 tables, or the bundle path.
3. Each prespecified rule outcome as one line.
4. Anything that failed or looked wrong, with the log path.

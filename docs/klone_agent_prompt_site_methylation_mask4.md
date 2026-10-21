# Klone agent instructions: mask-4 follow-up of the site-methylation run

Hand this to an agent on Klone in the same checkout used for job 42090780 (`/mmfs1/gscratch/stergachislab/bmallo/large_home/code/mine/FS_MPRA_smk`). Use the same settings as last time unless told otherwise: `--account=stergachislab --partition=cpu-g2`, PUSH = yes.

---

You are rerunning only the analysis stage of a prespecified analysis. It reuses the existing extraction `results/rep3_site_methylation_20261021/extract` and needs no BAM pass. Read `docs/rep3_site_methylation_20261021.md`, especially "Follow-up run, prespecified before it is executed". The **hard rules** in `docs/klone_agent_prompt_site_methylation.md` all apply:
- compute only in Slurm jobs;
- never edit code, tests, audits or source data to make something pass;
- never overwrite or delete `results/` directories;
- no changes after seeing results;
- report numbers, not conclusions.

## Step 0 — preflight (login node)

```bash
cd "$REPO"; git branch --show-current; git status --short
git fetch origin marginal-null-model && git merge --ff-only origin/marginal-null-model
git log --oneline -3      # must include "Interpret site-methylation run; add mask-4 follow-up"
git diff fc4731a -- config/ldlr_site_methylation.json fsmpra/site_extract.py fsmpra/io.py fsmpra/processing.py   # must be empty
PY=/mmfs1/home/bmallo/miniconda3/bin/python
$PY -c "import json; a=json.load(open('results/rep3_site_methylation_20261021/extract/audit.json')); print(a['complete'], a['molecules'])"   # True 2659690
ls -d results/rep3_site_methylation_20261021/mask4 2>/dev/null && echo "OUTPUT EXISTS: stop and report"
```

The `git diff` must be empty because the extraction audit pins those files. A non-empty diff means stop and report.

## Step 1 — one Slurm job (~1 h)

Write `results/site_methylation_jobs/mask4.sbatch`, submit it from `$REPO`, and wait (poll every few minutes):

```bash
#!/usr/bin/env bash
#SBATCH --job-name=site_mask4
#SBATCH --account=stergachislab
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=results/site_methylation_jobs/mask4_%j.log
set -euo pipefail
cd "<REPO>"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
export MPLCONFIGDIR=/tmp/fsmpra-matplotlib-$SLURM_JOB_ID
PY=/mmfs1/home/bmallo/miniconda3/bin/python
CFG=config/ldlr_site_methylation_mask4.json
EX=results/rep3_site_methylation_20261021/extract
OUT=results/rep3_site_methylation_20261021/mask4
echo "host=$(hostname) job=$SLURM_JOB_ID start=$(date)"
$PY -m unittest discover -s tests -p 'test_*.py'
[ -d $OUT/effects ] || /usr/bin/time -v $PY -m fsmpra.site_effects --config $CFG --extract $EX --output $OUT/effects
for pair in "0.02 0.7" "0.005 0.5"; do
  set -- $pair
  [ -d $OUT/sensitivity/u$1_m$2 ] || /usr/bin/time -v $PY -m fsmpra.site_effects --config $CFG --extract $EX \
      --output $OUT/sensitivity/u$1_m$2 --u-max $1 --m-min $2
done
[ -d $OUT/effects_supported_barcodes ] || /usr/bin/time -v $PY -m fsmpra.site_effects --config $CFG --extract $EX \
    --output $OUT/effects_supported_barcodes --cohort eligible_with_barcode
echo "MASK4 DONE $(date)"
```

Expect 133 tests passing, with one optional skip. Report each run's wall time and max RSS.

## Step 2 — prespecified readouts

Use the four `mask4` directories. Values are fractions, so ×100 gives percentage points, and blank cells are unavailable measurements, not zeros. `effects.tsv.gz` now has a `context` column (`homopolymer` / `other`), and `kernel.tsv` has `context` too.

1. **Calibration (rule 1):** `calibration_summary.tsv` `gate_pass` counts for non-`:U` endpoints, the four primary endpoints in full, and every failing non-U endpoint.
2. **Homopolymer reach (rule 2):** from the primary `kernel.tsv`, for `measure` in {M, C_std}, `at_change` in {gain, loss, retained_count} and `context` in {homopolymer, other}: `mean_delta` (pp), `se` and z for offsets −14…+14. State mechanically whether any homopolymer-context M row has |z| > 2 at |offset| ≥ 5; if so, local endpoints for homopolymer-context SNVs are "unresolved".
3. **Residual artifact (rule 3):** for `local:A_norm`, give the mean Δ over SNVs with `at_change == gain` and with `loss`. Their difference has SE sqrt(var_gain/n_gain + var_loss/n_loss), using per-SNV deltas, with only SNVs with ≥ 300 molecules. Give both the difference and its z, plus the gate status. `local:A_norm` is usable only if the gate passes and |z| < 3. Repeat this for `context == other` SNVs only.
4. **Position 125 (rule 4):** the same table as the first run (`WT_P3:M`, `WT_P3:C_std`, `WT_P3:A_norm`, `local:*`, under all three threshold pairs). Apply the position-125 rule from `docs/site_methylation_usage.md` (its rule 3) unchanged, both the observability clause and the chromatin clause. Also give `sites` for each endpoint, since the wider mask leaves fewer.
5. **Primary family:** counts of `q_bh_primary_family` < 0.05 per primary endpoint. For each endpoint's hits, give:
   - the number up and down;
   - the fraction whose sign agrees with `activity_LDLR` (opposite sign counts as agreement for `nuc_85_250_C_std`);
   - counts by element: P1 core 54–62, Repeat 1 89–104, SRE-1 124–139, Repeat 3 140–155, other (oligo positions);
   - counts by `context`.

   Then the top 20 rows.
6. **Library association:** `activity_association.tsv` rows for the four primary endpoints plus `msp_ALL`, `nuc_85_250_ALL`, `element:C_std` and `local:C_std`, for both activity sources, primary and supported-barcode runs.
7. **Regions:** as in the first run: mean Δ and z over in-region SNVs for `<region>:A_norm`, `<region>:M`, `msp_C_std` and `nuc_85_250_C_std`, with `sites`. Note when the ±4 mask leaves no sites.

## Step 3 — bundle and report

Create `docs/evidence/site_methylation_mask4_20261021/`, stopping if it exists, with the same file selection as before for each of the four directories. Use no read IDs, `.npy`/`.npz` or `calibration_draws`. Add a `README.md` with provenance (commit, job ID, config) and the step-2 tables verbatim. Commit only that path and push to `marginal-null-model`.

End with a short message to the user:
- did it run (job, times, memory);
- each rule outcome in one line;
- anything odd, with log paths.

Mention that `results/site_smoke_42090480/` can be deleted by the user; do not delete it yourself.

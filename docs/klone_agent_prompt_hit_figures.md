# Klone agent instructions: per-base profiles for the top-hit figures

Hand this to an agent on Klone in the same checkout used for jobs 42090780 and 42093363 (`/mmfs1/gscratch/stergachislab/bmallo/large_home/code/mine/FS_MPRA_smk`). Use the same settings as before: `--account=stergachislab --partition=cpu-g2`, PUSH = yes.

---

You are running one descriptive step that reads the existing site-methylation extraction `results/rep3_site_methylation_20261021/extract`. It needs no BAM pass and computes no new statistics. It writes per-base WT and variant occupancy and m6A profiles for all 634 alleles, plus small molecule rasters, so that the WT-vs-variant figures in `tools/variant_hit_figures.py` can be drawn. The hard rules from `docs/klone_agent_prompt_site_methylation.md` all apply:
- compute only in Slurm jobs;
- never edit code, tests, audits or source data to make something pass;
- never overwrite or delete `results/` directories;
- report numbers, not conclusions.

## Step 0 — preflight (login node)

```bash
cd "$REPO"; git branch --show-current; git status --short
git fetch origin marginal-null-model && git merge --ff-only origin/marginal-null-model
git log --oneline -3      # must include "Add top-hit evidence figures and footprint-motif analysis"
git diff fc4731a -- config/ldlr_site_methylation.json fsmpra/site_extract.py fsmpra/io.py fsmpra/processing.py   # must be empty
PY=/mmfs1/home/bmallo/miniconda3/bin/python
$PY -c "import json; a=json.load(open('results/rep3_site_methylation_20261021/extract/audit.json')); print(a['complete'], a['molecules'])"   # True 2659690
for d in results/rep3_site_methylation_20261021/hit_profiles results/rep3_site_methylation_20261021/hit_figures docs/evidence/top_hits_20261021/klone_profiles; do
  ls -d $d 2>/dev/null && echo "OUTPUT EXISTS: stop and report"
done
```

Stop and report if any check fails.

## Step 1 — one Slurm job (expected well under 1 h)

Write `results/site_methylation_jobs/hit_profiles.sbatch`, submit it from `$REPO`, and wait (poll every few minutes):

```bash
#!/usr/bin/env bash
#SBATCH --job-name=hit_profiles
#SBATCH --account=stergachislab
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=results/site_methylation_jobs/hit_profiles_%j.log
set -euo pipefail
cd "<REPO>"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
export MPLCONFIGDIR=/tmp/fsmpra-matplotlib-$SLURM_JOB_ID
PY=/mmfs1/home/bmallo/miniconda3/bin/python
EX=results/rep3_site_methylation_20261021/extract
OUT=results/rep3_site_methylation_20261021
echo "host=$(hostname) job=$SLURM_JOB_ID start=$(date)"
$PY -m unittest discover -s tests -p 'test_*.py'
/usr/bin/time -v $PY -m tools.variant_hit_figures profiles --extract $EX --output $OUT/hit_profiles
/usr/bin/time -v $PY -m tools.variant_hit_figures figures --profiles $OUT/hit_profiles --output $OUT/hit_figures
ls -la $OUT/hit_profiles $OUT/hit_figures
echo "HIT PROFILES DONE $(date)"
```

Expect 135 tests passing, with one optional skip. Report wall time and max RSS for the two steps.

## Step 2 — consistency readout (login node is fine; it reads two small files)

```bash
$PY - <<'EOF'
import csv, gzip, json
import numpy as np
from tools.variant_hit_figures import TOP_HITS
p = np.load('results/rep3_site_methylation_20261021/hit_profiles/profiles.npz')
names = list(p['source_variant']); wt = names.index('WT')
eff = {}
with gzip.open('docs/evidence/site_methylation_mask4_20261021/effects/effects.tsv.gz', 'rt') as h:
    for r in csv.DictReader(h, delimiter='\t'):
        if r['endpoint'] in ('msp_C_std', 'nuc_85_250_C_std'):
            eff[(r['source_variant'], r['endpoint'])] = r['delta']
print('allele\tmolecules\tU\tM\tC\tprofile_dMSP_pp\teffects_dMSP_pp\tprofile_dNUC_pp\teffects_dNUC_pp\tfinite_bases_msp_C_std')
for h in ['WT'] + TOP_HITS:
    i = names.index(h); m = p['molecules'][i]
    d = lambda k: 100 * np.nanmean(p[k][i] - p[k][wt])
    e = lambda k: 100 * float(eff[(h, k)]) if eff.get((h, k)) else float('nan')
    print(f"{h}\t{m.sum()}\t{m[0]}\t{m[1]}\t{m[2:].sum()}\t{d('msp_C_std'):+.2f}\t{e('msp_C_std'):+.2f}\t"
          f"{d('nuc_C_std'):+.2f}\t{e('nuc_85_250_C_std'):+.2f}\t{np.isfinite(p['msp_C_std'][i]).sum()}")
print('alleles with finite msp_C_std:', int(np.isfinite(p['msp_C_std']).all(1).sum()), 'of', len(names))
r = np.load('results/rep3_site_methylation_20261021/hit_profiles/rasters.npz')
print('raster alleles:', len(r['names']), 'molecules per allele:', sorted({len(r[f'state_{k}']) for k in range(len(r['names']))}))
print(json.dumps(json.load(open('results/rep3_site_methylation_20261021/hit_profiles/audit.json')), indent=1))
EOF
```

Report the table verbatim. The profile column averages per-base fractions over the 211 bases, while the effects column uses each molecule's element MSP fraction from the extraction. Both use the same strata, so they should be close, but they are not the same quantity. Report any difference larger than 1 pp; do not try to fix it.

## Step 3 — bundle and report

Create `docs/evidence/top_hits_20261021/klone_profiles/`, stopping if it exists, containing only:
- `profiles.npz`, `rasters.npz` and `audit.json` from `hit_profiles` (aggregate per-base fractions and 150 anonymous molecule rasters per allele; no read IDs);
- the job log, as `hit_profiles_<jobid>.log`;
- a `README.md` with provenance (commit, job ID, command lines) and the step-2 output verbatim.

If the bundle exceeds 30 MB, stop and report the sizes instead of committing. Do not commit `hit_figures/`; it stays in `results/` for the user to browse. Commit only the bundle path and push to `marginal-null-model`.

End with a short message to the user:
- did it run (job, times, memory);
- the step-2 table;
- anything odd, with log paths;
- the path of `hit_figures/` on Klone.

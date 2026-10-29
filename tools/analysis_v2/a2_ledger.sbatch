#!/usr/bin/env bash
#SBATCH --account=stergachislab
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --job-name=a2_ledger
# Step A2: the run-3 ledger workflow (workflow/ledger.smk) on a realigned HepG2 run.
# Usage: sbatch --output=<log> tools/analysis_v2/a2_ledger.sbatch <config.json> <aligned.bam> <run_dir>
set -euo pipefail
cd /mmfs1/gscratch/stergachislab/bmallo/large_home/code/mine/FS_MPRA_smk
export XDG_CACHE_HOME=${TMPDIR:-/tmp}/fsmpra-smk-cache OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
SMK=/mmfs1/home/bmallo/miniconda3/envs/snakemake/bin/snakemake
$SMK -s workflow/ledger.smk --cores ${SLURM_CPUS_PER_TASK:-1} --resources mem_mb=12000 \
  --config processing_config=$1 processing_python=/mmfs1/home/bmallo/miniconda3/bin/python \
  ledger_run_dir=$3 ledger_bam=$2

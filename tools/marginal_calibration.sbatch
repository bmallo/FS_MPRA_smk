#!/bin/bash
# Multi-seed calibration of --null-model marginal.
#
# SLURM array: each task runs one permutation (one Monte-Carlo replicate)
# of tools/marginal_permutation_test.py with a distinct seed. Aggregate
# the per-seed JSONs afterwards with tools/aggregate_marginal_calib.py.
#
#   sbatch tools/marginal_calibration.sbatch
#
#SBATCH --job-name=marginal_calib
#SBATCH --partition=ckpt-g2
#SBATCH --account=stergachislab
#SBATCH --qos=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --array=1-8
#SBATCH --output=/mmfs1/gscratch/stergachislab/bmallo/large_home/projects/Mouse_MPRA/logs/marginal_calib_%A_%a.out
#SBATCH --error=/mmfs1/gscratch/stergachislab/bmallo/large_home/projects/Mouse_MPRA/logs/marginal_calib_%A_%a.err

set -euo pipefail

PYTHON=/mmfs1/home/bmallo/miniconda3/envs/general/bin/python3
REPO=/gscratch/stergachislab/bmallo/large_home/git_repos/FS_MPRA_smk
PROJECT=/mmfs1/gscratch/stergachislab/bmallo/large_home/projects/Mouse_MPRA
SEED=$SLURM_ARRAY_TASK_ID

echo "[$(date)] marginal calibration seed=$SEED node=$(hostname)"

"$PYTHON" "$REPO/tools/marginal_permutation_test.py" \
    --bam "$PROJECT/variant_calling_merged/Gata4/PS01760_merged_Gata4.tagged.bam" \
    --target-region Gata4 \
    --promoter-region 1734-2033 \
    --output-dir "$PROJECT/marginal_validation/calib_Gata4/seed_${SEED}" \
    --seed "$SEED" \
    --n-null-iterations 2000 \
    --min-reads 50 \
    --threads 16

echo "[$(date)] seed=$SEED done"

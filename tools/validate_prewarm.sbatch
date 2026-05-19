#!/usr/bin/env bash
#SBATCH --job-name=fsmpra_warmchk
#SBATCH --account=stergachislab-ckpt
#SBATCH --partition=ckpt-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --requeue
#SBATCH --output=results/phase0_hmm_full/analysis_motif/warmchk_%j.log

# Validate the worker-pool pre-warm on a real 64-core ckpt-g2 node
# (where the ~233s variant-1 cold penalty appeared). profile_fullpath
# now exercises run_variant_testing_parallel, which pre-warms the pool;
# expect a "Worker pool pre-warmed (...)" line then variant 1 ~ later
# variants (~3.5s) instead of ~233s.
set -euo pipefail
cd /mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/FS_MPRA_smk
export PATH="/mmfs1/home/bmallo/.pixi/bin:$PATH"
echo "host=$(hostname) cpus=$SLURM_CPUS_PER_TASK start=$(date)"
pixi run python tools/profile_fullpath_prod.py
echo "validate_prewarm DONE $(date)"

#!/usr/bin/env bash
#SBATCH --job-name=fsmpra_coocc_calib
#SBATCH --account=stergachislab-ckpt
#SBATCH --partition=ckpt-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --output=results/phase0_hmm_full/analysis_motif/coocc_calib_%j.log

# P3.6 BLOCKING gate: WT-vs-WT co-occupancy calibration on the real
# LDLR production data. Must PASS (no false dependencies under H0)
# before trusting any biological co-occupancy call. ckpt-g2 + shared-fs
# log + sbatch (the Phase-2 SLURM gotchas). Single-process numpy.
set -euo pipefail
cd /mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/FS_MPRA_smk
export PATH="/mmfs1/home/bmallo/.pixi/bin:$PATH"
echo "host=$(hostname) cpus=$SLURM_CPUS_PER_TASK start=$(date)"
pixi run python tools/cooccupancy_wt_vs_wt_prod.py \
  --bam results/phase0_hmm_full/variants_snv/LDLR_hmm_full.tagged.bam \
  --h5 results/phase0_hmm_full/analysis_motif/LDLR_hmm_full.h5 \
  --promoter-region 3184-3501
echo "coocc calibration DONE $(date)"

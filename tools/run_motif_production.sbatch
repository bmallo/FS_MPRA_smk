#!/usr/bin/env bash
#SBATCH --job-name=fsmpra_motif
#SBATCH --account=stergachislab-ckpt
#SBATCH --partition=ckpt-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --output=results/phase0_hmm_full/analysis_motif/slurm_%j.log

# Phase 2 production motif run on the full LDLR_hmm SNV-only library.
# Validated config: per-variant nulls, 10000 iterations, min-reads 50,
# stratify OFF (calibrated). Single-node job — the null build
# parallelizes across cpus-per-task within one node (ProcessPool +
# shared memory), so cores-per-node is the lever. ckpt-g2 = the same
# Genoa Gen2 hardware as cpu-g2 (128-192 cores, 1.5 TB) but on the
# scavenger account (no saturated lab QOS cap, idle nodes available).
# Preemptible -> --requeue: a preempted job reruns cleanly (Stage 3
# overwrites its output dir; ~30-45 min wall so preemption is unlikely).
set -euo pipefail
cd /mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/FS_MPRA_smk

export PATH="/mmfs1/home/bmallo/.pixi/bin:$PATH"

BAM=results/phase0_hmm_full/variants_snv/LDLR_hmm_full.tagged.bam
REF=/gscratch/stergachislab/bmallo/large_home/ref_seqs/LDLR_MPRA/LDLR.fa
PROMOTER=3184-3501
OUTDIR=results/phase0_hmm_full/analysis_motif
mkdir -p "$OUTDIR"

echo "host=$(hostname) cpus=$SLURM_CPUS_PER_TASK mem=${SLURM_MEM_PER_NODE:-?}MB start=$(date)"
pixi run python workflow/scripts/03_analyze_library.py \
  --bam "$BAM" \
  --output-dir "$OUTDIR" \
  --sample-name LDLR_hmm_full \
  --promoter-region "$PROMOTER" \
  --reference "$REF" \
  --n-null-iterations 10000 \
  --min-reads 50 \
  --threads "$SLURM_CPUS_PER_TASK" \
  -v
echo "motif production DONE $(date)"

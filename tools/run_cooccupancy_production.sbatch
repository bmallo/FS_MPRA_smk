#!/usr/bin/env bash
#SBATCH --job-name=fsmpra_coocc_prod
#SBATCH --account=stergachislab-ckpt
#SBATCH --partition=ckpt-g2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --output=results/phase0_hmm_full/analysis_cooc/run_%j.log

# Production co-occupancy biological run: full Stage 3 (Phase 1+2) plus
# the Phase-3 co-occupancy module on the real LDLR library. Validated
# config (10000 null iters, min-reads 50, stratify OFF, --reference);
# co-occupancy WT-vs-WT was certified by P3.6 (exact FDR). Separate
# output dir so the motif-only production results are preserved.
# ckpt-g2 + shared-fs log + sbatch (learned SLURM conventions).
set -euo pipefail
cd /mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/FS_MPRA_smk
export PATH="/mmfs1/home/bmallo/.pixi/bin:$PATH"
echo "host=$(hostname) cpus=$SLURM_CPUS_PER_TASK start=$(date)"
pixi run python workflow/scripts/03_analyze_library.py \
  --bam results/phase0_hmm_full/variants_snv/LDLR_hmm_full.tagged.bam \
  --output-dir results/phase0_hmm_full/analysis_cooc \
  --sample-name LDLR_hmm_full \
  --promoter-region 3184-3501 \
  --reference /gscratch/stergachislab/bmallo/large_home/ref_seqs/LDLR_MPRA/LDLR.fa \
  --n-null-iterations 10000 \
  --min-reads 50 \
  --enable-co-occupancy \
  --cooccupancy-n-iter 10000 \
  --threads "$SLURM_CPUS_PER_TASK" \
  -v
echo "co-occupancy production DONE $(date)"

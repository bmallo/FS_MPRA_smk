#!/usr/bin/env bash
#SBATCH --job-name=fsmpra_sweep
#SBATCH --account=stergachislab-ckpt
#SBATCH --partition=ckpt
#SBATCH --array=1-42
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --requeue
#SBATCH --output=results/phase0_hmm_full/sweep/cell_%a.log

# Best-powered WT-vs-WT calibration sweep (disjoint pseudo-variants).
# One array task per cell in tools/calibration_sweep/cells.tsv
# (idx arm N stratify seed). ckpt + --requeue: a preempted cell just
# reruns; cells are independent so partial completion is safe.
#
# Edit BAM / region for a different dataset; run artifacts go to a
# gitignored results/ dir. Aggregate with aggregate_sweep.py.
set -euo pipefail
cd /mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/FS_MPRA_smk

BAM=results/phase0_hmm_full/variants_snv/LDLR_hmm_full.tagged.bam
TARGET=LDLR:3184-3501
PROMOTER=3184-3501
CELLS=tools/calibration_sweep/cells.tsv
OUTDIR=results/phase0_hmm_full/sweep
mkdir -p "$OUTDIR"

ROW=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'NR==i+1{print}' "$CELLS")
IDX=$(echo "$ROW" | cut -f1); ARM=$(echo "$ROW" | cut -f2)
N=$(echo "$ROW" | cut -f3);  STRAT=$(echo "$ROW" | cut -f4)
SEED=$(echo "$ROW" | cut -f5)
STRAT_FLAG=""; [ "$STRAT" = "on" ] && STRAT_FLAG="--stratify"

echo "cell $IDX arm=$ARM N=$N stratify=$STRAT seed=$SEED cpus=$SLURM_CPUS_PER_TASK"
pixi run python tools/wt_vs_wt_calibration.py \
  --bam "$BAM" --target-region "$TARGET" --promoter-region "$PROMOTER" \
  --pseudo-n "$N" --n-pseudo 1000000 \
  --n-null-iterations 10000 --threads "$SLURM_CPUS_PER_TASK" \
  --seed "$SEED" --disjoint-pseudo $STRAT_FLAG \
  --out "$OUTDIR/cell_${IDX}_${ARM}_seed${SEED}.json"
echo "cell $IDX DONE"

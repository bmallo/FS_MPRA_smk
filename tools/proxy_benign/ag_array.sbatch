#!/bin/bash
# AlphaGenome scoring, one array task per element of tools/proxy_benign/elements.tsv (rules.json "alphagenome_detail").
# Each task scores into results/proxy_benign_20261106/ag/<element>/ (own windows file, links to the shared gnomAD and
# reference pulls), so tasks never write the same file. Submit from the repo root:
#   sbatch --array=0-15%8 tools/proxy_benign/ag_array.sbatch
#SBATCH --job-name=ag-proxy
#SBATCH --account=stergachislab
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=results/proxy_benign_20261106/logs/ag_%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
O=results/proxy_benign_20261106
E=tools/proxy_benign/elements.tsv
W=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==i{print $1}' $E)
D=$O/ag/$W
mkdir -p $D/gnomad $D/reference
awk -F'\t' -v w=$W 'NR==1 || $1==w' $E > $D/windows.tsv
ln -sf $(realpath $O/gnomad/$W.json) $D/gnomad/$W.json
ln -sf $(realpath $O/reference/$W.txt) $D/reference/$W.txt
echo "$W on $(hostname -s), start $(date)"
/gscratch/stergachislab/bmallo/large_home/envs/alphagenome/bin/python tools/promoter_variant_census/ag_score.py \
    --windows-tsv $D/windows.tsv --output $D
echo "end $(date)"

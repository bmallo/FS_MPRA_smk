#!/bin/bash
# Figure 1 (WASF1): FIRE elements as in the LRSS slides (05_fire_fiberhmm.sh: fibertools-rs 0.13 add-nucleosomes ->
# fire, no --fire-filter) on the reads fig1_fiberhmm.sh staged, then ft extract --all (fig1_rules.json: wasf1.draft1).
# Usage: bash tools/r01_figure/fig1_fire.sh <work dir> <threads> <dataset ...>
set -euo pipefail
W=$1; T=$2; shift 2
FT=/mmfs1/gscratch/stergachislab/bin/fibertools-rs/ft-wrapper.sh
mkdir -p $W/fire
for b in "$@"; do
  $FT add-nucleosomes -t $T $W/bams/$b.sorted.bam - | $FT fire -t $T - - | samtools sort -o $W/fire/$b.fire.bam -
  samtools index $W/fire/$b.fire.bam
  $FT extract --all - $W/fire/$b.fire.bam | gzip > $W/fire/$b.extract.tsv.gz
  echo "$b: $(samtools view -c $W/fire/$b.fire.bam) reads"
done

#!/bin/bash
# Figure 1 (EF1a): stage pEF-GFP replicate 3 (subsample) and the endogenous EEF1A1 reads, call footprints with fiberhmm
# and FIRE (fig1_rules.json: ef1a). Usage: bash tools/r01_figure/fig1_ef1a_stage.sh <out_dir> <cores>
set -euo pipefail
OUT=$1; CORES=${2:-1}
LRSS=/gscratch/stergachislab/bmallo/large_home/resources/presentations/lrss2026_slides56
D=/gscratch/stergachislab/bmallo/large_home/data/kopah_restore
mkdir -p $OUT/bams $OUT/fiberhmm
export FIBERHMM_NO_UPDATE_CHECK=1
samtools view -b -F 0x900 -s 42.03 $LRSS/pefgfp/pefgfp_r3.span.bam > $OUT/bams/pEFGFP.bam
samtools view -b -F 0x900 $D/Hek293_Whole_Genome/target_bams/EF1a_all_reads.bam > $OUT/bams/endogenous.bam
for b in pEFGFP endogenous; do
  samtools sort -o $OUT/bams/$b.sorted.bam $OUT/bams/$b.bam && samtools index $OUT/bams/$b.sorted.bam
  echo "$b reads: $(samtools view -c $OUT/bams/$b.sorted.bam)"
  fiberhmm-call -i $OUT/bams/$b.sorted.bam -o $OUT/fiberhmm/$b.fp.bam --enzyme hia5 --seq pacbio --min-mapq 0 \
      -c $CORES --no-qc 2> $OUT/fiberhmm/$b.log
  samtools index $OUT/fiberhmm/$b.fp.bam
done
bash $(dirname $0)/fig1_fire.sh $OUT $CORES pEFGFP endogenous

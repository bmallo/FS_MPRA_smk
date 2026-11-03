#!/bin/bash
# Figure 1: subsample DPAGT1 datasets and call footprints with fiberhmm (settings: tools/r01_figure/fig1_rules.json).
# Usage: bash tools/r01_figure/fig1_fiberhmm.sh <out_dir> <cores>
set -euo pipefail
OUT=$1; CORES=${2:-1}
D=/gscratch/stergachislab/bmallo/large_home/data/kopah_restore
mkdir -p $OUT/bams $OUT/fiberhmm
export FIBERHMM_NO_UPDATE_CHECK=1

# 300-bp constructs: ~6,000 reads per construct per replicate (samtools view -s, seed 42)
for rep in 1 2; do
  B=$D/IGVF_construct_mix/source_bams/IGVF_construct_mix_rep_${rep}.ft.aligned.bam
  for c in DPAGT11 DPAGT13; do
    n=$(samtools view -c -F 0x900 $B $c)
    frac=$(python3 -c "print(min(0.999, 6000/$n))")
    samtools view -b -F 0x900 -s 42${frac#0} $B $c > $OUT/bams/${c}_rep${rep}.bam
  done
done
samtools merge -f $OUT/bams/DPAGT11.bam $OUT/bams/DPAGT11_rep1.bam $OUT/bams/DPAGT11_rep2.bam
samtools merge -f $OUT/bams/DPAGT13.bam $OUT/bams/DPAGT13_rep1.bam $OUT/bams/DPAGT13_rep2.bam
# full-length 4.1-kb construct (all reads on the DPAGT1 contig) and endogenous reads
samtools view -b -F 0x900 $D/genomic_promoters/merged_bams/DPAGT1/DPAGT1.ft.m6a.merged.bam DPAGT1 > $OUT/bams/DPAGT1_full.bam
samtools view -b -F 0x900 $D/Hek293_Whole_Genome/target_bams/DPATG1_all_reads.bam > $OUT/bams/endogenous.bam
for b in DPAGT11 DPAGT13 DPAGT1_full endogenous; do
  samtools sort -o $OUT/bams/$b.sorted.bam $OUT/bams/$b.bam && samtools index $OUT/bams/$b.sorted.bam
  echo "$b reads: $(samtools view -c $OUT/bams/$b.sorted.bam)"
  fiberhmm-call -i $OUT/bams/$b.sorted.bam -o $OUT/fiberhmm/$b.fp.bam --enzyme hia5 --seq pacbio --min-mapq 0 \
      -c $CORES --no-qc 2> $OUT/fiberhmm/$b.log
  samtools index $OUT/fiberhmm/$b.fp.bam
done

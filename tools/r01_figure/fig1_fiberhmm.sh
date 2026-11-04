#!/bin/bash
# Figure 1: subsample the gene's datasets and call footprints with fiberhmm (settings: tools/r01_figure/fig1_rules.json).
# Usage: bash tools/r01_figure/fig1_fiberhmm.sh <out_dir> <cores> [DPAGT1|WASF1]
set -euo pipefail
OUT=$1; CORES=${2:-1}; GENE=${3:-DPAGT1}
D=/gscratch/stergachislab/bmallo/large_home/data/kopah_restore
case $GENE in
  DPAGT1) C300="DPAGT11 DPAGT13"; ENDO=DPATG1_all_reads.bam ;;
  WASF1)  C300="WASF11 WASF12";   ENDO=WASF1_all_reads.bam ;;
  *) echo "unknown gene $GENE" >&2; exit 1 ;;
esac
mkdir -p $OUT/bams $OUT/fiberhmm
export FIBERHMM_NO_UPDATE_CHECK=1

# 300-bp constructs: ~6,000 reads per construct per replicate (samtools view -s, seed 42)
for rep in 1 2; do
  B=$D/IGVF_construct_mix/source_bams/IGVF_construct_mix_rep_${rep}.ft.aligned.bam
  for c in $C300; do
    n=$(samtools view -c -F 0x900 $B $c)
    frac=$(python3 -c "print(min(0.999, 6000/$n))")
    samtools view -b -F 0x900 -s 42${frac#0} $B $c > $OUT/bams/${c}_rep${rep}.bam
  done
done
for c in $C300; do
  samtools merge -f $OUT/bams/$c.bam $OUT/bams/${c}_rep1.bam $OUT/bams/${c}_rep2.bam
done
# full-length construct (all reads on the gene's contig) and endogenous reads
samtools view -b -F 0x900 $D/genomic_promoters/merged_bams/$GENE/$GENE.ft.m6a.merged.bam $GENE > $OUT/bams/${GENE}_full.bam
samtools view -b -F 0x900 $D/Hek293_Whole_Genome/target_bams/$ENDO > $OUT/bams/endogenous.bam
for b in $C300 ${GENE}_full endogenous; do
  samtools sort -o $OUT/bams/$b.sorted.bam $OUT/bams/$b.bam && samtools index $OUT/bams/$b.sorted.bam
  echo "$b reads: $(samtools view -c $OUT/bams/$b.sorted.bam)"
  fiberhmm-call -i $OUT/bams/$b.sorted.bam -o $OUT/fiberhmm/$b.fp.bam --enzyme hia5 --seq pacbio --min-mapq 0 \
      -c $CORES --no-qc 2> $OUT/fiberhmm/$b.log
  samtools index $OUT/fiberhmm/$b.fp.bam
done

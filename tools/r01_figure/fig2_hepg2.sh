#!/bin/bash
# Figure 2 draft 5: HepG2 SLC39A4 WT / c.-169A>G plasmid Fiber-seq through the same steps as the HEK293T data
# (projects/SLC39A4_variant_slide/scripts/run_fiberhmm.sh, then occupancy.py unchanged). Rules: fig2_rules.json
# revision_draft5. Usage: tools/r01_figure/fig2_hepg2.sh <work dir> <threads>
set -euo pipefail
W=$1; T=${2:-1}
SRC=/gscratch/stergachislab/bmallo/large_home/data/kopah_restore/SLC39A4_HepG2
FH=/mmfs1/home/bmallo/miniconda3/envs/general/bin/fiberhmm-call
ST=/mmfs1/home/bmallo/miniconda3/envs/MPRA/bin/samtools
PY=/mmfs1/home/bmallo/miniconda3/envs/general/bin/python
OCC=/gscratch/stergachislab/bmallo/large_home/projects/SLC39A4_variant_slide/scripts/occupancy.py
mkdir -p $W/bams $W/footprints $W/logs
for r in 1 2 3 4; do for a in WT MUT; do
  ln -sf $SRC/rep_$r/plasmid_bams_processed/SLC39A4$a/SLC39A4${a}_50_range_exact_reads.bam $W/bams/rep_$r.SLC39A4$a.bam
  ln -sf $SRC/rep_$r/plasmid_bams_processed/SLC39A4$a/SLC39A4${a}_50_range_exact_reads.bam.bai $W/bams/rep_$r.SLC39A4$a.bam.bai
done; done
for r in 1 2 3 4; do for a in WT MUT; do
  o=$W/footprints/rep_$r.SLC39A4$a.fp.bam
  [[ -s $o.bai ]] && continue
  $FH -i $W/bams/rep_$r.SLC39A4$a.bam -o - --enzyme hia5 --seq pacbio --circular --min-mapq 0 --primary -c $T \
     2> $W/logs/rep_$r.$a.fh.log | $ST sort -@ $T -m 2G -T $W/footprints/tmp_$r$a -o $o -
  $ST index $o
  echo "$o $($ST view -c $o)"
done; done
cd $W && rm -f analysis/reads.tsv.gz && OMP_NUM_THREADS=1 $PY $OCC footprints/rep_*.fp.bam

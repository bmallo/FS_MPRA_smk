# Replicate-3 m6A and nucleosome-score provenance — 2026-10-18

The proposed ML≥128 versus ML≥204 comparison cannot test probability
sensitivity on the available replicate-3 input. In the first 1,000 records of
each of the raw ft, aligned, and footprint BAMs, every target m6A ML byte is
240–255. This restriction is present **before alignment and FiberHMM**. Other
modification channels in those records cover the full 0–255 range, so this is
specific to the retained m6A calls, not a blanket compression of every ML
value. These are bounded convenience samples from different BAM orders, not
matched molecules or a distribution estimate for the full run.

The [reproducible audit](../tools/audit_rep3_caller_provenance.py) and its
[JSON output](../results/rep3_caller_provenance_20261018/audit_v2.json) examine at
most 1,000 records from each BAM. The parser consumes every MM group with its
corresponding ML values, including 5mC/5hmC channels, and does not use
`pysam.modified_bases` (which can segfault for these mixed tags). Sample counts:

| Stage | MM/ML records | Target m6A ML bytes | Target range | Nucleosome AQ |
|---|---:|---:|---:|---:|
| Raw ft BAM | 999/1,000 | 221,571 | 240–255 | no MA |
| Aligned ft BAM | 998/1,000 | 222,449 | 240–255 | no MA |
| Footprint BAM | 997/1,000 | 153,375 | 240–255 | 13,112 zeros |

The raw ft BAM header records PacBio `jasmine` 26.1.3, followed by `lima` and
fibertools-rs 0.13.0. The [historical ft script](/mmfs1/gscratch/stergachislab/bmallo/large_home/data/pacbio/r84046_20261002_184948_1_A01/logs/ft_r84046_20261002_184948_1_A01.sh)
runs `ft add-nucleosomes` then `ft clear-kinetics`; it does not run
`ft predict-m6a`. The existing Rep3 processing configuration also specifies
`input.modifications: existing`, so this workflow imports those MM/ML calls.
The original barcoded HiFi BAM named in the ft script is unavailable at that
path, and the PacBio run's Jasmine settings or earlier raw probabilities have
not been recovered. The evidence identifies the earliest *available* stage
with restricted m6A ML values; it does not establish Jasmine's internal
threshold, calibration, or what omitted A/T positions would have scored.

The prior WT sample contained 331,304 nucleosome-class calls, all with first
AQ byte zero. The independent first-1,000-record footprint sample here has
13,112 nucleosome calls, again all zero. The historical circular and linear
FiberHMM logs identify fused apply plus TF recall, model path, `min_llr=5`,
`min_opps=3`, `unify_threshold=90`, and uplift 1.0, but do not report a
confidence-score flag. The **currently installed** FiberHMM code makes
`--with-scores` optional and writes `nq=0` as an “unverified” sentinel when
nucleosome scores are not supplied. This is a plausible explanation of the
historical zeros, not proof that the same code version or option state was
used. Zero AQ must therefore be treated as an unavailable nucleosome quality
measure, not a posterior probability of zero or biological absence.

The tool records BAM identities, full @PG chains, SHA256 of small source
files and caller logs, separate non-m6A ML histograms, and the historical caller
log banners. BAM device/inode/size/mtime/ctime and small-file SHA256 are checked
before and after sampling; full BAM content is not hashed. The tool records its
own source hash, writes atomically, and refuses to replace an existing audit.
To reproduce in an allocated Slurm job, choose a new output path:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /mmfs1/home/bmallo/miniconda3/bin/python \
  tools/audit_rep3_caller_provenance.py --limit 1000 \
  --output results/rep3_caller_provenance_20261018/audit_recheck.json
```

No caller settings, scientific filters, or original annotations were changed.
The next measurement check needs an earlier Jasmine output or primary kinetics
with a documented m6A caller setting, and a separate test of how sequence
context affects both m6A detection and TF recall.

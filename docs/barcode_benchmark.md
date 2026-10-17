# Barcode benchmark

The benchmark compares the production exact-barcode baseline, the optional
`quality_directional` method, and the **unmodified Pacybara clustering core**.
It runs separately from production and does not change processing defaults.
It preserves individual molecules and writes barcode assignments for joining
back to chromatin measurements; it never collapses footprint measurements.

## Reproduce

Run inside an allocated compute job. These commands use bounded, single-threaded
processing and no nested scheduler submissions. The source checkout is small;
the indexed BAM traversal and extraction are the substantial work. The current
benchmark reads the existing replicate-3 alignment, rather than realigning it.

```bash
mkdir -p results/barcode_benchmark_20261016/vendor
git clone https://github.com/rothlab/pacybara.git \
  results/barcode_benchmark_20261016/vendor/pacybara
git -C results/barcode_benchmark_20261016/vendor/pacybara checkout \
  b5e39bd0e2731283eee2ee88c0e0852ddf83fe33

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m fsmpra.barcode_benchmark sample \
  --config config/ldlr_processing.json \
  --bam /gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/fiberseq/HepG2/replicate_3/aligned/PS01898.hg38_LDLR.ft.aligned.bam \
  --fraction 0.01 --output results/barcode_benchmark_20261016/sample

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m fsmpra.barcode_benchmark run \
  --config config/ldlr_processing.json \
  --sample results/barcode_benchmark_20261016/sample \
  --checkout results/barcode_benchmark_20261016/vendor/pacybara \
  --output results/barcode_benchmark_20261016/comparison
```

Use a Python environment with `pysam`. Existing output directories are protected;
choose a new directory for a new run. The Pacybara checkout must be clean and at
the pinned commit. Its GPL source remains external and is invoked as a subprocess.
No downloaded source is incorporated into this repository's implementation.

The 1% random sample proved too sparse for family-level validation. A second,
separately labeled benchmark collects full-depth neighborhoods around 200
repeated and 200 singleton discovery barcodes:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m fsmpra.barcode_benchmark families \
  --config config/ldlr_processing.json \
  --bam /gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/fiberseq/HepG2/replicate_3/aligned/PS01898.hg38_LDLR.ft.aligned.bam \
  --discovery results/barcode_benchmark_20261016/sample \
  --families-per-stratum 200 --output results/barcode_benchmark_20261016/families

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m fsmpra.barcode_benchmark run \
  --config config/ldlr_processing.json \
  --sample results/barcode_benchmark_20261016/families \
  --checkout results/barcode_benchmark_20261016/vendor/pacybara \
  --output results/barcode_benchmark_20261016/family_comparison_final
```

Seeds require Q30 barcode support, sequence-qualified WT/SNV genotypes, and no
observed conflicting qualified genotype in the discovery stratum. Selection
within each stratum uses a fixed barcode hash. Collection includes every primary
plasmid with anchored barcode edit distance <=1 from any seed, regardless of
genotype, quality, length eligibility, or methylation. The full extraction audit
then records those properties. `families/seeds.tsv` records the selection.
These neighborhoods exclude more distant errors and unreadable flanks. They
are enriched for discovered barcodes and cannot estimate library-wide yield.

## What is being tested

Sampling uses a fixed SHA256 function of the read ID across **all primary plasmid
alignments**, with 1% inclusion probability. It does not take the first records
or select by footprint, methylation, genotype, or barcode. Extraction reuses the
production sequence/QC code, including barcode qualities and indels between
anchored flanks. The sampled BAM and rejected observations are retained.

The comparison cohort has confidently observed WT/SNV mutable sequences, passes
sequence QC, and has anchored ACGT barcodes of length 14–16. Barcode quality can
be low; each method decides how to use that evidence. Designed insertion/deletion
**genotypes** and unassigned sequences remain in the extraction audit but are
outside this first clustering comparison. Barcode indels are included.

Pacybara's protein-oriented extraction and translation wrappers are not run.
Our adapter creates exact-barcode preclusters, an ID alias table, nucleotide
genotype tokens, and exhaustive Levenshtein distances up to two. Radius two is
needed for the core's complete-link check with maximum merge distance one.
The adapter fails before exceeding three million expanded read pairs. Raising
the sample fraction may require a different execution strategy.

Core parameters are `maxDiff=1`, `minJaccard=0.2`, `minMatches=1`, `minQual=19`.
Trusted genotype tokens receive Q20, the conservative lower bound guaranteed by
the shared cohort filter. The core's strict `> minQual` comparison retains these
tokens. Barcode qualities are not consumed by this core; its FASTQ quality lines
are placeholders, not calibrated consensus qualities. This is a specified core
benchmark, **not a validation of the full upstream workflow or its defaults**.
Upstream may invoke MUSCLE for a tied consensus; a missing dependency fails the
run rather than silently substituting an algorithm.

All three methods require at least two supporting molecules and an unambiguous
N15 barcode/variant assignment for export. Ambiguous Pacybara consensus barcodes
are rejected even if the core produces separate genotype clusters with the same
barcode. Molecules disagreeing with the inferred genotype remain unassigned.
Consequently training genotype agreement is not an independent accuracy metric.

## Evaluation and outputs

* `sample/sampling.json`: seed, sampling definition, denominator, input provenance.
* `sample/extracted/observations.tsv.gz`: all sampled primary plasmids and QC.
* `comparison/report.json`: cohort sizes and method comparisons.
* `family_comparison_final/report.json`: full-depth targeted-family comparison;
  `heldout_predictions.tsv.gz` records each frozen-map prediction and disagreement.
* `comparison/{real,training,synthetic}/variant_coverage.tsv`: assigned molecules,
  unique/effective barcodes, and largest-barcode fraction for the common cohort.
* `comparison/{real,training,synthetic}/corrections.tsv`: individual correction
  transitions, including Q30 support of the original and destination barcode.
* Each method directory: molecule assignments and headerless gzip TSV
  `barcode → element_id` export. Pacybara additionally retains its inputs, raw
  clusters, command log, source hash, and ID aliases.

The held-out test learns maps from a deterministic 80% split, freezes the export,
then looks up the remaining 20% by barcode alone. Test genotypes are used only
after prediction to count disagreements. This tests barcode/variant linkage
consistency for exact matches to exported representatives. It does not estimate
correction precision, biological reproducibility, or accuracy for unseen barcodes.
For the targeted-family follow-up, all molecules from the discovery sample are
forced into training and excluded from held-out validation, because their
genotypes contributed to seed selection. The remaining molecules use the 80/20
hash split; the overall training fraction therefore need not equal exactly 80%.

Synthetic tests cover low/high-quality substitution errors, barcode insertions
and deletions, genuine neighboring barcodes with identical/different genotypes,
barcode collisions, and a child compatible with two possible parents. Their
known truth exposes false merges; their deliberately adversarial frequencies
must not be treated as empirical error rates in replicate 3.

Real-data corrections between two barcodes each supported by at least three Q30
reads are flagged for inspection. They are **potential** false merges, not proven
errors: genuine library sequences are not independently known. A read subsample
reduces family depth and changes abundance ratios. Full-depth confirmation is
needed before making production corrections or interpreting library complexity.

Run tests with:

```bash
PACYBARA_CHECKOUT=results/barcode_benchmark_20261016/vendor/pacybara \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python -m unittest discover -s tests -p 'test_*.py' -v
```

Without `PACYBARA_CHECKOUT`, the external-core integration test is skipped;
candidate completeness and held-out leakage checks still run.

Atomic publication means some subordinate audit files record temporary execution
paths. Their relative locations and content hashes remain valid after publication;
the sampling manifest identifies the final run, and its source BAM is unchanged.

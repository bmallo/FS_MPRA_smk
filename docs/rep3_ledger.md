# Complete replicate-3 barcode/genotype ledger

This entrypoint reuses the existing indexed HepG2 replicate-3 alignment to produce
a complete primary-plasmid sequence ledger. It applies the same sequence catalog,
QC definitions, and exact-barcode assignment policy used in the benchmark.
It does not realign the raw reads or generate new footprint calls.

## Run

Use a compute allocation, with at least two CPUs. The recorded run used four
allocated CPUs: one reader plus three extraction workers. Heavy intermediate BAMs
and SQLite databases use node-local scratch and are removed after completion.

```bash
XDG_CACHE_HOME=/tmp/fsmpra-smk-cache \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
snakemake -s workflow/ledger.smk --cores 4 --resources mem_mb=12000 \
  --config processing_config=results/rep3_ledger_20261016/config.json \
  processing_python=/mmfs1/home/bmallo/miniconda3/bin/python \
  ledger_run_dir=results/rep3_ledger_20261016 \
  ledger_bam=/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/fiberseq/HepG2/replicate_3/aligned/PS01898.hg38_LDLR.ft.aligned.bam
```

For a new run, copy `config/ldlr_processing.json` into a new run directory first
and update the two run-directory arguments above. `ledger_run_dir` overrides
the general processing config's output directory. Existing completed stages are
protected from direct overwrite. Snakemake can resume from completed stages;
an interrupted extraction stage restarts its traversal. No nested Slurm jobs
are submitted. The supplied paths are specific to this installation.

## Outputs

All paths below are relative to `results/rep3_ledger_20261016/`.

| Output | Contents |
|---|---|
| `ledger/molecule_ledger.tsv.gz` | One row per primary plasmid molecule: observed sequence/genotype, barcode and qualities, QC, labeling summary, accepted barcode assignment and rejection status |
| `ledger/validation.json` | Row conservation, read-ID join validation, assignment/export checks and scope |
| `extracted/observations.tsv.gz` | All sequence observations, including unresolved and off-design genotypes |
| `extracted/qc.json` | Full-contig traversal counts, QC counts, chunk accounting, source provenance and alignment program history |
| `barcodes/barcode_map.tsv.gz` | Observed barcode, representative, genotype, status, total reads and seed flag |
| `barcodes/barcode_variant_support.tsv.gz` | Qualified genotype evidence for each barcode, including competing genotypes at conflicting barcodes |
| `barcodes/molecule_assignments.tsv.gz` | Barcode assignments/rejections keyed by original read ID |
| `barcodes/variant_barcode_coverage.tsv` | Per-variant molecule coverage, observed and supported barcode counts, concentration and effective barcode count |
| `barcodes/mprasnakeflow_assignments.tsv.gz` | Headerless two-column barcode → stable element ID export |
| `catalog/design.fa` and `catalog/oligo_aliases.tsv` | Matching element sequences and original oligo names |
| `code_snapshot/`, `config.json`, `run_manifest.json` | Code/configuration snapshot and execution provenance |

The final ledger separates `observed_variant_id` from `assigned_variant_id`.
An observed ID means the mutable sequence exactly matches a designed element;
`genotype_status` separately records its quality. An assigned ID means the barcode
has an accepted map consistent with that observed sequence. Unresolved observations
are not replaced with WT or rescued from barcode identity.

Exact mapping needs at least two sequence-qualified, barcode-qualified observations
and no competing qualified designed genotype. Conflicting barcodes and singletons
remain in the evidence tables. Physical molecules are never deduplicated by barcode.
Assignments do not use methylation or footprint outcomes.

Coverage cohorts have explicit meanings:

* `sequence_assigned`: all exact catalog sequence matches, including low-quality
  observations. This is a descriptive upper envelope, not a confidence filter.
* `sequence_qualified`: exact-design genotypes with the configured base, pass,
  read-quality, and mapping-quality criteria; includes designed indels.
* `analysis_eligible`: the current WT/SNV, sequence, full-plasmid and labeling
  eligibility criteria. It does not require an accepted barcode assignment.

`observed_barcodes` counts distinct N15 ACGT sequences in each cohort, including
unsupported or conflicting observations. `unique_barcodes` counts only accepted
barcode assignments consistent with the observed variant. Effective barcode count
is the inverse concentration of assigned molecule counts, not a biological
replicate count. The minimum seed-quality setting defines
`seed_quality_molecules`; it is Q30 in this run.

The ledger includes the existing m6A summary at ML≥128 and preserves missing or
unknown modification data. The inferred high-methylation population is not
labeled as an independently prepared control. No occupancy significance or
activator/repressor classifications are produced by this stage.

## Scaling and validation

The reader traverses the indexed plasmid contig once and dispatches bounded
100,000-molecule BAM chunks. Workers call the same production extractor. Completed
compressed tables are concatenated as standard gzip members, with exactly one
TSV header. Global duplicate IDs are rejected during barcode ingestion, which
uses a primary key across the entire table. The final join checks every read ID,
observed barcode, accepted genotype, and total molecule count before publication.

The modification-summary optimization validates MM skip bounds and ML lengths
without creating a position object for every modification. Repeated target groups
fall back to the reference decoder so overlapping positions are counted once.
Barcode support is aggregated with grouped SQL instead of millions of separate
barcode queries. Assignment semantics are unchanged.

Before the full run, optimized extraction reproduced every observation field for
all 87,618 molecules in the random benchmark. Barcode assignments, maps, exports,
and existing coverage fields also matched the full-depth benchmark baseline.
The 38-test suite includes sparse/multiple modification groups, reverse reads,
CIGAR indels, chunk boundaries, gzip concatenation and duplicate IDs across chunks.

This run's completeness refers to primary alignments on the plasmid contig.
Genomic and unmapped reads are not part of the plasmid ledger. Source BAM metadata
and index counts are recorded; the large BAM is not content-hashed. Future
MPRAsnakeflow data must come from the same physical library, and barcode orientation
must be confirmed against its RNA/DNA sequencing design.

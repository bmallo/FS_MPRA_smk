# Raw-BAM processing foundation

This workflow develops the redesigned analysis on HepG2 replicate 3. It starts
with unaligned HiFi reads and produces molecule-level features, a sequence
catalog, barcode coverage, and an MPRAsnakeflow-compatible assignment export.
It does not yet fit statistical models or replace the existing dashboard data.

To generate a complete barcode/genotype ledger from an existing indexed
alignment, use the separate [ledger entrypoint](rep3_ledger.md). It retains the
same extraction and barcode policies while avoiding unnecessary realignment or
footprint calling for this task.

## Run

Use an allocated compute node. The initial reference alignment requires roughly
40 GB requested RAM; the other rules request 0.5–16 GB. The workflow's maximum
per-rule CPU request is four. Supply at least two allocated cores because the
footprint caller uses a worker and driver. Snakemake owns scheduling; no rule
submits nested Slurm jobs. Select resources from the actual allocation.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
snakemake -s workflow/process.smk --cores 4 --resources mem_mb=48000 \
  --config processing_config=config/ldlr_processing.json \
  processing_python=/mmfs1/home/bmallo/miniconda3/bin/python
```

Add `-n` to inspect the DAG without processing data. `processing_run_dir` can
override the configured output directory. Paths in the supplied config are
specific to this installation. Set `scratch_dir` to node-local storage; SQLite
and alignment/caller temporary files are cleaned after normal completion or
exceptions. A hard process termination can leave hidden staging directories.
Choose a new run directory for scientific revisions; preserve its configuration.

The tools currently pinned by version are pbmm2 26.1.99 and FiberHMM 2.12.0.
The Hia5 PacBio model is additionally pinned by SHA256. A version mismatch fails
before that stage executes. These are the available local versions, not claims
that they are optimal or identical to historical processing. Large BAM/reference
provenance records path, size and mtime, not full-file content hashes. Small
configuration/catalog/model files are hashed. Reused `.mmi` indexes are checked
against the reference dictionary; this cannot prove their nucleotide identity.
Set `alignment.index` to null to align directly against the configured FASTA.

## Stages and retained evidence

1. **Preflight:** inspect 1,000 unaligned records and their program history.
   Existing MM/ML m6A calls are used by default. The alternative `predict` mode
   requires fi/fp/ri/rp kinetics and calls `ft predict-m6a`; it is not exercised
   by the replicate-3 smoke test. Missing calls on individual reads are allowed
   and represented explicitly downstream. An inspected input with no m6A calls
   fails rather than silently proceeding.
2. **Catalog:** validate the plasmid reference and design; export one element
   per unique mutable sequence and a separate table of original oligo aliases.
   `design.fa` IDs equal the assignment export's element IDs. Multiple source
   oligos encoding the same sequence cannot be distinguished from that sequence.
3. **Alignment:** retain unmapped reads, align to plasmid plus genomic background,
   then sort and index. Original read IDs and modification tags are retained.
4. **Extraction/QC:** audit all alignment records, count genomic/unmapped records,
   and retain observations for every primary plasmid molecule. Duplicate primary
   IDs fail. Extract the mutable sequence and barcode between aligned flanks,
   preserving inserted and deleted bases. Exact design assignment and quality
   flags are explicit; novel/off-design sequences are retained without forced
   assignment. No genotype consensus rescue is performed yet.
5. **Caller partition:** independently of sequence and methylation eligibility,
   use circular calling for whole-plasmid reads with acceptable clipping and no
   supplementary-alignment tag; use linear calling for other plasmid molecules.
   Missing/unknown modification observations bypass the caller and retain missing
   footprint measurements. Old caller tags are removed before recalling.
6. **Barcodes:** reconstruct an auditable assignment map and report barcode
   support before/after analysis eligibility. This stage does not depend on
   footprint results and does not select a map by the observed chromatin effect.
7. **Footprints/features:** call with the pinned model, require molecule-count
   conservation, then join by read ID. Export reference-coordinate footprint
   intervals, aligned positions, and promoter MSP coverage. All primary plasmid
   observations remain present, including exclusions and missing measurements.

Changing barcode settings reruns barcode assignment and features. Changing
sequencing/methylation eligibility thresholds reruns extraction, barcode coverage,
and features, while retaining alignment and calls. Changing structural boundary
rules or the modification-probability threshold also reruns caller partitioning
and footprint calling. Code changes intentionally invalidate affected workflow
inputs broadly. No p-values or confidence intervals are generated at this stage.

## Initial cohort and measurement definitions

Defaults require exact WT/SNV design sequence, minimum base quality 20 across the
mutable element, at least 10 passes, read quality >=0.99, MAPQ >=20, and a
whole-plasmid structure. These are initial technical criteria, not optimized
scientific cutoffs. Barcode readability is separate from chromatin eligibility.

Methylation bounds default to [0,1]: neither labeling extreme is excluded.
Missing or unknown modification observations fail methylation eligibility.
The summary is called m6A at ML>=128 divided by all A/T opportunities across
the read. Sparse absent calls within declared channels are treated as unmodified;
MM groups marked `?` are not interpreted as fully observed. This is a labeling
summary, not a direct estimate of biological accessibility. Partial-read summaries
cover a different sequence span from full plasmids and must not be pooled blindly.

The primary element remains reference `[1234,1445)` (211 bp). Small TF-sized
footprints are 10–60 bp and nucleosome-sized footprints are >=85 bp. A present
MA tag omitting a track means evaluated zero under the pinned caller's sparse
format. An absent MA tag means missing, never zero. MSP fraction is reported only
when the full element has aligned bases. Footprint intervals and eligible positions
use **reference 0-based, half-open coordinates**. These are descriptive caller
measurements; sequence-dependent detectability still needs calibration.

The local older pysam/htslib crashes when decoding some mixed modification tags
through `modified_bases`. This pipeline uses an explicit MM/ML parser instead,
including multi-modification groups and reverse-read orientation; malformed tag
lengths fail with an exception.

## Barcode baseline and correction benchmark

The default is **exact barcode grouping**, requiring >=2 sequence-qualified reads
and an unambiguous designed genotype for export. One qualified observation of a
different designed genotype marks a collision. Singleton/low-support and conflicting
barcodes remain in audit tables. This conservative baseline is not Pacybara and
is not yet the final barcode reconstruction method.

An optional `barcodes.method: quality_directional` benchmark only corrects one
substitution, requires a unique fixed high-quality seed, compatible high-quality
genotype evidence, low quality at the mismatching base, adequate quality elsewhere,
and a >=5-fold abundance difference. Seeds need >=3 reads with barcode minimum
quality >=30; error positions must be <=19. Child groups are capped at 10 reads.
No transitive merges, indel rescue, or merging of well-supported neighbouring
barcodes is performed. The separate [barcode benchmark](barcode_benchmark.md)
compares this mode with exact grouping and the pinned Pacybara clustering core,
using known-truth synthetic cases and sampled replicate-3 molecules. Independent
library truth and full-depth validation remain necessary before production use.

One barcode can represent many different physical molecules. Their methylation
and footprints are **never collapsed or deduplicated by barcode**. Assignment maps
are inferred using sequence evidence independently of length/methylation cohort
selection. Unique-barcode counts describe observed confident assignments, not
the complete library complexity. Effective barcode count is inverse concentration
of read support, not a number of biological replicates.

## Output interfaces

| Output | Meaning |
|---|---|
| `catalog/design.fa` | Unique mutable sequences with stable element IDs |
| `catalog/oligo_aliases.tsv` | Original design names linked to element IDs |
| `extracted/observations.tsv.gz` | Observed sequences, barcode qualities, genotype/QC flags |
| `barcodes/barcode_map.tsv.gz` | Observed/corrected barcode map with support and ambiguity |
| `barcodes/barcode_variant_support.tsv.gz` | Qualified genotype evidence, including competing variants at conflicting barcodes |
| `barcodes/molecule_assignments.tsv.gz` | Read-level assignment join; one row per plasmid molecule |
| `barcodes/variant_barcode_coverage.tsv` | Molecule and barcode support, imbalance, and cohort |
| `barcodes/mprasnakeflow_assignments.tsv.gz` | Headerless barcode sequence → element ID TSV |
| `footprints/footprints.bam` | Original molecules with recalled footprint annotations |
| `features/molecule_features.tsv.gz` | QC, genotype, barcode, and footprint features joined per molecule |

For the complete replicate-3 ledger with existing caller output, use the
[separate footprint join](footprint_join.md). Its `molecule_footprints.tsv.gz`
also retains individual region-overlapping calls, original lengths and raw AQ
quality values. It explicitly retains ledger rows missing from the annotation
BAM. This entrypoint reuses historical calls and does not run the recall stage
listed above.

The `sequence_assigned` coverage cohort includes observed exact design matches,
including those subsequently rejected for sequencing quality; inspect genotype
status for confidence. `sequence_qualified` applies genotype and sequencing QC
without restricting the designed variant class or chromatin eligibility.
`analysis_eligible` applies all configured primary criteria. `observed_barcodes`
counts observed N15 ACGT sequences, while `unique_barcodes` counts supported
assignments; these differ for singleton and conflicting evidence.

MPRAsnakeflow export uses only supported representative barcode sequences;
sequencing-error aliases are not exported as extra real barcodes. It follows the
documented external assignment interface (`type: file`, `assignment_file`).
The default orientation is reference-forward. Confirm the future RNA/DNA read
orientation before use; do not independently canonicalize reverse complements.
Mappings apply only to the physical library from which they were inferred. A
barcode absent here is not evidence of absence from a future RNA/DNA dataset.

MPRAsnakeflow itself is not run by this workflow. Import of its barcode-level
counts, replicate metadata, QC and activity summaries is a subsequent milestone.
Reference: https://mprasnakeflow.readthedocs.io/stable/1_getting_started/config.html

## Validation

Run `python -m unittest discover -s tests -p 'test_*.py' -v`. New tests cover barcode
indels and orientation, collision handling, singleton retention, quality-aware
correction, missing modification data, mixed MM groups, duplicate IDs, atomic
failure, cohort coverage, assignment export, and missing versus empty footprints.

`results/processing_smoke_20261016` contains a real integration test using the first
10,000 raw replicate-3 records. This is deliberately an engineering smoke test,
not a representative sampling experiment or a full-library barcode analysis.
The original input and previous pilot results remain unchanged.

Completed validation on 2026-10-16: 29 Python tests passed; the raw-data workflow
produced 219 primary plasmid records (218 with evaluated footprints and one missing),
1,417 catalog elements, and 23 supported barcode assignments. Read IDs reconcile
across observation, assignment, and feature tables. The export has unique 15-bp
ACGT barcodes and valid catalog element IDs. A no-change dry run does no work;
a methylation-filter change schedules only extraction, barcode coverage, and
feature export. The saved `validation.json` and logs record these checks.

# Replicate-3 footprint join — 2026-10-17

The complete barcode/genotype ledger is now joined to the existing FiberHMM
annotations. All **8,715,827 primary plasmid molecules matched by read ID**;
there were no missing BAM records, duplicate primary IDs, or mutable-sequence,
barcode, read-length or strand disagreements. Every ledger row is retained.

| Measurement | Molecules |
|---|---:|
| Ledger rows, matched BAM records and output rows | 8,715,827 |
| Present, valid MA annotation | 8,712,120 |
| Missing MA annotation, retained as missing | 3,707 |
| All 211 reference bases aligned | 6,591,700 |
| Current analysis-eligible WT/SNV molecules with usable annotations | 2,594,368 |
| Analysis-eligible molecules with supported barcode assignments | 2,102,906 |

All current analysis-eligible molecules have evaluated annotations, complete
211-bp reference coverage and exact projected design sequence. Thus this join
causes **no additional loss from the current primary analysis cohort**. The
3,707 unannotated molecules lie outside that cohort. Missing annotations are
not converted to unoccupied or inaccessible bases.

## Variant coverage

| Class | Eligible annotated molecules | With supported barcode assignment |
|---|---:|---:|
| WT | 1,243,891 | 985,264 |
| SNVs | 1,350,477 | 1,117,642 |

All 633 designed SNVs have eligible annotated molecules; 632 have eligible
annotated molecules with supported barcode assignments. Median SNV coverage is
1,955 eligible annotated molecules, or 1,619 when a supported barcode is required.
These are molecule counts, not independent barcode or biological replicate counts.

The three low-coverage SNVs identified in the barcode ledger remain limited:

| Variant (oligo 1-based position) | Eligible annotated molecules | With supported barcode assignment |
|---|---:|---:|
| 160A>G | 1 | 0 |
| 222T>C | 2 | 2 |
| 176G>A | 24 | 8 |

Insertions, deletions, unresolved genotypes and other excluded molecules remain
in the joined table, with their annotation and QC fields. They are outside the
current WT/SNV primary cohort. Deletion gaps are absent reference coverage,
not evidence for accessibility or protection.

## What is retained

Each row preserves the original ledger columns and adds annotation status,
aligned reference coverage, individual TF/nucleosome/MSP calls overlapping
the 211-bp element, raw caller quality bytes and optional annotation names.
Call lengths are retained before clipping. Convenience coverage unions and
fractions use 10–60-bp TF calls, ≥85-bp nucleosome calls and all MSPs.
The detailed [field definitions and workflow instructions](footprint_join.md)
describe coordinates, missingness, size thresholds and circular-call limitations.

Under the explicit sparse-output contract, omitted tracks in a present MA tag
are evaluated empty: 417,069 reads omit MSP, 720,670 omit TF and 1,106,721 omit
nucleosome sections. These counts overlap. Another audit field records 3,518
evaluated reads without an AQ tag; absent per-call quality is retained as null.

## Validation and artifacts

The run completed in about 29 minutes on one CPU within the existing Slurm
allocation. A bounded SQLite cache and node-local scratch supported the ID join;
scratch was disposable and the final output was published atomically. No
footprints were recalled. The historical circular and linear caller logs and
sparse-output contract are recorded in the audit.

Five targeted tests passed. Projection also matched the existing per-base
implementation across 10,000 real reads for all three tracks. During the full
run, the BAM index, read identities, sequence checks, source metadata and final
counts reconciled. The per-variant coverage totals independently sum to every
corresponding audit count.

- [Joined molecule table](../results/rep3_footprints_20261017/molecule_footprints.tsv.gz) (1.47 GB compressed)
- [Per-variant annotation coverage](../results/rep3_footprints_20261017/variant_footprint_coverage.tsv)
- [Machine-readable summary](../results/rep3_footprints_20261017/summary.json)
- [Full audit and input provenance](../results/rep3_footprints_20261017/audit.json)
- [Run manifest and validation record](../results/rep3_footprints_20261017_run/run_manifest.json)

The next analysis can aggregate WT/SNV occupancy and MSP-accessibility profiles
while retaining barcode grouping. Confidence intervals and tests still need
to account for barcode dependence, uneven coverage and sequence-dependent
observability. No significance tests or protein identities are inferred by
this join; MSP fraction remains a caller-derived accessibility proxy.

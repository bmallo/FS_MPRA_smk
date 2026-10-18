# Joining existing footprint annotations to the complete ledger

This stage joins the existing replicate-3 FiberHMM BAM to the validated
barcode/genotype ledger by original read ID. It retains every ledger row,
including unassigned genotypes, indels, rejected barcodes and missing annotations.
Neither barcode assignment nor sequence QC is recomputed from footprint outcomes.

## Reproduce

Run inside a compute allocation. This stage uses one CPU, a 128-MiB SQLite cache,
and a 40-GiB database cap on the configured local scratch directory. The recorded
run requests 2 GiB of workflow memory; the database cap excludes its rollback
journal. No Slurm jobs are submitted by the workflow.

```bash
XDG_CACHE_HOME=/tmp/fsmpra-smk-cache \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
snakemake -s workflow/footprint_join.smk --cores 1 --resources mem_mb=2048 \
  --config processing_python=/mmfs1/home/bmallo/miniconda3/bin/python \
  ledger_run_dir=results/rep3_ledger_20261016 \
  footprint_contract=config/rep3_existing_footprints.json \
  footprint_output=results/rep3_footprints_20261017
```

The separate entrypoint consumes the completed ledger without rebuilding it.
Choose a new output directory for a changed contract or implementation. Results
are published atomically after the complete join; scratch is removed on exit.
An interrupted join restarts from the ledger rather than a partial database.

## Outputs

| File | Contents |
|---|---|
| `molecule_footprints.tsv.gz` | All original ledger columns, followed by annotation status, region coverage, individual calls, coverage unions and fractions |
| `variant_footprint_coverage.tsv` | Counts for every designed variant plus an unresolved-genotype row; includes missingness and eligible measured counts with/without barcode assignment |
| `audit.json` | Conserved molecule counts, input metadata and small-file hashes, caller contract, thresholds, coordinate conventions and limitations |

Run logs, source snapshots, pre-run validation results and execution metadata
are stored separately in `results/rep3_footprints_20261017_run/`.

## Coordinates, calls and missingness

The primary interval is reference **[1234, 1445)**: the user-confirmed 211-bp
mutable LDLR element, corresponding to oligo positions 21–231. All exported
reference intervals are zero-based, half-open, and clipped to this region.
`eligible_ref_intervals0` records aligned reference coverage. It is an alignment
denominator, not a statement that every base has equal methylation evidence.
Insertions have no reference positions; deletions and reference skips remain
gaps. Reverse alignments use stored BAM SEQ coordinates, as does the MA tag.

The three `*_ref_calls0` columns contain JSON lists for TF, nucleosome and MSP
calls intersecting the region. Each call preserves:

- `query_start0` and `query_length`: original stored-SEQ start and length before clipping;
- `ref_intervals0`: one or more projected reference intervals, retaining deletion gaps;
- `quality_spec` and `quality`: the caller's quality specification and corresponding raw AQ bytes;
- `name`: the AN annotation name, when provided.

Calls of all lengths remain available. The convenience unions use TF lengths
10–60 bp, nucleosome lengths ≥85 bp, and all MSPs. Thresholds apply before
clipping. Overlapping calls are unioned when calculating covered bases.
Circularly split calls retain their piece lengths and AN names when available;
this stage does not reconstruct their combined lengths.

The existing caller logs record `min_llr=5`, `min_opps=3` and
`unify_threshold=90`, with separate circular and linear processing. These
historical call settings differ from the convenience ≥85-bp size filter.
No new caller threshold is applied to the retained AQ values. Their raw values
remain available for subsequent sensitivity analysis; exact LLRs cannot be
recovered above a saturated quality byte.

`footprint_status` distinguishes:

- `evaluated`: present, valid MA tag. An omitted track is empty under the explicit FiberHMM sparse-output contract.
- `missing_MA`: the molecule exists in the annotation BAM but has no MA tag.
- `missing_bam`: no primary alignment for that ledger ID in the configured plasmid contig.

Missing annotations produce blank call/fraction fields. Evaluated empty tracks
produce `[]` and, when the whole interval is aligned, fraction 0. Missing AQ
quality is `null`, not zero. Malformed MA/AQ/AN relationships fail the run.
The sparse contract is backed by historical caller logs and the installed
writer source; that source is not asserted to be the historical binary.

Fractions are the union of called bases divided by 211 and are blank unless all
211 reference bases align. `footprint_analysis_eligible` additionally requires
the ledger's WT/SNV analysis eligibility, an evaluated MA tag, and exact
agreement of the projected region sequence with the design. The accompanying
`footprint_analysis_eligible_with_barcode` also requires the ledger's supported
barcode assignment. Neither flag changes which rows are retained.

## Identity and validation

The join enforces unique ledger IDs and rejects duplicate or unknown primary
footprint IDs. Read lengths, strands, extracted mutable sequences and extracted
barcodes must agree with the ledger. It also checks catalog/config provenance,
reference sequence identity, BAM index counts and unchanged input metadata.
Large BAM and ledger files are checked by size and modification time, not
content-hashed. Missing BAM rows are retained and counted explicitly.

Five targeted tests cover missing versus empty calls, indel projection,
reverse alignments, quality/name retention, lengths before clipping, overlaps,
malformed annotations, reordered joins, missing reads, duplicate IDs and changed
sequences. Projection also matched the existing implementation for all three
tracks across 10,000 real reads (9,983 evaluated and 17 without MA).

The output supports subsequent barcode-aware occupancy summaries. It does not
yet perform hypothesis tests, identify binding proteins, or update the dashboard.
MSP fraction remains a caller-derived accessibility proxy. Sequence-dependent
observability and within-barcode dependence still need statistical treatment.

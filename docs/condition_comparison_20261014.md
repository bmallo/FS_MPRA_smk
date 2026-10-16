# Reaction-aware LDLR comparison — 2026-10-14

Open `results/conditions_20261014/dashboard/index.html` for the three-run
dashboard. It adds condition labels, a QC table, WT profile overlays, and
allele-preserving switches between runs. The primary interval is now
user-confirmed: the 211-bp mutable LDLR element (oligo positions 21–231).

## Experimental clarification

| Runs | Duration | Hia5 volume | Cells per reaction |
|---|---:|---:|---:|
| HepG2 1 and 2 | 10 min | 1.5 µL | 2 million |
| HepG2 3 | 5 min | 3 µL | 2 million |

These are independent biological preparations under **two reaction
conditions**, not three exchangeable replicates of one condition. Enzyme stock
activity and other preparation details remain incompletely documented.
Volume multiplied by time is not used as a claim of equivalent treatment.

## Processing evidence

The October reanalysis scripts for runs 1–2 explicitly use the replicate-3
NotI-rotated reference and genotype script. All three footprint BAMs contain
the expected 4,647-bp LDLR contig. Circular and linear caller logs document
FiberHMM's fused apply/TF-recall path with Hia5/PacBio models, k=3, min_llr=5,
min_opps=3, unify_threshold=90, and uplift=1. The same sparse MA annotation
contract is now declared for each sample: omitted evaluated tracks are empty;
absent annotations remain missing.

Upstream processing differs:

| Program | Runs 1–2 | Run 3 |
|---|---|---|
| CCS | 8.7.0 | 26.1.0 |
| Jasmine modification caller | 2.7.99 | 26.1.3 |
| fibertools-rs | 0.8.0 | 0.13.0 |
| pbmm2 | 26.1.99 | 26.1.99 |

These versions come from BAM program headers. Logged model paths/settings do
not prove identical historical caller code or model-file contents. Reaction,
sequencing run, software, and biological effects cannot be disentangled using
one preparation under the newer condition.

## Matched measurement rules, unequal coverage

Each run uses a deterministic 1% sample of its complete genotype ledger, the
same seed, full-plasmid classification, at least 10 polymerase passes, exact
agreement with the assigned WT/SNV mutable sequence, and available MA calls.
There is no methylation filter. Existing replicate-3 measurements were reused;
updated experimental annotations are recorded in the new collection without
rewriting its original measurement bundle.

| Run | Imported molecules | Full plasmids¹ | Measured WT/SNV | WT | SNVs with ≥10 molecules | Median SNV n |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 18,950 | 34.1% | 2,180 | 1,026 | 2/633 | 1 |
| 2 | 28,475 | 2.9% | 214 | 107 | 0/633 | 0 |
| 3 | 87,385 | 95.2% | 28,415 | 13,549 | 603/633 | 22 |

¹ Among metrics-matched sampled genotype-ledger molecules, **not** all
sequenced reads or all transfected plasmids. All selected reads matched a
metrics row. Full-plasmid classification is an imported upstream annotation.

Only two SNVs meet the ten-molecule display minimum in both runs 1 and 3;
none meet it in either pair involving run 2. The comparison reports null
correlations for these insufficient intersections. This does not establish
non-reproducibility or a lack of biological effects. These small samples are
unsuitable for variant-effect concordance testing.

WT mean MSP coverage is 0.615, 0.481, and 0.537 for runs 1–3. These are
descriptive caller-derived fractions over the mutable interval, not calibrated
accessibility estimates or evidence that one Hia5 condition is biologically
more accessible. The within-run WT comparator remains separate for every
allele estimate.

## Implementation and validation

- `fsmpra compare` writes a versioned multi-sample collection, a provenance/QC
  report, and sample/allele TSVs. It rejects incompatible references, intervals,
  allele sets, measurement parameters, missingness contracts, duplicate sample
  IDs/BAM sources, and conflicting condition annotations.
- `workflow/compare.smk` packages comparison and dashboard export from existing
  measurements; its default input paths are in `config/ldlr_comparison.json`.
- Fourteen Python tests pass, including metadata-overlay preservation,
  condition assignment, coverage intersections, and incompatible-input failures.
  JavaScript checks validate all three real datasets. The comparison Snakemake
  workflow passes a dry run.
- Headless Firefox interaction checks pass, including reaction labels, the
  cross-run allele table, preservation of the selected allele while switching
  samples, missing-data handling, and the existing profile/coverage controls.
- Source experimental data remain read-only. Full input identities, hashed
  logs, model paths, program versions, and current metadata provenance are
  available in `comparison/comparison.json` and the dashboard's evidence panel.

## Next analysis step

Build a bounded, disk-backed import and measurement path to use substantially
more of runs 1–2, rather than raising a memory cap and loading the entire ledger.
Audit why full-plasmid recovery differs, including read structures and genotype
confidence. Then assess small-footprint observability and null behavior within
each run/condition, preserving all methylation levels. Use the two preparations
under the older reaction for within-condition agreement; treat comparison to
run 3 as a separate reaction/run comparison. No generalized variant discoveries
or treatment-effect estimates are justified by this pilot alone.

# FS_MPRA_smk

**Start here:** [Scientific and implementation handoff — 20 October 2026](docs/PROJECT_HANDOFF_20261020.md).
It summarizes the project from the initial review through the raw-methylation audit,
including unsuccessful statistical approaches, biological candidates, reproduction
instructions and untested next steps. [Compact evidence tables and provenance](docs/evidence/20261020/README.md)
are included in Git; full BAMs, ledgers and generated dashboards remain on Klone.
**No replicate-3 Fiber-seq variant significance has been released.**

**Then read:** [Independent review and recommended direction — 20 October 2026](docs/INDEPENDENT_REVIEW_20261020.md).
It re-checks the evidence snapshot, finds a robust library-level link between promoter
architecture and published activity, and recommends a library-wide naked-DNA control
plus a same-condition replicate before any per-variant claims or construct panel.

> **Raw-BAM redesign:** start with the new [processing workflow](docs/processing_workflow.md)
> for replicate 3. It produces molecule QC/features, barcode assignments and
> coverage, and an MPRAsnakeflow-compatible assignment export. Barcode correction
> remains conservative. The [barcode benchmark](docs/barcode_benchmark.md)
> compares exact grouping, quality-based correction, and Pacybara's clustering core.
> See the [replicate-3 benchmark findings](docs/barcode_benchmark_20261016.md).
> The [complete-ledger workflow](docs/rep3_ledger.md) reuses an indexed alignment
> to generate all primary-plasmid genotype/barcode records and coverage summaries.
> See the [completed replicate-3 ledger and coverage report](docs/rep3_ledger_20261016.md).
> The [footprint join](docs/footprint_join.md) links existing FiberHMM calls to
> this ledger, preserving molecule IDs, barcode assignments and annotation missingness.
> See the [completed replicate-3 footprint join](docs/rep3_footprints_20261017.md).
> The new [statistical specification](docs/rep3_statistics_design.md) and
> [statistics workflow](docs/statistics_usage.md) use barcode-level comparisons
> in fixed windows, with WT A/A and simulation checks before releasing effects.
> These replace the legacy statistical assumptions for the replicate-3 redesign.
> See the [first full-data statistical report](docs/rep3_statistics_20261017.md):
> effect estimates are available, but the initial calibration failed and
> significance labels are withheld.
> The [heterogeneity design](docs/rep3_heterogeneity_design.md) and
> [heterogeneity workflow](docs/heterogeneity_usage.md) characterize observed
> whole-plasmid, distal and element methylation/footprint variation before
> choosing filters or a replacement statistical test.
> See the [completed replicate-3 heterogeneity report](docs/rep3_heterogeneity_20261017.md).
> The [WT footprint-quality design](docs/rep3_footprint_quality_design.md) and
> [audit workflow](docs/footprint_quality_usage.md) inspect original call lengths,
> sequence observability and barcode outcomes, then benchmark Welch and
> studentized permutation tests without releasing variant significance.
> See the [completed WT audit and benchmark](docs/rep3_footprint_quality_20261018.md)
> for caller-score limitations, reproducibility and candidate-test results.
> The [expanded calibration design](docs/rep3_heldout_calibration_design.md) and
> [workflow](docs/heldout_calibration_usage.md) evaluate all fixed endpoints
> against a much larger WT reference, including sparse unequal-variance stress
> tests, smaller effects and interval coverage. The
> [caller provenance audit](docs/rep3_caller_provenance_20261018.md) traces
> retained m6A scores and the uninformative nucleosome-score field.
> See the [completed expanded calibration](docs/rep3_heldout_calibration_20261018.md):
> permutation testing improves false-positive control, but sparse effects and
> confidence intervals still need calibration before variant significance is released.
> The [WT landscape and variant biology report](docs/rep3_landscape_20261019.md)
> now categorizes protection by original call length and compares all designed
> SNVs with WT. Open the [interactive landscape atlas](results/rep3_landscape_20261018/atlas_biology_20261020/index.html)
> for profiles, mutation maps, molecule combinations and published functional context.
> The [variant browser update](docs/landscape_navigation_20261019.md) adds arrow
> navigation and synchronized functional-effect lollipops above each comparison.
> The [genomic orientation audit](docs/genomic_coordinate_audit_20261020.md)
> verifies the forward-strand mapping against GRCh38 and UCSC; the latest atlas
> uses genomic coordinates and preserves the original measurements.
> The [candidate biology report](docs/rep3_candidate_biology_20261020.md) adds
> sequence/A/T tracks, JASPAR/HOCOMOCO candidates, and barcode sensitivity checks.
> See [annotated atlas usage and validation](docs/biology_atlas_usage.md).
> The [raw-call and filtering audit](docs/rep3_raw_candidates_20261020.md)
> examines 27 candidate substitutions at unchanged A/T sites, including
> barcode-independent genotype comparisons and 18 sensitivity settings.
> Open its [interactive audit](results/rep3_raw_candidates_20261020/viewer/index.html)
> for raw-call profiles, footprint profiles, support counts and exports.
> These are descriptive candidate patterns; additional confidence audits are deferred.
> See [landscape workflow usage](docs/landscape_usage.md) to reproduce the analysis.

> **Redesign pilot:** the new, replicate-aware LDLR implementation is available
> through `python -m fsmpra`. See [pilot usage and limitations](docs/pilot_usage.md)
> and the [scientific design plan](docs/redesign_scope_and_plan.md). It imports
> the newer HepG2 library and produces descriptive small-footprint and
> accessibility measurements with a self-contained interactive HTML viewer.
> The workflow described below is the legacy analysis; the pilot does not yet
> replace its statistical analysis.
>
> **Condition comparison:** the pilot now distinguishes the two Hia5 reaction
> conditions and compares runs without pooling. See the
> [comparison report](docs/condition_comparison_20261014.md).
>
> **Full-data processing:** version 0.2 adds bounded SQLite processing and read
> recovery QC. See the [complete-ledger report](docs/full_data_20261016.md).

A Snakemake pipeline for analyzing **Fiber-seq MPRA** (Massively Parallel Reporter Assay) data. It processes PacBio Fiber-seq reads from a saturation mutagenesis plasmid library to measure how single nucleotide variants affect chromatin architecture — protein footprints, nucleosome positioning, and transcription factor binding.

## Overview

The pipeline takes aligned Fiber-seq BAM files through three stages:

```
Raw Fiber-seq BAM
  -> 01_filter_reads.py    -> {sample}.plasmid.bam   (QC, separate genomic/plasmid)
  -> 02_call_variants.py   -> {sample}.tagged.bam     (variant calls + barcode clusters)
  -> 03_analyze_library.py -> {sample}.h5 + .tsv + .pdf (null calibration + testing)
```

1. **Filter reads** — Quality filtering, methylation checks, genomic/plasmid separation, nucleosome tag annotation
2. **Call variants** — Single-pass CIGAR walk for SNV/indel calling, two-pass barcode clustering by Levenshtein distance
3. **Analyze library** — Empirical null calibration from WT reads, coverage-matched variant testing across footprint size bins

## Installation

Dependencies are managed with [pixi](https://pixi.sh):

```bash
pixi install
```

Or install manually:

```bash
pip install pysam numpy scipy h5py matplotlib seaborn tqdm snakemake
```

## Quick Start

### Via Snakemake (recommended)

Edit `config/config.yaml` to define your samples and parameters, then:

```bash
snakemake --cores 8        # run pipeline
snakemake --cores 8 -n     # dry run
```

### Standalone scripts

Each pipeline stage can be run independently:

```bash
# Stage 1 — Filter reads
python workflow/scripts/01_filter_reads.py input.bam \
    --reference ref.fa --output-dir out/ --sample-name MySample

# Stage 2 — Call variants and cluster barcodes
python workflow/scripts/02_call_variants.py filtered.plasmid.bam ref.fa \
    --promoter-region PlasmidName:3184-3501 \
    --barcode-region 3502-3516 \
    --output-dir out/ --sample-name MySample

# Stage 3 — Library-scale analysis
python workflow/scripts/03_analyze_library.py \
    --bam tagged.bam --output-dir out/ --sample-name MySample \
    --min-reads 50 --threads 8
```

### Interactive browser

After the pipeline completes, explore results with the Dash browser:

```bash
python extras/fs_mpra_browser.py --h5 results/MySample/analysis/MySample.h5 --port 8050
```

## Input Requirements

Input BAM files must be PacBio HiFi reads aligned with pbmm2, with Fiber-seq HMM-called footprint tags:

| Tag | Description |
|-----|-------------|
| `ns` | Footprint starts (array, query coordinates) |
| `nl` | Footprint lengths (array) |
| `as` | Nucleosome starts (array) |
| `al` | Nucleosome lengths (array) |
| MM/ML | m6A base modifications |

## Footprint Size Bins

Footprints are classified into biologically meaningful categories:

| Bin | Range | Interpretation |
|-----|-------|----------------|
| sub_TF | 10-19 bp | Sub-transcription factor |
| TF | 20-40 bp | Transcription factor |
| PIC | 41-80 bp | Pre-initiation complex |
| NUC | 81+ bp | Nucleosome |

## Statistical Approach

Variant testing uses **empirical null calibration** rather than parametric assumptions. WT reads are subsampled at a grid of coverage depths to build null distributions of footprint occupancy. Each variant is tested against the null matched to its read depth, with FDR correction across all positions and size bins.

## Directory Structure

```
workflow/
  Snakefile                  # Pipeline DAG
  scripts/
    common.py                # Shared utilities, constants, analysis core
    01_filter_reads.py       # Stage 1: QC and filtering
    02_call_variants.py      # Stage 2: Variant calling and barcode clustering
    03_analyze_library.py    # Stage 3: Null calibration and variant testing
config/
  config.yaml                # Sample sheet and pipeline parameters
extras/
  fs_mpra_browser.py         # Interactive Dash browser
  mpra_analysis_qc.py        # Single-variant deep-dive analysis
  MPRA_qc_visualization.py   # PDF report generation
  mpra_library_browser.ipynb # Jupyter notebook for library exploration
```

> Earlier versions of the pipeline (package prototype, v1 variant caller) are preserved in the initial commit history.

## License

MIT

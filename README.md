# FS_MPRA_smk

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

## License

MIT

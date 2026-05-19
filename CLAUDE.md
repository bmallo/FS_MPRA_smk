# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Repository Is

A **Fiber-seq MPRA (Massively Parallel Reporter Assay)** Snakemake pipeline. It processes PacBio Fiber-seq reads from a saturation mutagenesis plasmid library to measure how single nucleotide variants affect chromatin architecture (protein footprints, nucleosome positioning).

## Directory Structure

```
workflow/
  Snakefile                  # Pipeline DAG (filter -> variants -> analysis)
  scripts/
    common.py                # Shared utilities, constants, analysis core
    01_filter_reads.py       # Filter raw BAM, separate genomic/plasmid, add nc/bn tags
    02_call_variants.py      # Tag reads with variant calls and barcode clusters
    03_analyze_library.py    # Library-scale null calibration and footprint analysis
config/
  config.yaml                # Sample sheet and pipeline parameters
extras/
  fs_mpra_browser.py         # Interactive Dash browser (reads HDF5 from stage 3)
  mpra_analysis_qc.py        # Single-variant deep-dive analysis
  MPRA_qc_visualization.py   # PDF report from QC HDF5
  mpra_library_browser.ipynb # Notebook for library exploration
  fs_mpra_tf_annotate.py     # P4.1: motif -> candidate-TF (JASPAR+HOCOMOCO)
  data/                      # Bundled PWM DBs (JASPAR 2024, HOCOMOCO v12)
results/{sample}/            # Created by Snakemake at runtime
  filtered/                  # Stage 1 output
  variants/                  # Stage 2 output
  analysis/                  # Stage 3 output
```

## Running the Pipeline

### Via Snakemake (recommended)

Edit `config/config.yaml` to define samples and parameters, then:

```bash
snakemake --cores 8
snakemake --cores 8 -n          # dry run
```

### Standalone (for development/debugging)

All pipeline scripts are independently runnable. They share utilities via `common.py` (imported via `sys.path`).

```bash
# Stage 1 — Filter reads, output single merged plasmid BAM
python workflow/scripts/01_filter_reads.py input.bam \
    --reference ref.fa --output-dir out/ --sample-name MySample

# Stage 2 — Tag with variant calls and barcode clusters
python workflow/scripts/02_call_variants.py filtered.plasmid.bam ref.fa \
    --promoter-region PlasmidName:3184-3501 \
    --barcode-region 3502-3516 \
    --output-dir out/ --sample-name MySample

# Stage 3 — Library-scale analysis
python workflow/scripts/03_analyze_library.py \
    --bam tagged.bam --output-dir out/ --sample-name MySample \
    --min-reads 50 --threads 8

# Interactive browser (post-pipeline)
python extras/fs_mpra_browser.py --h5 results.h5 --port 8050
```

### CLI conventions

- All flags use hyphens: `--output-dir`, `--sample-name`
- Every pipeline script accepts `-o/--output-dir`, `-n/--sample-name`, `-v/--verbose`, `-q/--quiet`
- Regions: `CHROM:START-END` (e.g., `PlasmidName:3184-3501`) or `START-END` when chrom is implied

## Pipeline Architecture

### Data flow

```
Raw Fiber-seq BAM
    -> 01_filter_reads.py    -> {sample}.plasmid.bam  (adds nc, bn tags)
    -> 02_call_variants.py   -> {sample}.tagged.bam   (adds PV, VC, PR, BK, CS tags)
    -> 03_analyze_library.py -> {sample}.h5 + _summary.tsv + _motifs.tsv
                                + _cooccupancy.tsv + .pdf
```

Stage 3 also runs the **TF binding-motif layer** (Phase 2): with an
optional `--reference` FASTA it annotates significant clusters with
DNA + sign-consistency + causal-variant distance, then aggregates
across independent SNVs into motif calls (`motifs` HDF5 group +
`_motifs.tsv`). Gated on the calibrated cross-variant FDR.

Optionally (`--enable-co-occupancy`, default off) Stage 3 runs the
**protein co-occupancy module** (Phase 3 / §2.4–2.5): for each
site 1 = a disrupted motif, test whether distal site-2 footprint
changes exceed the WT `P(O2|O1)` channel across the many independent
SNVs disrupting site 1, with a directional-consistency call gate and
the §2.5 secondary-mutation control; `cooccupancy` HDF5 group +
`_cooccupancy.tsv`. WT-vs-WT calibrated (P3.6, exact FDR) but
conservative power on sparse data — read the `underpowered`/`mde`
columns (a null ≠ no dependency). Calibrate per dataset with
`tools/cooccupancy_calibration.sbatch` before trusting calls.

### BAM tags that flow between stages

| Tag | Set by | Meaning |
|-----|--------|---------|
| `nc` | 01_filter_reads | Nucleosome count per read |
| `bn` | 01_filter_reads | Basepairs per nucleosome |
| `ns`/`nl`/`nq` | upstream HMM caller | Footprint starts/lengths/quality (array, query coords) |
| `PV` | 02_call_variants | Variant ID JSON array or `"WT"` |
| `VC` | 02_call_variants | Variant count (0=WT, 1=single, 2+=multi) |
| `PR` | 02_call_variants | Raw (pre-cluster) variant call |
| `BK` | 02_call_variants | Barcode cluster centroid |
| `CS` | 02_call_variants | Cluster size |

### Footprint size bins (defined in `common.py`)

| Bin | Range | Biological interpretation |
|-----|-------|--------------------------|
| `sub_TF` | 10-19 bp | Sub-TF |
| `TF` | 20-40 bp | Transcription factor |
| `PIC` | 41-80 bp | Pre-initiation complex |
| `NUC` | 81+ bp | Nucleosome |

## Key Design Patterns

**Shared utilities in `common.py`**: Constants, tag parsing, `ReadData` class, BAM parsing, NC-matched subsampling, shared memory helpers, cluster detection, and null calibration are all centralized. Pipeline scripts import from it; `extras/mpra_analysis_qc.py` should also import from it (migration pending).

**Empirical null calibration**: WT reads are subsampled at a grid of coverage depths to build a null distribution. Each variant is tested against the nearest null by coverage. This avoids parametric assumptions.

**Two-pass barcode clustering** (02_call_variants.py): Pass 1 collects all barcodes; pass 2 clusters by Levenshtein edit distance using neighborhood hashing, then calls consensus variants per cluster.

**HDF5 interchange format**: Stage 3 writes HDF5 storing null calibration data and per-variant results. The browser and visualization scripts read from these files.

## Dependencies

**Pipeline**: `pysam`, `numpy`, `scipy`, `h5py`, `matplotlib`, `seaborn`, `tqdm`, `snakemake`

**Browser (extras)**: `plotly`, `dash`, `pandas`

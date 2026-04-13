# Fiberseq MPRA Analysis Pipeline

A Python pipeline for analyzing Fiber-seq footprint data from Massively Parallel Reporter Assays (MPRA) to identify how single nucleotide variants affect chromatin architecture.

## Overview

This pipeline analyzes BAM files containing PacBio Fiber-seq data with HMM-called footprints to:
1. Compare footprint landscapes between wild-type and variant sequences
2. Identify statistically significant differential footprints
3. Generate interactive HTML reports for data exploration

## Installation

```bash
# Clone the repository
git clone https://github.com/your-org/fiberseq-mpra.git
cd fiberseq-mpra

# Install with pip
pip install -e .

# Or install dependencies directly
pip install -r requirements.txt
```

## Quick Start

```bash
# Basic analysis
fiberseq-mpra --bam sample.bam --output results/

# With nucleosome filtering (only analyze chromatinized plasmids)
fiberseq-mpra --bam sample.bam --output results/ \
    --nuc-min 10 --nuc-max 25

# With custom configuration
fiberseq-mpra --bam sample.bam --output results/ \
    --config config.yaml

# Generate a default configuration file
fiberseq-mpra --generate-config config.yaml
```

## Input Format

### BAM File Requirements

The input BAM file must contain the following custom tags:

| Tag | Description | Format |
|-----|-------------|--------|
| `ns` | Footprint starts | Array of integers (0-based, query coordinates) |
| `nl` | Footprint lengths | Array of integers |
| `nc` | Nucleosome count | Integer |
| `PV` | Promoter variant | JSON array (e.g., `["245:A>G"]`) or `"WT"` |
| `VC` | Variant count | Integer (0=WT, 1=single, 2+=multi) |

### Example BAM record tags:
```
ns:B:S,0,143,340,524,674    # Footprint starts
nl:B:C,113,178,148,140,13    # Footprint lengths
nc:i:23                      # 23 nucleosomes
PV:Z:WT                      # Wild-type sequence
VC:i:0                       # 0 variants
```

## Configuration

### Default Size Bins

Footprints are categorized into biologically meaningful size bins:

| Bin Name | Size Range | Interpretation |
|----------|------------|----------------|
| TF | 20-49 bp | Transcription factor-sized |
| mid | 50-79 bp | Medium (large TF or small complex) |
| nucleosome | 80-200 bp | Nucleosome-sized |

### Configuration File (YAML)

```yaml
# Filtering parameters
min_variant_reads: 500
nucleosome_range: [10, 25]  # or null for no filtering
require_single_variant: true

# Footprint size bins
size_bins:
  TF: [20, 49]
  mid: [50, 79]
  nucleosome: [80, 200]

# Statistical parameters
fdr_threshold: 0.05
min_effect_size: 0.5  # Minimum |log2FC|

# Output settings
output_format: [html, tsv]
per_variant_output: false
```

## Output Files

```
output/
├── report.html           # Interactive HTML report
├── all_results.tsv       # Complete results table
├── significant_results.tsv  # Significant hits only
├── wt_baseline.tsv       # WT footprint landscape
├── config.yaml           # Configuration used
└── variants/             # (if per_variant_output=true)
    ├── 245_A_G.tsv
    └── ...
```

## Statistical Methods

### Differential Testing

For each position and footprint size bin, we perform a **Fisher's exact test** comparing:
- Number of footprints at that position in WT reads
- Number of footprints at that position in variant reads

### Multiple Testing Correction

**Benjamini-Hochberg FDR correction** is applied to control the false discovery rate across all tests within each variant.

### Effect Size

**Log2 fold change** is calculated as:
```
log2FC = log2(variant_rate / wt_rate)
```

Where rates are the proportion of reads with a footprint at each position.

## Interactive Report Features

The HTML report includes:

1. **Summary Dashboard**: Overview statistics, top hits
2. **Variant Explorer**: Interactive heatmaps for each variant
3. **Results Table**: Sortable, filterable table of all results
4. **WT Baseline**: Visualization of wild-type footprint landscape

## Python API

```python
from fiberseq_mpra.cli.main import run_analysis
from fiberseq_mpra.config import Config

# Create configuration
config = Config(
    min_variant_reads=500,
    nucleosome_range=(10, 25),
    fdr_threshold=0.05
)

# Run analysis
results = run_analysis(
    bam_path="sample.bam",
    output_dir="results/",
    config=config
)
```

## Snakemake Integration

Example Snakemake rule:

```python
rule fiberseq_mpra_analysis:
    input:
        bam="data/{sample}.bam"
    output:
        html="results/{sample}/report.html",
        tsv="results/{sample}/all_results.tsv"
    params:
        outdir="results/{sample}",
        nuc_min=10,
        nuc_max=25
    shell:
        """
        fiberseq-mpra --bam {input.bam} --output {params.outdir} \
            --nuc-min {params.nuc_min} --nuc-max {params.nuc_max}
        """
```

## License

MIT License

## Citation

If you use this tool in your research, please cite:
[Citation information to be added]

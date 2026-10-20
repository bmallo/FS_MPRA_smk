# Descriptive landscape workflow

The latest [variant browser](../results/rep3_landscape_20261018/atlas_navigation_20261019/index.html)
adds previous/next arrows and synchronized published-functional-effect lollipops.
See [navigation and read-count normalization](landscape_navigation_20261019.md).
The original analysis and atlas below remain archived as the initial landscape.

The [design](rep3_landscape_design.md) and
[completed findings](rep3_landscape_20261019.md) explain the analysis. Open
`results/rep3_landscape_20261018/atlas/index.html` in a browser for the local
interactive atlas. It uses embedded data and no external services. Keep the
atlas directory together for its table/figure download links; the viewer itself
also works as a standalone HTML file. About 12 MB of embedded data can take a
moment to load. No server or installation is needed.

The source cohort is the existing replicate-3 supported-barcode WT/SNV set.
Indels and unsupported barcodes are outside this stage. This is a run-specific
analysis with strict conservation checks, not a generic new-sample processor.
No additional inference or confidence audit is performed.

On an allocated Slurm compute node:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
export MPLCONFIGDIR=/tmp/fsmpra-matplotlib XDG_CACHE_HOME=/tmp/fsmpra-smk-cache
/mmfs1/home/bmallo/miniconda3/envs/snakemake/bin/snakemake \
  -s workflow/landscape.smk --cores 1 --resources mem_mb=4096
```

Each stage refuses to overwrite an existing output directory. For a clean
reproduction, copy `config/ldlr_landscape.json`, change its `output_dir` to a new
workspace path, and pass `--config landscape_config=/absolute/path/to/copy.json`
to Snakemake. Keep the original configuration intact for its recorded audits.
The equivalent stages are:

```bash
/mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.landscape \
  --config config/ldlr_landscape.json
/mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.landscape_atlas \
  --config config/ldlr_landscape.json
/mmfs1/home/bmallo/miniconda3/bin/python -m tools.summarize_landscape_patterns \
  --atlas results/rep3_landscape_20261018/atlas \
  --output results/rep3_landscape_20261018/patterns
```

`landscape/` exports the full 634-allele universe (WT plus 633 SNVs), with
per-base profiles for eight annotation/length categories, category coverage and
call-presence summaries, eight molecule-combination fractions, WT call-length
and boundary tables, an NPZ array bundle, and an audit. The variant without
eligible data has NaN arrays and blank table measurements, not zero occupancy.
Projected WT call boundaries at element edges are clipped/censored and must not
be interpreted as complete physical footprint boundaries. Circular splitting
can affect original annotation segment lengths.

`atlas/` exports the local HTML viewer, JSON data, all variant summaries,
25 region comparisons per variant (four WT peak neighborhoods and 21 fixed
windows), separate descriptive correlations with each published MPRA source,
and PNG/SVG figures. The 50-barcode/100-molecule toggle only filters the variant
table; the mutational heatmap deliberately includes all variants with data.
Low-coverage extremes can dominate that heatmap's color scale. The table ranks
absolute observed differences for navigation, not confidence.

`patterns/` exports two expressly post hoc examples selected after inspecting
the atlas: the three substitutions at position 125 and all 24 substitutions
at positions 143–150. These are fixed exports of exploratory selections, not
independent motif validation or statistical discoveries.

The primary effect is the difference between means of per-barcode molecule
fractions. Covered bases are not independent trials. Call-presence fractions
refer to any overlapping call anywhere in the element; they differ from mean
covered-base fractions. The three state bits are small TF, nucleosome-sized
(85–250 bp), and long nucleosome-class (>250 bp) protection. State zero does not
exclude other TF or short nucleosome-class calls and does not mean naked DNA.

Published MPRA records are imported from the previously sequence-verified
catalog. Current original-file identities and variant sequences are checked;
the two source files remain separate and are not assumed independent biological
replicates. No caller filters are chosen to improve MPRA concordance.

Validation includes nine Python unit checks, an actual-atlas JavaScript DOM/
canvas smoke test, full-depth conservation, independent TF/MSP reconciliation,
and comparison of nucleosome-length categories to the prior union summary.
The JavaScript check is not a visual browser test. Static scientific figures
were visually inspected. The completed run used direct CLI stages; Snakemake
was dry-run checked and consequently reports missing execution metadata for
the three rules. Per-stage source audits remain available.

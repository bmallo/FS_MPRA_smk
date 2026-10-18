# Replicate-3 Fiber-seq heterogeneity extraction

The extraction streams the frozen existing FiberHMM BAM in primary alignment order alongside `results/rep3_footprints_20261017/molecule_footprints.tsv.gz`. It validates each read ID, stored query length, strand, and the BAM's recorded path/size/mtime before writing. Inputs remain unchanged. The selected feature table includes every exact-design, sequence-QC-passing WT or SNV molecule, regardless of its older length, methylation, footprint, or barcode flags. `retention.tsv` counts all primary molecules, including unresolved genotypes and overlapping QC failures.

Use the existing allocated compute job, one CPU, and the configured Python with NumPy and pysam:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 \
  /mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.heterogeneity \
  --config config/ldlr_heterogeneity.json \
  --output results/rep3_heterogeneity_20261017/features
```

The Snakemake entrypoint runs extraction before `fsmpra.heterogeneity_summary`:

```bash
snakemake -s workflow/heterogeneity.smk --configfile config/ldlr_heterogeneity.json --cores 1
```

Both rules request one thread and 8 GiB. For a bounded smoke run, choose a separate output directory and add `--limit-primary 1000` to the direct Python command. Its `audit.json` has `complete: false`, and full-cohort conservation checks apply only to an unrestricted run. Output directories publish atomically after all checks pass and cannot overwrite an existing run.

`molecule_heterogeneity.tsv.gz` retains the joined genotype, barcode, quality, and prior eligibility fields. Each of `whole`, `element`, and `distal` has aligned bases, coverage fraction, aligned A/T count, m6A calls and fractions at ML 128 and 204, and union fractions for small TF, nucleosome, and MSP annotations. The whole interval is `[0,4647)`, the element `[1234,1445)`, and distal is `[0,1134) ∪ [1545,3269) ∪ [3334,4647)`. All coordinates are zero-based and half-open. Regional fractions use only aligned bases or aligned A/T opportunities; coverage is reported even for partial alignments. Sparse MA with an absent track means evaluated zero, while missing MA yields blank footprint fractions. Missing MM/ML yields blank methylation outcomes. A target MM group marked `?` retains observed positive counts but has blank fractions because omitted opportunities are unknown.

The table also reports projected whole-plasmid MA nucleosome **segment** count, median original call length, and median spacing between projected call midpoints. Each MA call counts once even if an indel splits its projection into several spans; its midpoint uses the outer projected bounds and therefore approximates location across the gap. Circular split MA calls are not reunited; spacing excludes the wrap. TF counts distinguish all projected calls from those with original query length 10–60. Distal tile quantiles summarize observed ML-128 methylation fractions in fixed 100-bp reference bins wholly within distal intervals and containing at least 20 aligned A/T sites. These are descriptive call-pattern summaries, not labeling-efficiency estimates.

`audit.json` records source identities, full or partial traversal counts, selected counts, coordinates, thresholds, and the older eligible-cohort conservation checks. On a full run, the extractor requires 8,715,827 matched primary reads, reconciles 2,660,013 selected WT/SNV molecules to the existing variant-barcode coverage table, conserves the prior 2,594,368 eligible and 2,102,906 barcoded footprint-eligible molecules, and checks complete-element footprint fractions against the joined table. No new genotype or eligibility class is assigned.

## Descriptive summaries

After extraction, run the summary alone with:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  MPLCONFIGDIR=/tmp/fsmpra-matplotlib \
  /mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.heterogeneity_summary \
  --config config/ldlr_heterogeneity.json \
  --features results/rep3_heterogeneity_20261017/features \
  --output results/rep3_heterogeneity_20261017/summary
```

The configuration and catalog must match the extraction audit. A partial
extraction requires the explicit `--allow-partial` flag, and its figures are
marked as partial. The complete run checks the expected cohort and barcode
counts. Both stages refuse an existing output directory; use a new destination
for a revised analysis. Snakemake accepts `heterogeneity_config` to select a
different JSON contract and `output_override` to select a new output root.

The summary contains:

- `barcode_heterogeneity.tsv.gz`: depth, per-feature means and observed
  within-barcode population SD. SD is blank when fewer than two valid molecules
  are available; it is not an estimated biological variance component.
- `variant_heterogeneity.tsv.gz`: molecule and equal-barcode means, descriptive
  differences from WT, coverage and missingness for each variant and cohort.
- `histograms.tsv.gz`, `histogram_boundaries.tsv.gz`, `count_histograms.tsv.gz`
  and `joint_histograms.tsv.gz`: fractional and discrete distributions, exact
  zero/one masses, joint distributions and missing weight.
- `variant_state_summaries.tsv.gz`: fixed distal methylation/nucleosome-coverage
  bins, their proportions, unique barcodes, and conditional element TF/MSP
  coverage. These are numerical display bins, not inferred biological states.
- `summary.json`: provenance, conservation checks, WT summaries, barcode depth,
  missingness, regional completeness and within-barcode spread.
- `wt_heterogeneity.png/.svg` and `variant_heterogeneity.png/.svg`: observed WT
  variation and descriptive variant differences. Sparse display bins and SNVs
  with fewer than 50 barcodes are omitted from the respective plots, but remain
  in the tables.

Molecule weighting treats every eligible molecule equally. Equal-barcode
histograms give every barcode total weight one, divided over all its eligible
molecules. They retain the spread among molecules rather than plotting only
barcode means. Within a numerical bin, these original weights are normalized
over the molecules in that bin: a barcode with half its molecules in the bin
contributes half a unit. A barcode can appear in several bins, so unique-barcode
counts cannot be added across bins. Missing values retain their own weight.

Fractional quantiles in JSON are approximated by 40-bin histogram midpoints;
consult the separate exact zero/one proportions for point masses. Discrete
call-count quantiles use the weighted empirical CDF. The TF heatmap uses a
data-dependent scale beginning at zero; the MSP scale is zero to one.

No output assigns a cellular compartment, labels a molecule as naked DNA,
estimates true reaction efficiency or TF residence time, or supplies new
p-values or significance calls. Observed m6A and inferred footprint features
share the same underlying assay, so correlations between them are not
independent confirmation of chromatin state. WT identifies the mutable-element
sequence; it does not establish that every base of the plasmid is identical.

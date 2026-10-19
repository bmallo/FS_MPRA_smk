# Replicate-3 WT footprint-quality audit

This audit inspects the original FiberHMM `MA` calls and raw `AQ` bytes on a deterministic sample of 10,000 supported WT barcodes. It retains every molecule in each sampled barcode's previously defined `footprint_analysis_eligible_with_barcode` cohort. The existing genotype, barcode, and primary QC calls are frozen. This is descriptive inspection of one WT sequence, with no new significance calls or plasmid-compartment assignments.

Run inside a Slurm compute allocation with one CPU. The configuration is `config/ldlr_footprint_quality.json`; its source paths, seed, barcode count, and band boundaries are checked by the extractor. The complete workflow produces `sample`, `summary`, and `benchmark` under `results/rep3_footprint_quality_20261018`:

```bash
snakemake -s workflow/footprint_quality.smk \
  --configfile config/ldlr_footprint_quality.json --cores 1
```

The three stages can also run independently, in order:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/fsmpra-matplotlib
/mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.footprint_quality \
  --config config/ldlr_footprint_quality.json \
  --output results/rep3_footprint_quality_20261018/sample
/mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.footprint_quality_summary \
  --config config/ldlr_footprint_quality.json \
  --sample results/rep3_footprint_quality_20261018/sample \
  --output results/rep3_footprint_quality_20261018/summary
/mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.footprint_quality_benchmark \
  --config config/ldlr_footprint_quality.json \
  --summary results/rep3_footprint_quality_20261018/summary \
  --output results/rep3_footprint_quality_20261018/benchmark
```

The independent validator reconciles sampled barcode selection and depths with the earlier aggregates and checks the endpoint matrix against the sampled molecule bitmaps. The plotting command reads the completed benchmark table and writes separate WT A/A and synthetic-null false-positive figures (PNG and SVG), a known-effect power table, and its own provenance record:

```bash
/mmfs1/home/bmallo/miniconda3/bin/python tools/validate_footprint_quality_results.py \
  --sample results/rep3_footprint_quality_20261018/sample \
  --summary results/rep3_footprint_quality_20261018/summary \
  --aggregates results/rep3_statistics_20261017/aggregates \
  --output results/rep3_footprint_quality_20261018/independent_validation.json
/mmfs1/home/bmallo/miniconda3/bin/python tools/plot_footprint_quality_benchmark.py \
  --benchmark results/rep3_footprint_quality_20261018/benchmark \
  --output results/rep3_footprint_quality_20261018/benchmark_plots
```

The benchmark figures show rejection counts out of 200 contrasts per row, exact binomial 95% intervals, and the nominal 0.05 line for Welch and studentized permutation. WT A/A contrasts are repeated barcode draws conditional on this WT sample; synthetic equal-mean nulls are a separate simulation. The figures make no pass/fail or variant-discovery claim. Each output command refuses to replace an existing path.

`sample/audit.json` records the SHA256 barcode rule (`seed|barcode`), barcode-disjoint A/B split (`split|seed|barcode`), source identities, expected depths, matched molecule and call counts, and completion status. Every sampled molecule has weight `1 / eligible_barcode_depth`, so each barcode has total weight one. Extraction verifies that all sampled molecule depths agree with the frozen barcode metadata and that every sampled primary BAM read was found. It retains original call query lengths, projected reference intervals, raw AQ quality bytes, and query-level A/T and ML-thresholded m6A counts. Nucleosome bands are `<85`, `85–250`, `251–500`, `501–1000`, and `>1000` bp; TF bands are `<10`, `10–60`, and `>60` bp. Circular split calls remain separate.

`summary/element_profiles.tsv.gz` gives 211 base-position rows for each split, observed nucleosome-call group, and weighting. The columns distinguish aligned positions, observed A/T opportunities, ML-128 m6A calls, opportunities without called m6A, and TF/nucleosome/MSP callability and covered weight. An A/T site without a called m6A is an observed opportunity, not proof of unmethylated DNA. Empty evaluated calls are zero; blank feature values remain missing. `sequence_context_opportunities.tsv.gz` adds centered trinucleotides and distance to the nearest reference A/T base. The profiles and `element_quality_profiles.png/.svg` show the A/B recurrence, numerical nucleosome-call groups, and methylation opportunities. A/B groups have no shared barcodes but share the same WT element sequence, so agreement is not independent TF validation.

The AQ sensitivity applies to scored TF and nucleosome calls when **their first recorded AQ byte is ≥100**. Later bytes in a multi-Q specification measure other call properties and are not part of this threshold. A molecule with any relevant call lacking a usable first AQ byte is excluded from that track's AQ-qualified per-base denominator; the corresponding missing-molecule counts are in `summary.json`. MSP calls have no defined first-byte score and are marked unscored, with no AQ sensitivity profile. In this sample, all 331,304 nucleosome calls have first AQ byte zero; the ≥100 rule therefore retains none and cannot rank nucleosome-call quality. AQ values have no cross-track calibration, so a TF threshold must not be interpreted as the same evidence level for nucleosomes or MSPs. `call_length_bands.tsv` separates original, AQ-qualified, missing-AQ, and unscored calls, along with query A/T and m6A totals and zero-A/T call counts. `call_observability_histograms.tsv.gz` preserves the per-call count distributions by track and length band. These raw scores are a sensitivity threshold, not calibrated binding probabilities.

`barcode_endpoint_matrix.npz` contains ordered barcode IDs, A/B splits, eligible depths, 22 fixed endpoints, and float64 `presence` and `coverage` matrices. For each molecule with aligned bases in an endpoint, presence is one if a selected original call covers any aligned base and zero otherwise; coverage is the fraction of its aligned endpoint bases covered by the call union. Both values are averaged over valid molecules within each barcode. Missing endpoints are `NaN`, with denominators in `valid_molecules`. The matrix uses 21 overlapping 20-bp TF windows stepped by 10 bp, including the terminal window, and one whole-element MSP endpoint. `barcode_endpoints.tsv.gz` exposes the same values in long form. This per-molecule presence cannot be reconstructed from a pooled per-base profile.

`group_feature_summaries.tsv.gz` records valid and missing weight for original and length-band features; `nucleosome_length_bands.png/.svg` shows original nucleosome-call lengths. `representatives/` contains whole-plasmid and element plots for deterministic reads chosen within frozen-feature groups before BAM inspection. `summary.json` records split and group conservation, source identities, and A/B TF-profile correlation. Figures and tables are descriptive and do not imply a naked-DNA class, methylase efficiency, TF residence time, or new primary QC thresholds.

Every stage publishes an output directory atomically and refuses to replace an existing directory. A bounded extraction pilot can be made with `--barcode-limit N` and a separate output path. Its audit has `complete: false`; the summary requires `--allow-pilot` and labels the outputs as a pilot. The benchmark requires a complete 10,000-barcode summary. Snakemake accepts `footprint_quality_config` and `output_override` to select a separate contract and new output root.

# Replicate-3 WT footprint quality and candidate statistics — 2026-10-18

**The audit and fixed benchmark are complete. Studentized permutation is a
promising candidate for further calibration, but variant significance remains
withheld.** The observed footprint profile is reproducible across WT barcode
partitions, yet substantially sensitive to caller score. These are separate
measurement and statistical issues; improving the test does not validate the
underlying calls.

The [design](rep3_footprint_quality_design.md) fixes sampling, endpoints and
benchmark settings. The [usage guide](footprint_quality_usage.md) documents
the three-stage Snakemake workflow, independent validator and plot commands.
Astra reviewed the design and results; Sol implemented the stages and validators.

## Sample and measurement checks

The deterministic sample contains **10,000 WT barcodes and all 37,075 eligible
molecules belonging to them**, drawn from 268,283 supported WT barcodes. The
barcode-disjoint partitions contain 5,018 and 4,982 barcodes. The BAM scan
examined 8,715,827 primary records, matched every expected molecule, and retained
1,031,171 original annotations. These are sampled WT results, not a new full
variant analysis. No endpoint was missing in the benchmark pools.

All 37,075 molecules cover the full 211-bp mutable element. An independent
implementation recomputed the hash selection, depths, splits, molecule-level
presence and coverage, and all 22 barcode coverage endpoints. Comparison with
the previous full-data per-base arrays gave a maximum absolute discrepancy of
**2.22 × 10⁻¹⁶**; comparison with the sampled molecule bitmaps was exact.

With each barcode contributing equal total weight:

| Measurement | Sample estimate |
|---|---:|
| Mean element coverage by 10–60-bp TF calls | 2.61% |
| Mean element MSP coverage | 52.38% |
| Molecules with any element 10–60-bp TF call | 22.02% |
| Molecules with no whole-plasmid nucleosome call ≥85 bp | 14.59% |
| Local TF-call presence within that no-nucleosome group | 30.56% |
| Exactly one ≥85-bp nucleosome-class segment, longer than 1,000 bp | 5.99% |

These estimates closely resemble the preceding full-cohort heterogeneity
summary. The last row denotes an annotation pattern, not a single physical
nucleosome. Among all 331,304 original nucleosome-class segments, 72.30% are
85–250 bp, 12.85% are 251–500 bp, 4.32% are 501–1000 bp, 2.74% exceed 1000 bp,
and 7.78% are shorter than 85 bp. Circular splitting can affect segment lengths;
these call-count proportions differ from molecule prevalence.

## Reproducibility, score sensitivity and observable sequence

The A/B element TF profiles have Pearson r=0.9987 and Spearman r=0.9935.
Both partitions share the same WT sequence and caller, so this demonstrates
repeatability within this preparation rather than independent evidence for
TF binding. The [profile figure](../results/rep3_footprint_quality_20261018/summary/element_quality_profiles.png)
shows the partition comparison, observed nucleosome-call groups, score
sensitivity and per-base A/T methylation opportunities.

Keeping TF calls with first AQ byte ≥100 retains 93,021 of 198,329 whole-plasmid
10–60-bp calls (46.90%). Mean element coverage falls from **2.61% to 1.66%**.
Some peaks, particularly near element offset 100, shrink more than others.
This sensitivity is informative, but AQ=100 is not established here as the
correct biological filter. Scores depend on the caller's evidence model and
are not calibrated binding probabilities. The additional TF quality bytes
describe boundary uncertainty and were not included in this threshold.

**Every sampled nucleosome-class call has first AQ byte zero.** Thus the same
threshold would discard all 331,304 nucleosome-class segments and provides no
ranking within that track. Do not interpret the zeros as demonstrated absence
of nucleosomes or use this cross-track filter. MSP calls have no score at all.
The historical caller logs and current installed writer documentation are
archived in the run's `caller_provenance.json`; the current code is supporting
evidence, not a byte-identical historical caller snapshot. The full sampled
first-byte distributions are in `raw_quality_histograms.json`.

The 10–60-bp TF calls contain a median of 10 A/T bases (10th/90th percentiles
5/18); none contains fewer than three. This is consistent with the documented
minimum-opportunity setting, not an independent validation. About 74.6% have
no retained m6A call inside their original query span. The ML≥128 and ML≥204
element positions are identical for every sampled molecule. Consequently,
that threshold pair still cannot test meaningful probability sensitivity on
this input. Raw-probability retention and sequence-dependent detection need
separate investigation.

Twenty-four deterministic representative plots show the extremes and mixed
patterns: [long protection](../results/rep3_footprint_quality_20261018/summary/representatives/00_exactly_one_nuc_longer1000.png),
[no ≥85-bp nucleosome call with local TF](../results/rep3_footprint_quality_20261018/summary/representatives/05_no_nuc_with_localTF.png),
and [mixed plasmid annotations](../results/rep3_footprint_quality_20261018/summary/representatives/07_other.png).
They illustrate why compartment, true labeling efficiency and protein identity
cannot be read directly from a single annotation class.

## Candidate-test results

The fixed benchmark completed **7,200 contrasts**, evaluating both methods on
each contrast: 4,000 WT-versus-WT comparisons, 2,400 simulated equal-mean nulls,
and 800 known-effect simulations. Each individual scenario has 200 iterations.
WT comparisons use 100 or 400 pseudo-variant barcodes versus 1,000 disjoint
WT controls. The new Welch diagnostic includes sparse outcomes without the
previous sparse-support gate, so this is not an exact rerun of the original
release-gated calibration. The permutation test uses the same studentized
mean-difference statistic and 499 label permutations.

Selected fixed comparisons below report false-positive rejection rates at a
nominal 5% threshold. All scenarios, counts and exact 95% intervals are in the
[WT figure](../results/rep3_footprint_quality_20261018/benchmark_plots/wt_aa_false_positive.png),
[simulation figure](../results/rep3_footprint_quality_20261018/benchmark_plots/synthetic_null_false_positive.png)
and `benchmark/summary.tsv`.

| Null comparison | Welch, A / B | Permutation, A / B |
|---|---:|---:|
| TF [120,140) coverage, 100 vs 1,000 barcodes | 11.5% / 9.0% | 6.0% / 5.5% |
| TF [140,160) coverage, 400 vs 1,000 | 7.0% / 6.0% | 6.0% / 6.0% |
| Element MSP coverage, 400 vs 1,000 | 6.5% / 8.0% | 4.5% / 6.0% |
| Synthetic mean .02, unequal dispersion/depth, 100 vs 1,000 | 14.0% / 11.5% | 0.5% / 1.5% |

Permutation improves several problematic comparisons, but not every result is
resolved. Its largest WT rejection rate is **8.5% (17/200, 95% interval
5.03–13.26%)** for [120,140) presence in split B at 100 barcodes. This is one
of many diagnostic scenarios, not an independently selected formal failure
test; nevertheless it prevents claiming uniformly demonstrated 5% control.
The strongly conservative sparse simulation also matters: reducing false
positives can come with reduced sensitivity to small real effects.

TF rejections remain predominantly apparent losses: 186 negative versus 39
positive for Welch, and 148 versus 35 for permutation. These totals combine
correlated endpoints and must not be treated as independent evidence counts.
Unequal allocation of rejection probability between directions is possible
with a skewed statistic; it warrants direction-specific calibration and power
checks rather than assuming a valid test must have a 50:50 split.

Both methods detected every simulated ±0.10 effect in these 200-iteration
scenarios, including complete losses. This establishes neither power for small
effects nor confidence-interval coverage. The permutation p-values cannot be
smaller than .002, which is inadequate to validate the small tails required
for the previous 13,926-hypothesis family. The 1,000-barcode control arm also
does not reproduce the full WT reference size. No new variant p/q-values,
confidence intervals or significance labels have been released.

## Recommended next stage

Keep current primary eligibility and original calls frozen while developing
inference. Retain TF score sensitivity and observed chromatin groups as
separate descriptive views; do not choose a filter because it produces more
variant hits. Preserve both barcode-level presence and mean coverage.

The next bounded statistical stage should test studentized permutation on
previously unused WT barcodes, across all fixed endpoints and the actual
variant/control barcode-size imbalance. Include rarer outcomes, reversed
variance differences, barcode depth–outcome dependence, realistic small gains
and losses, and explicit confidence-interval validation. Define the scope and
Monte Carlo precision before running; assess multiple-testing tails separately.
The current results nominate this method for that work, not for release.

In parallel, trace upstream m6A probability filtering and the historical
nucleosome-score encoding. Develop a sequence-aware caller check that preserves
A/T and local sequence context; A/B recurrence on identical WT sequence is
insufficient. Before motif inference, verify that variant effects are not
primarily changes in methylation opportunities or caller sensitivity.

## Artifacts and validation

Outputs are under `results/rep3_footprint_quality_20261018/`: `sample/` retains
molecules, original calls and representatives; `summary/` contains profiles,
length/observability distributions and float64 barcode endpoint matrices;
`benchmark/` contains every contrast and summary; `benchmark_plots/` contains
PNG/SVG figures and the known-effect table. Audit records preserve input
identities and completion flags; source mutation checks bracket computation.

Nineteen focused tests passed across extraction, summarization, inference and
independent validation. All endpoint checks passed, and the principal figures
and representative examples were visually inspected. Each processing stage
used one CPU in Slurm job 41830288: extraction took 6m38s (323 MiB peak), summary
1m45s (450 MiB), and benchmark 1m11s (92 MiB). Commands, logs, snapshots and
independent check reports are in `results/rep3_footprint_quality_20261018_run/`.
Execution used the direct CLIs; the final Snakemake dry run checks the completed
workflow state and reports no work remaining. Snakemake has no execution
metadata for these direct runs; use the stage audits and archived logs for their
provenance. Original ledgers, annotations and release decisions are intact.

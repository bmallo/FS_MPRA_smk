# Expanded replicate-3 calibration: fixed design

This stage follows the 2026-10-18 WT footprint-quality audit. Preserve the
original annotations, barcode assignments, eligibility and endpoint grid.
Use the fixed configuration `config/ldlr_heldout_calibration.json`; do not tune
the sample, seeds or filters after seeing results. The purpose is to assess
false-positive rates, effect sensitivity and interval coverage before any
variant significance release.

## Cohort and outcomes

Exclude exactly the prior audit's 10,000 WT barcodes from the supported primary
WT cohort. The remaining 258,283 barcodes contain 948,189 eligible molecules.
They are disjoint from the recent audit sample, but **not a globally untouched
validation population**: the earlier full-data calibration sampled from all WT
barcodes. They also come from the same biological preparation.

Reconstruct per-molecule presence and covered-base fraction from the existing
joined footprint table, not by recalling the BAM. Average within barcode,
retaining one barcode as one observation. Use 21 overlapping 20-bp TF windows
and whole-element MSP coverage. Evaluate 43 measurements: TF presence and
coverage for each window, plus MSP coverage. Export both metrics for all 22
endpoints for reuse, although MSP presence is not a benchmark target.
Check all barcode depths, full 211-bp alignment and empty-versus-missing
annotations. Independently reconcile coverage with the prior per-base arrays.

## WT A/A comparison and computational strategy

For each pseudo-variant size 50, 100, 400 and 1,000 barcodes, draw 1,000 random
subsets without replacement. Compare each subset with its entire complement
in the 258,283-barcode pool. Arms are disjoint within a contrast; barcodes may
recur across contrasts. The reference is approximately 96% of the full WT
cohort, approaching the actual imbalance rather than using only 1,000 controls.

Compare Welch inference with studentized permutation inference. For each size,
draw a separate bank of 9,999 random subsets from this same fixed pool,
recomputing the Welch denominator for each subset and its complement. The
observed and reference-bank random streams are independent. Share each bank
across endpoints and observed contrasts to avoid repeatedly processing the
large reference. Draw only the smaller arm; obtain complement moments by
subtraction from centered pool moments. Do not allocate a permutations-by-pool
array. Validate this optimization against direct complements and independent
test calculations on bounded fixtures.

Two-sided permutation p-values are `(1 + exceedances) / 10000`, counting
ties in absolute statistic. Handle pooled constants, zero standard errors and
infinite statistics explicitly. Reusing a bank induces dependence among
p-values. Exact binomial intervals for repeated observed rejection rates are
conditional on **both the WT pool and this bank** and omit reference-bank
Monte Carlo uncertainty. They do not establish biological reproducibility or
family-wide error control. Report .05 and .01 rejection rates and directional
rejections; do not assume a skewed two-sided test must split errors 50:50.

Permutation exactness requires exchangeability of labels; studentization does
not guarantee finite-sample validity under equal means with unequal
distributions. The [methodological basis](https://arxiv.org/abs/1304.5939)
therefore motivates the separate stress tests below. Minimum p=.0001 is still
insufficient to validate the tails needed for the prior 13,926-test family.

## Stress tests and realistic effects

Each simulation uses 10,000 control barcodes, 200 iterations, and 999 independent
label permutations per contrast. This remains a bounded approximation to the
full reference size. Generate a latent barcode probability from a beta
distribution with specified mean and intrabarcode dispersion; generate a
binomial molecule count and divide by barcode depth.

The 32 null scenarios cross means .002, .02, .1 and .5; variant sizes 100 and
400; dispersion pairs (.02,.2) and (.2,.02); and independent versus associated
barcode depths. Independent depths are geometric with success probability .3
for variant and .5 for control. Associated depths are
`1 + Poisson(1 + 20 * latent_probability)` for both arms. Conditional expected
barcode fractions remain the latent probabilities, so the equal-mean null is
preserved despite depth association. These are stress tests, not claims that
one latent model explains the experiment.

The 18 alternative scenarios cross sizes 100, 400 and 1,000 with six
(control mean, variant mean) pairs: (.002,0), (.02,0), (.02,.015), (.02,.025),
(.1,.08) and (.1,.12). Use dispersion .02 versus .2 and independent depths.
These include complete losses and changes of 0.5 or 2 percentage points,
which are smaller than the prior ±10-point check. Record effect direction,
nonzero support, rejection at each nominal level and interval coverage.

## Confidence intervals are a separate question

Evaluate ordinary 95% Welch intervals, clipped to the feasible difference
range [-1,1], including their zero-standard-error degeneracies. Do not label
these as permutation intervals or assume that permutation p-values validate
them. Report unconditional truth coverage, width, directional noncoverage and
uninformative outcomes.

Include a deliberately conservative bounded-outcome reference interval:

`delta ± [sqrt(log(4/.05)/(2*n_variant)) + sqrt(log(4/.05)/(2*n_control))]`,

clipped to [-1,1]. Apply a two-sided Hoeffding bound to each [0,1] barcode mean
at error .025, then a union bound for the difference. This supplies a transparent
coverage baseline under independent barcode sampling; the finite-population
bound also applies to uniform sampling without replacement. See
[Hoeffding (1963)](https://www.tandfonline.com/doi/abs/10.1080/01621459.1963.10500830).
The interval is expected to be wide; coverage alone is not a useful-method
criterion. Correlated barcodes across cells or other experimental units can
invalidate the independence premise. No simultaneous intervals or automatic
interval-method adoption are claimed.

## Provenance and completion

Separately inspect historical caller commands, score encoding, BAM headers and
bounded raw/processed MM/ML samples. Distinguish documentary evidence from
inference based on currently installed code. Do not restore absent probability
values by treating missing calls as measured zeros, and do not reinterpret
nucleosome score zero as evidence against nucleosomes.

Use source checks before and after computation, atomic output directories,
explicit partial-run status, one CPU per stage and the existing Slurm allocation.
Run bounded numerical/unit tests and a timing pilot before the complete
benchmark. Report limitations and failures alongside improvements. Original
variant results and their withheld-significance status remain unchanged.

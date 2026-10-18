# Replicate-3 statistical specification — 2026-10-17

This specification supersedes the legacy inference assumptions for the new
replicate-3 pipeline. The dataset has already been explored; these choices are
fixed before the new variant-effect calculation, not a claim of preregistration.

## Question, population and limits

Estimate how a designed SNV changes small-footprint coverage and promoter MSP
coverage relative to WT in this one Fiber-seq preparation. The immediate target
is **caller output in the supported-barcode, analysis-eligible cohort**. This
is narrower than protein occupancy in all transfected cells. We have one
biological preparation under this reaction condition; thousands of barcodes
cannot replace biological replicates or identify unrecorded cell-level effects.

The primary cohort is the 2,102,906 joined molecules with
`footprint_analysis_eligible_with_barcode=1`. The 2,594,368-molecule cohort
without the barcode requirement is a descriptive selection-sensitivity cohort.
Do not discard the high-methylation population, subtract it as a naked-DNA
control, or adjust the primary contrast for methylation or nucleosome occupancy.
Those variables may be consequences of the regulatory change itself.

## Estimands and units

For barcode b and endpoint e, let y_be be the mean measured fraction across
all eligible molecules carrying that barcode. The primary genotype mean is the
unweighted mean of y_be across its observed barcodes. Delta is variant minus WT;
multiply by 100 for percentage points. Thus one abundant barcode does not count
as many independent construct backgrounds. Retain every molecule in the source
table and use all of them to estimate its barcode's mean.

The target is the average **observed supported barcode**, not the complete
original construct library. Barcode inclusion and exact-map filtering may
select constructs. Report molecule-weighted effects on the same cohort, and
unbarcoded-cohort molecule-weighted effects, as distinct sensitivities. They
are different estimands when barcode abundance and phenotype are related.

Treat different barcodes as approximately independent construct clusters for
within-sample inference. This is an assumption supported by the construction
description, not proven independence of cells, transfections, or every possible
barcode error. The conservative exact barcode map is frozen before outcomes.

## Fixed endpoints and multiplicity

The primary interval is reference [1234,1445), 211 bp. Preserve base-resolution
TF (10–60-bp caller intervals) and MSP profiles for visualization.

Inferential endpoints are fixed without variant-specific peak selection:

1. TF coverage fraction in 20-bp windows at relative starts 0,10,...,190, plus
   the terminal window at 191. This gives 21 overlapping windows covering the
   element. Each molecule contributes its union-covered bases divided by 20;
   positions within a molecule are never treated as independent trials.
2. Promoter MSP fraction across all 211 bases, one endpoint.

The family is all 633 designed SNVs × 22 endpoints = **13,926 hypotheses**.
Use one joint Benjamini–Yekutieli adjustment (BY, q=0.05), including untestable
members internally as p=1; export their p/q as missing with an explicit reason.
BY accommodates arbitrary dependence of valid marginal p-values, including
shared WT and overlapping windows. It does not repair invalid marginal tests.
The inferential claim concerns a window's mean coverage, not a precise footprint
boundary, the mutated nucleotide itself, or an identified DNA-binding protein.
Base-resolution profiles and subsequently selected motifs are descriptive;
they do not inherit window-level error control.

## Initial uncertainty model: auditable cluster-summary inference

The first candidate uses unequal-variance Welch contrasts on barcode means.
For G variant and H WT barcodes, SE² = s_v²/G + s_w²/H. Use Welch–Satterthwaite
degrees of freedom, two-sided t p-values, and pointwise 95% t intervals for delta
(clipped to the parameter range [-1,1]). Include WT uncertainty. The common WT
is reused; it is not downsampled and its dependence is handled in the joint BY
family. No molecule-level binomial or Fisher test may pretend reads from a
barcode are independent. No base-level binomial likelihood may treat covered
bases in one footprint as independent.

This is an approximate, transparent baseline, not a claim that Welch is optimal
for every endpoint. Defaults require ≥50 barcodes and ≥100 molecules in each
arm. For an endpoint, each arm must have ≥10 nonzero barcode means and ≥10
barcode means below one, and strictly positive empirical variance. Otherwise
retain point estimates/coverage but omit CIs and p-values with a sparse or
degenerate status. In particular, a complete apparent loss can be biologically
interesting while remaining outside this initial inferential method; it must
not be labeled 'no effect.' Do not invent zero-width intervals at boundaries.

Compare this baseline against hierarchical beta-binomial or cluster bootstrap
methods when richer site-level occupancy definitions are established. A
one-random-molecule-per-barcode exact binary test can be a conservative sparse
endpoint benchmark; it is not silently substituted here. Avoid fitting a
random-effect model merely because it is more complex: its estimand, boundary
behavior and calibration still need validation.

## Observability and interpretation

An SNV can change A/T opportunities, local sequence context, methylation
calling, and footprint detection. Store ref/alt bases, A/T gain/loss status,
distance from mutation to window and whether the window overlaps the mutation.
Calculate a secondary effect after omitting the mutation ±2 reference bases
from both variant and WT denominators. Keep primary effects unmasked so the
central scientific question remains visible. This masking is a diagnostic;
it does not prove that sequence bias has been removed, because caller boundaries
and evidence can extend beyond the omitted positions.

No primary adjustment for global methylation, nucleosome state, MPRA activity,
or motif annotations is allowed. Those are secondary mechanistic analyses with
separate targets and multiplicity. Fixed-window results may support within-run
caller associations; binding changes require sequence-aware validation.

## Calibration and release criteria

Freeze the seed, endpoint grid, thresholds and family before computing effects.
Produce candidate p/q values separately from a supported-effect flag. The latter
must remain false unless the statistical calibration stage passes; a failed
calibration must be reported, not fixed by trying seeds or loosening thresholds.

The root review will implement these checks independently of the Sol engine:

- WT barcode A/A contrasts: sample 50,100,400,1000 complete barcode families
  without replacement, compare each pseudo-variant with its disjoint WT
  complement, and repeat 200 times per size. Use the same endpoint and testing
  functions. Never split molecules of one barcode across arms. Evaluate each
  endpoint separately, and report testability, nominal false-positive rate and
  exact binomial uncertainty. There is no variant peak selection to repeat.
- Hierarchical simulations: barcode heterogeneity, unequal read depth,
  informative depth, rare and common outcomes, unequal group sizes, and known
  increases/decreases. Report unconditional nominal rejection rates, fractions
  not testable, and CI coverage among testable contrasts. Nominal power is not
  an FDR-adjusted minimum detectable effect.
- A calibration stratum with ≥100 testable iterations must have the upper
  exact 95% bound for unconditional type-I error ≤0.10 at nominal alpha=0.05.
  Simulation CI coverage among testable iterations must have lower exact 95%
  bound ≥0.90. Insufficiently testable strata remain explicitly unresolved;
  they cannot validate sparse inference. Report the number of assessed strata.
- For this first release, any assessed failure blocks all supported-effect
  flags. Any unresolved null simulation also blocks a global pass: assigning
  real comparisons to a validated outcome regime has not yet been implemented.
  An endpoint additionally needs an assessed WT A/A size and no assessed A/A
  failures. This deliberately conservative release rule can withhold inference
  even when many individual strata behave well; report the reasons separately.
- Test zero-variance handling, rare outcomes, barcode imbalance, sign reversal,
  fixed family size with missing members, BY against an independent reference,
  row-order invariance within floating-point tolerance, and count conservation.

These are operational development gates, not mathematical proof of calibration
or validation of tiny p-value tails. We will also report real WT A/A p-value
quantiles and any failed strata. If the gate fails, publish descriptive effects
and diagnostics with candidate inference clearly withheld.

Sequence-aware caller calibration remains a separate gate for biological
interpretation: simulate the same protected states on WT and alternate
sequences, generate labeling/modification observations using documented model
assumptions, and rerun the caller; benchmark against actual independent control
data when available. WT A/A alone cannot test this bias. Motif matching and
published MPRA agreement are downstream context, not calibration targets.

## Implementation contract

Implement new `fsmpra/statistics.py`, `workflow/statistics.smk`, and focused tests.
Root owns `fsmpra/statistics_calibration.py`, this design and the run report.
Configuration: `config/ldlr_statistics.json`. Do not alter the existing joined
table, catalog, barcode map, or prior pipelines. Use bounded, explicit memory
capacity (default maximum 1,000,000 barcode groups; two 211-position uint32
coverage arrays plus counts), one CPU per stage, and atomic outputs. Read the
joined gzip once, parse only the selected coverage intervals for eligible rows,
and reconcile total/cohort counts to the join audit. Fail on unsupported barcode
identity, genotype inconsistency, malformed intervals or resource-cap overflow.

Persist barcode IDs, variant IDs, molecule counts and per-position TF/MSP covered
counts in reusable arrays. Also retain all-eligible per-variant molecule totals
and coverage sums for selection sensitivity. Use these arrays for all subsequent
estimation/calibration; do not repeatedly reread 8.7 million joined rows.

Pure engine API for independent calibration:
`welch_contrast(x, y, min_barcodes=50, min_nonzero=10, min_nonone=10)` accepts
2D arrays (barcodes × endpoints), variant first and WT second, and returns arrays
`mean_variant`, `mean_wt`, `delta`, `se`, `df`, `ci_low`, `ci_high`, `p_value`,
`status`. Ineligible uncertainty/p values are NaN. The molecule-count gate is
applied by the inference stage, separately from this API.

Outputs: barcode arrays and metadata/audit; base-resolution profile table;
fixed-endpoint effect table with coverage, effect, CI, candidate p/q, statuses
and observability annotations; masked and weighting sensitivity effects;
calibration reports; a final release status. No dashboard, motif calling or
functional classifier is implemented in this milestone.

## Methodological sources

- [Kahan et al., cluster estimands (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9908044/)
  distinguishes equal-cluster and equal-individual weighting when cluster size
  is informative. We apply that distinction to barcodes and molecules; this
  does not turn this experiment into a randomized cluster trial.
- [Ashuach et al., MPRAnalyze (2019)](https://link.springer.com/article/10.1186/s13059-019-1787-z)
  models barcode structure in RNA/DNA MPRA measurements. It motivates keeping
  that structure for future functional integration; its count model is not
  automatically appropriate for these bounded chromatin summaries.
- [SciPy's false-discovery-control documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html)
  states the dependence distinction between BH and BY and cites the original
  Benjamini–Yekutieli (2001) theorem. Our software checks BY against that
  independent implementation; marginal Welch approximation remains a separate
  assumption.
- [Stergachis et al., Fiber-seq (2020)](https://doi.org/10.1126/science.aaz1646)
  establishes the single-molecule chromatin-stenciling measurement underlying
  this experiment. The sequence-observability sensitivity and calibration
  requirements here are project-specific design decisions, not an assertion
  that the original paper validates this variant analysis.

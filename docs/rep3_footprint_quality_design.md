# Replicate-3 footprint-quality audit and candidate-test benchmark

This diagnostic follows the completed heterogeneity characterization. It asks
whether existing small-footprint annotations are interpretable, and whether a
candidate test improves false-positive control for sparse barcode outcomes.
It does not select biological filters or release variant significance.

## Frozen sample and measurements

The configuration is `config/ldlr_footprint_quality.json`. Select the 10,000 WT
barcodes with the lowest SHA256 of `20261018|barcode` from the existing supported
primary cohort. Retain every eligible molecule belonging to those barcodes.
Selection precedes inspecting their footprint phenotypes. Split barcodes into
A/B by an independently prefixed hash; these are disjoint technical partitions,
not biological replicates. Summaries give both molecule weights and weights of
1 / eligible barcode depth, so each barcode contributes equal total weight.
Conditional groups retain those original weights; barcodes can contribute to
more than one observed group.

Inspect original query-coordinate call lengths, including tails: TF <10,
10–60 and >60 bp; nucleosome <85, 85–250, 251–500, 501–1000 and >1000 bp.
These are descriptive annotation bands, not assignments of protein identity.
Circular calls may have been split at the read origin. Project intervals to the
reference using CIGAR blocks, retaining observed/aligned denominators and
missingness. Reconcile the existing TF, nucleosome, MSP and methylation
measurements molecule by molecule against the previous feature extraction.

For each original call, retain raw quality bytes and count A/T opportunities
and retained m6A marks within its stored-query span. Plot whole-plasmid calls
and the 211-bp mutable element, including A/T opportunities and m6A positions.
Choose up to six representative reads per group by a deterministic read-ID
hash: no nucleosome call with local TF, no nucleosome call without local TF,
exactly one >1000-bp nucleosome-class call, and other. These examples are
illustrative; their balanced display is not an estimate of group prevalence.

Compare A/B per-base profiles and observed nucleosome-call groups. Overlay
sequence context and observable A/T positions. Shared sequence-dependent
artifacts can recur in both partitions: recurrence does not validate TF binding.
Raw first quality byte >=100 is a descriptive sensitivity analysis, with
missing scores kept distinct from failures. It is not a calibrated probability.
The installed FiberHMM writer documents a scaled LLR for TF quality, but is not
a frozen historical software snapshot. The historical logs report `min_llr=5`,
`min_opps=3`, `unify_threshold=90`, Hia5/PacBio mode and k=3. Archive this
supporting provenance separately from the immutable call input.

## Barcode outcomes and fixed comparison

For each molecule and each of the 21 fixed overlapping TF windows, calculate
(1) whether any 10–60-bp TF call overlaps the window and (2) the fraction of
aligned window bases covered by such calls. Average each measurement within
barcode. Presence must be calculated before averaging: thresholding a
barcode's mean per-base profile would answer a different question. Also retain
MSP coverage of the whole 211-bp element. Export float64 endpoint matrices,
barcode depths, callable molecule counts and partition labels.

The benchmark fixes five targets in advance: TF presence and coverage in
element offsets [120,140) and [140,160), and element MSP coverage. These TF
windows showed problems in the previous calibration. Compare ordinary Welch
inference with label permutation using the studentized Welch statistic. For
each of A/B, draw 200 independent random contrasts of 100 or 400 pseudo-variant
barcodes versus 1,000 disjoint control barcodes. Barcode sets may overlap across
iterations, but never between arms within a contrast. Both methods use the
same contrast; neither is tuned on A or B. Report endpoint missingness and the
common complete-case pool explicitly.

The 1,000-barcode control arm is a bounded diagnostic choice. It does not match
the much larger full WT reference cohort; passing here would not establish
calibration at that more extreme arm-size imbalance.

For permutation, shuffle whole barcode labels, recomputing the unequal-variance
standard error for every draw. Use 499 random draws, a two-sided statistic and
`p = (1 + exceedances) / 500`, counting ties. Pooled constant outcomes are
uninformative and receive p=1. The minimum possible p is 0.002. Random-label
permutation is finite-sample exact under exchangeable distributions, not
generally under equal means with unequal variances. Studentization provides
asymptotic robustness under suitable conditions; this is not a guarantee for
rare finite samples. See [Chung and Romano](https://arxiv.org/abs/1304.5939) and
[Phipson and Smyth](https://pubmed.ncbi.nlm.nih.gov/21044043/).

Additionally simulate barcode fractions with equal means of .02, .1 and .5,
unequal latent-propensity dispersions (.02 versus .2), and unequal independent
geometric depth distributions (.3 versus .5). Use the same arm sizes and 200
iterations. Simulate known losses (.1 to 0) and gains (.1 to .2) at 400 versus
1,000 barcodes to check direction and nominal power. These simple simulations
do not reproduce all biology or caller artifacts. Report null rejection rates
at .05, exact binomial 95% intervals, uninformative contrasts, nonzero support,
and known-effect detection. A/A intervals describe repeated random draws
conditional on the sampled WT pool, not uncertainty across experiments.

## Decision boundary

No real variant p/q-values or effect confidence intervals are produced here.
Two hundred null iterations and p-value resolution .002 cannot validate the
small tails needed for the earlier 13,926-hypothesis family. If a candidate
improves nominal control, the next stage must evaluate sparse support,
barcode-depth dependence, unequal-variance stress tests and tail resolution
before any release. If neither works, preserve the failure rather than tuning
seeds or filters until one passes. Keep the total variant-associated effect
separate from future comparisons conditioned on distal chromatin measurements.

All stages use one CPU within the existing Slurm allocation, bounded memory,
atomic outputs, pre/post source checks and explicit partial-run labels.
Validation includes independent comparison of sampled barcode mean coverage
against the previous full-data per-base count arrays. No source annotations,
ledger assignments, primary eligibility rules or existing release decisions
are modified.

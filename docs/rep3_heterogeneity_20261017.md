# Replicate-3 heterogeneity — 2026-10-17

The full characterization is complete. WT molecules show broad variation in
observed methylation and protection, including variation within the same
supported barcode groups. Two findings matter immediately: molecules without
nucleosome calls can still have small footprint calls, and a single very long
protected annotation cannot be interpreted as one nucleosome. Simple filtering
on nucleosome count would therefore mix very different molecular patterns.

Astra specified and reviewed the analysis; Sol implemented the extraction,
summary modules, workflow and focused tests. Astra performed independent
reconciliation and the supplemental annotation diagnostics. See the
[design](rep3_heterogeneity_design.md) and [usage](heterogeneity_usage.md).
This is descriptive analysis of the five-minute HepG2 replicate-3 preparation.
It does not pool the other reaction conditions or release statistical discoveries.

## Data and weighting

| Population | Molecules | Supported barcodes |
|---|---:|---:|
| All primary plasmid alignments checked | 8,715,827 | — |
| Sequence-qualified exact WT/SNV designs | 2,660,013 | Not required |
| Previous footprint-analysis cohort | 2,594,368 | Not required |
| Previous cohort with supported barcodes | 2,102,906 | 568,354 |
| WT in the supported-barcode cohort | 985,264 | 268,283 |
| SNVs in the supported-barcode cohort | 1,117,642 | 300,071 |

The broader cohort includes 65,645 molecules excluded by the previous filters
(2.47%). Length/alignment QC excludes 65,345; methylation QC excludes 323;
23 fail both. Missing MA footprint annotations occur in 337 broader-cohort
molecules. Missing measurements are preserved as missing, not scored as zero.

Unless otherwise stated, the findings below use the previous supported-barcode
cohort and give each barcode total weight one, divided among its molecules.
This retains molecular variation while limiting the influence of abundant
barcodes. There are 633 SNVs; 630 have at least 50 supported barcodes and appear
in the variant scatterplots. This is a display threshold, not a significance gate.

## What the WT data show

The element is the complete 211-bp mutable interval. Distal measurements exclude
that element, 100 bp on either side, and the barcode with 25-bp flanks.

| WT measurement | Result |
|---|---:|
| No small footprint call overlapping the element | 78.16% |
| Mean element coverage by 10–60-bp footprint calls | 2.60% |
| Mean element MSP coverage | 52.59% |
| No element MSP coverage | 20.97% |
| Mean distal observed m6A fraction | 34.10% |
| Mean distal nucleosome-call coverage | 47.11% |
| No nucleosome calls anywhere on the aligned plasmid | 14.57% |
| Nucleosome-call count: 10th percentile / median / 90th percentile | 0 / 8 / 17 |

These are annotation measurements. A small footprint call is not proof of a
particular TF, and its covered fraction is not instantaneous binding occupancy.
MSPs can contain small footprints; these tracks are not mutually exclusive
fractions that must sum to one. The five-minute reaction can erase transient
protection, but this dataset alone cannot quantify that loss.

The [WT figure](../results/rep3_heterogeneity_20261017/summary/wt_heterogeneity.png)
shows a broad continuum of methylation/protection patterns with concentrations
near its extremes. The strong relationship between observed methylation and
called protection is partly built into footprint inference, so it is not an
independent validation of biological chromatin states.

### No nucleosome calls does not establish protein-free DNA

The 14.57% of WT barcode weight with no nucleosome calls has mean distal m6A
fraction **66.91%** and mean element MSP coverage **99.98%**. Nevertheless,
**29.79% of that group has a small footprint call in the element**.

Those calls might include real protein protection and calling/labeling
artifacts. Their presence makes it inappropriate to automatically label the
whole group as naked or extranuclear DNA. Cellular origin remains unmeasured.
Also, highly accessible molecules do not need an observed m6A fraction near one.

### Long protected segments complicate nucleosome counting

**5.93%** of WT barcode weight has exactly one nucleosome-class annotation
longer than 1,000 bp. In this group, mean distal m6A fraction is **0.27%** and
mean distal nucleosome-call coverage is **99.75%**. That is extensive unresolved
protection, not evidence for one histone octamer.

Across all WT molecules, **11.74%** of barcode weight has a median
nucleosome-call length greater than 500 bp; **8.74%** exceeds 1,000 bp. Low
labeling, extensive protein protection and merged annotations can contribute
to these patterns. The present data do not identify their relative contributions.
Call count, call length and coverage must therefore be considered together.
These length cutoffs describe the outputs; they were not used to filter reads.

### Variation remains within barcode groups

WT barcodes have median depth three, with 10th–90th percentiles of two to six
molecules. Of 268,283 WT barcodes, 266,752 have at least two eligible molecules.
Among those repeated barcode groups:

| Feature | Groups with observed within-barcode variation |
|---|---:|
| Nucleosome-call count | 260,192 / 266,752 = 97.54% |
| Distal observed methylation fraction | 266,601 / 266,752 = 99.94% |
| Element small-footprint coverage | 146,977 / 266,752 = 55.10% |
| Element MSP coverage | 243,504 / 266,752 = 91.28% |

Variation is scored when observed population SD exceeds 10⁻⁶, avoiding numerical
round-off. The mean within-barcode SD is 17.83 percentage points for distal
methylation and 30.39 points for element MSP coverage. These are descriptive
spreads among observed molecules, not estimates separating technical and
biological variance. Barcode groups are not independent biological replicates
or established cell lineages. WT also fixes the mutable sequence, not every
base of the plasmid backbone.

## Filtering and measurement checks

All sequence-qualified molecules completely cover the mutable element.
Whole/distal coverage needs a less literal completeness criterion: the median
whole-reference coverage is 99.914%, corresponding to four uncovered reference
bases. Only 31 primary barcoded molecules cover every reference base, yet
**99.886% cover at least 99%** of the reference. The cause of the small gaps is
not established here. A requirement for exactly 100% whole-reference coverage
would discard almost the entire cohort.

The broader sequence-qualified cohort includes substantially partial plasmids:
97.58% cover at least 99% of the reference, and its first coverage percentile is
31.59%. Whole-plasmid call counts on these fragments are not directly comparable
to counts on nearly complete molecules. Keep the existing alignment/length QC
as the primary baseline and use the broader cohort as an explicit sensitivity
analysis. Regional fractions use actual aligned bases or aligned A/T opportunities.

WT population means are fairly stable across the existing filters. Element
small-footprint coverage is 2.621% in the broader cohort and 2.600% in the prior
all-eligible cohort; MSP coverage is 52.774% versus 52.581%. Within the supported
WT cohort, molecule weighting versus equal-barcode weighting gives MSP means of
52.622% versus 52.591%. This does not establish negligible filtering or weighting
effects for every individual variant.

The proposed m6A cutoffs, ML ≥128 and ≥204, give **identical reported regional
means for every variant/cohort**. A bounded check of the first 1,000 annotation
BAM records found retained target m6A ML bytes only from 240 to 255. Consequently,
these two cutoffs do not provide an informative probability-threshold sensitivity
analysis on this input. Future raw-BAM processing should audit retained
probabilities and upstream thresholding before selecting such comparisons.
The bounded check does not establish the byte distribution of the entire BAM.

## Consequences for variant analysis

The [variant figure](../results/rep3_heterogeneity_20261017/summary/variant_heterogeneity.png)
shows descriptive differences across the 630 displayed SNVs. Equal-barcode mean
element TF-coverage differences span −0.715 to +0.946 percentage points; MSP
differences span −5.51 to +6.36 points. These ranges do not establish significance.
Across these observed means, element MSP differences correlate with distal m6A
differences (Pearson r = 0.694). Shared biological variation and correlated
measurement/sampling noise can both produce this pattern.

The next statistical development should preserve two questions: the total
variant-associated change in observed footprint/accessibility signal, and
whether a local difference remains among molecules with comparable distal
patterns. The current fixed-bin tables support descriptive investigation of
the second question. Conditioning on distal measurements can remove a real
variant effect on chromatin, so it should not silently replace the total effect.

For sparse TF endpoints, benchmark a two-part description—whether a small
footprint is detected, and its extent when detected—alongside the existing mean
coverage endpoint. Calibrate candidate barcode-level inference using the observed
zero masses, barcode depths and heterogeneity, with held-out WT A/A checks and
simulation of unequal variances and apparent losses. Sequence-dependent A/T
observability and caller artifacts still need separate checks. Neither a
different test nor stratification can recover erased transient footprints or
identify compartment and true labeling efficiency from these data alone.

Before introducing new biological filters, separate nucleosome-sized annotations
from long unresolved protection in subsequent processing summaries and investigate
the long-protection/near-unlabeled group. Do not select filters by the number of
significant variant hits. The previous failed calibration and withheld
significance labels remain unchanged.

## Outputs and verification

`results/rep3_heterogeneity_20261017/features/` contains the molecule features,
all-primary retention table and extraction audit. `summary/` contains barcode
and variant summaries, marginal/count/joint distributions, exact zero/one
masses, fixed-bin conditional comparisons, JSON audit and PNG/SVG figures.
See [usage](heterogeneity_usage.md) for schemas, weighting and commands.

Twelve focused tests passed, including mixed/reverse MM decoding, indel
projection, missing-versus-zero handling, ordered stream checks, cohort counts,
barcode weighting, source identity and bin boundaries. An independent bounded
real-read check validated regional methylation on 505 molecules. Full-data
validation compared TF/MSP means for all **568,354 barcodes** against the earlier
per-base arrays; maximum discrepancy was **2.22×10⁻¹⁶**. All 152 marginal, 32
joint and 16 discrete-count histogram groups conserve their expected mass.
Both figures were inspected, and the final Snakemake dry run reports no work
remaining. Direct CLI execution supplies stage audits rather than Snakemake's
own execution metadata.

Extraction took 50 minutes 31 seconds with about 115 MiB peak resident memory;
summarization took 10 minutes 49 seconds with about 754 MiB. Each used one CPU
in existing Slurm allocation 41830288. Logs, code/configuration snapshots,
independent checks and supplemental diagnostics are in
`results/rep3_heterogeneity_20261017_run/`. The supplemental annotation audit is
reproducible with `tools/inspect_heterogeneity_annotations.py`; its final result
is `annotation_diagnostics_v2.json`. No original ledger, barcode assignment,
footprint annotation or statistical release decision was changed.

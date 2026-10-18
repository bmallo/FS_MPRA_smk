# Replicate-3 heterogeneity characterization — 2026-10-17

This stage measures observable variation before selecting filters or the next
inferential model. It follows the user's four concerns: protection during a
five-minute reaction, uneven labeling, heterogeneous chromatinization, and a
possible contribution from extranuclear plasmids. It does not estimate true
labeling efficiency, instantaneous TF occupancy, residence time or compartment.
The prior failed Welch calibration remains unchanged.

## Populations and preserved evidence

Traverse all 8,715,827 primary plasmid molecules in the existing annotation BAM
and joined table. Count retention and overlapping QC failures for every row.
Extract new features for **sequence-qualified WT/SNV molecules**: exact designed
genotype, sufficient mutable-element base quality, and sequence QC passing the
existing passes/read-quality/MAPQ requirements. Do not require length QC,
methylation QC, barcode assignment or footprint availability for this broader
feature cohort. Missing measurements remain missing. Unselected molecules and
all original evidence remain in the unchanged source ledger and joined table.

Summarize three explicitly named cohorts: sequence-qualified, prior
footprint-analysis eligible, and prior eligible with supported barcode. The
last cohort must reconcile to 2,102,906 molecules and 568,354 barcodes; the
prior all-eligible cohort must reconcile to 2,594,368 molecules. Do not imply
that the broader cohort includes low-sequence-quality or off-design molecules.

## Regions and measurements

Reference length is 4,647 bp. Use zero-based half-open coordinates:

- Whole reference: [0,4647).
- Mutable element: [1234,1445), the user-confirmed 211 bp.
- Distal intervals: [0,1134), [1545,3269), [3334,4647). Exclude the element
  plus 100 bp on each side and the N15 barcode plus 25 bp on each side.

Distal regions reduce direct reuse of the target outcome, but are not assumed
invariant to genotype or free of regulatory biology. No normalization to these
regions is imposed.

For each region export aligned bases and coverage completeness, aligned A/T
opportunities, observed m6A counts/fractions at ML thresholds 128 and 204, and
union-covered fractions for 10–60-bp TF calls, ≥85-bp nucleosome calls and all
MSPs. Fraction denominators use actual aligned bases or aligned A/T opportunities;
they are distinct and explicitly stored. Missing annotation is not zero. MM
groups with unknown omitted observations cannot establish an unmethylated
denominator. Apply call-size thresholds to original query length before clipping.

Summarize distal methylation variation across fixed 100-bp reference tiles
wholly contained in a distal interval, requiring ≥20 aligned A/T opportunities
per tile. Export the number of informative tiles and rate q10/median/q90.
These describe observed variation; they are not independent measurements of
methylase efficiency. Tile qualification and incomplete coverage can themselves
select observations and are visible in the exported denominators.

Whole-plasmid nucleosome call-segment counts, median call lengths and internal
midpoint spacing describe annotation structure. Circularly split calls are not
reconstructed, and wraparound spacing is excluded. Thus segment count is not a
literal count of bound histone octamers. Keep TF counts both before and after
the 10–60-bp restriction to expose potential cutoff effects. Raw AQ presence
can be recorded without treating byte scores as calibrated probabilities.

## Summaries and weighting

Produce molecule-weighted distributions for the broader and previous cohorts.
For the primary supported-barcode cohort also give each barcode total weight
one, split equally among its eligible molecules. A histogram with these weights
retains within-barcode molecular heterogeneity; a histogram of barcode means
would obscure it. Separately export per-barcode means and within-barcode spread
for comparisons across construct backgrounds.

Output per-variant means using equal-barcode and molecule weighting, coverage,
missingness, and QC retention. WT is the reference for descriptive differences.
SNV pooling is labeled explicitly and does not represent an average variant;
plots of individual variant effects retain coverage requirements only for
readability, with all measurements retained in tables.

Joint plots examine distal methylation versus distal nucleosome coverage and
local TF/MSP coverage. Fixed display bins have edges 0,.1,.25,.5,.75,1 on both
distal axes. They are numerical bins, not inferred biological states or quality
classes. Describe each variant's bin proportions and local TF/MSP coverage
within bins. Molecule weights remain 1/(eligible molecules in its barcode),
then normalize within a bin; report unique barcodes and effective weight there.
Within-bin effects are conditional descriptive associations, not direct causal
effects. Do not call high methylation/no protection 'naked', 'extranuclear', or
'technical failure'. Do not select bins or filters by significant variant hits.

## Validation and execution

The existing join was written in primary annotation-BAM order. A streaming
lockstep traversal is permitted only after verifying the exact frozen BAM source
and matching read ID, length and orientation for every primary record. Reject
trailing, reordered or mismatched rows. Respect missing/evaluated annotation
semantics and verify source metadata before/after execution. Reconcile decoded
whole-read m6A counts to the ledger and element coverage to the previous join.

Never use pysam.modified_bases for this dataset: earlier work found a crash on
mixed modification groups. Test explicit MM/ML decoding for forward/reverse
orientation, multiple modifications, repeated groups, unknown omissions and
malformed records. Projection must preserve insertion/deletion gaps. Baseline
and stricter ML thresholds must be monotone. Include missing/empty MA tests,
count conservation, barcode weighting and state-bin boundary tests.

Use one CPU per stage in the existing Slurm allocation, explicit bounded
arrays, and 8 GiB workflow memory. Benchmark a labeled partial run before the
full traversal; incomplete runs must never masquerade as full datasets. Publish
stage directories atomically, preserving input files, configurations, hashes,
logs and code snapshots. Root/Astra owns this scientific design, independent
validation, interpretation and review. The explicitly requested Sol agents
implement extraction, descriptive summaries, their tests and the workflow.

Deliverables: molecule features, retention table, barcode and variant summaries,
weighted histograms, within-bin descriptive comparisons, static scientific
figures, and a readable report. No new statistical discoveries, automated
filter choice, live server or dashboard deployment is part of this stage.

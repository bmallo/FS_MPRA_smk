# LDLR plasmid Fiber-seq redesign: scope and proposed plan

Updated 2026-10-14. This plan now has a descriptive pilot implementation;
calibrated inference remains planned. Existing analyses are exploratory and
subject to reassessment.

## Confirmed experimental context

- Project data: `/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA`.
- Initial development dataset: `fiberseq/HepG2/replicate_3` (PS01898), identified
  by the user as the newest and best dataset.
- HepG2 runs 1–3 are independent biological preparations. Runs 1–2 used
  10-minute reactions with 1.5 µL Hia5 per 2 million cells; run 3 used 5 minutes
  and 3 µL per 2 million cells. Treat these as two reaction conditions, with
  two and one biological preparations respectively. Other protocol details,
  enzyme stock activity, transfection conditions, and harvest timing remain
  incompletely documented.
- The user confirmed the 211-bp mutable element as the primary accessibility
  interval. Do not silently expand it to the rest of the plasmid promoter.
- HEK293T uses an older library; HepG2 uses a newer library with a smaller LDLR
  element. Cell type and library version must remain separate metadata fields.
  Their effects cannot be isolated by a simple comparison between these groups.
- Functional data are existing published MPRA measurements, not newly collected
  RNA/DNA data matched to the Fiber-seq biological replicates.
- Reads called "naked" in the October scripts are high-methylation plasmids
  from the same Fiber-seq reaction. They are not separately prepared controls.
  Use "high-methylation population" as the observational label; unchromatinized
  plasmids are its biological interpretation, not an independently known label.
- All October filtering, normalization, and analysis choices are exploratory.
- The user's top scientific priority is how individual base changes alter
  chromatin architecture at both the small, TF-sized footprint scale and the
  overall promoter-accessibility scale. Published MPRA activity provides
  functional context; motif identity and dependencies between sites are later
  extensions.

## Existing assets and checks

The old repository contains filtering, barcode-consensus variant assignment,
empirical occupancy testing, motif/co-occupancy analyses, calibration tools,
and a Dash browser. Its documentation and interfaces have diverged from the
implementation.

The newer replicate-3 directory contains sequence/design-based genotyping,
per-read ledgers, filter and purity reports, FiberHMM outputs, likelihood-based
footprint analyses, published-MPRA comparisons, and downstream scripts reused
for the other HepG2 replicates. These are candidates for reuse, not defaults
to adopt without review.

Small reference/table checks performed during the initial review:

- The rotated HepG2 plasmid reference has 4,647 bases.
- The mutable interval is 211 bases: oligo positions 21–231, corresponding to
  zero-based half-open rotated-reference interval [1234, 1445).
- Both `mpra_published/GRCh38_LDLR.tsv` and `GRCh38_LDLR.2.csv` contain all 633
  SNVs in that interval, with reference-base agreement for every SNV.
- The mapping verified for these files is genomic position = 11,089,264 +
  oligo position, on chr19 in GRCh38.
- All 954 score/p-value pairs in `extras/data/LDLR_mpra_scores.tsv` match
  `GRCh38_LDLR.2.csv` after complementing alleles and mapping genomic position
  = 11,092,732 minus the old table's position. Do not generalize this offset
  to other constructs without sequence validation.

The existing replicate-3 final summary reports 3,280,497 retained reads,
779,516 WT reads, and median coverage of 1,229 reads across the 633 SNVs.
Those are saved report values, not a fresh BAM recount. That final BAM is
already selected by methylation and other criteria; upstream outputs are
needed to reassess filtering and measure the high-methylation population.

## Scientific questions and proposed comparisons

The organizing objective is to connect sequence variation, single-molecule
chromatin organization, and regulatory activity. Keep the following questions
distinct in both analysis outputs and the dashboard.

### Primary readouts: two complementary spatial scales

**Small footprints:** for each base substitution, measure changes in occupancy,
location, and length distribution of small protected intervals across the
promoter, including effects away from the mutated base. Describe these as
TF-sized footprints; size alone does not identify the bound protein. Preserve
continuous footprint lengths alongside configurable size ranges, and evaluate
caller agreement and sequence-dependent observability. Report signed effect
sizes, uncertainty, molecule support, and biological-replicate consistency.

**Promoter accessibility:** measure both the spatial accessibility profile and
the distribution of promoter-wide accessibility across individual molecules.
Evaluate a continuous accessibility summary accounting for observable sites;
if an open/closed classification is useful, define and validate its rule
explicitly and retain the continuous measurements. Distinguish shifts in the
fraction of broadly accessible promoters from changes within more protected
promoters. Nucleosome occupancy and organization help interpret these effects.
The promoter interval and any subregions must be reference-defined and stored
as metadata rather than selected separately to maximize each variant's effect.

For both scales, compare each designed single-SNV genotype with WT within each
biological replicate. Preserve the high-methylation population in the source
data and quantify how selection changes the result. Separate missing or
uninformative sequence positions from measured absence of protection. An A/T
gain or loss changes the observation opportunities, so naive methylation-count
comparisons are not sufficient. Any normalization or masking strategy must be
evaluated with sequence-aware checks and sensitivity analyses.

Small-footprint changes and promoter-accessibility changes are parallel primary
outcomes; neither is required to be significant for the other to be reported.
Do not require MPRA significance to test or display a chromatin effect.

1. **Overall variant effect:** within an assay- and genotype-quality-qualified
   population, does the variant change methylation, accessibility, footprint
   occupancy, nucleosome organization, or the distribution of molecule states?
2. **Within-state effect:** among molecules with comparable global chromatin
   properties, does the variant change local or distal organization?
3. **Functional association:** do those effects agree with published activity
   measurements for the same genomic allele and regulatory sequence context?
4. **Mechanistic hypotheses:** do multiple perturbations and biological
   replicates support a particular binding element or relationship between
   sites? Distinguish such evidence from identification of a bound protein or
   proof of a causal pathway to expression.

### Inferred high-methylation population

Proposed policy: retain this population in the data model and study its
abundance as a possible biological outcome. Do not silently discard it or use
it as a known protein-free reference.

An effect on entry into this population and an effect within the remaining
population answer different questions. Conditioning on methylation or
nucleosome count could remove part of a variant-associated change; these
adjustments therefore require explicitly labeled comparisons.

The high-methylation reads can inform exploratory background and caller
sensitivity analyses. Their interpretation as controls, and the October
difference-in-differences subtraction, remain hypotheses to validate. Preserve
raw contrasts alongside any adjusted contrast. Do not tune adjustments merely
to increase agreement with MPRA scores.

### Replicates and uncertainty

Estimate effects separately in each biological replicate before combining
evidence. Molecules provide within-replicate precision; they are not additional
biological replicates. Keep replicate effects, uncertainty, coverage, and
heterogeneity visible. Sparse variants can remain descriptive in individual
replicates rather than being forced into a pooled significance claim.

## Proposed implementation sequence

### 1. Experimental manifest and coordinate contract

Record library version, construct/reference sequence and checksum, biological
replicate, cell type, preparation/batch, assay and caller versions, and source
paths. Define stable sequence/haplotype IDs and explicit genomic, oligo, and
rotated-plasmid coordinate mappings, including strand and indel equivalence.

Acceptance: all designed alleles are accounted for; reference and allele joins
are verified; unmatched or ambiguous mappings are reported explicitly.

### 2. Auditable molecule-level ingestion and QC

Adapt existing upstream outputs and genotyping code behind a common interface.
Preserve read IDs, sequence/genotype confidence, raw and assigned variants,
barcodes, methylation observations, footprint calls, and exclusion flags.
Separate technical/genotype QC from biological state selection. Keep
homopolymer indels and multi-mutant molecules identifiable rather than silently
assigning their effects to a single SNV.

Acceptance: every analysis cohort can be reconstructed from recorded criteria;
counts reconcile with source data; missing observations are distinguishable
from observed absence of a footprint. Existing final BAMs are useful comparison
artifacts, not the sole source for this stage.

### 3. Replicate-3 measurement and analysis pilot

Begin with WT and designed single-SNV molecules; support indels and complete
haplotypes in the schema from the outset. Summarize overall methylation/state
distributions, accessible regions, nucleosome organization, and smaller
footprints. Compare candidate callers and inspect molecules at representative
sites. Evaluate effects of A/T changes, sequence context, genotype uncertainty,
and methylation selection before choosing primary statistical comparisons.

Acceptance: effects and uncertainty have clearly defined populations and
baselines; apparent losses of observability are not labeled binding changes;
technical and scientific assumptions are explicit.

### 4. Calibration and biological replication

Use small synthetic fixtures for coordinates, reverse strands, indels, missing
observations, and format conversion. Assess null behavior with appropriate
within-replicate splits/permutations and assess power with injected effects.
When data drive site discovery or normalization, include those steps in the
validation or evaluate on held-out data. Define multiplicity across the actual
reported variants, sites, readouts, and comparisons.

Bring in HepG2 replicates 1 and 2 after documenting protocol differences.
Evaluate reproducibility without selecting procedures on the same data later
presented as independent validation. Keep sensitivity analyses distinct from
the designated primary results.

Acceptance: representative end-to-end runs are reproducible, null and power
limitations are reported, and biological replicate agreement is assessable.

### 5. Functional integration

Import both published MPRA tables with their original values, provenance, and
identities preserved. Resolve what LDLR and LDLR.2 represent before choosing a
combined estimate. Keep assay context differences explicit. Summarized DNA/RNA
columns do not substitute for barcode-by-biological-replicate count data.

Acceptance: allele-specific joins retain every alternate allele, report
coverage/missingness, and make no claim of matched experimental measurement.

### 6. Dashboard and workflow packaging

First views: experiment/QC overview; allele-level activity versus chromatin;
and variant detail with sequence, molecule tracks, observed and adjusted
effects, uncertainty, and replicate comparison. Show molecule-state abundance
and genotype confidence. Add motif and co-occupancy exploration after the
underlying measurements and inference are validated.

Prioritize a sequence-position-by-alternate-base overview linked to two views:
(1) small-footprint occupancy profiles and a variant-by-footprint-position
effect matrix, and (2) promoter-accessibility profiles and per-molecule
distributions. Selecting an allele should synchronize these views with WT,
replicate estimates, and sampled molecule tracks. Published activity is an
optional linked annotation/comparison. Show uncertainty and coverage alongside
effects; thresholding the display must not redefine the statistical test family.

Use a versioned results contract with separate tables for samples, constructs,
variants, effects, and annotations, plus indexed molecule/track data. The
dashboard reads precomputed results and lazily loads detail. Retain Snakemake,
declare resources and dependencies explicitly, and use allocation-bounded
thread counts and Slurm execution for substantive computation. Hosting and
frontend selection follow the access and usage requirements.

Acceptance: a scientist can trace a displayed effect to its allele, sample,
cohort, method, and molecule evidence; export a view and its data; and
distinguish unmeasured, uncertain, and statistically supported effects.

## Remaining experimental and product information

- Other preparation/transfection conditions, harvest timing, controls, and
  enzyme stock activity; interpretation of methylation extremes.
- Library construction and barcode use, and the intended role of designed
  indels and non-design multi-mutants in the scientific questions.
- Any pre-specified subregions within the confirmed 211-bp primary interval,
  and the first decisions the dashboard should support.
- Intended users, expected data scale, and local versus shared dashboard access.

## First milestone

A reproducible HepG2 replicate-3 pilot with validated sequence mappings,
auditable molecule selection, per-SNV small-footprint effects and promoter-wide
accessibility measurements, and a linked variant-detail view. Include published
activity as context and extend to biological-replicate validation before making
generalized claims. Establish credible measurement and comparison choices
before scaling discovery analyses.

### Implementation status — 2026-10-14

The first descriptive implementation is working: a validated sequence catalog,
deterministic whole-ledger molecule sampling, WT/SNV measurements from upstream
FiberHMM annotations, an independent Snakemake workflow, and a portable local
dashboard. See [usage](pilot_usage.md) and the [pilot run report](pilot_run_20261014.md).
The corrected replicate-3 run measures 28,415 molecules from a 1% ledger sample.

This is an initial measurement and visualization milestone. Full-scale storage,
genotype reassessment, sequence-observability calibration, uncertainty and
multiple-testing procedures, biological replicate
comparison, and a shared production dashboard remain open work. The current
viewer keeps replicate datasets separate and makes no discovery calls.

The next implementation adds explicit reaction metadata and a three-run
comparison dashboard. The primary interval is now confirmed. The 1% samples
from runs 1–2 have insufficient shared SNV coverage for effect concordance;
larger imports and within-condition calibration are prerequisites. See
[condition comparison report](condition_comparison_20261014.md).

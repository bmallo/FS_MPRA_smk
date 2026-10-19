# Replicate-3 descriptive landscape

This stage follows the user's decision to prioritize biological description and
candidate generation over further significance calibration. It uses the full
existing supported WT/SNV cohort, with no new biological filters or caller
thresholds. The primary interval remains the user-confirmed 211-bp element,
reference [1234,1445), oligo positions 21–231. Existing whole-plasmid and distal
summaries supply context; the joined call intervals cover the element only.

Every eligible molecule contributes 1 / eligible depth of its barcode. Dividing
by the number of barcodes gives each barcode equal total weight. Molecule-weighted
profiles are also exported. WT and 633 SNVs remain separate, including variants
with limited coverage. No p-values, significance labels, or validated effect
intervals are produced. Barcodes are not biological replicates.

Categorize calls using their original query length, before element clipping:
TF 10–60 bp; other TF calls; nucleosome-class calls <85, 85–250, 251–500,
501–1,000, and >1,000 bp; and MSPs. These are annotation/length categories,
not protein identities. Union overlapping intervals per molecule and category.
Export per-base coverage, mean element coverage, and the fraction of molecules
with any overlapping call. Categories can overlap, especially TF and MSP, so
their fractions need not sum to one. Record WT original call lengths and
TF boundary distributions separately from molecule coverage.

Additionally classify molecules into eight mutually exclusive combinations of
small-TF, 85–250-bp nucleosome-class, and >250-bp nucleosome-class overlap with
the element. The all-absent state means none of these three calls, not naked,
unbound, or extranuclear DNA. Other-TF and <85-bp nucleosome-class calls remain
visible in their own tracks. State fractions describe calls anywhere within the
element; co-occurrence does not establish interaction or simultaneous in vivo
binding.

For descriptive local comparisons, use the pre-existing fixed 20-bp windows
(step 10, plus terminal [191,211)). Also outline WT small-footprint peaks using
a five-base moving average, minimum prominence 0.005, and minimum separation
10 bp; display 20-bp neighborhoods around peak centers. Peak neighborhoods may
overlap and are not inferred physical binding-site boundaries. Parameters are
specified before inspecting variant profiles; sites are selected from WT only.
Keep the full unsmoothed profile and fixed-window effects available to avoid
hiding signals outside these neighborhoods. These selections are exploratory
and must not be reused as independent confirmation in later audits.

The atlas should show WT coverage and length categories, molecular state
fractions, all-variant effects across position, selected variant-versus-WT curves,
coverage, A/T gain/loss status, and total/distal heterogeneity context. Rank raw
effect magnitudes only for navigation, with support displayed; such ranks are
subject to sampling noise and winner's bias. Adjacent-base patterns are candidate
regions, not motif identifications. Published functional measurements, if
included, retain their original assay context and sequence-verified mapping.

Mandatory verification: all supported barcode depths and variant identities
reconcile, eligible read IDs are unique, all eligible molecules span the element,
empty calls mean evaluated zero while missing calls fail, category coverage is
bounded by call presence, molecular state fractions sum to one, and full-cohort
TF/MSP profiles independently reconcile with the existing barcode count arrays.
Record immutable source identities, atomic outputs, and one-CPU execution.
Further sequence observability, caller, and biological confidence audits are
explicitly deferred, as requested.

# HepG2 replicate-3 barcode benchmark — 2026-10-16

Keep exact barcode grouping as the production baseline. This benchmark compares
it with quality-based correction and the pinned Pacybara clustering core. It
tests correction yield and risks to the number of independently constructed
barcodes; it does not validate footprint statistics or infer biological effects.

## Random-sample QC

The indexed scan covered 8,715,827 primary plasmid molecules. A deterministic 1%
read sample retained 87,618 molecules. Of these, 25,972 met the common WT/SNV
sequence/QC and anchored barcode criteria used for the comparison. No methylation
or footprint criterion was used to select this comparison cohort.

There were 25,521 observed barcodes: 25,085 singletons, 421 doubletons and 15
triplets. This sample is useful for QC but is too sparse to optimize correction.

| Method | Assigned molecules | Exported barcodes | Changed barcode sequences on molecules |
|---|---:|---:|---:|
| Exact | 872 | 429 | 0 |
| Quality directional | 872 | 429 | 0 |
| Pacybara core, edit distance 1 | 874 | 430 | 0 |

The two additional Pacybara assignments came from one barcode with a Q40 read
and a read containing a Q10 barcode base. They reflect different barcode-quality
acceptance, **not demonstrated error correction**. The held-out test had only six
matches per method, with zero genotype disagreements; that denominator cannot
establish reliable accuracy.

## Full-depth barcode neighborhoods

A second complete contig traversal collected all anchored edit-one neighbors of
400 discovery barcodes: 200 repeated barcodes and 200 singletons, selected by a
fixed hash within those strata. It recovered 3,738 molecules; 2,620 passed the
comparison criteria, spanning 449 observed barcodes. This is a targeted follow-up,
not an unbiased library sample. Discovery molecules were excluded from validation:
608 were forced into training, leaving 397 held-out molecules from the remaining
80/20 hash split.

| Method | Assigned molecules | Exported barcodes | Molecules with changed barcodes | Held-out export matches | Genotype disagreements |
|---|---:|---:|---:|---:|---:|
| Exact | 2,533 | 392 | 0 | 372 / 397 | 1 |
| Quality directional | 2,533 | 392 | 0 | 372 / 397 | 1 |
| Pacybara core, edit distance 1 | 2,556 | 390 | 18 | 370 / 397 | 1 |

Pacybara added 23 assignments (0.88% of the comparison cohort). Of its 18 barcode
changes, eight affected molecules already assigned by exact matching. Three
source barcodes with 2, 4, and 2 Q30-supported reads were merged into other
barcodes with the same genotype. The strongest example is
`GGATGAATGAGTGGA` (four Q30 WT reads) → `GAATGAATGAGTGGA` (eight Q30 WT reads).
These are questionable merges, not independently proven errors. The remaining
ten changed barcodes were previously unassigned, including four barcode indels.

The same held-out disagreement occurred in all three methods:
`m84046_261003_041040_s2/230167521/ccs`, barcode `CCCACAATTAATGAG`, observed WT
versus predicted `195C>A` (catalog oligo-position label). Across the complete
qualified cohort, this barcode has 13 `195C>A` observations and one WT observation.
Exact grouping quarantines it as conflicting. Pacybara moves the WT observation
to neighboring WT barcode `CCCACAATTAACGAG` and thereby also accepts the 13
substitution-bearing reads. This is a plausible correction, but barcode error,
genotype error, and mixed construct origins cannot be distinguished conclusively
with this evidence alone. Held-out matching against the frozen barcode export
appropriately exposes the conflict instead of using the test genotype to repair it.

The final run is `family_comparison_final/report.json`, with per-read held-out
predictions in `heldout_predictions.tsv.gz` and all proposed transitions in
`real/corrections.tsv`. Per-variant barcode counts and effective counts are in
`real/variant_coverage.tsv`.

## Known-truth stress tests

The synthetic set contains 2,980 molecules across eight deliberately difficult
scenarios, with 20 independently generated instances of each scenario. These
frequencies do not represent measured error frequencies in replicate 3.

| Method | Assigned molecules | Incorrect assignments | Merged genuine neighboring barcode families |
|---|---:|---:|---:|
| Exact | 2,520 | 0 | 0 |
| Quality directional | 2,540 | 0 | 0 |
| Pacybara core, edit distance 1 | 2,620 | 72 | 20 |

Quality directional corrected all 20 low-quality substitution-error molecules
and abstained on high-quality errors, indels, and ambiguous parents. Pacybara
corrected all 80 introduced substitution/indel error molecules. It also merged
all 20 genuine neighboring families with the same genotype (60 affected reads),
and chose the wrong parent for 12 of 20 ambiguous children. Separate neighboring
genotypes remained distinct; ambiguous barcode-to-genotype exports were rejected
by the shared export policy.

This demonstrates a risk relevant to this library: genotype agreement and unequal
read abundance do not establish that neighboring barcodes represent the same
constructed plasmid. Neither method has been proven free of false merges on real
data, whose barcode truth is not independently known.

## Decision

Keep `barcodes.method: exact` as the production default, retain conflicting and
singleton evidence, and keep correction proposals separate for inspection.
Quality directional supplied no extra assignments in the tested real cohorts;
Pacybara's modest gain came with mergers of supported barcodes and no improvement
in the held-out linkage check. This does not establish that Pacybara is unsuitable
in general: only its clustering core, these settings, and these cohorts were tested.

The next processing milestone should generate the complete replicate-3 barcode
and genotype ledger with these audit fields, then assess barcode diversity and
concentration per variant before choosing barcode-aware occupancy statistics.
Any broader correction or genotype rescue should be validated using independent
construct evidence or calibration data, including error simulations matched to
the observed qualities and family depths. Barcode absence in Fiber-seq remains
distinct from absence in future RNA/DNA data.

## Reproducibility and scope

See [commands and method details](barcode_benchmark.md). Results are in
`results/barcode_benchmark_20261016/`; the random-sample comparison is
`comparison/report.json`. Pacybara commit is
`b5e39bd0e2731283eee2ee88c0e0852ddf83fe33`. Its clustering core is unmodified;
the shared extraction, exact distance calculation, and export rejection policy
are adapters specific to this benchmark. The full protein-oriented Pacybara
workflow was not run.

The barcode extractor now traverses CIGAR blocks to locate flanking anchors.
It produced identical barcode sequences, qualities, and extraction statuses to
the prior per-base projection on 10,000 real molecules. In that measured check,
the extraction calls took 4.898 versus 0.091 seconds (about 54-fold faster).
This timing excludes BAM decoding and the rest of the pipeline. Randomized CIGAR
tests cover clips, substitutions, insertions, deletions, and reference skips.
All 36 unit/integration tests pass, including the unmodified upstream clustering
core, complete edit-distance candidate checks, and discovery/validation separation.

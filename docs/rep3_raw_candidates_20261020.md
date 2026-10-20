# Replicate-3 raw-call and filtering audit — 2026-10-20

The position-125 candidate survives comparison at unchanged methylatable bases and grouping by directly observed genotype without barcode support. The adjacent 143–150 region retains a broader accessibility/protection pattern, with weaker consistency for individual variants. These are descriptive findings from one preparation, not newly significant results or identified protein-binding events.

Open `results/rep3_raw_candidates_20261020/viewer/index.html` for an offline interactive audit with 27 variants, four regions and 18 comparison settings. The earlier functional/footprint atlas remains unchanged at `results/rep3_landscape_20261018/atlas_biology_20261020/index.html`.

## Position 125: a localized raw-call change supports the footprint observation

Oligo125 is GRCh38 **chr19:11089389**, two bases upstream of the mapped SRE1 core. Within WT_P3 (chr19:11089378–11089397), comparisons exclude the mutated position **and two flanking bases on either side**, leaving nine unchanged A/T positions: five reference A and four reference T. Reference A measures adenine on the forward genomic strand, and reference T measures adenine on the reverse strand.

| Substitution | Shared-site retained m6A change | Small-footprint coverage change | Whole-element shared-site m6A change | Primary molecules / barcodes |
|---|---:|---:|---:|---:|
| 125A>C | +4.12 pp | −3.68 pp | +0.93 pp | 1,360 / 380 |
| 125A>G | +3.65 pp | −3.82 pp | +0.34 pp | 1,786 / 480 |
| 125A>T | +7.00 pp | −3.66 pp | +2.33 pp | 1,602 / 438 |

All three have positive raw-call changes on both genomic adenine strands. A>T retains an A/T opportunity while switching its strand; its effect at unchanged sites is especially useful evidence against simple deletion of the informative base as the entire explanation. One A>T molecule fails the complete-window raw-call requirement, leaving 1,601 for its combined raw endpoint.

These patterns persist with equal molecule weights and with direct-genotype grouping independent of barcode support. The barcode-independent circular cohort gives raw-call increases of +4.31, +4.18 and +5.92 pp and footprint changes of −3.64, −3.64 and −3.34 pp (C, G, T). Both alignment orientations retain the direction. In all 50 selections of one molecule per barcode, the raw changes remain positive and footprint changes negative for each substitution.

The strict ≥40-pass subset attenuates footprint changes to −1.42, −1.44 and −1.91 pp, with only 158, 234 and 205 variant molecules. This is a real sensitivity to report; it does not establish whether read quality, molecule selection or sampling caused the attenuation. The same filters applied to WT also change the reference landscape.

The local raw effects exceed the corresponding whole-element effects. That favors a localized alteration over a purely uniform change in methylation. However, A>T also has a broader accessibility increase and less long protection. The local signal is therefore not automatically an isolated TF-binding event. Shared-site comparisons remove the direct opportunity-count artifact, but do not eliminate sequence-dependent methylation kinetics, basecalling context effects, or footprint-boundary changes extending beyond ±2 bp. They also do not assign a protein or an activator/repressor role.

## Positions 143–150: broad architecture is more informative than a tiny local footprint

Across 24 substitutions, the primary barcode-weighted whole-element comparison has:

- Lower retained methylation at shared A/T sites for **21/24**, median −1.27 pp.
- Lower MSP accessibility for **20/24**, median −1.79 pp.
- Higher 85–250-bp protection coverage for **21/24**, median +2.09 pp.

With direct-genotype grouping and no barcode-support requirement, the corresponding counts are **23/24**, **22/24** and **22/24**, with medians −1.21, −1.67 and +2.10 pp. This supports a regional trend toward more protection/lower accessibility that is not solely a barcode-assignment artifact. The prior functional data provide complementary context, not independent validation of these molecular measurements.

The primary counts differ slightly from the previous full-WT atlas (for example, 20 rather than 22 variants with lower MSP) because this raw audit uses a **fixed 10,000-barcode WT sample**. Differences near zero can change sign with the WT reference sample; the group trend is more stable than those individual signs. The original full-WT results have not been overwritten.

The eight-base local interval has only two WT A/T positions. Excluding the mutation ±2 bp leaves **no local raw measurement for all substitutions at positions 145 and 146**. Missing values remain missing. Use the flanking interval and whole-element measurements instead of interpreting unavailable local measurements as no effect.

Noncircular molecules do not show a consistent adjacent-region direction, but they are a small, different population (18–93 molecules per variant; WT 927). Adding them to the full-element genotype cohort barely changes the broad trend. This is insufficient evidence to claim replication across plasmid forms.

## What was measured and normalized

The extraction scanned all **8,715,827 primary alignments** and recovered **134,529 selected molecules**, verifying full mutable-element sequence identity, coordinate projection and agreement between two retained-call decoders. It retained all qualified candidate molecules, plus the union of two distinct WT samples. Those WT samples are kept separate in analysis:

- Supported-barcode comparisons: the fixed 10,000-barcode WT sample and all its 37,075 eligible molecules.
- Barcode-independent comparisons: a deterministic 1/32 sample of sequence-qualified WT molecules, with 39,022 circular/eligible molecules or 39,949 when circularity is not required.

**Every plotted outcome is normalized for its usable support.** Equal-barcode analyses first average molecules within each barcode and then average barcodes equally. Equal-molecule analyses average molecules directly. Raw outcomes are the fraction of unchanged A/T positions with a retained call. Footprint outcomes are the fraction of region bases covered by the specified original calls. Neither analysis compares unnormalized read counts; equal weighting is not a claim of independent biological replication.

Per-site profiles use the available denominator at each site. Shared-window endpoints require all selected sites to be known and pass the chosen base-quality threshold in a molecule. A barcode with no available outcome is missing for that endpoint. The comparison export includes callable molecule/barcode counts, and the retention export gives selected counts before endpoint-specific missingness.

All extracted positive m6A ML scores are **240–255**. “Raw” means stored MM/ML calls before interpretation as footprints, not original kinetic signals or the probability distribution before upstream filtering. Unlisted known sites mean no retained call, not a measured probability of zero. Missing channels and unknown skipped sites are never converted to zero. Decoding follows the [SAM modification-tag specification](https://github.com/samtools/hts-specs/blob/master/SAMtags.tex). Threshold 250 is a sensitivity analysis of retained calls only. No footprints were re-called, including under the base-quality or ML-threshold scenarios.

## Sensitivity comparisons and their limits

The 18 settings cover equal barcode/molecule weights; direct genotype with and without circularity; noncircular molecules alone; read quality ≥0.9999; ≥20 or ≥40 passes; mapQ ≥60; read length within 1% of 4,647 bp; both alignment orientations separately; barcodes with ≥3 or ≥5 primary molecules; site baseQ ≥30; retained ML ≥250; whole-plasmid m6A <0.90; and distal m6A 0.10–0.80.

The existing sequence qualification remains in force even when barcode support is ignored: exact designed genotype, passes ≥10, read quality ≥0.99, mapQ ≥20 and genotype baseQ ≥20. Thus “barcode-independent” does not mean unfiltered or newly re-genotyped. The barcode-depth requirement uses original eligible cohort depth; other filters recompute weights over remaining molecules.

MapQ ≥60 and whole-plasmid m6A <0.90 retain the entire primary cohort here and supply **no additional validation**. BaseQ/ML changes affect raw outcomes only, so unchanged footprint results in these settings are not independent confirmations. Whole-plasmid and distal methylation restrictions can select on the biological outcome itself; they are diagnostics, not recommended default filters. Alignment orientation and genomic adenine strand are distinct concepts.

The export also contains 50 one-molecule-per-barcode repeats and 100 same-size WT barcode subsets for four candidate endpoints. The WT subset ranges describe a sampling scale; they are **not calibrated p-values or confidence intervals**. Joint molecule-state tables cross four retained-call fraction bins with presence of nucleosome-sized and/or long protection in each region. Their fractions use equal barcode weights and sum to one within each available comparison; they describe observed mixtures without identifying nuclear location or protein-free molecules.

## Next experimental checks

The proposed ten-construct panel in `validation_panel.tsv` includes WT; all three position-125 substitutions; 143C>G, 146C>G and 150C>G as GC-preserving adjacent-region changes; 144T>A and 147T>A as adenine-strand switches; and 147T>G as a functional-discordance example. These are hypotheses and design contrasts, not a ranked list of significant variants. The exact 211-bp mutable sequences and stable variant IDs are exported; no barcodes or constructs have been ordered.

Prioritize independent biological preparations of this small panel with balanced multiple barcodes. In a dedicated, separately prepared protein-free reaction, compare the same sequences at the same unchanged sites to test sequence-dependent labeling/calling. A modest reaction-time/enzyme series can test whether local footprint contrast improves before saturation; it must preserve comparable preparation conditions and have explicit controls rather than treating high-methylation molecules as a naked-DNA control. Recover original modification probabilities/kinetics if available in new sequencing. Pair new functional measurements with these same constructs and use an orthogonal binding assay or candidate-factor perturbation before assigning protein identity or regulatory sign.

The practical next step is this focused experimental panel, while retaining the present observations as prioritized biological hypotheses. Additional filter tuning on replicate 3 should not be used to manufacture significance.

## Reproduce and inspect

Run on an allocated compute node with BLAS/OpenMP threads set to one. Existing result directories are immutable; choose fresh output directories when rerunning.

```bash
python -m fsmpra.raw_candidate_extract --output results/rep3_raw_candidates_20261020/extracted
python -m fsmpra.raw_candidate_summary --output results/rep3_raw_candidates_20261020/summary
python -m tools.render_raw_candidates --output results/rep3_raw_candidates_20261020/viewer
python -m unittest tests.test_raw_candidate
python -m tools.validate_raw_candidates --viewer results/rep3_raw_candidates_20261020/viewer --output results/rep3_raw_candidates_20261020_run/validation.json
node tools/validate_raw_candidates.js results/rep3_raw_candidates_20261020/viewer/index.html
```

Audits record source fingerprints and definitions. The dashboard embeds its data and needs no server. Validation checks scientific invariants and the DOM/canvas interface with actual data; it does not replace a full browser visual check.

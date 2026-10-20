# Independent review and recommended direction — 20 October 2026

This review was requested after the [project handoff](PROJECT_HANDOFF_20261020.md). It asks whether to re-run and validate the existing tests, start over, or build on what exists. Everything below was computed from the committed [evidence snapshot](evidence/20261020/README.md); no BAMs or `results/` directories were available to the reviewer. New tables are in [`evidence/independent_review_20261020/`](evidence/independent_review_20261020/), reproducible with `python -m tools.independent_review_20261020 --output <new dir>`.

## Verdict

**Build on the data foundation, but change the analysis question and the next experiment.** Do not re-run the existing tests and do not start over.

1. **The numbers are sound; re-running them would reproduce the same answers.** Evidence checksums verify (34/34). The full test suite passes in a fresh Python 3.13 / current SciPy environment (127 tests including three new ones, one optional skip). Every reported library correlation recomputed exactly from the evidence tables. The difficulties are in framing and design, not arithmetic.
2. **The existing data already answer the central question at the library level.** Changes in promoter accessibility and nucleosome-sized protection track published activity across all 630 measured SNVs. The association survives a spatial null, every A/T-change class, independent barcode halves, and removal of the strongest region. The handoff left this descriptive. It can be stated with valid inference now.
3. **The per-variant × 20-bp-window significance program should stop in its current form.** The design has 13,926 hypotheses, sparse caller-derived endpoints, BY correction, and a release gate that requires every stratum to pass. It cannot release anything. It also answers a narrower question than the one the project is asking.
4. **The handoff's lead candidate (oligo 125, chr19:11089389) is the one most exposed to the main confounder.** The handoff's proposed ten-construct panel is centred on it. Better-supported leads are the two Sp1 repeats (accessibility plus activity) and the activity-critical core of WT_P1 (local footprint plus activity).
5. **The next experiment should be a library-wide, sequence-matched naked-DNA control plus a same-condition biological replicate**, not a construct panel first.

## What the new checks show

### Which per-variant effects are reproducible?

Barcodes were split by SHA256 parity into two disjoint halves for both WT and every variant. The table compares effects across the 630 SNVs with ≥50 barcodes ([table](evidence/independent_review_20261020/split_half_reliability.tsv)).

| Endpoint | Half-A vs half-B r | Reliability of full estimate | Signal SD / noise SD (pp) |
|---|---:|---:|---:|
| WT_P1 small footprint, mutation inside window | 0.76 | 0.86 | 1.61 / 0.65 |
| WT_P2, inside | 0.63 | 0.78 | 1.05 / 0.58 |
| WT_P3, inside | 0.87 | 0.93 | 1.63 / 0.44 |
| WT_P4, inside | 0.61 | 0.76 | 0.56 / 0.32 |
| Any peak, mutation **outside** window | 0.01–0.24 | 0.01–0.38 | — |
| Whole-element MSP | 0.19 | 0.32 | 0.86 / 1.26 |
| Whole-element small footprint | 0.22 | 0.36 | 0.13 / 0.18 |

Local footprint changes caused by mutations inside a window are highly reproducible within this preparation. Distal footprint effects and per-variant whole-element MSP are mostly noise. Per-variant MSP is therefore underpowered by design. Power has to come from aggregation (positions, elements, the library) and from removing nuisance variance. Within-preparation reproducibility includes any **deterministic** sequence-dependent measurement bias. It cannot by itself separate biology from observability.

### How much reproducible local signal follows A/T-opportunity change?

The question here is how much of the reproducible in-window variance is explained by the A/T-change class (gain/loss/retained), using cross-half class means ([table](evidence/independent_review_20261020/at_class_local_effects.tsv)):

| Window | A/T-class share | Pattern |
|---|---:|---|
| WT_P1 (oligo 52–71) | 2% | Gains, losses and retained-count changes all reduce the footprint (78–86% negative) |
| WT_P2 (74–93) | 1% | Same pattern as P1, weaker |
| WT_P3 (114–133) | **44%** | Gains *increase* coverage (median +0.68 pp), losses reduce it (−2.01 pp); unchanged after ±2-bp masking |
| WT_P4 (149–168) | 32% | Losses all negative; gains near zero |

The A/T-class pattern in P3 is the expected signature of observability, not protein protection. Adding an A/T opportunity should not create a protein footprint. In P3, making flanking C/G bases A/T raises the call (116G>T +3.1, 132C>T +2.2 pp), while removing A/T in the core lowers it. Masking ±2 bp does not remove this, so the dependence extends beyond the mutation, as expected from kinetic m6A-calling context and HMM boundaries.

The strongest P3 losses sit over the A/T-rich stretch at oligo 117–127 (`ACATTTGAAAA`, a poly-A run at 124–127). Activity there is largely unaffected: median +0.04 log2, 18/33 positive, and only 126A>G and 127A>T below −1. The activity-critical positions 128–133, at the start of the SRE-1 core, have strong losses (median −1.8 log2; 15/18 below −1), but their footprint changes follow A/T class in sign. Position 125 may still carry real signal: A>T, which keeps the A/T count, also loses −3.6 pp. But it is the candidate whose interpretation depends most on a naked-DNA control.

In P1, substitutions that add an A/T *reduce* the footprint (e.g., 57G>A −2.7, 56G>T −1.9, 58G>A −2.2 pp), which is the opposite of the artifact direction. In the core at oligo 54–62 (chr19:11089318–11089326, `ACGGGTTAA`), all 24 substitutions outside oligo 55 reduce published activity (element mean −1.57 log2), while MSP is unchanged. Oligo 55 is functionally neutral, and 55C>A raises the footprint by +6.0 pp, the artifact direction. That is the expected signature of a factor binding inside an already-accessible region. P1 was not among the handoff's candidates. Its identity is unknown; motif matches at threshold do not settle it.

### Chromatin change tracks activity across the library

Pearson correlations with published activity (source 1; source 2 agrees within ±0.03) are below ([table](evidence/independent_review_20261020/activity_association.tsv)). The *p* values come from an exact circular position-shift test. Position means keep a site's three substitutions together, and rotating one profile against the other preserves both profiles' spatial autocorrelation. With 211 positions only 192 rotations exist, so *p* = 0.005 is the floor. The last column shows how far the observed position-level *r* sits from the largest rotated value.

| Chromatin change | Variant r | Position r | Rotation-null max \|r\| | p |
|---|---:|---:|---:|---:|
| Whole-element MSP | +0.29 | +0.38 | 0.20 | 0.005 |
| MSP standardized for distal m6A composition | **+0.43** | +0.51 | 0.28 | 0.005 |
| MSP standardized for distal protection composition | +0.45 | +0.52 | 0.29 | 0.005 |
| Nucleosome-sized (85–250 bp) protection | **−0.41** | −0.51 | 0.28 | 0.005 |
| Element-wide small footprint | +0.17 | +0.22 | 0.24 | 0.02 |

Robustness checks:

- **Within each A/T-change class**, standardized MSP gives r = 0.41 (gain), 0.43 (loss) and 0.45 (retained), so opportunity change does not drive the association.
- **Excluding the Repeat-3 neighbourhood (oligo 133–160)**: 0.32. **Excluding |activity| > 1 log2**: 0.30. **Each barcode half alone (raw MSP)**: 0.22 and 0.23.
- **Correcting raw MSP for its measured reliability (0.32)** gives r ≈ 0.52.
- **Barcode count** is unrelated to activity (r = −0.01) and to MSP change (r = −0.04).
- **Element-wide small-footprint change** is weak and absent within A/T gains (p = 0.52).

Two consequences follow. The architecture–function link the project set out to measure is present and robust in replicate 3. And the plasmid-state mixture behaves mainly as **noise** for MSP: standardizing for distal composition cut the spread of variant MSP changes by 25% (SD 1.52 → 1.15 pp) and *raised* the functional association. The handoff treated that standardization as a sensitivity check that might remove real effects. Here it removes more noise than signal.

### Literature-defined promoter elements

Classical LDLR repeats located by exact sequence in the WT element ([table](evidence/independent_review_20261020/promoter_elements.tsv)). "Rank" is where the element's mean standardized MSP change falls among all non-overlapping windows of the same width (0 = most negative).

| Element | Oligo / chr19 | Activity (log2) | MSP Δ, half A / B (pp) | Standardized MSP rank | 85–250-bp Δ (pp) |
|---|---|---:|---:|---:|---:|
| Repeat 1 (Sp1) | 89–104 / 11089353–368 | −1.09 | −0.75 / −0.46 | 0.10 | +0.51 |
| Repeat 2 (SRE-1) | 124–139 / 11089388–403 | −1.02 | +0.05 / +0.19 | 0.84 | −0.01 |
| Repeat 3 (Sp1) | 140–155 / 11089404–419 | −1.33 | −1.18 / −1.27 | **0.006** | +1.34 |
| P1 core (data-selected) | 54–62 / 11089318–326 | −1.57 | +0.28 / −0.43 | 0.39 | +0.17 |

Each element is similarly important for activity, but they differ in architecture. Sp1 repeats maintain accessibility, and their loss is accompanied by nucleosome-sized protection. Repeat 3 is the strongest accessibility element in the promoter; Repeat 1's CTCCTCC core (oligo 92–98) shows the same pattern more weakly. SRE-1 changes activity without changing accessibility. P1 changes a local footprint and activity without changing accessibility. This is the kind of mechanistic distinction the project set out to make. The repeats were defined in advance by the literature, but this review looked at the data first, so these remain hypotheses to confirm on new data.

### A calibration detail hidden by aggregate rates

The reports summarize WT A/A rejection rates over both directions. Split by direction, TF-coverage null rejections are almost all apparent losses. Welch at 50 barcodes gives 2,875 losses vs 74 gains in 21,000 contrasts. This is the expected behaviour of t-statistics on right-skewed, zero-inflated barcode means, where the mean and SD shrink together. The studentized permutation test is exact overall in WT A/A, because labels are exchangeable, but it uses |t|. At 100 barcodes its 5.4% overall rate splits into 5.2% for the loss tail and 0.2% for the gain tail, against 2.5% each nominally. A "loss" call from that test therefore carries roughly twice the nominal directional error, and gains are nearly undetectable. Any reuse should report equal-tailed or directional permutation *p* values.

## Assessment of the existing work

**Strengths worth keeping:** coordinate and sequence contracts; the full molecule ledger; conservative barcode grouping with a benchmark behind it; strict missing-versus-zero handling; the custom MM/ML parser; reusable barcode arrays; the WT A/A harness; and unusually honest reporting.

**What went wrong:**

- **Inference was aimed at the hardest unit.** Single variants × 20-bp windows of a caller-derived endpoint with a ~2–6% baseline, where per-variant noise is ~0.5 pp, have little information. The all-strata release gate then guaranteed that nothing would be released. Better tests cannot fix that.
- **The strongest evidence was left descriptive.** Library-level associations were reported without any inferential framework, even though spatially aware nulls are straightforward.
- **The largest nuisance variance was treated as a possible effect.** Distal/plasmid-state composition was relegated to sensitivity checks. It can be handled as a precision covariate once variants are shown not to shift it.
- **Candidates were prioritized by visual inspection rather than confounder-aware criteria.** Each audit layer (18 settings, resampling, joint states, viewers) increased within-preparation confidence in an observation, but none can address sequence-matched observability or biological replication.
- **Process cost is high relative to scientific decisions.** Immutable outputs, source guards and dashboards for every exploratory step slowed iteration. A legacy pipeline also still sits alongside `fsmpra/`, and `CLAUDE.md` still describes the legacy pipeline as primary.

## Recommended plan

### Stop

- Per-variant × window significance gating on caller coverage, and further single-preparation sensitivity layers for the same candidates.
- Building the experimental panel around oligo 125.
- New viewer features until endpoint definitions settle.

### Analysis on existing replicate-3 data (no new experiments)

1. **Freeze a small confirmatory family now, before new data.**
   - Elements: Repeat 1, Repeat 2/SRE-1 core, Repeat 3, P1 core.
   - Endpoints: distal-adjusted MSP, 85–250-bp protection, and local per-site m6A at unchanged A/T sites.
   - Predicted directions, as given in this review.
   - The library-level chromatin–activity association, with its circular-shift test.

   Replicate 3 is the discovery set; new preparations confirm.
2. **Make covariate-adjusted MSP/nucleosome endpoints primary**, adjusting for distal m6A at molecule or barcode level, and keep the unadjusted total effect alongside. First confirm across the library that variants do not shift distal m6A beyond noise, using split-half reliability of distal-m6A changes. If they do, keep total effects primary.
3. **Replace sparse small-footprint coverage with per-A/T-site retained m6A as the local readout, for all 633 SNVs.** Extend `fsmpra.raw_candidate_extract` beyond 27 candidates; it is the same full-BAM scan. The WT shared-site value in P3 is ~29%, against a ~5% footprint baseline. The outcome is dense, so the sparse-endpoint problem that consumed most of the statistical effort disappears. A starting model:
   - Per site, logistic on molecule-level retained calls.
   - Covariate: molecule labeling (distal m6A).
   - Uncertainty: barcode cluster bootstrap, with the existing WT A/A harness for calibration.
4. **Run a cheap internal observability diagnostic first.** Compute the oligo-125 and P1-core effects separately in the most highly methylated molecules and in chromatinized molecules. This is not a naked control. But if the 125 effect is just as large in near-saturated molecules, deprioritize it before ordering anything.
5. **If per-variant calls are still wanted**, use empirical-Bayes shrinkage on dense endpoints with bootstrap SEs, and report local false sign rates. Calibrate them with WT pseudo-variants, and use equal-tailed or directional p-values.
6. **Optional:** check directions at the element level in HepG2 run 1 (older condition, median ~150 molecules/SNV), comparing rather than pooling.

### Experiments, in order of information per cost

1. **Library-wide naked-DNA control.** Methylate the same purified plasmid library *in vitro* with the same Hia5 lot and buffer, at two or three doses that bracket the in-cell accessible-region labeling. Sequence and call it through the identical pipeline. This gives sequence-matched observability for all 633 SNVs at once and directly tests the P3/125 question. Protection can then be estimated as the in-cell minus naked per-site effect (difference on the logit scale).
2. **One or two biological replicates under the replicate-3 condition** (5 min, 3 µL per 2 M cells; same library). Collect RNA/DNA from the same transfections for MPRAsnakeflow. Activity and chromatin then come from the same cells, and step 1 of the analysis plan becomes a true out-of-sample confirmation.
3. **Keep kinetics, or the pre-`ft clear-kinetics` BAM, and the full ML range** for all new runs, so m6A probabilities can be recalibrated against the naked control.
4. **Only then a focused construct panel**, chosen from effects that replicate after naked normalization. Weight it toward the P1 core, Sp1 repeats and SRE core; keep oligo 125 only as an observability test case.

## Limitations of this review

- **One preparation.** Split-half reproducibility cannot distinguish deterministic sequence effects from biology.
- **The activity data are not an independent replicate or matched to these cells.** They come from a different published experiment, and the two activity sources correlate at r = 0.97, so they are not independent of each other either. Associations are not mediation.
- **The circular-shift null has limits.** It assumes the two profiles have no shared alignment beyond their own autocorrelation, and its p floor is ~0.005; within-class tests have fewer positions.
- **Standardization depends on the existing table.** It uses the evidence table's quartile standardization, not a new molecule-level model.
- **The candidate regions were inspected after seeing the data.** The P1 core is data-selected; the repeats are literature-defined. All of them require confirmation on new data.
- **The joint-state analysis covers only the 27 extracted candidates.**
- **No BAM-level stage was re-executed.** The handoff's earlier independent reconciliations were relied on for those.

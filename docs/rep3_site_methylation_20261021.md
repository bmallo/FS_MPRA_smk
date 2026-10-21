# Replicate-3 per-site methylation with proxy controls — results, 21 October 2026

This report interprets the [site-methylation run](evidence/site_methylation_20261021/README.md) (commit `fc4731a`, Klone job 42090780) and its [mask-4 follow-up](evidence/site_methylation_mask4_20261021/README.md) (commit `f8840e8`, job 42093363) against rules fixed before each run, in the [usage doc](site_methylation_usage.md) and the follow-up section below. Read the [independent review](INDEPENDENT_REVIEW_20261020.md) first. All numbers here were recomputed from the bundle's tables; corrections to the run summary are listed at the end.

## Bottom line

1. **Promoter architecture tracks function across the library, now with a calibrated per-variant readout.**
   - In chromatinized molecules standardized for labeling, MSP accessibility change correlates with published activity at r = 0.50. Nucleosome-sized protection gives r = −0.48. Raw element-wide m6A at unchanged sites, which involves no footprint caller, gives r = 0.49. All three reach the floor of the spatial null.
   - In the calibrated primary endpoints, 80 SNVs change accessibility and 37 change nucleosome-sized protection (q < 0.05), selected on chromatin data alone. Their direction matches the published activity change for 86% and 100% of them.
2. **Sequence-dependent observability is large and local, and it explains the local "footprint" candidates.**
   - Changing an A/T base shifts retained calls at neighbouring adenines by about ±30 pp in heavily methylated plasmids, about ±15 pp in chromatinized ones, and nothing in unmethylated ones. The effect reaches ±4 bp, so the earlier ±2-bp mask was too narrow.
   - Position 125 meets the prespecified "predominantly observability" criterion. The P1-core and SRE-1 local signals show the same signature; all three sit in A runs.
3. **The Sp1 repeats are the strongest, cleanest architectural elements.**
   - Repeat 3: accessibility −2.0 pp (z −20) and nucleosome-sized protection +1.9 pp (z +16), with little change in the methylation proxy.
   - Repeat 1 shows the same pattern, weaker.
   - SRE-1 lowers activity just as much without changing accessibility.
4. **After a ±4-bp mask, a calibrated local readout exists outside A runs, and it yields three new local candidates** (mask-4 follow-up; all 73 non-U endpoints pass calibration).
   - Oligo 115–117 (chr19:11089379–81): substitutions open the site and the whole element, reduce nucleosome-sized protection and *raise* activity. This is a repressive or nucleosome-positioning candidate at the left edge of the old P3 "footprint".
   - Oligo 56 in the P1 core (chr19:11089320): local protection loss with activity loss. This is an activator-footprint candidate, at a position free of the A-run artifact.
   - Oligo 165–166 and 190–191 (`TGACGTGG`): local protection loss without any activity change.
   - Oligo 125 stays unsupported: its context is homopolymer, which the prespecified rule marks unresolved, and the proxy class moves as much as chromatin.

## The measurement: an adjacent-adenine observability kernel

Mean change in retained m6A at unchanged A/T sites against WT, by offset from the substitution (pp; 633 SNVs; `effects/kernel.tsv`):

| Substitution | Class | ±1 bp | ±2 bp | ±3–4 bp | +6 / −8 bp | other 5–14 bp |
|---|---|---:|---:|---:|---:|---:|
| Adds an A/T | M (≥ 0.60 distal m6A) | −32, −33 | −10 | −3 | +0.4 / +1.9 | mostly < 1 |
| Adds an A/T | C (standardized) | −16 | −5 to −7 | −2 to −3 | 0 / +0.9 | mostly < 1 |
| Removes an A/T | M | +31, +25 | +6 to +8 | +2 to +4 | −3.4 / −3.1 | mostly < 1 |
| Removes an A/T | C | +17, +14 | +5 to +7 | +2 to +3 | −1.1 / −0.8 | mostly < 1 |
| Keeps A/T count | M / C | ±1–2 | ±1 | ±1–2 | ≈ 0 | ≈ 0 |
| Any | U (≤ 0.01) | ≤ 0.4 | ≈ 0 | ≈ 0 | ≈ 0 | ≈ 0 |

The ±1–4-bp core has z = 4–15 per offset. A smaller lobe of opposite sign follows at +6 and −8 bp (M: z −3.5 and −2.2 for losses).

Neighbouring adenines compete for retained calls. Adding one suppresses calls at its neighbours, and removing one raises them, in proportion to how heavily the molecule is labeled. U molecules show nothing, which locates the effect in the labeling/calling of methylated DNA rather than in false positives. This is a property of the measurement, not of chromatin.

Two consequences follow:
- Any footprint or methylation window that includes sites within ~4 bp of an A/T-changing substitution is strongly biased, and weakly biased out to ~8 bp. FiberHMM's small-footprint calls, which depend on runs of unmethylated A/T, inherit the same bias.
- Pooling A/T gains with losses cancels the kernel. The run summary evaluated rule 2 on pooled substitutions, which is how it reported only scattered offsets. By class, rule 2 is met decisively.

## Local candidates, against the prespecified rules

**Position 125 (chr19:11089389) — rule 3: predominantly observability.**
- In the WT_P3 window, M-class methylation rises with the in-cell signal for all three substitutions: M +4.4 / +6.1 / +6.2 pp against C +4.1 / +4.4 / +5.5 pp (A>C/G/T).
- The chromatin criterion (normalized z > 3 for all three under all threshold pairs) fails: z = 2.82 and 2.56 at 0.005/0.50.
- Position 125 lies in `AAAA` (oligo 124–127). The library-average kernel predicts only +0.4 pp of M-class change at this window's sites; the observed +4.4 to +6.2 pp points to a much stronger homopolymer-specific effect, or to protein contacts in some M molecules, which this proxy cannot distinguish.
- Normalized accessibility remains +2.4 to +4.6 pp, but that residual is not robust. It should not be presented as a protein footprint; this matches the review's concern.

**P1 core (oligo 54–62) — rule 4: prediction "no M change" fails.**
- The region mean is M +2.8 pp (z 16).
- It is concentrated in 61A>T and 62A>C/G/T, which break the `AAAAA` run at 61–65: M +9 to +11 pp, C +10 to +12 pp.
- Elsewhere in the core (54–60), M and C changes are small (mostly 0–3.5 pp), while activity falls for every substitution at 54 and 56–60.
- P1 remains a functional element, but its chromatin footprint is unresolved.

**SRE-1 (Repeat 2).** Local normalized +0.9 pp (z 5.3) with M +1.9 pp (z 16): M-dominated, and the region contains the 124–127 `AAAA` run. Accessibility change is +0.2 pp (z 2.4), which fails "≈ 0" only on the |z| ≤ 2 criterion; the magnitude is a tenth of Repeat 3's. The data support the review's distinction: SRE-1 affects activity, not accessibility.

## Promoter architecture: calibrated per-variant and region results

The calibration gate passed for `flank:A_norm`, `msp_C_std` and `nuc_85_250_C_std`. `local:A_norm` failed narrowly (15/200, upper bound 0.121). It also keeps an A/T-class dependence (gain −0.63, loss +0.67 pp), the kernel's signature, so it stays excluded.

Discoveries use BH over those three endpoints (`q_bh_primary_family` over four endpoints gives 78/80/39):

| Endpoint | SNVs at q < 0.05 | Direction | Agrees with activity sign | Concentration | Median activity of hits (log2) |
|---|---:|---|---:|---|---:|
| Accessibility (`msp_C_std`) | 80 | 62 down / 18 up | 86% | Repeat 3: 25, Repeat 1: 11, SRE-1: 7 | −0.76 |
| Nucleosome-sized protection (`nuc_85_250_C_std`) | 37 | 30 up / 7 down | 100% (more protection, less activity) | Repeat 3: 21, Repeat 1: 4 | −2.36 |
| Methylation 11–30 bp away (`flank:A_norm`) | 74 | 50 down / 24 up | 77% | Repeat 3: 20, Repeat 1: 7 | −0.68 |

Sign agreement would be ~50% if chromatin and activity were unrelated. Selection used only chromatin data.

Region means over SNVs inside each element are below (primary run). Across the three threshold pairs and the supported-barcode cohort, accessibility and nucleosome-sized changes agree within 0.2 pp, and signs agree everywhere. Normalized and M-proxy values vary by up to ~1 pp; for example, Repeat 3's M change ranges from +0.1 to +1.1 pp.

| Element | Normalized local | M proxy | Accessibility | Nucleosome-sized | Activity |
|---|---:|---:|---:|---:|---:|
| Repeat 1 (Sp1) | −1.42 (z −7.6) | −0.64 | −1.00 (z −9.4) | +0.85 (z +6.9) | −1.09 |
| SRE-1 | +0.92 (z +5.3) | +1.85 (z +15.8) | +0.24 (z +2.4) | −0.18 | −1.02 |
| Repeat 3 (Sp1) | −2.24 (z −11.3) | +0.47 | −2.02 (z −20.6) | +1.92 (z +16.3) | −1.33 |
| P1 core | +2.96 (z +10.7) | +2.82 (z +15.9) | −0.39 | +0.25 | −1.57 |

## Library association (one common SNV set, ≥ 300 molecules)

| Change in | Variant r | Position r | Rotation-null max \|r\| | p |
|---|---:|---:|---:|---:|
| Accessibility, all molecules | +0.34 | +0.42 | 0.22 | 0.005 |
| Accessibility, chromatinized + standardized | **+0.50** | +0.56 | 0.31 | 0.005 |
| Nucleosome-sized protection, all | −0.46 | −0.55 | 0.25 | 0.005 |
| Nucleosome-sized protection, standardized | **−0.48** | −0.57 | 0.30 | 0.005 |
| Raw m6A, whole element, standardized | **+0.49** | +0.53 | 0.34 | 0.005 |
| Raw m6A, 11–30 bp, normalized | +0.35 | +0.37 | 0.31 | 0.005 |
| Small-footprint coverage (caller) | +0.19 | +0.23 | 0.25 | 0.03 |

- **630 SNVs; LDLR activity source.** The supported-barcode cohort agrees within ±0.01–0.02.
- **Removing proxy-class molecules and standardizing for labeling raises the accessibility association from 0.34 to 0.50.** Nucleosome-sized protection changes little (−0.46 → −0.48).
- **The raw-methylation result needs no footprint caller.** It is an independent measurement of the same association.

**What the M class contains.** Beyond the kernel's reach (11–30 bp), M-class changes have only modest excess variance across SNVs (z SD 1.42 against 2.22 in C) and are unrelated to activity (r = +0.06, p = 0.42). Heavily methylated plasmids therefore behave as a non-functional proxy, as intended. They still share some non-functional, sequence-linked variation with C at a distance (r = 0.33), so they are not a perfect naked control. **Class membership is independent of genotype:** standardized class-fraction deviations have SD 1.01 (U) and 1.04 (M), and neither fraction tracks activity.

## Corrections to the run summary

1. **`msp_ALL` r = 0.151 compared different SNV sets.** It included three nearly unmeasured SNVs (160A>G has one molecule and a −53 pp "effect") that the standardized endpoint excludes. Like for like, the gain is 0.34 → 0.50. Rule 5 still holds. `site_effects` now uses one common SNV set (≥ `min_association_molecules`, default 300).
2. **The class-dispersion ratios (1.62 and 0.60) were a non-robust statistic** distorted by those tiny cohorts. They are replaced by a Pearson dispersion: about 1.0 for both classes here.
3. **Rule 2 should be read per A/T-change class.** Pooling cancels the kernel's opposite signs. The kernel is now also stratified by homopolymer context.
4. **With 200 draws, the gate fails an exactly calibrated endpoint 30% of the time.** Its 17/73 failures (23%) are consistent with chance; the median z SD is 0.99 and the median rejection rate 0.050. The follow-up uses 1,000 draws, where a false failure is very unlikely. No exclusion here is reversed after the fact.

## Updated priorities

- **Lead mechanistic result:** Sp1 repeats maintain accessibility, and their disruption brings nucleosome-sized protection and lower activity. Repeat 3 is strongest. SRE-1 changes activity without accessibility.
- **Deprioritize oligo 125** as a protein footprint, as rule 3 dictates. In any construct panel, keep it only as an observability test case. Weight the panel toward Repeat 3/Repeat 1 substitutions spanning A/T gain, loss and retention, plus SRE-1 core substitutions as the activity-without-accessibility contrast.
- **Local footprints:** unresolved wherever substitutions touch A runs (154 of 630 SNVs). Elsewhere, oligo 117, 56 and 191 are the best-supported local candidates; see the follow-up results. A separately prepared, in-vitro-methylated library is still the decisive control, most of all for homopolymer contexts.
- **Revised construct panel:** [`evidence/validation_panel_v2_20261021.tsv`](evidence/validation_panel_v2_20261021.tsv) has 11 constructs with exact mutable sequences and supporting numbers. It replaces the 20 October panel's emphasis on oligo 125.

## Follow-up run, prespecified before it is executed

Config `config/ldlr_site_methylation_mask4.json` reuses the existing extraction:
- mask ±4 bp;
- `local` = 5–10 bp, `flank` = 11–30 bp, `near` = 1–4 bp (diagnostic only);
- 1,000 WT pseudo-variants.

Mask and window widths come only from the library kernel above. The mask covers the dominant ±1–4-bp core; the smaller +6/−8 lobe stays inside `local` on purpose (a wider mask would leave no local window), and rule 3 tests whether normalization removes it. Instructions: [`klone_agent_prompt_site_methylation_mask4.md`](klone_agent_prompt_site_methylation_mask4.md).

1. **Calibration:** `gate_pass` per endpoint, as before.
2. **Homopolymer reach:** if the M kernel for homopolymer-context SNVs (`context == 'homopolymer'`, per A/T class) has |mean| > 2 SE at any |offset| ≥ 5, local endpoints for homopolymer-context SNVs are reported as unresolved.
3. **Residual artifact:** `local:A_norm` is usable only if it passes the gate and its mean Δ for A/T gains and losses differ by less than 3 SE (SE from per-SNV deltas within each class).
4. **Position 125:** rule 3 is reapplied unchanged to the new windows.
5. **Library association:** the same table on the common SNV set; no new endpoints are added to the primary family.

## Mask-4 follow-up: results

Klone job 42093363, 6:39 for all four runs; the extraction was reused unchanged, and the pinned files are identical to `fc4731a`. Numbers below were recomputed from the [bundle](evidence/site_methylation_mask4_20261021/README.md).

1. **Calibration:** 73/73 non-U endpoints pass in the primary run (z SD 0.97–1.08, median rejection rate 0.053, 1,000 draws), including all four primary endpoints. One M-proxy endpoint fails at thresholds 0.02/0.70 (`adjacent_143_150:M`, z SD 1.24); none fails in the other runs.
2. **Homopolymer reach: rule triggered.** The M kernel for homopolymer-context A/T gains exceeds 2 SE at offsets −14, −8 and −5 (−5: −6.1 pp, z −3.1), so local endpoints for the 154 homopolymer-context SNVs are unresolved. Their core artifact is also stronger: about −46/+44 pp at ±1 bp and ±10 pp at ±3 bp, against about −30 pp (gains) and +6 to +18 pp (losses) at ±1 bp, and ≤ 4 pp at ±3 bp, elsewhere.
3. **Residual artifact:** `local:A_norm` is usable. The gain-minus-loss difference is −0.59 pp (z −2.7) overall and −0.24 pp (z −0.9) outside A runs.
4. **Position 125:** the observability clause is met for 3/3 substitutions (M +5.1/+3.8/+6.3 pp against C +4.5/+3.9/+5.9 pp), and the chromatin clause is not. All three are homopolymer context. Not a supported protein footprint.
5. **Library association:** unchanged for MSP and nucleosome endpoints (r = 0.50 and −0.48). The local 5–10-bp endpoint gives r = 0.24, just beyond its rotation null (0.26 against 0.23).

**Primary-family calls (q < 0.05):** local 93, flank 75, accessibility 80, nucleosome-sized 38, with sign agreement with activity of 74%, 78%, 87% and 100%. Of the 93 local calls, 28 are in homopolymer context and therefore unresolved. Of the 65 interpretable ones, 38 also have a flank call and 30 an accessibility call, so most local calls are the broad Sp1-repeat architecture seen at shorter range.

**Position-consistent local calls outside A runs.** These are positions where at least two substitutions are called with the same sign; 13 of 159 testable positions qualify:

| Oligo (chr19) | WT context | Called | Local Δ | M proxy | Accessibility | Activity | Reading |
|---|---|---|---:|---:|---:|---:|---|
| 146, 148–150 (11089410–14) | Repeat 3 `CTCCTCCCCC` | 2–3/3 each | −2.0 to −4.4 | −0.8 to +6.4 | −3.1 to −4.4 | −2.1 to −2.7 | Sp1-repeat architecture at short range |
| 115, 117 (11089379, 81) | `GTGAAG[A]CATTTG` | 2/3, 3/3 | +3.7, +6.2 | ≈ 0 | +2.6, +3.7 | +1.1, +0.9 | Protected element whose loss opens chromatin and raises activity |
| 56 (11089320) | `CTTCAC[G]GGTTAA` | 2/3 | +4.2 | +1.0 | +0.3 | −1.25 | Local protection loss with activity loss |
| 165, 166 (11089429–30) | `AAACCT[C]ACATTG` | 3/3, 2/3 | +2.8 | +1.8, −0.3 | +0.2, +0.8 | +0.1, +0.3 | Protection without function |
| 190, 191 (11089454–55) | `TGACGT[G]GGCCCC` | 2/3 each | +3.3, +3.8 | −3.3, +1.1 | +1.2, +0.9 | +0.2, +0.1 | Protection without function (`TGACG`, CRE-like) |
| 198, 202 (11089462, 66) | `GCCCCG[A]GTGCAA`, `CGAGTG[C]AATCGC` | 2/3 each | −3.2, −3.5 | small | −1.4, −2.0 | ≈ 0 | Small local losses, no function |

The figures are mean percentage-point changes over the position's substitutions. The 117, 56 and 191 candidates hold in all four runs; their local z ranges are 3.1–4.3 (117A>T), 3.2–4.7 (56G>C) and 6.4–7.4 (191G>T). Twenty-four local calls have no flank or accessibility counterpart; several of those carry large M-proxy changes (78T>C, 113G>C, 129C>G), and normalization can turn such changes into apparent effects, so single-substitution local-only calls are not prioritized. Motif identity is not assigned here.

## Transcription-factor footprint evidence (post hoc)

These analyses were chosen after seeing the results above. They describe one preparation and are not prespecified tests. They reuse only committed evidence (mask-4 effects, published activity, JASPAR 2024 and HOCOMOCO v12 matrices). `python -m tools.variant_hit_figures figures --output <dir>` reproduces the tables and the [overview](figures/top_hits_20261021/top_hits_overview.png) and [motif-site](figures/top_hits_20261021/motif_sites.png) figures. Per-variant WT-vs-variant occupancy figures need one Klone step; see [the Klone prompt](klone_agent_prompt_hit_figures.md).

**Proxy-dominated check.** A local change counts as proxy-dominated when the heavily methylated proxy moves at the same site (|z| > 3) by at least as much as chromatin. Normalization can then turn an observability change into apparent opening or closing. This check flags 202 of 630 SNVs.

**1. Library level: a modest activator-footprint signature.** The question is whether local chromatin change at 5–10 bp tracks activity once broad accessibility and nucleosome-sized change are accounted for. Local change and activity were each regressed on the MSP and nucleosome-sized changes, and the residuals correlated. The test uses SNVs outside A runs with at least 300 molecules, and the circular position-shift null ([table](figures/top_hits_20261021/footprint_association.tsv)).

| SNVs | Local endpoint | n | Partial r | Position-level partial r | Null max | p (shift) |
|---|---|---:|---:|---:|---:|---:|
| Outside A runs | chromatin (`local:C_std`) | 473 | −0.23 | −0.26 | 0.33 | 0.021 |
| Outside A runs | normalized (`local:A_norm`) | 472 | −0.12 | −0.17 | 0.22 | 0.036 |
| … and not proxy-dominated | chromatin | 306 | −0.24 | −0.29 | 0.24 | 0.008 (floor) |
| … and not proxy-dominated | normalized | 305 | −0.13 | −0.15 | 0.26 | 0.14 |

The sign is the one expected for activator footprints: beyond broad opening, SNVs that unprotect their immediate neighbourhood lower activity. The effect is small (about 6% of residual variance) and clearest in the chromatin class. The normalized version does not survive the proxy filter.

**2. Element level: motif disruption against chromatin and function.** For each named element, every substitution inside the best WT motif match is ranked by how much it weakens that match. Spearman ρ is given with p and n; local chromatin excludes A-run and proxy-dominated points.

| Element (oligo; chr19) | Best WT motif model (p) | Activity | Local chromatin | Accessibility | Nucleosome-sized |
|---|---|---|---|---|---|
| P1 core 52–62 `TCACGGGTTAA` (11089316–26) | ZN652, HOCOMOCO v12 (7e-5) | **−0.87** (6e-11, 33) | **+0.46** (0.05, 19); +0.49 (0.004, 33) unfiltered | −0.19 (0.3) | +0.23 (0.2) |
| Repeat 1 82–102 (11089346–66) | SP4, HOCOMOCO v12 (3e-5) | **−0.64** (2e-8, 63) | −0.10 (0.6, 26) | **−0.42** (5e-4) | **+0.40** (0.001) |
| SRE-1 127–136 (11089391–400) | SREBF1, JASPAR MA0595 (1e-6) | **−0.46** (0.01, 30) | +0.19 (0.6, 12) | −0.30 (0.1) | +0.19 (0.3) |
| Repeat 3 143–152 (11089407–16) | WT1, HOCOMOCO v12 (1e-6) | **−0.40** (0.03, 30) | +0.02 (1.0, 16) | **−0.46** (0.01) | **+0.54** (0.002) |

The four elements behave in three different ways:
- **P1 core: motif, local footprint and activity line up, with no broad accessibility change.**
  - Motif strength nearly ranks activity (ρ −0.87). Substitutions that strengthen the match (55C>A, 55C>G) are neutral (+0.06, +0.08 log2).
  - 56G>C and 56G>T open the 5–10-bp neighbourhood in chromatin by +3.1 and +2.9 pp (z 6.2, 7.1), against +1.6 and +0.8 pp in the proxy. Activity falls by 1.2–1.6 log2, and MSP does not change.
  - At 54, 58 and 59, chromatin opens at 8 of 9 substitutions (+0.9 to +2.4 pp), but the proxy rises nearly as much, so these sites support the pattern only weakly.
  - The motif model is not an identification. ZNF652 has mostly been described as a repressor, so the activator acting here may be another C2H2 zinc-finger or a factor with similar specificity.
- **Sp1 repeats: motif loss acts through the open region, not through a resolved footprint.** Disruption lowers accessibility, raises nucleosome-sized protection and lowers activity. The 5–10-bp readout does not follow the motif.
- **SRE-1: function without a visible chromatin change.** Motif loss lowers activity without changing accessibility. Most local readouts here are unresolved (A run), so a footprint can be neither seen nor excluded.

**3. Oligo 115–117: opening, and creating an AP-1/CRE-like match adds activity.**
- Every substitution at 115 and 117 opens chromatin (MSP +1.3 to +4.6 pp; local chromatin +1.3 to +3.5 pp, proxy −1.5 to +2.0 pp) and raises activity.
- The alleles that create or strengthen AP-1/CRE-like matches gain the most:
  - 117A>T (CREM p 3e-5 → 6e-7, JDP2): +1.45 log2;
  - 115A>C (CREB5, ATF1, ATF4): +1.31;
  - 115A>T (FOS, FOSB): +1.22.
- The others gain less (+0.54 to +0.84), although 117A>C opens just as much (+4.6 pp) as 117A>T.
- A reading consistent with this is loss of a repressive or nucleosome-positioning protection, plus a new activator site for some alleles. These are small numbers, and the analysis is post hoc.

**4. Not supported: oligo 73–74.**
- All three 73 substitutions lower the proxy by 21–23 pp in the local window (z ≈ −30), more than chromatin (−8 pp). The normalized readout shows apparent opening (+3.7 to +5.7 pp) only because of that. The window covers the A run at oligo 61–65.
- Two of three 74 substitutions are proxy-dominated.
- The activity losses at 73 (T>C, T>G) and 74 are not explained by any resolved chromatin change here.

## Limitations

- **One preparation.** All results come from a single HepG2 preparation.
- **The controls are proxies.** U and M come from the same reaction, and M carries some non-functional sequence-linked variation.
- **Retained calls only.** Calls were stored at ML ≥ 240, so sub-threshold probabilities are unobserved; the Globus kinetics would allow re-calling.
- **Activity data are external.** They come from separate published experiments.
- **Correlation is not mediation.** Agreement between chromatin and activity does not show that one causes the other.

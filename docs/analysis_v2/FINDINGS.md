# Analysis v2: approach and findings (Phase A + follow-ups, 2026-10-28 to 2026-11-02)

Each step below gives the question, what was done, the numbers, and what carries forward. Numbers come
from `config/analysis_v2/rules.json` (each step's `results` / `decisions`) and the result files named.
"Post hoc" marks anything chosen after seeing results. Plain-language slides for every step are in the deck
(`deck/`, published at https://claude.ai/artifact/NYSJ7cHat24e2SRj4SErwG).

Vocabulary: an adenine with an m6A call is **marked = accessible**; one without is **unmarked =
inaccessible**. Unmarked can also mean the adenine is hard for Hia5 to methylate (A runs), so "inaccessible"
is a measurement, not proof of protein. **Labeling** = a molecule's m6A fraction on the plasmid backbone
("distal"); classes U ≤ 0.01, C 0.01–0.6, M ≥ 0.6. **Oligo** coordinates follow Kircher et al. 2019: the 211-bp
element spans oligo 21–231, and SNV names use them (`79C>A` = oligo 79).

## Data in one paragraph

LDLR promoter saturation-mutagenesis MPRA library (633 designed SNVs plus 1-bp indels), transfected into
HepG2 (runs 1–3) and HEK293T (1 run), Hia5 Fiber-seq in cells, PacBio HiFi. Each HiFi read carries its
own genotype and an N15 barcode (clone). The main dataset is **HepG2 run 3**: 10.5 M HiFi reads →
2.59 M analysis-eligible WT/single-SNV molecules; 633 SNVs, median 775 clones and ~2,100 molecules per
SNV. **HEK293T** is a different construct (element reversed), pipeline and m6A caller: 197 SNVs observed
with enough molecules, median 4 clones each. Activity comes from the published MPRA (Kircher 2019, HepG2,
two elements "LDLR" and "LDLR.2"), a separate experiment.

## Phase A steps

### A1 Data inventory (`a1_inventory/`, `config/analysis_v2/runs.json`)
- HepG2 runs 1 and 2 differ only in sequencing date: the one true replicate pair (but shallow).
- Run 3 differs from them in 7 ways (5 vs 10 min reaction, 3 vs 1.5 µL Hia5, chemistry, CCS, Jasmine
  26.1.3 vs 2.7.99, fibertools, date), so run 1/2 vs 3 is not a clean reaction test.
- HEK293T differs from every HepG2 run in 8–9 ways, including the m6A calling tool, so cells are never
  pooled or compared directly. HEK Hia5 volume and cells per reaction are not recorded.

### A2 Read accounting (`a2_accounting/`, `a2_ledgers/`)
- Run 3 funnel (first filter failed): not plasmid 17.1%; not full length 4.0%; <10 passes or low rq 17.5%;
  element matches no design 13.8%; exact design but a base < Q20 4.5%; **designed 1-bp indel 18.0%
  (1.89 M, never analysed)**; eligible 24.7% (2.59 M).
- HEK293T: 3.06 M plasmid reads, 0.31 M eligible (145 k SNV, 163 k WT).
- Designed indels (146 deletions, 637 insertions, ~1,200 molecules each) are the largest unused resource.

### A3 Genotype and barcode audit (`a3_genotypes/`)
- Per-read genotypes reproduce across runs (99.4–99.98%, same barcode in independent runs) and across two
  independently written genotypers (99.99%).
- 1.6% of WT-clone reads are misread as a variant, mostly homopolymer indels.
- Rules fixed in rules.json (`molecules`, `clones`): clone = barcode + genotype; a read is dropped when its
  barcode has ≥ 3 full-length reads with a ≥ 90% majority genotype that the read contradicts; no
  readable barcode = its own clone; partial reads stay out; Q20 and ≥ 10 passes kept.

### A4 Informative bases (`a4_informative/`)
- The element has 101 A/Ts (one methylatable adenine per A/T pair). Spacing never limits: a footprint
  ≥ 11 bp is detectable at every position.
- Adenines inside A runs are nearly invisible on naked DNA too (Hia5 sequence bias, not protein): oligo
  59–65 (A run), the TTT/AAAA runs near 114–130, and others.
- **The plasmid-LLR caller drops the edited base for 303 SNVs** (every A/T→G/C and every A↔T swap: it only
  reads reference-matching A/Ts) and ignores gained adenines. 472 of 1,416 designs change a neighbour's
  run context. Every SNV carries three annotations (A/T change, caller drop, neighbour change):
  `a4_informative/visibility/variants_visibility.tsv`.

### A5 m6A call QC (`a5_m6a_qc/`)
- About 0.3% calls in U molecules (upper bound on false positives) vs 67–69% in labeled ones; up to 3% at a
  few positions (oligo 189, 198, 135, 78).
- No dam methylation at the 35 backbone GATC sites (prep is dam−/dcm−).
- Calls hold up at 10–14 passes, but unlabeled molecules have more passes (U share 3% → 12–30%), and read
  orientation shifts some adenines by up to 28 points: passes and strand are covariates (checked in A9).

### A6 M and U as stand-ins for controls (`a6_proxies/`)
- Labeling is a continuum (no separate naked population). Run 3: U 7.5%, C 70.9%, M 21.6%.
- **M is not protein-free**: element spans ≥ 80 bp unmarked in 2.0% of M (naked plasmid 0.01%); at oligo
  60–86, 5.6–10.9% of M molecules have no m6A (0.1% when shuffled).
- So the LLR emissions (built from WT M = accessible, U = protected) are "relative to in-cell M", and absolute
  occupancy is not reportable before naked controls (B1).

### A7 How SNVs alone change calls: in-silico edits of WT molecules (`a7_snv_edits/`)
- Dropping one A/T changes 20–80 bp call coverage by a median 0.9% of C molecules (max 4–5%); calls are
  reshaped more than added or removed.
- Masks: the **A7 primary mask** (the edited base when the caller loses it + adenines whose run context
  changes; 93 distinct sets) is applied symmetrically to variant and WT. The **±5 bp mask** (all A/Ts within
  5 bp) is the sensitivity analysis; it costs about twice as much signal.
- The plan's check failed: the caller artifact explains only ~15% (29% with a context model) of the v1
  "A/T-loss" offset. A/T-changing SNVs therefore stay flagged and are tested apart from the rest.

### A8 Caller confidence: spike-ins and reproducibility (`a8_confidence/`)
- A perfectly protected 10-bp stretch inserted into M molecules is found in 95–98% of them; 20 bp
  everywhere. Harder at oligo 56–65, 78–90, 120–129, 135–141, 149–152, 188–199, 205–211.
- **82–86% of low-labeled C molecules (0.01–0.2) already have a call at any 20-bp spot**, so call-based
  endpoints use labeling ≥ 0.2; low-labeled molecules enter only the caller-free endpoint.
- Called lengths are 0–4 bp longer than spikes up to 80 bp; neighbouring size classes are not different
  proteins (report length as continuous).
- Clone halves of one run reproduce (r 0.956–1.000). HepG2 runs 1 vs 3 do not (|z| > 1.96 at 45–87% of
  positions): compare variant vs WT within a run only.

### A9 Statistical test (`a9_statistics/`)
- Endpoints: **E1** = adenine unmarked (caller-free; labeling > 0.01); **E2** = position covered by an LLR
  call of a size bin (labeling ≥ 0.2).
- Test: D = mean over the variant's molecules of (y − p_WT(labeling stratum)), strata 0.05 wide.
  Variance = binomial (both sides) + τ² · Σ_c (n_c/n)², with τ² the between-clone variance from WT clones
  (**random effects; amendment 1, post hoc**: the prespecified clone sandwich failed, HEK 57% p < 0.05).
  A feature is testable only with ≥ 10 molecules, ≥ 3 clones and ≥ 5 expected minority events
  (**amendment 2, post hoc**). Calibrated on WT-only pseudo-variants: HepG2 E1 4.97% p < 0.05, λ 1.00.
- Covariate balance (labeling, passes, strand, read quality) between variant and WT molecules is fine.
- Implementation: `tools/analysis_v2/a9_framework.py` (`Reference`, `contrast`); reusable through
  `tools/analysis_v2/access.py`.

### A10 Variant effects (`a10_effects/`)
- HepG2, A7 mask: **602 of 633 SNVs change the inaccessible fraction somewhere** (5,633 hit adenines);
  ±5 bp mask 467; ±5 bp + empirical null (post hoc) 432. E2 (calls): 209 SNVs.
- **Changes fall off with distance**: hit share 59.8% at 0–5 bp, 35.1% at 6–10, 15.8% at 11–20, 8.6% at
  21–50, 3.3% beyond 50 bp. **All 224 changes above 15 points lie within 10 bp of the SNV** → labelled
  "possible Hia5 sequence effect until B1".
- **17 SNVs at oligo 143–151 (Repeat 3 Sp1 site) raise protection across the whole element** (all three
  alleles at 146, 148, 149). 12 other SNVs have element-wide effects.
- Global test per SNV = Simes over its features (amendment 3, before results), Benjamini–Bogomolov
  hierarchical FDR q = 0.05, separately for "drops its base" SNVs and the rest.
- HEK293T: the prespecified pseudo-pool null was invalid (each SNV ≈ one clone); a clone-mean motif test
  (post hoc) gives 7 of 9 E1 and 6 of 9 E2 motif pools at BH; "drops its base" pools are negative in all five
  motifs and are not interpreted before B1.

### A11 Motifs and function (`a11_motifs_function/`)
- Motif windows: literature repeats R1 Sp1 89–104, R2 SRE-1 124–139, R3 Sp1 140–155 (primary; scored with
  JASPAR MA0079.5 SP1 and MA0596.1 SREBF2); v1 FIMO windows ZNF652-like 52–62, SP4 82–102, SREBF1/2
  127–136, SP1/SP3/KLF 140–156, CREB1/CREM 182–192 (secondary).
- **Only Repeat 3 shows a dose-response**: breaking the Sp1 site lowers 10–40 bp call coverage there
  (+0.26 points per log2-odds, q 0.017). No motif dose-response on the caller-free endpoint.
- **Element-wide shift** = mean change in inaccessible fraction at adenines > 50 bp from the SNV (±5 mask
  run). It tracks nucleosome-size call coverage (r 0.87) and, within position, lower published expression:
  −0.35 log2 per point (CI −0.50 to −0.20; same in LDLR.2). Robust to removing the Repeat 3 SNVs.
- The 17 Repeat 3 SNVs have median activity −2.59 log2 (about 6-fold lower) vs −0.10 for others.
- All function links are correlations across SNVs from separate experiments, not mechanisms.

## Follow-ups after Phase A (exploratory)

### Footprint at oligo 60–86 (`fp60_86/`, `tools/analysis_v2/fp60_86.py`)
Real protection over ~66–81 is supported: coverage rises with labeling (oligo 73: 2.3% → 6.0–6.3%) while
simulated naked DNA falls to 0; variants in the A/T-free gaps (66–69, 79–81) remove it and raise marking at
core adenines 5–10 bp away; alleles at one position go opposite ways (71T>A −3.5 vs 71T>C +2.9). The
20–40 bp size is stretched by the unmarkable A run 59–65 and the TT at 85–86. The 86T>X spike seen earlier
was an unmasked-display artifact. Link to activity is weak (core variants −0.5 to −1.3 log2). A Hia5
sequence effect is not excluded until B1. The 182 (CREB1/CREM window) call falls with labeling like the
naked simulation: a candidate artifact.

### Footprint at oligo 114–138 (`fp114_138/`, `tools/analysis_v2/fp_region.py`)
Mostly a 12-bp call at 119–130 resting on four weakly marked adenines (120, 127, 128, 130) next to the
TTT/AAAA runs. Coverage falls with labeling like simulated naked DNA (artifact signature); only 0.6–2.2%
remains at labeling ≥ 0.6. No sequence-specific variant evidence. The SRE-1 (A/T-free CCCC at 131–134) is
functionally strong but nearly invisible to Fiber-seq.

### Data-defined footprint map (`fp_features/`, `tools/analysis_v2/fp_features.py`)
Footprints are peaks in the WT (midpoint × length) density of LLR calls (C labeling 0.2–0.6, lengths
10–90), found with v1's `fsmpra.footprint_sizes.find_features`: **24 footprints F01–F24**
(`footprints.tsv`: box in midpoint/length, span, WT occupancy, confident share, occupancy by labeling).
Examples: F04 (span 60–86) 86% confident calls; F11 (span 117–129, the 114–138 call) 3% confident.
Per-SNV occupancy tests (A9 test, BH across footprints within an SNV) for masks none / A7 / ±5 and all vs
confident calls: significant SNV × footprint pairs A7 166 (confident 89), none 210 (87), ±5 87 (56).
With the A7 mask and confident calls, 79 of the 89 pairs are at F04 (60–86); then F18 5, F07 3.

### Confidence of single calls (`docs/footprint_confidence_design.md`)
log10 P_chance = dose(L) · Σ log(1 − p_a) / ln 10 over the call's unmarked adenines, with
dose(L) = ln(1 − L) / ln(1 − 0.70) and p_a the run's accessible rates. Confident if P_chance < 10⁻³.
Display approximation only (explorer); the full model (per-position local FDR) waits for B1/B2.

### Promoter explorer (https://claude.ai/artifact/JTJkTUZ5NtDK6WSVt4BYsH, Version 7)
Interactive page over HepG2 run 3 and HEK293T (separate tabs): lollipop overview of every SNV on several
measures, then a per-variant track view (motifs, WT footprints, sequence, activity, WT and variant
molecule pileups, footprint occupancy and change, inaccessible fraction and change), mask (A7 / ±5 /
none) and call (confident / all) controls. How to rebuild, render and republish: `DATA_GUIDE.md`.

## Phase B (waiting for data the user will generate)
- **B1 naked controls on this variant library**: unmethylated (PCR, no-Hia5), saturated Hia5, and a
  matched-dose Hia5 series on bare DNA. Include the variant sequences with the largest near-SNV changes
  (25C>A, 211G>T, 113G>T, 73T>A/C/G) and the 60–86 set (68C>A/G/T, 80G>A/C/T, 71T>A/C, 77A>G/T).
- **B2 WT time course** in cells and on naked DNA.
- Then: rebuild emissions per variant sequence, rerun A6–A10, build the confidence model and the
  sequence-information QC tool, and revamp the pipeline.

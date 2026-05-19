# FS-MPRA statistical framework (plain-language reference)

What the pipeline does, end to end, and why — every term defined.
Status: Phase 1 complete and validated (exact FDR control at
production scale with the default config).

## The question being answered

For each single-nucleotide variant in a saturation-mutagenesis
promoter library, does that mutation **change protein/nucleosome
footprint occupancy** along the promoter, compared to wild-type — and
is that change bigger than sequencing/sampling noise?

## Glossary (terms used below)

- **Read / molecule**: one sequenced plasmid. Single-molecule, no PCR,
  so every read is an independent physical molecule.
- **WT (wild-type)**: reads of the unmutated construct — the baseline.
- **Variant**: reads whose *only* promoter change is one specific SNV
  (multi-mutation reads are excluded by default).
- **Footprint track** (FiberHMM `MA` tag): `tf` = transcription-factor
  footprints, `nuc` = nucleosomes, `msp` = accessible patches. The
  analysis runs per track.
- **Occupancy**: at one position, the fraction of molecules whose
  footprint (of a given track) covers that position — a number 0–1.
- **Δ (delta)**: variant occupancy − WT reference occupancy, per
  position. The effect of the mutation.
- **NC (nucleosome count)**: nucleosomes on a molecule. A strong
  confounder of footprint occupancy.
- **Null distribution**: the spread of Δ expected from noise alone
  (no real effect), built empirically by resampling WT.
- **p-value**: chance of seeing an effect this big if nothing is
  going on. **FDR / BH q-value**: false-discovery-rate-controlled
  p-value across many variants (Benjamini–Hochberg).
- **N**: a variant's read count. **W**: the WT read count.
- **B**: number of null resamples (default 10,000).

## The pipeline, stage by stage

**Stage 1 — filter.** Keep full-length, methylated plasmid reads;
attach nucleosome count (`nc`) from the FiberHMM `nuc` track.

**Stage 2 — call variants.** One CIGAR walk per read calls promoter
SNVs and extracts the barcode. Barcodes are clustered by edit distance
and a per-cluster consensus variant is assigned (`--snv-only` drops
indel-bearing reads as error noise).

**Stage 3 — test each variant.** For every variant with ≥ `min_reads`
(default 50) reads that is the *sole* promoter mutation:

1. **Occupancy curves.** Build per-track occupancy along the promoter
   for the variant's reads and for WT.

2. **Option-A reference (the nucleosome-confounder fix).** Don't
   compare to the plain WT average. Recompute the WT occupancy *as if
   WT had the variant's nucleosome-count mix* (age-adjustment analogy:
   standardize WT to the variant's NC profile). Δ = variant occupancy
   − this NC-reweighted WT reference. Now a non-zero Δ reflects the
   mutation, not a difference in nucleosome composition.

3. **Per-variant empirical null (the sampling-scheme fix).** For this
   variant's exact N, draw N WT reads *with replacement*, NC-matched
   to *this variant*, B times; each draw's Δ vs the same reference is
   one null sample. Under "no effect" the variant's real Δ is just one
   more such draw. The variant's full read set is used as-is (no
   sub-sampling) so there is no null-vs-observed scheme mismatch.

4. **Family-wide max-statistic (the multiple-testing fix).** A real
   footprint change is a contiguous **cluster** of same-direction Δ.
   The test statistic is the single largest cluster signal across
   *all tracks* for the variant. Its significance is read against the
   null distribution of *the same "largest-anywhere" quantity* —
   best-of-many vs best-of-many, so taking the maximum doesn't inflate
   significance. One detection rule is used identically for the
   observed data and the null.

5. **Cross-variant FDR.** BH across all variants' p-values → q-values
   (`variant_fdr_q`). A q<0.05 hit means ≤5% expected false discoveries.

**Extra per-variant readouts:**
- **NC-shift**: signed mean-NC change + Wasserstein distance with its
  own empirical p — detects variants that *reposition nucleosomes*
  (the effect Option-A deliberately conditions out of the footprint
  test).
- **MDE (minimum detectable effect)**: the smallest Δ this variant's N
  could have detected — distinguishes a true "no effect" from
  "underpowered".

## The TF binding-motif layer (the biological readout)

Significance tells you *that* a variant changed footprints; the motif
layer turns many such hits into *where the protein-binding DNA elements
are*. Terms:

- **Motif cluster**: a significant footprint-loss (or gain) region that
  is contiguous, **sign-consistent** (≥90% of its positions move the
  same way — not a noisy mix), the right width (5–25 bp, TF-sized), and
  in a TF-sized bin (`TF`/`sub_TF`). One per affected spot in a
  significant variant.
- **Causal-variant cross-check**: for each motif cluster we record how
  far it sits from the SNV that (presumably) caused it (0 = the cluster
  covers the SNV). A local effect is the expected signature of
  disrupting a binding site you sit in.
- **Reference DNA**: with an optional `--reference` FASTA, every cluster
  carries the underlying genomic sequence (a built-in self-check
  confirms the SNV's reference base matches the FASTA, locking the
  coordinate convention).
- **Cross-variant aggregation** (the deliverable): a protein-binding
  element should be hit by **many independent SNVs** tiling it, not one.
  Per reference position we count *distinct* variants (that pass the
  calibrated cross-variant FDR) whose motif cluster covers it — a
  "disruption-density" track. Where density ≥ a threshold (default 2
  independent variants) we call a **motif**: its reference sequence
  plus a per-position/per-base **sensitivity profile** (which base
  substitutions at which positions drive the disruption). Gating on the
  validated FDR means the motif layer inherits Phase 1's exact
  error control — it never invents calls the statistics don't support.

## What was fixed vs the original pipeline

The original test was ~93%/55% false-positive under a no-effect
control. Three flaws were corrected: (3) null vs observed used
mismatched sampling and an unconditioned WT baseline; (1,2) the
"best of 4 bins × many clusters" was reported uncorrected and the
null/observed cluster rules differed with a selection-biased
reference. After the fixes, a properly-powered WT-vs-WT control
(pooled n=616 independent no-effect pseudo-variants, B=10,000) gives
**0 false discoveries at q<0.05 and q<0.10, p-values ≈ Uniform** —
exact FDR control.

## What is and isn't trustworthy

- **Trustworthy (default config):** effect sizes (Δ, occupancy,
  cluster boundaries), per-variant significance and cross-variant FDR,
  the NC-shift and MDE readouts.
- **Opt-in / not yet calibrated:** `--null-stratify` (a speed
  optimization) is ~2× anti-conservative and OFF by default; do not
  use it for reported significance until its per-member-reference fix
  lands.
- **Detection floor is physical:** the smallest detectable Δ ≈
  `z·√(p(1−p)/N)` — set by molecules sequenced, not computation. See
  `docs/sensitivity_and_power.md`.
- **Per dataset:** re-run the disjoint WT-vs-WT sweep
  (`tools/calibration_sweep/`) to certify FDR for that dataset's
  occupancy structure.

# What this dataset can and cannot tell you

Read this before starting a new question. Numbers are from `config/analysis_v2/rules.json` and the results
named in `FINDINGS.md`.

## Power

| | HepG2 run 3 (main) | HEK293T run 1 |
|---|---|---|
| SNVs analysable | 633 (all designed) | 197 in the analysis table (119 with ≥ 3 clones); motif-level only |
| Clones per SNV | median 775 | median 4 (78 SNVs with 1–2) |
| Molecules | 1.35 M SNV; 1.24 M WT (311 k in the 25% WT clone sample used for tests) | 145 k SNV, 163 k WT |
| MDE, one adenine's inaccessible fraction (E1) | 2.5 points (10–90%: 2.0–3.1; Bonferroni 3.9) on a WT level of ~68% | 3.6 points (2.4–5.7) |
| MDE, 20–40 bp call coverage at one position | 0.72 points on a WT level of ~0.75% | 2.3 points |
| MDE, 10–19 / 41–80 / 81+ bp coverage | 0.60 / 0.65 / 3.0 points (WT ~0.50 / 0.64 / 55%) | – |

MDE = minimum detectable effect at α 0.05 per feature, 80% power (`a9_statistics/*/mde.parquet`).

What that means in practice:
- **Per-adenine accessibility (E1) is the strong endpoint.** A few-point change at a single adenine is
  detectable for most HepG2 SNVs, and averages over many adenines (the element-wide shift) are far more
  precise.
- **Per-position TF-size call coverage is weak.** The MDE is about the WT level itself: a single SNV can
  only show that a TF-size footprint is gone, not that it is 30% weaker. Footprint-level occupancy (F01–F24,
  pooling the calls in a footprint's box) is somewhat better.
- **Pooling SNVs helps**: motif dose-response (A11) and within-position regressions pool many SNVs.
- **HEK293T supports WT landscape and motif-pooled tests only.** Its library is bottlenecked (one clone
  per SNV is typical), so per-SNV claims are not possible.
- The detection floor (A8 spike-ins): a perfectly protected 10-bp stretch is found in ~95% of M molecules
  at 160 of 211 positions; 15–20 bp at the rest (harder at oligo 56–65, 78–90, 120–129, 135–141, 149–152,
  188–199, 205–211).

## Limitations

### No true controls yet (the biggest one)
- **No naked-DNA or no-Hia5 controls on this library** (Phase B1). The caller's emissions come from in-cell
  M (high labeling) and U (no labeling) molecules, and M itself carries footprints (A6). Therefore:
  - footprint calls are "relative to in-cell M", never "above naked DNA";
  - **no absolute occupancy** can be reported;
  - **Hia5 sequence preference cannot be separated from protein** near an SNV. Every effect within 10 bp
    of its SNV (all changes > 15 points are) is "possible Hia5 sequence effect until B1".
- **No time course** (B2): occupancy vs kinetics (how fast a site fills in) cannot be separated.

### Sequence sets what you can see
- Each A/T pair carries one methylatable adenine; G/C pairs carry none. A runs are nearly invisible even on
  naked DNA (oligo 59–65, 114–130 region). "Inaccessible" there mostly means "unmethylatable".
- The SRE-1 core (CCCC at 131–134) has no A/T: it is functionally strong but cannot be footprinted directly.
- Footprint lengths are stretched by neighbouring unmarkable stretches (60–86 is really ~66–81).

### Variants change the evidence, not only the biology
- The plasmid-LLR caller drops the edited base for 303 SNVs (all A/T→G/C and A↔T) and ignores gained
  adenines. Always compare a variant with WT masked the same way (A7 mask primary; ±5 bp sensitivity).
- The A7 check failed: the caller artifact explains only ~15–29% of the old "A/T-loss" offset, so
  "drops its base" SNVs are analysed as their own group and not interpreted alone.

### Caller and labeling
- Low-labeled molecules (0.01–0.2) are called almost everywhere (82–86% at any spot): call endpoints use
  labeling ≥ 0.2. Labeling is a continuum and must be standardised (strata 0.05 wide) in every comparison.
- Size bins 10–19/20–40/41–80/81+ are arbitrary; called lengths are 0–4 bp longer than true. Use the
  data-defined footprints F01–F24 or continuous length.
- Calls at weakly marked positions can be chance: use the confident-call flag (P_chance < 10⁻³) or the
  footprint's `confident_share`. F11 (oligo 117–129) is only 3% confident.

### Design and replication
- One deep run per cell. HepG2 runs 1–2 (the replicate pair) are shallow and differ from run 3 in 7 ways;
  runs 1 vs 3 differ well beyond sampling. Compare within a run; combine runs only at the effect level.
- HepG2 and HEK293T differ in 8–9 ways (construct orientation, m6A caller, Hia5 dose...): never pool or
  overlay them; compare conclusions, not numbers.
- Activity is from Kircher 2019, a different experiment: chromatin–function links are correlations across
  SNVs. Use within-position models (fixed effect per position) to avoid motif clustering.
- Designed 1-bp indels (1.89 M molecules) and double mutants were never analysed.
- HEK293T genotype/barcode checks (A3) were not done.

### Statistics already used: be careful with re-use
- Several choices were made after seeing data and are labelled post hoc in rules.json: the random-effects
  test and the ≥ 5 expected-events rule (after WT-only calibration), the HEK clone-mean test, the
  empirical-null sensitivity, the element-wide shift robustness checks, and all footprint follow-ups
  (60–86, 114–138, fp_features, confidence). Treat new hypotheses suggested by these as exploratory.
- 633 SNVs × 101 adenines × several endpoints have been looked at: new tests on the same data need their
  own multiplicity control and, ideally, a prespecified setting written to `rules.json` first.

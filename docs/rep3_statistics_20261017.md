# Replicate-3 barcode statistics — 2026-10-17

The first statistical implementation completed on the full joined dataset.
**The candidate Welch method did not pass its sampling-calibration gate.
No variant effect is released as statistically supported.** The barcode arrays,
effect estimates, profiles and calibration diagnostics are ready for subsequent
method comparison and descriptive visualization.

Astra specified the design and reviewed the implementation. A `gpt-6-sol`
agent implemented the aggregation/inference engine, workflow and focused tests;
Astra implemented independent calibration and full-data validation. The
[specification](rep3_statistics_design.md) fixes the estimand, endpoints,
thresholds and interpretation; [usage](statistics_usage.md) documents execution.

## Data and outputs

The aggregation scanned all 8,715,827 joined molecules once. Counts reconcile
exactly with the previous footprint-join audit.

| Quantity | Count |
|---|---:|
| Analysis-eligible molecules, including unsupported barcodes | 2,594,368 |
| Eligible molecules with supported barcodes | 2,102,906 |
| Supported barcode groups | 568,354 |
| WT barcode groups / molecules | 268,283 / 985,264 |
| SNV barcode groups / molecules | 300,071 / 1,117,642 |
| SNVs passing 50-barcode and 100-molecule coverage gates | 630 of 633 |
| Fixed SNV × endpoint hypotheses | 13,926 |
| Comparisons meeting the candidate test's data gates | 13,745 |
| Sparse endpoint comparisons | 115 |
| Comparisons with insufficient molecule coverage | 66 |
| Descriptive base-profile rows | 267,548 |
| Sensitivity rows | 41,778 |

Primary effects give each observed supported barcode equal weight, averaging
all its eligible molecules first. There are 21 fixed 20-bp TF-coverage windows
and one MSP-coverage endpoint across the complete 211-bp mutable element.
These are covered fractions, not binary occupancy of an identified binding
site. Confidence intervals are pointwise; candidate BY adjustment includes
all 13,926 members, including untestable members internally as p=1.

The candidate calculations yield 330 rows with q ≤ 0.05 across 174 variants
(325 TF-window rows and five promoter-MSP rows). **These are diagnostic counts,
not discoveries:** every `caller_effect_supported` and `statistical_supported`
field is zero, and release status explicitly withholds candidate inference.
The existence of untestable rows likewise does not establish absence of an
effect, particularly for sparse or apparently complete losses.

All outputs are under `results/rep3_statistics_20261017/`:

- `aggregates/`: reusable barcode arrays, barcode/genotype metadata, molecule
  totals, all-eligible coverage arrays and provenance audit.
- `inference/effects.tsv.gz`: fixed-endpoint estimates, uncertainty, candidate
  p/q values, coverage gates, mutation overlap and A/T opportunity-change flags.
- `inference/profiles.tsv.gz`: TF/MSP curves with equal-barcode, supported-cohort
  molecule-weighted and all-eligible molecule-weighted means.
- `inference/sensitivities.tsv.gz`: weighting comparisons and mutation ±2-bp
  masking, applied to the same reference positions in WT and variant.
- `calibration/report.json` and `calibration/iterations.tsv.gz`: every scenario,
  split, testability status and operational gate result.
- `inference/release_status.json` and `summary.json`: machine-readable release
  decision, independently checked counts and diagnostics.

The existing barcode ledger and joined footprint table were not changed.

## Calibration findings

Each check used 200 iterations, with a seed fixed before computing variant
effects. WT A/A draws compare complete barcode families with their disjoint WT
complement. Simulations cover barcode heterogeneity, unequal/informative read
depth, rare/common outcomes and known positive/negative effects.

| Check family | Passed strata | Failed strata | Insufficiently testable |
|---|---:|---:|---:|
| WT A/A: four barcode counts × 22 endpoints | 29 | 40 | 19 |
| Hierarchical simulations: 24 null + three alternative scenarios | 19 | 4 | 4 |

For an assessed null stratum, the preset gate requires the upper exact 95%
binomial bound on its unconditional rejection rate to be ≤10% at nominal
alpha 5%. Simulation CI coverage additionally requires its lower exact 95%
bound to be ≥90%. Any assessed failure blocks release. Unresolved null
simulation regimes also prevent a global pass; the current engine cannot
assign real comparisons to separately validated regimes.

The strongest WT A/A failure was the TF window at relative [140,160), reference
[1374,1394), with 400 pseudo-variant barcodes: 26/200 null rejections, **13%**,
exact 95% interval **8.7–18.5%**. Another TF window, relative [120,140), rejected
23/200 at 100 barcodes. Across all WT A/A TF endpoints and sizes there were
650 negative versus 150 positive rejections. These counts are correlated
across overlapping windows and are not independent evidence counts, but the
directional imbalance warrants investigation of the approximation's lower tail.

Some failures are less decisive. For example, promoter MSP at 400 barcodes had
13/200 rejections (6.5%, interval 3.5–10.9%), which fails the conservative gate
without establishing that its true error rate exceeds 5%. Three null simulation
failures had empirical rejection rates of 6–7%. The alternative simulation
failure covered the true +0.10 effect in 183/200 intervals (91.5%). Thus the
failed release combines real warning signals with limited Monte Carlo precision
and an intentionally strict all-strata gate. It does not mean every endpoint
has a demonstrated inflated error rate.

These checks assess sampling behavior within an empirical WT population or a
specified simulation model. They do not validate barcode independence across
cells, biological reproducibility, sequence-dependent footprint detectability,
or the very small p-value tails needed for this large hypothesis family.

## Verification and resources

Nine focused tests passed, including independent SciPy Welch/BY comparisons,
barcode/complement moment equivalence, sparse and constant outcomes, count
conservation, masking/export behavior, calibration gating and stale provenance
rejection. The Snakemake dry run and complete calibration→inference workflow
passed. A final dry run confirms no work remains for these outputs.

Full-data validation independently checked count bounds, cohort/barcode totals,
WT base means, unique table keys, fixed family size, missing uncertainty for
untestable comparisons, and BY values against SciPy. Ten additional real-WT
contrasts across 22 endpoints agreed with direct SciPy unequal-variance tests;
the maximum absolute p-value difference was 3.90×10⁻¹⁴. This supports the
numerical implementation, not the scientific validity of the approximation.

Each stage used one CPU within existing Slurm job 41830288. Aggregation took
3 minutes 3 seconds with approximately 1.15 GiB maximum resident memory;
calibration plus inference took 39 seconds, with approximately 1.30 GiB peak
resident memory for the timed workflow. Rules declare 8 GiB. Logs, version
information, source/config snapshots and checks are preserved in
`results/rep3_statistics_20261017_run/`.

## Next statistical decision

Keep these frozen results as the baseline. Benchmark barcode-level resampling
methods, including studentized bootstrap/permutation contrasts, against the
same difficult outcome regimes, with separate development and validation WT
splits. Explicitly test apparent losses, rare calls and unequal variances;
examine whether a change in the estimator or a stricter testability domain is
needed. Increase simulation precision under a new recorded protocol rather
than retrying this seed or weakening its thresholds until it passes. Assess
small-p behavior before treating adjusted values as reliable discoveries.

Descriptive profiles and effect estimates can already feed the dashboard, with
the release-status flag controlling significance labels. Binding-site occupancy,
motif inference and functional activator/repressor interpretations remain later
steps. Sequence-aware caller calibration and biological replication remain
separate requirements even after a sampling method passes.

# Expanded replicate-3 calibration — 2026-10-18

The expanded benchmark is complete. Ordinary Welch testing remains unsuitable
for sparse footprint endpoints. Studentized permutation testing substantially
reduces that problem, but this run does **not** establish a generally calibrated
variant test with useful, reliable effect intervals. Small changes remain hard
to detect. No variant significance labels were released.

## What was run

The frozen [design](rep3_heldout_calibration_design.md) and
[workflow instructions](heldout_calibration_usage.md) describe the implementation.
The analysis retained **258,283 WT barcodes and 948,189 molecules**, excluding
the 10,000 barcodes from the recent footprint-quality audit. These are disjoint
from that audit sample, but not globally unused: previous full-data calibration
drew from the entire WT pool. They all come from one biological preparation.

The builder streamed 8,715,827 joined rows and reproduced every retained
barcode's existing coverage endpoints to a maximum absolute error of
2.22e-16. Each molecule contributes a presence indicator and covered-base
fraction before averaging within barcode; barcodes then receive equal weight.
Only molecules meeting the existing eligibility criteria and spanning the
entire 211-bp element enter this table. This is not a new filtering experiment.

There were 4,000 random WT splits: 1,000 each with 50, 100, 400, or 1,000
pseudo-variant barcodes, compared with the entire remaining WT pool. Each
split was evaluated at 43 outcomes: 21 TF-window presence outcomes, 21
TF-window coverage outcomes, and whole-element MSP coverage. The resulting
172,000 endpoint contrasts are correlated across windows; they are not 172,000
independent biological experiments. Another 10,000 simulated contrasts covered
32 equal-mean scenarios and 18 alternatives, with 200 draws per scenario.

Results: [complete summary table](../results/rep3_heldout_calibration_20261018/calibration/summary.tsv),
[calibration audit](../results/rep3_heldout_calibration_20261018/calibration/audit.json),
and [power/interval table](../results/rep3_heldout_calibration_20261018/plots/effect_and_interval_summary.tsv).

## WT against WT: sparse footprints remain the main difficulty

At a nominal 5% threshold, the observed rejection-rate ranges across the 43
outcomes were:

| Pseudo-variant barcodes | Welch | Studentized permutation | Welch 95% interval coverage |
|---:|---:|---:|---:|
| 50 | 4.6–26.4% | 0.0–5.9% | 73.6–95.4% |
| 100 | 4.9–19.6% | 3.4–7.3% | 80.4–95.1% |
| 400 | 4.1–10.2% | 3.9–8.0% | 89.8–95.9% |
| 1,000 | 4.3–7.7% | 3.6–6.7% | 92.3–95.7% |

The worst Welch case was coverage in element-relative window [180,200) with
50 barcodes: 264/1,000 rejections, or 26.4% (binomial 95% interval 23.7–29.2%).
Of those, 263 were apparent losses. A typical draw had only two barcodes with
nonzero footprint coverage. The small group can miss most rare footprints and
estimate very little variation, while the huge WT reference has a precise
nonzero mean. That makes an ordinary standardized difference look much more
certain than it is. More WT reads do not fix sparse support in the variant arm.

Permutation testing improved this greatly, but the most elevated observed
rate was still 8.0% (80/1,000; 95% interval 6.4–9.9%) for presence in
[100,120) at 400 barcodes. Conversely, several sparse endpoints at 50 barcodes
had zero rejections. These extremes must not be treated as a blanket pass or
as proof of universal failure: endpoints are correlated, many rates were
examined, and each size uses a shared finite reference bank. The reported
binomial intervals condition on both the WT pool and that bank and omit bank
Monte Carlo uncertainty. A fresh-bank check is needed before interpreting an
individual outlier as persistent miscalibration.

Whole-element accessibility was less problematic: permutation rejection rates
were 4.5%, 5.7%, 4.9%, and 5.3% as barcode count increased. Welch interval
coverage for that endpoint was 94.1–95.4%. This is encouraging for the broad
accessibility outcome, but does not validate sparse TF endpoints or biological
generalization.

![WT coverage calibration](../results/rep3_heldout_calibration_20261018/plots/wt_coverage_rates.png)

[Presence calibration figure](../results/rep3_heldout_calibration_20261018/plots/wt_presence_rates.png).
Cells show percentages; the colorbar uses proportions.

## Unequal barcode variability and depth

The synthetic nulls give both groups the same mean, while changing barcode
variability and the relationship between depth and occupancy. They test a
harder question than splitting an exchangeable WT pool.

At a mean of 0.2% and 100 variant barcodes, Welch falsely rejected 54–68.5%
of simulations. The worst case was 137/200 (95% interval 61.6–74.9%). The
median number of nonzero variant barcodes was zero. Permutation testing
rejected none of these four scenarios, demonstrating its conservatism in this
very sparse regime, rather than demonstrating adequate sensitivity.

Across all 32 null scenarios, the highest permutation rejection rate was
9.0% (18/200; 95% interval 5.4–13.9%), in two scenarios with mean 10%,
400 variant barcodes, independent depth, and opposite dispersion contrasts.
Two hundred draws per scenario cannot establish accurate tail probabilities,
and these intervals do not correct for inspecting many scenarios. Studentizing
a permutation statistic does not make an unequal-distribution mean test
finite-sample exact; this limitation is part of the frozen design.

![Synthetic null calibration](../results/rep3_heldout_calibration_20261018/plots/synthetic_null_rates.png)

## How large a change can be detected?

These are observed permutation detection percentages at an **uncorrected 5%**
threshold, with 200 simulations per cell and 10,000 control barcodes. They are
conditional on the specified beta/binomial depth model, not measured power
for an actual variant or for the full multiple-testing family.

| Control → variant mean | 100 variant barcodes | 400 | 1,000 |
|---|---:|---:|---:|
| 0.2% → 0% | 0% | 0% | 62% |
| 2% → 0% | 96% | 100% | 100% |
| 2% → 1.5% | 3.5% | 22.5% | 40% |
| 2% → 2.5% | 1% | 5.5% | 21.5% |
| 10% → 8% | 17.5% | 47% | 81.5% |
| 10% → 12% | 6.5% | 29% | 69% |

A complete loss can be difficult to establish when the original footprint is
extremely rare. The low detection rates for a 0.5-percentage-point change also
show why a nonsignificant result cannot be interpreted as no regulatory
effect. Welch's apparent higher sensitivity is not a solution given its null
failures. Gains and losses need separate power checks, as these results show.

## Effect intervals are still unresolved

The nominal 95% Welch intervals covered truth only 73.6% of the time in the
worst WT case and 31.5% in the worst synthetic null. Even in the alternative
scenarios, coverage ranged from 86.5% to 96%. A better permutation p-value
does not repair those intervals: they are Welch intervals, not intervals
obtained by inverting a permutation test.

The conservative bounded-data Hoeffding intervals covered truth in every draw
here, but were too wide for small effects. Their total widths were about
42.4, 30.2, 15.4, and 10.0 percentage points in the four WT group sizes.
With 10,000 simulated controls, widths were 32.6, 17.8, and 12.3 percentage
points for 100, 400, and 1,000 variant barcodes. This is a useful conservative
reference, not an automatically adopted analysis. Its independence assumptions
also do not establish independence between plasmids sharing a cell.

## Caller provenance

The accompanying [caller audit](rep3_caller_provenance_20261018.md) found
retained m6A probabilities restricted to ML bytes 240–255 in the first 1,000
records of each available raw-ft, aligned, and footprint BAM. The restriction
therefore precedes alignment and FiberHMM in the available files. Changing a
downstream cutoff from 128 to 204 cannot assess lower-confidence calls in
these sampled records. The original pre-ft BAM path is unavailable; the exact
upstream filtering setting remains unresolved.

The installed FiberHMM code can use zero nucleosome AQ as an unverified-score
sentinel. That is a plausible explanation for the observed zero scores, but
does not prove the historical caller version or flags. Those zeros cannot be
used as quantitative evidence of low nucleosome confidence.

## Decision and next implementation

Keep descriptive barcode-weighted effects available and continue withholding
variant significance labels. Do not optimize methylation or chromatinization
filters to make this benchmark pass; that would change the biological target
and select on the observed result.

The next development target is **calibrated effect estimation for sparse
barcode outcomes**, with explicit support limits:

1. Summarize actual per-variant barcode depth and nonzero footprint support,
   without looking at significance, to map variants to the regimes tested here.
2. Freeze the next model and validation criteria before fitting variants.
   Compare a barcode-level presence/count model that allows unequal dispersion
   with a conservative bounded-mean interval baseline. Retain coverage as a
   distinct fractional outcome: covered bases are correlated and must not be
   treated as independent binomial trials. A model-based interval must earn its
   use through coverage and misspecification checks; a model alone does not
   overcome missing information in a rare endpoint.
3. Recheck the elevated WT permutation cases with independent reference banks
   and extend null simulation precision before declaring calibration. Include
   all-zero groups, depth associated with occupancy, and unequal dispersion.
   Once calibration holds, assess multiplicity and attainable power for the
   prespecified testing family. Current p-value resolutions (0.0001 for WT,
   0.001 for simulations) do not establish the required far-tail behavior.
4. Expose effects, support, uncertainty, and an explicit “insufficient
   information” state in the dashboard. Keep within-preparation evidence
   separate from later biological replication and functional validation.

This benchmark is now development data for that next method. Any revised
method needs fresh random draws and prespecified new stress scenarios; merely
fitting this table is not independent validation. The earlier probability/kinetic
data remain useful for the separate sequence-dependent measurement-bias check,
but their absence does not block statistical development on current outcomes.

## Verification and execution

All 15 targeted tests passed. An independent SciPy/direct-complement check
replayed eight pilot WT splits across all 43 outcomes and 76 reference draws.
Maximum absolute difference was 4.45e-12 for Welch p-values and zero for
permutation p-values; effect and interval differences were below 4.5e-14.
See the [validation audit](../results/rep3_heldout_calibration_20261018_run/independent_validation/audit.json).

The outcome build took 3:21 (peak resident memory about 1.45 GiB), and full
calibration took 6:48 (about 0.75 GiB), each using one CPU inside Slurm job
41830288. The three figures were inspected. A final Snakemake dry run found
all requested outputs up to date. These runs used the equivalent direct CLI
commands, so Snakemake correctly reports missing execution metadata for its
three rules; each stage instead retains its own source audit. No claim of an
end-to-end Snakemake execution is made.

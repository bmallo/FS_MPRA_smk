# Library-wide per-site m6A with same-reaction proxy controls

This stage implements analysis-plan steps 2–4 of the [independent review](INDEPENDENT_REVIEW_20261020.md) on the existing replicate-3 data:
- covariate-adjusted MSP;
- per-site retained m6A as a dense local readout for **all 633 SNVs**;
- an observability check using molecules that apparently never met Hia5, or were labeled close to saturation.

It needs only the inputs already on Klone (heterogeneity feature table, catalog and FiberHMM BAM). The raw-kinetics BAM is **not** needed for this stage.

## Design, fixed before running

**Molecule classes use distal methylation only.** Distal means the plasmid outside the element ±100 bp and the barcode ±25 bp, as in the heterogeneity stage. A substitution in the 211-bp element therefore cannot move its own molecules between classes. The run reports this check as `class_dispersion` in `summary.json` (between-SNV variance of class fractions over the binomial expectation).

| Class | Default rule | What it estimates |
|---|---|---|
| **U** | distal m6A ≤ 0.01 | Sequence-dependent *false-positive* calls (plasmids that apparently never met Hia5) |
| **M** | distal m6A ≥ 0.60 | Sequence-dependent *labeling/calling* near saturation (presumed unchromatinized plasmids) |
| **C** | in between | The chromatin signal, directly standardized over 10 WT distal-labeling strata |

U and M are proxies from the same reaction, not separately prepared controls. M may include open nuclear plasmids, and U may include plasmids that are protected rather than unexposed. Thresholds are set in `config/ldlr_site_methylation.json` before any variant result is viewed. If you used specific cutoffs before, put them there now and record why; two sensitivity pairs (0.02/0.70 and 0.005/0.50) run automatically.

**Readouts.** For every SNV, using only A/T sites that are unchanged between the SNV and WT (the substituted base is always excluded):

- Windows by distance from the substitution: `near` (1–2 bp), `local` (3–10), `flank` (11–30), `element` (all, excluding ±2).
- The fixed regions in the config: WT_P1–P4, adjacent 143–150, P1 core, and Repeats 1–3 (excluding ±2 bp).
- Measures per window:
  - `U` and `M`: retained-call fraction in each proxy class.
  - `C_raw`: pooled fraction in C.
  - `C_std`: C standardized over WT labeling strata.
  - `ALL`: all classes, comparable to earlier reports.
  - `A_norm = (C_std − U)/(M − U)`: normalized accessibility. It cancels any sequence effect that scales labeling/calling, as tested in `tests/test_site_methylation.py`.
- Element scalars: MSP, 85–250-bp protection and 10–60-bp TF coverage (C-standardized and pooled), plus U/M class fractions.

**Uncertainty.** A barcode-cluster Poisson bootstrap (200 replicates) of each SNV, with an added reference-size variance term (≈0.3% for real WT). The WT pseudo-variant check draws whole barcode clusters, matched to the median SNV cohort size, against the WT remainder. It uses the same windows at random substituted positions and records z-scores for every endpoint. Molecules without a barcode are their own cluster. The primary cohort is direct genotype (`footprint_analysis_eligible`, no barcode requirement); a supported-barcode run is included as a sensitivity.

**Library kernel.** `kernel.tsv` gives the mean per-site change against WT as a function of offset from the substitution (±40 bp), per class and per A/T-change class. It shows how far sequence-dependent observability reaches, and so whether ±2-bp masking was ever sufficient.

## Interpretation rules, fixed now

1. **An endpoint's z/p/q are usable** only if its WT pseudo-variant `sd_z` is within 0.85–1.20 and the upper exact 95% bound of `rate_abs_z_gt_1_96` is ≤ 0.10. `U` endpoints are a nuisance correction inside `A_norm`, not results; rare false-positive calls make their z-values slightly anti-conservative in simulation.

   Variant-level claims use `q_bh_primary_family`: one BH family across all SNVs for the four prespecified `primary_endpoints` (`local:A_norm`, `flank:A_norm`, `msp_C_std`, `nuc_85_250_C_std`). `q_bh_within_endpoint` is exploratory only. With ~90 endpoint families, about 5% of families show a false q < 0.05 under a global null; the null simulation below shows exactly that.
2. **Observability reach.** An M-class kernel with |mean| > 2 SE at offsets beyond 2 bp means ±2-bp masking was insufficient. `A_norm` (or wider masking) is then the primary local readout.
3. **Position 125 (chr19:11089389).** The signal is classed as predominantly observability if `WT_P3:M` changes in the same direction as `WT_P3:C_std` with at least half its magnitude for at least two of the three substitutions. It is classed as supporting a chromatin change only if `WT_P3:A_norm` is positive with z > 3 in all three, under all three threshold pairs.
4. **P1 core, Repeats 1–3.** The review's predictions:
   - P1 core: in-region SNVs lose protection, so local methylation rises (`P1_core:A_norm` > 0) without an M-class change.
   - Sp1 repeats: accessibility falls (`msp_C_std` < 0, `Repeat*:A_norm` < 0) and `nuc_85_250_C_std` > 0.
   - SRE-1: `msp_C_std` ≈ 0.

   Report agreement and disagreement equally.
5. **Library association.** `activity_association.tsv` gives Pearson r with published activity and the exact circular position-shift p for every endpoint. The prediction from the review is r(`msp_C_std`) > r(`msp_ALL`) > 0 and r(`nuc_85_250_C_std`) < 0.

Do not tune thresholds, windows or cohorts against these outcomes; changes go in a new config with a stated reason, and the earlier run is kept.

## Run on Klone

From the repository root inside a Slurm allocation, with one CPU:

- **Extraction:** about 1–2 h (≈0.5 ms of decoding per molecule plus the 8.7M-record BAM scan), ≤ 8 GB.
- **Each effects run:** minutes. It took 0.8 min with 5.6 GiB peak RSS (mapped pages included) on a 1.8M-molecule synthetic library at production settings; request 16 GB.

Start with a smoke test of the extraction:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
export MPLCONFIGDIR=/tmp/fsmpra-matplotlib
PY=/mmfs1/home/bmallo/miniconda3/bin/python
$PY -m unittest tests.test_site_methylation
$PY -m fsmpra.site_extract --output results/site_smoke/extract --limit-molecules 20000
$PY -m fsmpra.site_effects --extract results/site_smoke/extract --output results/site_smoke/effects --aa-draws 20
```

The smoke output is marked incomplete and covers only the first 20,000 molecules in BAM order, so few SNVs have enough data. Check that it finishes, then delete it. The full workflow is then:

```bash
/mmfs1/home/bmallo/miniconda3/envs/snakemake/bin/snakemake \
  -s workflow/site_methylation.smk --cores 1 --resources mem_mb=16384
```

or directly:

```bash
$PY -m fsmpra.site_extract --output results/rep3_site_methylation_20261021/extract
$PY -m fsmpra.site_effects --extract results/rep3_site_methylation_20261021/extract \
  --output results/rep3_site_methylation_20261021/effects
```

Outputs refuse to overwrite. The extraction (≈ 2.3 GB of `.npy` arrays plus metadata) is reused by every effects run, so threshold or cohort sensitivities take no further BAM passes.

## Outputs

`extract/`:
- `scores.npy`: retained ML byte; 0 = no retained call.
- `known.npy`: A/T opportunity.
- `qualities.npy`, `coverage.npy`: coverage bit flags 1 TF 10–60, 2 MSP, 4 nucleosome 85–250, 8 nucleosome >250.
- `meta.npz`: genotype, cluster, cohort flags, distal/whole methylation.
- `molecules.tsv.gz`: read IDs and row numbers.
- `audit.json`.

`effects/`:
- `effects.tsv.gz`: one row per SNV × endpoint, with molecule counts per class, variant/WT values, Δ, SE, z, p, BH q within endpoint, and published activity.
- `calibration_summary.tsv`, `calibration_draws.tsv.gz`: WT pseudo-variants.
- `kernel.tsv`, `kernel.png`: offset kernel.
- `activity_association.tsv`.
- `observability_regression.tsv`: across SNVs, Δ in C regressed on Δ in M per window.
- `distal_m6a_histogram.tsv/.png`: the thresholds on the WT/SNV distal distribution.
- `class_counts.tsv`, `observability_scatter.png`, `calibration_z.png`.
- `summary.json`.

## Validation before real data

`tests/test_site_methylation.py` covers:
- **Extraction** end to end on an indexed synthetic BAM: mixed MM channels, a reverse-strand read, MA tracks, a shared barcode, an unbarcoded singleton, a missing distal value and a skipped feature row.
- **Analysis** of a simulated library with one injected multiplicative observability artifact and one injected chromatin change.

In that simulation:

| Injected effect | Readout | Δ | z |
|---|---|---:|---:|
| Observability artifact | raw C | +2.7 pp | 8.4 |
| Observability artifact | `A_norm` | −0.3 pp | −0.7 |
| Chromatin change | `A_norm` | +8.5 pp | 12.5 |
| Chromatin change | M | −1.9 pp | −1.95 |

So a raw-call analysis would report the artifact as an accessibility gain, and normalization removes it. The chromatin change survives normalization and is essentially absent in M.

A production-scale null simulation (633 SNVs, 1.8M molecules, 200 bootstrap replicates, 200 WT pseudo-variants) gave:
- **WT pseudo-variants (non-U endpoints):** median z SD 1.007 (max 1.09); median |z| > 1.96 rate 5.0%.
- **SNVs:** median z SD 0.99.
- **Primary family:** no BH rejections among 1,899 tests (minimum q 0.45).
- **Class fractions:** variance ratios 0.94/1.01.

Real-data calibration still has to be read from each run's own `calibration_summary.tsv`.

## Limitations

- **Retained calls only.** These are stored MM/ML calls, all with ML 240–255 in this BAM, not kinetics. Calls below the upstream threshold are unobserved, not zero.
- **Proxies, not prepared controls.** The proxy controls cannot separate compartment from labeling and are not a substitute for an in-vitro-methylated library.
- **Normalization assumes a multiplicative effect.** `A_norm` assumes sequence effects scale labeling/calling. Additive or chromatin-state-specific caller effects are not fully removed.
- **WT is not resampled.** The WT reference contributes only the reference-size variance term.
- **Single preparation.** All results are within one preparation; replication still requires a new one.

## Raw kinetics (later)

Transfer the Globus raw-kinetics BAM to Klone for a later stage; it is not needed here. With kinetics, m6A can be re-called with full per-site probabilities instead of the retained ≥240 calls. The U and M molecules then become labeled training/validation sets for a sequence-context calling model, and the same comparisons can be repeated threshold-free.

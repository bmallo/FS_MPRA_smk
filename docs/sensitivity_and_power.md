# Detection sensitivity & power — FS-MPRA Stage 3

How much sequencing (physical) and how much null computation
(computational) you need to detect an occupancy change of a given
size. Use this for experimental design and for interpreting whether a
"no effect" call is real or underpowered.

> Scope note: the pipeline does **not** apply a z-threshold — it uses
> an empirical family-wide max-statistic null. The closed-form formulas
> here are a *planning approximation*. For the realized, dataset-exact
> sensitivity of a finished run, read the per-variant `mde` track /
> `mde_median` the pipeline emits, and trust the WT-vs-WT calibration
> over any back-of-envelope number.

## 1. The model

At each promoter position, **occupancy** = the fraction of molecules
(reads) whose footprint (in a given track: `tf`, `nuc`, …) covers that
position. For a variant with **N** reads it is a mean of N Bernoulli(p)
draws, so its standard error is

```
SE(occupancy) = sqrt( p (1 - p) / N )          p = local WT occupancy
```

The effect is **Δ = variant occupancy − Option-A WT reference**. The WT
reference is averaged over the WT pool of size **W**; when **W ≫ N** its
error is negligible, so

```
SE(Δ) ≈ sqrt( p (1 - p) / N )      (WT side negligible iff W ≫ N)
```

**This is the hard floor.** No amount of null computation lowers it —
it is set by how many *molecules* you sequenced. p(1−p) is maximal at
p=0.5 (SE = 0.5/√N) and smaller toward p→0 or 1.

## 2. Variant coverage (N) → smallest detectable Δ

Minimum detectable |Δ| ≈ `z · sqrt(p(1−p)/N)`, so

```
N  ≈  ( z · sqrt(p(1-p)) / Δ )²
```

`z` is the effective significance multiplier. Two references:
- **z ≈ 1.96** — nominal single-position (optimistic lower bound).
- **z ≈ 3.5** — practical *family-wise* value: the test maxes over
  4 tracks × ~promoter positions and then BH across variants, so the
  effective per-effect bar is ≈ p < 10⁻³–10⁻⁴ (z ≈ 3.3–3.9). Use 3.5
  as the planning rule of thumb.

**Reads needed, at the worst case p = 0.5 (SE = 0.5/√N):**

| Target Δ | N (z=1.96, nominal) | N (z=3.5, family-wise) |
|---|---|---|
| 0.5 % | ~38,000 | ~122,000 |
| 1 %   | ~9,600  | ~30,600 |
| 2 %   | ~2,400  | ~7,700 |
| 5 %   | ~390    | ~1,200 |
| 10 %  | ~100    | ~310 |
| 20 %  | ~25     | ~80 |
| 40 % (full footprint knockout) | ~6 | ~20 |

Lower baseline occupancy needs fewer reads — multiply N by
p(1−p)/0.25 (e.g. at p=0.1 or 0.9, ≈ 0.36× the reads above).

**For this LDLR dataset** (WT≈193k; single-variant N: p10≈200,
p50≈1248, p90≈2866), at p=0.5, z=3.5:

| Coverage | N | Detectable Δ |
|---|---|---|
| Low (p10)    | ~200  | ~12 % |
| Median (p50) | ~1250 | ~5 % |
| High (p90)   | ~2870 | ~3.3 % |

So at typical coverage you reliably catch ≥~5 % absolute occupancy
shifts and any full footprint gain/loss; sub-2 % effects need
thousands of reads per variant.

## 3. WT coverage (W) — two requirements

1. **W ≫ N.** The Option-A self-reference inflates effective variance
   by ≈ N/W. Keep **N/W ≲ 5 %** ⇒ W ≥ ~20 × the largest variant N you
   want to support. With W≈193k and max N≈5k that is ~2.6 % — fine.
   If W is small relative to N the test drifts mildly anti-conservative
   (and the closed form above no longer holds — add the WT term
   `sqrt(p(1−p)/W)` in quadrature).
2. **Per-NC-bin support.** The null is NC-matched to the variant, so
   every nucleosome-count bin a variant occupies needs enough WT reads
   (≳ a few hundred) to estimate the NC-conditional WT occupancy and to
   draw NC-matched samples. Rare-NC variants are the limiting case.

Practical target: **W ≥ 50k (≥100k ideal)**, spanning the variant NC
range. WT also gives you the NC-shift readout for free.

## 4. The null distribution — physical vs computational

| | Sets… | Lever |
|---|---|---|
| **Physical** (W, N molecules) | the *noise floor* — the smallest Δ that is in principle detectable | sequence more |
| **Computational** (B null subsamples, stratification) | the *resolution* of the p-value and the multiple-testing headroom — **not** the noise floor | bigger B / more compute |

- **B (null subsamples, default 10,000).** Empirical p resolves to
  ≈ 1/(B+1). For cross-variant FDR over ~150 variants at q=0.05 you
  need p resolvable to ≲ 1/150 ≈ 0.007 (B ≥ ~2,000 minimum). B=10,000
  gives a 10⁻⁴ floor — ample headroom and needed only to *assign* a
  small p, never to *create* sensitivity. Detecting a small Δ that is
  genuinely there still requires the molecules from §2.
- **Per-iteration cost ∝ N · L · (#tracks·bins)** (L = promoter
  length); memory ∝ B·L·bins (null Δ, float32) + N·L·bins (coverage).
  Stratification reuses one null across similar (N, NC) variants —
  pure compute savings, no statistical cost when strata are tight
  (validated by the WT-vs-WT sweep).

**Takeaway:** to detect smaller Δ, add reads (physical). Increasing B
only buys p-value resolution and multiplicity headroom.

## 5. Why cluster width does *not* lower the Δ floor

A footprint of width *w* is present/absent **coherently per molecule**
— the same Bernoulli at all w positions on a given read. The cluster
statistic Σ|Δ| then has signal ≈ w·δ *and* noise ≈ w·SE (the w
positions move together; no √w averaging). Signal/noise = δ/SE,
**independent of w**. So clustering improves *localization* and
*rejection of incoherent single-position noise*, not raw sensitivity:
the per-position `z·sqrt(p(1−p)/N)` is the footprint detection floor.
Don't expect wide footprints to be detectable at proportionally
smaller Δ. (Fuzzy footprint edges give a small partial-averaging
benefit only.)

## 6. Where you *do* gain power: saturation aggregation

Individually marginal variants combine at the **motif** level. If k
independent SNVs tile one binding site and each independently shifts
occupancy, the motif-level evidence aggregates roughly like k
independent instruments — a much lower effective Δ floor for the
*site* than for any single variant. This is the intended route to
weak/cooperative sites and is the basis of the Phase 2 motif layer:
design for saturation coverage even if per-variant N is modest.

## 7. Secondary / co-occupancy effects need more reads

A variant's *primary* footprint loss is usually large (often a near
knockout, Δ ≈ 0.3–0.5 → easy). A *dependent* protein's secondary loss
is **attenuated** (partial dependency → smaller Δ) and therefore needs
substantially higher N than the primary — budget per-variant coverage
for the secondary effect size, not the primary, if cooperativity is a
goal.

## 8. Design rules of thumb

- Decide the **smallest biologically meaningful Δ**; read N off the
  z=3.5 column of §2 (scale by p(1−p)/0.25 for non-0.5 baselines).
- Full footprint knockouts: N ≥ ~300–500 is plenty.
- Partial / cooperative / secondary effects (Δ ≈ 1–3 %): N ≥ a few
  thousand per variant.
- WT: ≥ 50k (≥100k ideal), spanning the variant NC range, W ≥ 20·N.
- Library: ~1 expected mutation per promoter (maximizes *single*-
  variant molecules, which is what Stage 3 tests); SNV-saturation;
  run Stage 2 `--snv-only` if indels are not designed.
- B = 10,000 (default) is sufficient for the Δ range above.
- Confirm per dataset with the disjoint WT-vs-WT sweep
  (`tools/calibration_sweep/`) — it certifies the realized FDR for
  *that* dataset's occupancy structure; trust it over these formulas.

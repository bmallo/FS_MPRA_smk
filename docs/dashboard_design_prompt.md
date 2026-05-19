# Prompt: FS-MPRA Fiber-seq results dashboard

> Paste everything below into Claude (design) as the brief. It is
> self-contained: it defines the science, the exact data contract, the
> baseline app to modernize, every required view, and the UX bar.

---

## Role & goal

You are designing and implementing a **dark-mode, highly interactive
single-page web dashboard** for exploring the results of a **Fiber-seq
MPRA saturation-mutagenesis pipeline**. A scientist uses it to go from
"which single-nucleotide variants in a promoter change chromatin
protein occupancy?" down to "which transcription factor, at which base,
and does that factor cooperate with another site?". Prioritize: (1)
**dark mode**, (2) **interactivity** (cross-linked, hover/click,
fast filtering), (3) **faithfully visualizing every analysis listed
below**. Aesthetics should be clean, dense-but-legible, scientific
(think a polished genomics browser, not a marketing page).

## Scientific context (plain language — keep tooltips this clear)

- A plasmid promoter is mutated to (nearly) every possible single
  nucleotide. Each sequenced molecule is one independent plasmid
  (single-molecule, no PCR). "Footprints" = stretches a DNA-binding
  protein protected, called per molecule and binned by size:
  **sub_TF (10–19 bp), TF (20–40), PIC (41–80), nuc (nucleosome 81+)**.
- **Occupancy** at a position = fraction of molecules footprinted
  there (0–1). **Δ (delta)** = a variant's occupancy minus the
  wild-type (WT) reference, per position, per bin. Negative Δ = the
  mutation *removed* protein binding (loss); positive = gain.
- The WT reference is **nucleosome-count-reweighted** (an
  "age-adjustment" so Δ reflects the mutation, not a difference in how
  many nucleosomes the molecules carried).
- Significance is a calibrated empirical FDR (`variant_fdr_q`).
  Significant contiguous Δ runs are **clusters**; clusters that are
  sign-consistent, TF-sized, and recur across many independent SNVs
  become **motif calls** (a putative protein-binding element). Each
  motif carries a **per-position × per-base sensitivity profile**:
  which substitution at which position most disrupts occupancy.
- **Co-occupancy**: when a variant disrupts site 1, does a *distal*
  site 2 also change beyond the WT `P(occupied@2 | occupied@1)`
  channel, consistently across many independent SNVs? (cooperativity).
  Calibrated but power-limited on sparse data — every pair has an
  `mde`/`underpowered` flag; **a non-call ≠ "no dependency"**.
- **TF candidates**: each motif's reference DNA scored against
  JASPAR/HOCOMOCO PWMs, *corroborated* by whether the PWM's
  information-rich positions coincide with the motif's
  mutation-sensitive positions (`ic_sensitivity_concordance`).

## Data contract (what the dashboard reads — be exact)

Inputs are produced per sample by the pipeline:

**`{sample}.h5`** (HDF5) groups:
- `metadata` (attrs): `reference_name`, `reference_length`,
  `promoter_start`, `promoter_end`, `analysis_start`, `analysis_end`,
  `analysis_length`, `bin_labels` (e.g. `["sub_TF","TF","PIC","nuc"]`),
  `n_variants_tested`, `n_null_iterations`.
- `wt_occupancy/{bin}`: float array length `analysis_length`
  (Option-A WT reference occupancy per position).
- `ground_truth_nc`: `nc_vals`, `nc_fracs`, `n_reads`.
- `summary` (parallel arrays, one row per tested variant):
  `variant_ids`, `positions`, `ref_bases`, `alt_bases`,
  `change_types`, `n_reads`, `best_cluster_p`, `variant_fdr_q`,
  `nc_delta`, `nc_wasserstein`, `nc_shift_p`, `nc_shift_q`,
  `mde_median`, and per bin `{bin}_max_abs_delta`,
  `{bin}_n_sig_positions`, `{bin}_n_sig_clusters`.
- `clusters` (flat table of all significant clusters):
  `variant_ids`, `bin_labels`, `abs_start`, `abs_end`, `width`,
  `sum_abs_delta`, `max_abs_delta`, `mean_signed_delta`, `direction`
  (`loss`/`gain`), `sign_consistency`, `is_motif` (bool),
  `variant_distance` (cluster–causal-SNV offset; 0 = overlaps),
  `ref_sequence`, `peak_position`.
- `variants/{key}`: per-variant attrs (`variant_id`, `position`,
  `ref_base`, `alt_base`, `n_reads`=`n_nc_matched`, `best_cluster_p`,
  `variant_fdr_q`) and per bin a subgroup with arrays length
  `analysis_length`: `delta_obs`, `variant_occ`, `empirical_p`,
  `q_values`, `z_scores`, `mde`; plus a `significant_clusters`
  subgroup (per-cluster attrs incl. `ref_sequence`,
  `sign_consistency`, `is_motif`, `variant_distance`).
- `motifs`: `density_tracks/{bin}__{loss|gain}` (int array length
  `analysis_length` = # distinct significant variants disrupting each
  position); `calls/m{i}` attrs (`bin`, `direction`, `abs_start`,
  `abs_end`, `width`, `ref_sequence`, `n_variants`, `peak_density`,
  `base_order`="ACGT") + datasets `sensitivity_count` [width×4],
  `sensitivity_mean_signed_delta` [width×4], `contributing_variant_ids`.
- `cooccupancy/pair{i}` attrs: `site1_id`, `site2_id`, `site1_bin`,
  `site2_bin`, `site1_abs`(2), `site2_abs`(2), `n_instruments`,
  `weighted_mean_excess`, `p_two_sided`, `fdr_q`,
  `frac_instruments_consistent`, `mde`, `underpowered` (bool),
  `wt_p2_given_1`, `wt_p2_given_0`, `is_call`, `validated_call`
  (the trustworthy flag), `secondary_artifact`, `sec_max_secondary_freq`;
  dataset `instrument_variant_ids`.

**Companion TSVs** (same fields, easier tabular ingest):
`{sample}_summary.tsv`, `{sample}_motifs.tsv`,
`{sample}_cooccupancy.tsv`, `{sample}_tf_candidates.tsv`
(columns: `motif_idx, bin, abs_start, abs_end, motif_seq, db, tf,
matrix_id, pwm_len, strand, offset, logodds, score_norm,
pwm_ic_mean, ic_sensitivity_concordance`).

**Optional** `--mpra` TSV: external functional scores keyed by
`position` (+ `ref`,`alt`,`pvalue`) to overlay.

> **Data-loading reality (important):** HDF5 cannot be read in a
> browser. Provide a thin Python preprocessing step that exports the
> needed slices to JSON/Parquet the front-end loads, **or** keep a
> Python backend (Dash/Plotly) that reads the HDF5 directly. Either is
> acceptable — state which you chose and why. Coordinates: variant IDs
> are 1-based; HDF5 cluster/site coords are 0-based; show both.

## Baseline to modernize (existing draft)

There is a working but light-themed Dash/Plotly app
(`extras/fs_mpra_browser.py`, ~1400 lines) with 5 tabs to preserve and
improve, NOT discard: **WT Landscape**, **Variant Browser**
(sortable/sliderable per-variant occupancy + Δ + significance),
**Library Heatmap** (variants × positions, selectable bin / metric /
diverging colorscale / clamp range), **MPRA Data** (functional-score
lollipop overlay), **Comparison** (two variants side by side). Keep
their good interactions (hover-linked positions, sort by
position/FDR/effect, colorscale + scale-clamp controls) but reskin to
dark mode and modern layout.

## Required views (organize as tabs or a left-nav)

1. **Library overview** — promoter-wide summary: per-bin
   disruption-density tracks (`motifs/density_tracks`), called motifs
   as highlighted intervals, a volcano (effect vs `variant_fdr_q`),
   counts (variants tested / significant / motifs / co-occ). Clicking
   anything cross-filters the other views.
2. **WT landscape** — Option-A WT occupancy per bin along the
   promoter; optional MPRA overlay.
3. **Variant browser** — step/search variants (sort by position, FDR,
   effect, reads); show WT vs variant occupancy, signed Δ with
   significant positions emphasized, `z`, `empirical_p/q`, the
   per-variant **MDE track** (so "flat" reads as *underpowered* vs
   *truly null*), the causal SNV marked, and that variant's
   significant clusters annotated with extracted `ref_sequence`.
4. **Motif explorer** (new, headline) — list/table of motif calls;
   for a selected motif: its reference sequence, the
   **per-position × per-base sensitivity heatmap** (4×width;
   `sensitivity_mean_signed_delta` colored diverging, size/alpha by
   `sensitivity_count`) as a sequence-logo-like panel, the
   contributing independent variants, and its **TF candidates** ranked
   by `score_norm` with `ic_sensitivity_concordance` shown as the
   discriminating signal (visually separate strong-concordance hits
   from promiscuous sequence-only hits).
5. **Co-occupancy** (new) — a site×site matrix / arc diagram over the
   promoter; cell/arc encodes signed `weighted_mean_excess`, opacity
   by significance, **explicit `underpowered` styling** (e.g.
   hatched/desaturated) and a `validated_call` badge; click a pair to
   show its instruments, `wt_p2_given_1` vs `wt_p2_given_0`, the
   secondary-mutation-screen status, and `mde`. Make the
   "absence ≠ no dependency" caveat unmissable.
6. **Library heatmap** — modernized existing heatmap (bin / metric:
   Δ, −log10 p, z / direction filter / diverging colorscale / clamp).
7. **Comparison** — two (or N) variants overlaid.

## Global UX requirements

- **Dark mode**: deep neutral background (~`#0d1117`/`#111317`),
  high-contrast text, a perceptually-uniform **diverging** scale for
  signed Δ (blue = loss, red = gain) and a sequential scale for
  density; colorblind-safe; never encode meaning by hue alone.
- **Interactive & cross-linked**: a position/variant selected in one
  view filters/scrolls the others; hover tooltips everywhere using the
  plain-language definitions above; debounced search; URL-encoded
  state so a view is shareable.
- **Honest stats surfacing**: always show `variant_fdr_q` (not raw p),
  and `mde`/`underpowered` wherever a null could be misread as
  "nothing there". Distinguish "calibrated significant",
  "underpowered", and "true null" visually.
- **Performance**: promoters are short (~300 bp) but variants reach
  thousands and per-variant arrays are length≈analysis_length; use
  virtualization/decimation so heatmaps and tables stay smooth.
- **Layout**: persistent sample/threshold controls; left-nav or tabs;
  responsive ≥1280px; export-PNG/SVG on every plot; a "what am I
  looking at?" info panel per view.

## Deliverable

A runnable dark-mode dashboard plus: the chosen stack & why, the
HDF5→front-end data step (or Dash backend), a small synthetic/sample
dataset or fixture so it runs without the real ~GB inputs, and a
README. Keep dependencies modest. Ask before inventing data fields not
in the contract above.

## Non-goals

Not re-running analysis or recomputing statistics in the dashboard
(it only visualizes pipeline outputs); not authentication/multi-user;
not editing data.

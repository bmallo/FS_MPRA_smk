# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Start here

1. `docs/PROJECT_HANDOFF_20261020.md` — the current scientific and implementation synthesis. Read it first.
2. `docs/INDEPENDENT_REVIEW_20261020.md` — a re-check of the evidence snapshot and the recommended direction. It says to build on the data foundation, but to stop the per-variant × 20-bp-window significance program in its current form.

Other dated reports in `docs/` record earlier decisions and may propose steps that were later completed. Where they disagree with the handoff, the handoff wins. **No replicate-3 Fiber-seq variant effect has been released as statistically significant.**

## What this repository is

Analysis code for an LDLR promoter saturation-mutagenesis plasmid library (plasmid MPRA) measured by PacBio Fiber-seq. The goal is to relate single-base substitutions to small (TF-sized) protection, promoter accessibility (MSPs) and published functional activity. The highest-priority data set is HepG2 replicate 3. Its data live on Klone; exact paths are in handoff §1/§8 and in `config/*.json`.

## Current code

```
fsmpra/          # Current package: python -m fsmpra <process|catalog|molecules|measure|plot|dashboard|compare>
                 #   or python -m fsmpra.<module> (ledger, footprint_join, statistics, heterogeneity, ...)
workflow/*.smk   # Current Snakemake workflows, one per stage: snakemake -s workflow/<name>.smk
config/*.json    # Per-workflow configs (ldlr_*.json, rep3_existing_footprints.json)
tools/           # Current: summarize_*, validate_*, render_*, plot_*, audit_rep3_caller_provenance,
                 #   inspect_heterogeneity_annotations, package_handoff_evidence, independent_review_*
dashboard/       # Viewer templates/JS; generated HTML goes under results/ (Git-ignored)
tests/           # unittest suite for fsmpra and current tools (+ Node checks for the viewers)
docs/evidence/   # Committed evidence snapshot (20261020) and independent-review tables
```

Entry points in dependency order (from handoff §8; usage docs are in `docs/`):

| Stage | Workflow / code | Instructions |
|---|---|---|
| Original pilot catalog / published function joins | `workflow/pilot.smk`, `fsmpra.catalog` | `pilot_usage.md` |
| New raw BAM through features | `workflow/process.smk`, `config/ldlr_processing.json` | `processing_workflow.md`; only a bounded smoke test was completed end-to-end |
| Existing-alignment complete ledger | `workflow/ledger.smk`, `fsmpra.ledger` | `rep3_ledger.md`; the saved run config is in the evidence bundle |
| Reuse historical calls | `workflow/footprint_join.smk` | `footprint_join.md` |
| Barcode arrays, diagnostic inference/calibration | `workflow/statistics.smk` | `statistics_usage.md` |
| Molecule heterogeneity and barcode summaries | `workflow/heterogeneity.smk` | `heterogeneity_usage.md` |
| WT sampled raw calls, quality, first benchmark | `workflow/footprint_quality.smk` | `footprint_quality_usage.md` |
| Large WT presence/coverage and calibration | `workflow/heldout_calibration.smk` | `heldout_calibration_usage.md` |
| Landscape and genomic atlas | `workflow/landscape.smk`, `fsmpra.genomic_coordinates` | `landscape_usage.md`, `genomic_coordinate_audit_20261020.md` |
| Motifs and barcode audits | `fsmpra.biology`, `tools.render_biology_atlas` | `rep3_candidate_biology_20261020.md`, `biology_atlas_usage.md` |
| Candidate raw extraction/sensitivities/viewer | `fsmpra.raw_candidate_extract`, `fsmpra.raw_candidate_summary`, `tools.render_raw_candidates` | `rep3_raw_candidates_20261020.md` |

`workflow/compare.smk` (`fsmpra compare`) compares existing pilot measurements across HepG2 runs with their reaction conditions kept explicit. It does not pool them (`condition_comparison_20261014.md`). Later stages read dated upstream paths under `results/`, so each workflow's config names the exact inputs it was run on.

**Coordinates:** reference `LDLR_MPRA_NotI` (4,647 bp). The mutable element is plasmid `[1234,1445)` = oligo 21–231 = GRCh38 chr19:11089285–11089495 (+), so `genomic_pos1 = 11089264 + oligo_pos1`. The barcode is `[3294,3309)`. Displayed positions are 1-based inclusive; arrays and BED files are 0-based half-open. Never give genomic coordinates to backbone or barcode bases.

## Legacy code (historical, not revalidated)

The following predate the barcode-aware cohort, and their calibration claims were **not revalidated** (handoff §2, §8):

- `workflow/Snakefile`, `config/config.yaml` and `workflow/scripts/`: `01_filter_reads` → `02_call_variants` → `03_analyze_library`, `common.py`, HDF5 output, the empirical-null, motif and co-occupancy layers.
- `extras/*.py` and the notebook (Dash/HDF5 browsers, QC and TF annotation).
- The legacy `tools/`: `*marginal*`, `*cooccupancy*` (including `tools/test_cooccupancy_*`), `wt_vs_wt_calibration.py`, `profile_*`, `phase0_benchmark.py`, `verify_null_equivalence.py`, `calibration_sweep/` and the `*.sbatch` files.
- The legacy docs: `statistical_framework.md`, `sensitivity_and_power.md`, `stage3_redesign_plan.md` and `dashboard_design_prompt.md`.

Keep them for history. Don't treat them as the pipeline, extend them by default, or present their output as an independent reproduction. The Pixi `run` task (`snakemake --cores 8`) and `dry-run` task target this legacy Snakefile; don't use them.

`extras/data/` is still used. `fsmpra.biology` reads the JASPAR 2024 and HOCOMOCO v12 files there. `fsmpra.genomic_coordinates` reconciles the legacy `LDLR_mpra_scores.tsv`, which needs `11092732 - legacy_position` and complemented alleles. Complementing alleles does not negate activity.

**Legacy definitions:** the BAM tags `nc`/`bn` (set by `01_filter_reads`) and `PV`/`VC`/`PR`/`BK`/`CS` (set by `02_call_variants`), and the `ns`/`nl`/`nq` caller tags, are used only by the legacy scripts. Current code reads FiberHMM annotations from `MA`/`AQ` (`fsmpra.footprint_join`). The `common.py` footprint bins (`sub_TF` 10–19, `TF` 20–40, `PIC` 41–80, `NUC` 81+ bp) are legacy too. Current analyses use small calls of 10–60 bp and nucleosome-sized protection of 85–250 bp, with lengths classified before clipping to the element.

## Operational rules

- Run substantive work only inside a verified Slurm allocation. Size threads and memory from the allocation and the per-stage resource declarations, never from visible node totals. Use node-local scratch for temporary SQLite/BAM work. Set threads to 1:
  ```bash
  export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
  export MPLCONFIGDIR=/tmp/fsmpra-matplotlib
  ```
- Outputs are atomic and refuse to overwrite (`fsmpra.io.atomic_output` / `new_output` raise `FileExistsError`). Use a new directory and a copied config. Changing `--output` alone does not retarget the downstream stages.
- Frozen-source guards check hashes and/or path/size/mtime, and a fresh checkout changes mtimes. Never disable a guard or edit a historical audit to force a pass; recreate the audit from verified inputs instead.
- Never use `pysam.modified_bases` on these mixed MM/ML tags, because it can segfault. Use the custom parsers in `fsmpra` (`processing`, `heterogeneity`, `footprint_quality`, `raw_candidate_extract`). They consume all MM groups with their matching ML columns, handle reverse orientation and preserve missing channels and unknown skipped sites.
- Missing measurements are not zeros. A present, evaluated-empty annotation track is zero. An absent annotation, a variant without supported molecules, or a known site without a retained call is missing.
- Tests: `python -m unittest discover -s tests -p 'test_*.py'` (with the thread settings above). Software agreement does not establish biological validity.

## Scientific guardrails

- Don't pool different libraries or reaction conditions. HEK293T uses an older library and HepG2 a shorter element. HepG2 runs 1–2 used 10 min with 1.5 µL Hia5, while run 3 used 5 min with 3 µL. Compare them; don't combine them.
- High-methylation molecules come from the same Fiber-seq reaction. They are not naked-DNA or protein-free controls.
- Long nucleosome-class annotations (e.g. >1,000 bp) are unresolved protection, not single nucleosomes.
- Don't tune filters against significance counts or agreement with functional data. Don't interpret the withheld candidate q-values as discoveries.
- Footprint coverage is an assay measurement, not occupancy by an identified protein. Published functional statistics do not validate footprint effects.

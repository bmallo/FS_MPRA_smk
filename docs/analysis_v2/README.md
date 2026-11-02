# LDLR Fiber-seq MPRA, analysis v2: start here

A step-by-step re-analysis of the LDLR promoter saturation-mutagenesis Fiber-seq MPRA (HepG2 and HEK293T),
run 2026-10-28 to 2026-11-02. Each confounder was tested before any biological claim, and every setting
was fixed in a config before its results were looked at. The results will feed a revamp of this
repository's MPRA Fiber-seq pipeline, aimed at confident footprint calls without inspecting each footprint
by hand.

**Status (2026-11-02):** Phase A (A1–A11) is complete, plus exploratory follow-ups: the 60–86 and
114–138 footprints, a data-defined footprint map and a confidence display. Phase B waits for data the user
will generate: naked-DNA controls on this library (B1) and a WT time course (B2). Nothing in Phase B has
been started.

## Read in this order

1. [POWER_AND_LIMITS.md](POWER_AND_LIMITS.md): what the data can and cannot support. Read it before
   proposing any analysis.
2. [FINDINGS.md](FINDINGS.md): approach and numbers for each step, with sources.
3. [DATA_GUIDE.md](DATA_GUIDE.md): files, schemas, coordinates, verified code recipes, and how to rebuild
   the explorer and deck.
4. `config/analysis_v2/rules.json`: the authoritative record of every setting, result, decision,
   amendment and post hoc label, one section per step.

## Living artifacts

| What | Where |
|---|---|
| Slide deck, plain-language slides per step | https://claude.ai/artifact/NYSJ7cHat24e2SRj4SErwG · source `docs/analysis_v2/deck/` |
| Promoter explorer, interactive tracks and lollipops | https://claude.ai/artifact/JTJkTUZ5NtDK6WSVt4BYsH · page `tools/analysis_v2/explorer.html`, data `results/analysis_v2_20261028/explorer_v6` |
| Original plan | `/mmfs1/home/bmallo/.claude/plans/fuzzy-toasting-sunbeam.md` |
| Confidence-call design note | `docs/footprint_confidence_design.md` |
| Code and tests | `tools/analysis_v2/*.py` (one script per step), `tests/test_analysis_v2_*.py`, `tests/test_explorer_data.py` |
| Data access helper | `tools/analysis_v2/access.py` (loads molecules, reruns the A9/A10 test) |

## Ground rules that carried through this work

- HepG2 and HEK293T are never pooled or overlaid. They differ in 8–9 ways, including the m6A caller.
- Compute runs on the standing node or in Slurm, with BLAS threads set to 1. Use the `envs/general` Python.
- New settings go into `rules.json` before their results are viewed; anything chosen later is labelled
  post hoc.
- Existing `results/` directories are never overwritten. Code, tests and data are never edited to make
  something pass.
- Results are reported as numbers, not conclusions; nulls are reported with their MDE. Effects within
  10 bp of an SNV are "possible Hia5 sequence effect until B1".
- Write "accessible (marked)" and "inaccessible (unmarked)". Use the data-defined footprints F01–F24 rather
  than fixed size bins.
- Commit only with the user's go-ahead; never commit `results/` or read IDs.

## Candidate sub-projects

Each is self-contained enough for its own session. For each one: where to start, and what is already known.

1. **Designed 1-bp indels**: 1.89 M unused molecules; 146 deletions and 637 insertions, ~1,200 molecules each.
   These test helical phasing and spacing. They need the same informative-base (A4) and masking (A7)
   treatment as SNVs. Start from A2 `a2_accounting/HepG2_run3`, A4
   `variants_visibility.tsv` (indels are classed there) and `a9_extract.py` (keeps WT and SNV only).
2. **Element-wide shift and the Repeat 3 Sp1 site**: 17 SNVs at 143–151 raise protection across the
   element, and the shift tracks nucleosome-size calls and lower expression. Open questions: nucleosome
   positioning versus occupancy; single-molecule co-occupancy of the Sp1 sites and nucleosomes. Start from
   `a11_motifs_function/function_snvs.tsv`, `a10_effects/followup`, and `access.py`.
3. **Footprint audit beyond 60–86 and 114–138**: apply the same tests (labeling dependence versus a naked
   simulation, gap variants, allele disagreement) to the other footprints F01–F24. The 20–40 bp coverage at
   oligo 182 (CREB1/CREM window) already falls with labeling like the naked simulation. Start from `fp_region.py`, `fp_features/`,
   and `fp60_86.py`.
4. **Pipeline revamp** (`workflow/`): fold in the A9 random-effects test, the A7 masks, the data-defined
   footprint map and the confidence flag, replacing fixed size bins and the old unmasked comparisons.
   Start from `a9_framework.py`, `a10_effects.py` and `fp_features.py`. Wait for B1 before settling the
   emission model.
5. **Phase B1/B2 analyses** (when data arrive): rebuild the emissions per variant sequence, rerun A6–A10,
   and fit the full confidence model with the sequence-information QC tool. Start from the plan's Phase B,
   `docs/footprint_confidence_design.md`, and the plasmid-llr skill (`floor` and `emissions` stages).
6. **HepG2 replicate runs 1–2**: the only true replicate pair, so far measured shallowly. Use them to check
   whether A10 effects replicate in direction (within-run tests, combined at the effect level). Start from
   `a2_ledgers/`, `a8_repro.py` and `a9_extract.py`.
7. **HEK293T**: the genotype/barcode audit was never done there (deferred in A3). The library is
   bottlenecked, so only motif-level and WT-landscape questions work. Start from
   `a10_effects/HEK293T_run1/motifs.tsv` and `a10_effects.clone_mean_test`.

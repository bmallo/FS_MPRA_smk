# Starting prompt: R01 figure from the LDLR plasmid Fiber-seq MPRA

Paste everything below the line into a new session in this repository.

---

I'm writing an R01 grant. Its goal is an assay plus a calibration method that classifies non-coding variants as
pathogenic or benign. Plasmid Fiber-seq on non-coding variants will give multi-modal data (transcriptional
output, chromatin architecture, cell-type specificity), and that data will provide "truth sets" of benign and
pathogenic variants for calibration. I want a preliminary-data figure built from our LDLR promoter saturation
mutagenesis data (HepG2 and HEK293T). It should show that this kind of multi-modal measurement of variant
effects is feasible.

The data does not need to be perfect or statistically significant. I want examples that are convincing and
move in the right direction, and are honest about what they are.

## Figure: three panels

1. **Activator disrupted.** A single variant or a cluster of variants that either lowers the occupancy of a
   TF-sized footprint and lowers promoter activity, or raises the occupancy of a TF-sized footprint and raises
   activity.
2. **Repressor disrupted.** Variants that remove a footprint and raise activity, or recruit a footprint and lower
   activity.
3. **Cell-type comparison.** Concordance or discordance of the panel 1–2 variants between HepG2 and HEK293T.

## Side question: clinically reported variants

I also want to look at clinically reported LDLR promoter variants. A handful of them are in the final figure of
our plasmid Fiber-seq paper. Do we have enough data to show these variants driving specific changes in chromatin
architecture and activity? If so, some could go into this figure.

## Before you start: read the existing analysis

This repository has a completed step-by-step re-analysis of these data ("analysis v2"). Read these, in order:
- `docs/analysis_v2/README.md`
- `docs/analysis_v2/POWER_AND_LIMITS.md`
- `docs/analysis_v2/FINDINGS.md`
- `docs/analysis_v2/DATA_GUIDE.md`

Everything you need is already computed there: per-variant tests, footprint maps, activity, and a data-access
helper. Build on those outputs and do not re-derive them. `tools/analysis_v2/access.py` loads per-molecule data
and reruns the variant-vs-WT test.

### Leads from that analysis (check them, don't assume them)

**Activator candidates**
- Repeat 3 Sp1 site (oligo 143–151). Variants that break the site lower 10–40 bp call coverage there (the one
  motif with a dose-response, q 0.017). The 17 SNVs at 143–151 have about 6-fold lower published activity.
- Footprint F04 at 60–86 (real protection over ~66–81; 86% confident calls). Variants in its core lower
  activity a little (−0.5 to −1.3 log2). Alleles at one position can move the footprint in opposite
  directions: 71T>A −3.5 vs 71T>C +2.9 points; 77A>T vs 77A>G. These are worth checking against activity
  direction.
- Repeat 1 (89–104) is a weak Sp1 match. Repeat 2, SRE-1 (124–139), has no A/T in its core (CCCC at
  131–134), so Fiber-seq barely sees it even though it matters functionally.

**Repressor candidates.** None has been looked for directly. Places to start:
- SNVs with the highest activity (published log2 values go up to about +1.7) and how their nearby footprints
  change.
- SNVs that lower protection across the element: 80G>C, 110G>C, 115A>C, 117A>C, 117A>T, 126A>T, 157A>G and
  209G>C (A10, ±5 bp mask).
- The "element-wide shift" in A11 links gained protection away from the variant (it tracks nucleosome-size
  calls) to lower activity within a position: −0.35 log2 per point. That is a "recruits protection, lowers
  activity" pattern, but at nucleosome scale and mostly from the Repeat 3 variants. Ask me before using it as
  the repressor panel.

**Places to avoid or caveat**
- F11 / oligo 114–138 is mostly chance calls (3% confident).
- The 20–40 bp coverage at oligo 182 falls with labeling, like naked DNA does: a likely artifact.

## Constraints that matter for this figure

- **Activity is HepG2 only.** It comes from the published MPRA (Kircher 2019; files in
  `/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/mpra_published/`). It is a separate experiment
  from the Fiber-seq. Ask me whether any HEK293T activity data exist before designing panel 3. Without them,
  panel 3 compares chromatin only.
- **HEK293T is bottlenecked.**
  - Of the 633 designed SNVs, 244 are observed. About 105 have ≥ 300 molecules, and most come from a single
    clone. Per-variant HEK statistics are not valid, so show molecule-level or clone-level views, labelled as
    such.
  - The construct is different: the element is reversed in a 318-bp mutable region, so convert to oligo
    coordinates.
  - The m6A caller is different too. Never pool or overlay HepG2 and HEK293T. Show them side by side, on the
    same oligo axis.
  - HEK footprints are named K01–K23 and HepG2 footprints F01–F24.
  - `access.load_element('HEK293T')` has not been tested yet. Check it against
    `a10_effects/HEK293T_run1/motifs.tsv` and the explorer's HEK data before relying on it.
- **No naked-DNA controls yet.**
  - Footprint occupancy is relative to in-cell, highly labeled molecules, not absolute.
  - A change within about 10 bp of a variant could be Hia5 sequence preference rather than protein. Prefer
    examples where the footprint change reaches beyond ±5 bp, or survives the ±5 bp mask.
  - "Drops its base" SNVs (any A/T→G/C or A↔T change) change the evidence the caller sees, so show their
    masked comparison.
- Always use the masked comparison (A7 mask by default) and labeling-matched WT. Prefer confident calls
  (P_chance < 10⁻³) for footprint occupancy.
- Use "accessible (marked)" and "inaccessible (unmarked)". Name every plotted metric by what it measures.

## Clinical variants

There is no ClinVar or paper-figure variant list on disk.
- Ask me for the variants in the paper's final figure.
- Offer to pull ClinVar / LOVD LDLR promoter variants for chr19:11,089,285–11,089,495 (hg38; oligo
  k = chr19:11,089,264 + k). Confirm the download with me first.
- Map each variant to the library. For each one, report:
  - whether it is designed (SNV or 1-bp indel);
  - molecules and clones in HepG2 run 3 and in HEK293T;
  - published activity;
  - its A10 and footprint (fp_features) results.

Designed 1-bp indels were never analysed (1.89 M molecules). If clinical indels matter, tell me what it would
take to bring them in.

## How to work

Go step by step. After each step, stop, explain it plainly (numbers, not conclusions), and wait for my OK. Use
the pub-analysis skill for the transparent, step-by-step workflow, and the dataviz skill for the figure.

Suggested steps:
1. **Re-orient.** Summarise what the existing results already say about each panel. List the candidate regions
   with their footprint, activity and HEK293T coverage.
2. **Shortlist.** Rank the candidates for panels 1 and 2 by:
   - effect direction agreeing between footprint and activity;
   - effect size;
   - confident-call support;
   - robustness to the ±5 bp mask;
   - distance of the footprint change from the variant;
   - number of independent variants (alleles or positions) agreeing.

   Show a table, plus one explorer-style view per candidate (WT vs variant molecules, footprint occupancy,
   activity).
3. **Clinical variants.** Map, check coverage, and decide whether any are figure-worthy.
4. **HEK293T comparison** for the shortlisted variants. Report concordant, discordant, or not measurable, and
   why.
5. **Draft panels.** Draft two or three layout options for the three-panel figure, then refine the one I pick.
   Include a short caption and a list of the caveats a reviewer would raise.

## Housekeeping

- **Where to run.** Compute on the standing node (`klone-node`) or in Slurm, never the login node. Set BLAS
  threads to 1. Use the Python at `/mmfs1/home/bmallo/miniconda3/envs/general/bin/python` with `PYTHONPATH=.`.
- **Where to write.**
  - New code goes in `tools/r01_figure/`, with tests.
  - Outputs go in a new `results/r01_figure_<YYYYMMDD>/`. Never overwrite existing `results/`.
  - Record any selection rule (for example, how candidates are ranked) in a small config before looking at the
    ranking, and label later choices as post hoc.
- **Checking visuals.** Before showing me an HTML page, look at it with
  `tools/analysis_v2/render_page.py` (headless Firefox).
- **Git.** Commit only when I say so. Never commit `results/` or read IDs.

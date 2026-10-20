# WT landscape and descriptive variant biology — 2026-10-19

The full supported WT/SNV landscape is now categorized, and all designed SNVs
can be explored against WT in the
[interactive atlas](../results/rep3_landscape_20261018/atlas/index.html).
The analysis supports concrete candidate patterns for follow-up, without new
significance testing or assigning protein identities. The run directory retains
the October 18 start-date identifier; analysis and reporting completed October 19.

The full joined stream contained 8,715,827 molecules. The analysis used
2,102,906 eligible molecules from 568,354 supported barcodes, including
985,264 WT molecules from 268,283 barcodes. All 633 designed SNVs are retained
in the atlas. Of these, 630 meet the existing display floor of 50 barcodes and
100 molecules; 160A>G has no eligible supported-barcode data and is explicitly
missing. The floor is not a significance or reliability threshold.

All coordinates below are **1-based oligo positions**, with the mutable element
at 21–231. The reference interval is [1234,1445), and genomic coordinates in
the tables are GRCh38. All primary summaries give each barcode equal total
weight. Results describe one preparation and inherited caller annotations.

## WT architecture

| Annotation category, using original call length | Mean element coverage | Molecules with any overlapping call |
|---|---:|---:|
| Small TF-class, 10–60 bp | 2.60% | 21.84% |
| Other TF-class lengths | 1.69% | 22.83% |
| Nucleosome-class <85 bp | 0.19% | 1.55% |
| Nucleosome-sized, 85–250 bp | 19.75% | 35.54% |
| Nucleosome-class 251–500 bp | 8.45% | 11.91% |
| Nucleosome-class 501–1,000 bp | 5.44% | 6.70% |
| Nucleosome-class >1,000 bp | 13.66% | 14.12% |
| MSP accessibility | 52.59% | 79.03% |

These categories can overlap and do not form a partition. In particular,
small footprints can occur within MSPs. Long protected calls are not single
nucleosomes. Lengths are classified before clipping to the element; an
element-wide protected interval can originate from a much longer call.

![WT landscape](../results/rep3_landscape_20261018/atlas/wt_landscape.png)

Four local peaks in the smoothed 10–60-bp profile define descriptive 20-bp
neighborhoods. The unsmoothed profiles and the original fixed windows remain
available. These windows are not estimates of protein contact boundaries.

| Neighborhood | Smoothed-profile peak position | Oligo interval | Mean small-footprint coverage |
|---|---:|---|---:|
| WT_P1 | 62 | 52–71 | 6.51% |
| WT_P2 | 84 | 74–93 | 5.05% |
| WT_P3 | 124 | 114–133 | 4.95% |
| WT_P4 | 159 | 149–168 | 2.26% |

Another prominent feature is visible in the **other-TF-class** track near
position 86 (13.42% per-base coverage). It is excluded from the 10–60-bp track
by definition. Across WT element-overlapping TF-class calls there are 244,650
calls shorter than 10 bp, 279,664 calls of 10–60 bp, and 25,478 longer than
60 bp. Those totals are counts of calls, not molecules. This makes the length
definition consequential; neither the other-length feature nor its disappearance
under a length filter should be assigned a protein identity without further work.

The WT molecule categories also show that a single average curve hides several
configurations. About 12.11% of barcode-weighted molecules have a small footprint
with neither an overlapping 85–250-bp nor >250-bp nucleosome-class call;
7.47% have small footprints plus nucleosome-sized calls, without long calls.
Overall, 31.85% have some >250-bp nucleosome-class protection overlapping the
element. These are measured combinations anywhere in the element, not evidence
for simultaneous TF–nucleosome contact at one site or a proven cellular compartment.

[Original length distributions and molecule combinations](../results/rep3_landscape_20261018/atlas/wt_lengths_states.png).

## Candidate pattern: a local change at position 125

All three substitutions at position 125 show lower coverage in WT_P3:

| Variant | Supported barcodes | WT_P3 coverage | Published LDLR log₂ activity effect |
|---|---:|---:|---:|
| WT | 268,283 | 4.95% | — |
| 125A>C | 380 | 1.37% | +0.08 |
| 125A>G | 480 | 1.23% | +0.47 |
| 125A>T | 438 | 1.39% | −0.07 |

![Position 125](../results/rep3_landscape_20261018/patterns/position_125.png)

The A→T substitution preserves the A/T count, so loss of an informative base
cannot be the sole count-based explanation shared by all three variants.
However, A→T changes strand orientation and local sequence context; this is
**not** a sequence-bias control. The three substitutions' published activity
effects are not uniformly reduced, despite their similar observed footprint
losses. This is a useful example of why footprint loss alone cannot identify
an activator or establish a direct occupancy–expression relationship.

The next audit for this pattern should compare unchanged informative positions,
call boundaries, and context-dependent detection. These checks have not yet
been performed on these three candidates.

## Candidate pattern: a region affecting broader architecture

Across all 24 substitutions at oligo positions **143–150**, WT sequence
`CTCCTCCC`:

- 22/24 have lower whole-element MSP coverage; the median difference is
  **−2.00 percentage points**.
- 22/24 have higher coverage by 85–250-bp nucleosome-class calls; the median
  difference is **+2.28 percentage points**.
- 23/24 have lower published activity in each separately retained source.
  The median log₂ effects are −2.50 and −3.08 for the two source files.

![Region 143–150](../results/rep3_landscape_20261018/patterns/region_143_150.png)

This is consistent with a candidate sequence region that helps maintain an
accessible, transcriptionally active promoter configuration. It does not identify
the bound protein, establish mediation by nucleosomes, or exclude labeling and
sampling effects. The region was chosen after inspecting the data; counts of
concordant substitutions are descriptive and are not binomial significance tests.

Examples include 145C>T (MSP −5.51 pp; nucleosome-sized coverage +3.00 pp;
published activity −2.87 log₂) and 146C>G (MSP −4.01 pp; nucleosome-sized
coverage +4.09 pp; published activity −2.70 log₂). The latter preserves A/T
count, but still changes context. The exception 147T>G has higher small-footprint
coverage (+0.67 pp across the element), little MSP change (−0.07 pp), and higher
published activity (+0.63 log₂); it is worth retaining as a contrasting candidate.

As another directional contrast, 157A>G has MSP +6.06 pp, nucleosome-sized
coverage −4.16 pp, and published activity +0.43 log₂. Its mean distal m6A
difference is approximately zero (−0.004 pp), which motivates examining local
versus global changes. This does not prove that other aspects of methylation
or chromatin heterogeneity are matched.

## Functional context across variants

Among the 630 SNVs meeting the coverage display floor, Pearson correlations
between whole-element coverage differences and published activity were:

| Coverage difference | Published LDLR | Published LDLR_2 |
|---|---:|---:|
| Small TF-class, 10–60 bp | +0.174 | +0.174 |
| MSP accessibility | +0.293 | +0.315 |
| Nucleosome-sized, 85–250 bp | −0.407 | −0.425 |

At this descriptive level, nucleosome-sized protection and broad accessibility
track published activity more closely than whole-element small-footprint
coverage. Whole-element averages can dilute local footprint effects. These
correlations have no significance claims, do not account for dependencies among
neighboring variants, and do not establish causality. The two published files
are kept separate and are not assumed independent biological replicates or
matched functional measurements from replicate 3.

## What is ready and what remains

The atlas contains all eight per-base categories, WT state summaries, all SNVs,
25 local comparison windows, coverage and A/T-count flags, distal context, and
separate functional-source values. Clicking a mutation cell or variant table
row displays its profile against WT. The underlying tables and SVG figures
are exportable. [Workflow usage](landscape_usage.md) describes reproduction.

The highlighted patterns are candidates for the next phase, not confidence
audits already passed. The most useful order is: check shared informative sites
and caller sensitivity around position 125; assess the spatial distribution and
barcode support of the 143–150 architectural pattern; then inspect contrasting
variants such as 147T>G and 157A>G. Motif assignment and activator/repressor
classification remain downstream of those checks.

## Verification

Nine targeted Python tests passed, and the actual-atlas JavaScript smoke test
passed initialization, selection, filtering, and missingness checks. All barcode
depths reconcile exactly. Independent TF/MSP aggregation agrees within floating
accumulation tolerance (maximum raw barcode-weight sum difference 2.24e-8).
All 267,548 prior TF/MSP per-base profile rows also reconcile, including 422
missing measurements, with maximum measured-fraction difference 2.00e-14.
Summed nucleosome-length category means agree with the previous independent
nucleosome-union means for all 633 measured alleles, with maximum difference
8.22e-15. Missing-data measurements remain missing in tables and the viewer.

The landscape build used one CPU for 3:46 and peaked at about 1.34 GiB of
resident memory inside the verified Slurm allocation. Atlas packaging took
about 20 seconds. Scientific figures were inspected. Browser interactions were
tested with a DOM/canvas harness, not a visual browser session. Snakemake was
dry-run checked; full execution used equivalent direct CLI commands with
per-stage source audits. No original BAM, barcode assignment, or footprint
annotation was changed.

# Replicate-3 biological candidates and sensitivity checks

The annotated atlas adds the genomic sequence, A/T opportunities, versioned
motif candidates, all local substitutions, and barcode-level sensitivity
comparisons. The original genomic atlas measurements are preserved unchanged.
These results are descriptive observations from one preparation, not new
statistical discoveries or protein assignments.

## What the sequence context adds

The published LDLR SRE-1 sequence, `ATCACCCCAC`, maps exactly to
**GRCh38 chr19:11089391–11089400** in this construct. This is a literature
sequence anchor; the cited experiments established its role in LDLR regulation,
not its occupancy in this Fiber-seq preparation.
[Original SREBP/LDLR study](https://pmc.ncbi.nlm.nih.gov/articles/PMC20191/).

The independent sequence scan identifies JASPAR SREBF1 and SREBF2 models at
that same interval, with HOCOMOCO SRBP models in overlapping intervals. Motif
orientations refer to the orientation of the stored matrix; they do not change
the genomic allele orientation or imply opposite regulatory activity.

Across the 30 substitutions in this 10-base interval, the SREBF1 motif-score
change has descriptive Spearman correlations of **0.76 and 0.72** with the two
published activity sources, but only **0.02** with small-footprint coverage
change within the motif interval. SREBF2 gives a similar pattern. Thus the
functional data support a known sequence preference much more clearly than
the local footprint readout does. This comparison is post hoc, and the two
functional sources are not assumed independent.

An overlapping HOCOMOCO SP1 model spans chr19:11089397–11089417. Across its
63 substitutions, motif-score changes correlate descriptively with promoter
accessibility (0.50) and published activity (0.69 / 0.75), more weakly with
local small footprints (0.21). Other SP/KLF and zinc-finger models overlap
the same sequence, so these observations do not identify SP1 specifically.
Published work provides a reason to investigate SREBP/Sp1 cooperation, but it
does not establish that mechanism here.
[Original cooperation study](https://www.sciencedirect.com/science/article/pii/S0021925818829817).

## Two priorities for follow-up

### chr19:11089389 — reproducible local footprint loss within this preparation

This is oligo position 125, **two bases upstream of the published 10-base
SRE-1 sequence**. It lies within some longer SREBP model windows, but mutation
at this position is not itself disruption of a base in the published core.

The WT_P3 neighborhood is chr19:11089378–11089397. Its WT mean small-footprint
coverage is 4.95%.

| Substitution | Barcodes | Variant coverage | Difference, pp | Barcode half A / B, pp | Shared A/T, excluding ±2 bp, pp | Published activity, source 1 / 2 |
| --- | ---: | ---: | ---: | --- | ---: | --- |
| A>C | 380 | 1.37% | −3.58 | −3.62 / −3.54 | −3.82 | +0.08 / +0.25 |
| A>G | 480 | 1.23% | −3.72 | −3.63 / −3.82 | −3.89 | +0.47 / +0.39 |
| A>T | 438 | 1.39% | −3.56 | −3.62 / −3.50 | −3.71 | −0.07 / 0.00 |

All three losses persist with molecule weighting, minimum-depth restrictions,
and separate standardization by distal methylation or distal protection.
Removing any single variant barcode changes the corresponding mean by at most
0.09–0.12 percentage points, much smaller than the observed loss. Only
29–37 variant barcodes have nonzero coverage in this window; the signal remains
sparse despite the consistent aggregate differences.

**Biological interpretation:** a useful local protection/detection phenotype
that is robust to these checks. The functional responses differ, so the data
do not support a simple claim that this footprint is an activating protein.
A>T preserves A/T count, but changes adenine strand and local sequence context.
Coverage on shared A/T positions reuses existing calls; it does not prove the
difference survives sequence-matched methylation calibration or re-calling.

**Discriminating follow-up:** test WT and all three substitutions in independent
preparations with multiple barcodes, a short reaction-time series, and separately
prepared sequence-matched protein-free controls. Compare raw methylation at
unchanged sites as well as re-called footprints. Pair those measurements with
functional activity from the same construct library and preparation conditions.

### chr19:11089407–11089414 — a coherent promoter-architecture candidate

This is oligo 143–150, sequence `CTCCTCCC`. All 24 substitutions are retained.
It is a post hoc selected region, not a prospectively defined test.

- 22/24 substitutions lower whole-element accessibility; median −2.00 pp.
- 22/24 increase nucleosome-sized coverage; median +2.28 pp.
- 23/24 lower activity in each published source; median log2 effects
  −2.50 and −3.08, respectively.
- Local small-footprint differences are less coherent: 16 decrease and 8
  increase. The WT local small-footprint coverage is only 1.11%, and this
  eight-base segment contains just two A/T opportunities.

| Accessibility comparison | Variants with lower accessibility | Median difference, pp |
| --- | ---: | ---: |
| Equal barcode weighting | 22/24 | −2.00 |
| Barcode half A | 21/24 | −1.84 |
| Barcode half B | 21/24 | −1.83 |
| Molecule weighting | 22/24 | −2.16 |
| Barcodes with at least five molecules | 21/24 | −2.79 |
| Distal methylation standardization | 22/24 | −2.62 |
| Distal protection standardization | 22/24 | −2.40 |
| Shared A/T, excluding mutation ±2 bp | 22/24 | −2.02 |

**Biological interpretation:** the region is a strong candidate for maintaining
an accessible, active promoter configuration. That is a statement about the
regulatory sequence, not an assignment of an activating protein. SP1/SP4,
KLF, WT1 and other models overlap here; the strongest sequence match alone
does not identify the protein. The T>G substitution at chr19:11089411 is an
informative functional exception and should be retained as a contrast.

**Discriminating follow-up:** select substitutions spanning the local sequence,
including C>G changes that preserve A/T count and the chr19:11089411:T>G
exception. Test whether accessibility and nucleosome-sized changes reproduce
alongside matched activity. Use candidate-specific binding or perturbation
experiments only after narrowing the factor hypotheses with expression and
binding evidence. Nucleosome-sized coverage is descriptive here; this stage's
barcode sensitivity audits cover small footprints and MSP accessibility.

## All four WT neighborhoods

The neighborhood windows are 20-base summaries, not inferred physical binding
site boundaries. Motif overlaps include partial overlaps at window edges.

| Region | GRCh38 chr19 interval | WT small-footprint coverage | Local substitutions with lower / higher footprint | Sequence candidates at the scan threshold |
| --- | --- | ---: | --- | --- |
| WT_P1 | 11089316–11089335 | 6.51% | 49 / 11 | HOCOMOCO ZN652; RARG is sensitive to background choice |
| WT_P2 | 11089338–11089357 | 5.05% | 44 / 16 | Several models, including HSFY1, ZNF684/ZN684, SP4; no unique assignment |
| WT_P3 | 11089378–11089397 | 4.95% | 42 / 18 | SREBF/SRBP at the right edge, overlapping SP/KLF candidates |
| WT_P4 | 11089413–11089432 | 2.26% | 44 / 15 | Many overlapping GC-rich models; one SNV has no eligible measurements |

WT_P3 is especially sensitive to A/T-changing substitutions: median local
differences are −2.01 pp for A/T loss, +0.68 pp for A/T gain, and −0.52 pp
for retained A/T count. Those groups alter different positions and sequences;
this is a diagnostic association, not a controlled estimate of assay bias.
The pattern strengthens the case for sequence-matched controls.

## Methods and interpretation limits

Motif inputs are the repository's **879 JASPAR 2024 CORE vertebrate matrices**
and **1,443 HOCOMOCO v12 CORE matrices**, with full input checksums. These are
fixed historical versions, not the latest releases. JASPAR includes vertebrate
profiles; this is not a human-only expression-filtered motif set.
[JASPAR 2024 documentation](https://jaspar2024.elixir.no/docs/),
[HOCOMOCO v12 publication](https://academic.oup.com/nar/article/52/D1/D154/7424436).

The local scanner uses log2 odds, one total background-distributed pseudocount
per motif column, uniform independent-base background, and both strands. A
dynamic program computes exact tail probabilities for contributions rounded
to a 0.01-bit grid. The per-window threshold is 0.0001: **114 matches across
923,098 motif/position/strand comparisons**; 110 also meet the threshold using
the element's strand-symmetric GC background. These are unadjusted sequence
matches, with extensive model redundancy. The probabilities are not calibrated
binding probabilities, variant-effect tests, or corrected discoveries. This is
a local scanner, not a run of FIMO. No variant-created sites are searched.

All overlapping SNVs are scored at each retained WT motif's fixed placement,
with correct reverse-strand complementation. Descriptive Spearman correlations
use every measured overlapping substitution, without choosing variants for a
favorable correlation. Functional measurements remain separate by source.

Barcode audits cover all 633 designed SNVs against the same WT cohort across
seven endpoints: four WT neighborhoods, the adjacent eight-base region, and
whole-element small-footprint/MSP coverage. They reconcile the original
568,354 barcodes and 2,102,906 eligible molecules. The barcode partition uses
SHA256 parity; it is one internal split, not independent replication. Minimum
depth comparisons restrict both WT and variant. Distal standardization uses
WT quartiles of barcode-mean features and WT bin weights, requiring at least
10 variant barcodes in every bin. Missing or insufficiently supported
comparisons remain missing. Distal features may be consequences of regulatory
activity, so standardization is a sensitivity check, not a causal correction.

Shared-site comparisons use the same unchanged reference A/T positions for
both genotypes, excluding the mutation and optionally two flanking bases on
each side. They summarize existing footprint calls and do not re-call them
or directly compare raw per-site m6A. No confidence intervals, new variant
significance, HepG2 expression evidence, or locus-specific ChIP evidence have
been added in this stage.

## Reproduce on an allocated compute node

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
export MPLCONFIGDIR=/tmp/fsmpra-matplotlib
/mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.biology \
  --output results/new_biology_run
/mmfs1/home/bmallo/miniconda3/bin/python -m tools.render_biology_atlas \
  --biology results/new_biology_run --output results/new_biology_atlas
```

The defaults use the frozen replicate-3 genomic atlas and supported-barcode
cohort. Outputs refuse overwrites. All files remain local; the HTML embeds its
data and works without a server. Keep the output directory together for linked
table and figure downloads. Tests check scoring against exhaustive enumeration,
reverse-strand allele changes, shared-site masks and standardization support.

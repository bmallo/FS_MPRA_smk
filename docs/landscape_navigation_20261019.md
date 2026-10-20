# Variant navigation and functional lollipops

This records the original navigation update. The
[current genomic atlas](../results/rep3_landscape_20261018/atlas_genomic_20261019/index.html)
retains these controls and uses verified genomic axes; see the
[coordinate audit](genomic_coordinate_audit_20261020.md). Reproduction with the
current code produces genomic coordinates, while the archived output below
retains its original oligo axes.

Open the [updated local atlas](../results/rep3_landscape_20261018/atlas_navigation_20261019/index.html)
and choose **Open variant-by-variant browser** near the top.

- Previous/Next steps through all 633 SNVs in oligo-position order, then
  alternate-base order A/C/G/T. The ends wrap; missing and low-coverage variants
  are retained. The variant dropdown and mutation-map/table selection also work.
- Left/right arrow keys navigate when the comparison panel has focus. Native
  keyboard behavior in inputs, dropdowns and buttons is preserved.
- The functional lollipop appears immediately above the footprint comparison
  on the same 1-based oligo-position axis (21–231). Gray marks are individual
  SNVs; the orange mark identifies the exact selected substitution, not an
  average across substitutions at that position. A dashed line locates the
  mutation in both plots. Its functional value is also printed above the plot.
- The source dropdown switches between the two separately retained published
  MPRA sources. The lollipop y-scale depends on the complete selected source,
  not the current variant. Missing functional values are omitted, not plotted
  as zero. Hover identifies a point and clicking selects that variant.
- Table search/coverage filters do not change the variant navigation order.

## Read-count normalization

The displayed profile at a base is

`mean across barcodes(covered eligible molecules / eligible molecules in that barcode)`.

WT and each variant are calculated separately. Every eligible molecule spans
the whole element. Each barcode therefore contributes equal total weight,
regardless of its read depth. A barcode with 20 eligible reads has the same
total weight as one with two. Mean element coverage averages this profile
across the 211 bases. The any-call fraction instead counts whether each
molecule has any overlapping call before averaging within and across barcodes.

These are fractions normalized for depth, not raw read counts. The analysis
does not downsample WT and variants to equal depths or equal barcode numbers.
More observations improve support; normalization does not make a sparsely
covered variant as informative as WT. No additional methylation/chromatin-state
normalization is applied. Functional effects are the published log₂ values as
imported; the viewer does not estimate new RNA/DNA ratios.

## Reproduction and verification

The underlying atlas data and numerical exports are unchanged. Only presentation
assets were updated. A fresh output directory preserves the earlier atlas:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
NUMEXPR_MAX_THREADS=1 MPLCONFIGDIR=/tmp/fsmpra-matplotlib \
  /mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.landscape_atlas \
  --config config/ldlr_landscape.json \
  --output results/rep3_landscape_20261018/atlas_navigation_20261019
```

Run inside an allocated Slurm compute job. Existing outputs cannot be replaced;
choose a fresh path to reproduce again. The original atlas audit refers to the
earlier presentation sources preserved in its run snapshots. The new atlas audit
records the current presentation sources.

The actual-atlas DOM/canvas test exercises all 633 navigation steps, end wrapping,
keyboard focus behavior, source switching, selected-point coordinate alignment,
point selection, and missing measurements. This is an interaction/controller
test, not a visual browser session. The source HTML/JS are embedded in the local
HTML; no service, installation, or hosted deployment is required.

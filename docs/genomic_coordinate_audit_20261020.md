# Genomic coordinate and orientation audit

Completed review: October 20, 2026. Run directories retain their October 19
creation names.

The current atlas has the correct orientation. The 211-base mutable construct
matches **GRCh38 chr19:11,089,285–11,089,495 on the forward strand** exactly.
No reversal, allele complementation, or functional-effect sign change is needed
for the current published-data joins.

Open the [genomic atlas](../results/rep3_landscape_20261018/atlas_genomic_20261019/index.html).
The older atlas directories remain archived. The new viewer and scientific
figures use genomic positions; original oligo labels remain in parentheses.

## Independent checks

- The physical plasmid FASTA slice `[1234,1445)` agrees with the catalog and
  all 211 genomic reference bases.
- A fresh [UCSC hg38 sequence response](https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr19;start=11089230;end=11089548)
  agrees with the local indexed reference across 318 bases, including both
  element flanks. The response, retrieval time, and source checksums are saved.
- In the element plus 10 kb on each side, the complete element has exactly
  one forward exact match and no reverse-complement exact matches. This is
  a local orientation check, not a whole-genome uniqueness claim.
- All 1,083 reference alleles in `GRCh38_LDLR.tsv` and all 1,093 in
  `GRCh38_LDLR.2.csv` match the local genome. Each source contains the expected
  633 element SNVs. Every one joins to the correct altered construct sequence,
  with its original effect and p-value unchanged.
- The legacy `extras/data/LDLR_mpra_scores.tsv` uses a different orientation:
  `genomic_pos1 = 11092732 - legacy_position`, with both alleles complemented.
  All 954 score/p-value pairs reconcile with the second published source.
  Complementing alleles does **not** negate activity effects. This legacy
  coordinate system is not used by the current atlas.

The [machine-readable audit](../results/rep3_landscape_20261018/genomic_coordinates_20261019/audit.json)
records the evidence and input identities.

## Coordinate conventions and exports

For this verified library only:

`genomic_pos1 = 11089264 + oligo_pos1`

`genomic_pos1 = 11089285 + element_relative_pos0`

| Feature | Original coordinates | GRCh38 genomic coordinates |
| --- | --- | --- |
| Mutable element | Oligo 21–231 | chr19:11089285–11089495 |
| Example substitution | 125A>T | chr19:11089389:A>T |
| Candidate adjacent-base region | Oligo 143–150 | chr19:11089407–11089414 |
| WT_P1 neighborhood | Oligo 52–71 | chr19:11089316–11089335 |
| WT_P2 neighborhood | Oligo 74–93 | chr19:11089338–11089357 |
| WT_P3 neighborhood | Oligo 114–133 | chr19:11089378–11089397 |
| WT_P4 neighborhood | Oligo 149–168 | chr19:11089413–11089432 |

Displayed positions and intervals are **1-based inclusive**. BED is
**0-based half-open**; the element BED interval is `chr19 11089284 11089495`.
The atlas exports both interval conventions explicitly. Its region fields
`start0`, `end0`, and `peak0` retain element-relative coordinates for slicing
the arrays; `genomic_*` fields describe genomic coordinates.

New exports include `base_coordinates.tsv`, `genomic_profiles.tsv.gz`,
`wt_regions.bed`, and `wt_tf_boundaries_genomic.tsv`. Region effects and variant
summaries also carry genomic coordinates, assembly, chromosome, and strand.
Selected-variant JSON exports include the genomic coordinate vector.
The [candidate-pattern figures](../results/rep3_landscape_20261018/patterns_genomic_20261019)
also use genomic coordinates.

The conversion covers the mutable element in the WT/SNV atlas. Plasmid
backbone and barcode coordinates remain plasmid coordinates; they are not
assigned invented genomic positions. Raw molecule ledgers are preserved.
Indel VCF normalization is outside this conversion. Exported WT call boundaries
remain clipped at element edges; the conversion does not recover full physical
footprint boundaries.

## Validation

Thirteen Python tests passed, including one-base and whole-element interval
conventions, out-of-range rejection, and rejection of reversed, shifted, or
ambiguous sequence matches. The actual-data JavaScript controller test passed
all 633 navigation steps and checked genomic axes, allele labels, genomic
search, lollipop selection, and missingness. This is not a visual browser test.
The WT landscape and both genomic candidate figures were visually inspected.

The full output comparison confirmed unchanged measurements in **1,070,192
profile rows, 15,850 region-effect rows, 634 variant rows, and 1,509 WT boundary
rows**. All original per-variant JSON fields and functional correlations are
unchanged. See the
[validation record](../results/rep3_genomic_20261019_run/validation/audit.json).
Read-count/barcode normalization is unchanged.

## Reproduction

Run in an allocated Slurm job with one thread. The saved independent reference
response avoids a network dependency during reproduction. Use fresh output
paths; the tools refuse to overwrite existing results.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
export MPLCONFIGDIR=/tmp/fsmpra-matplotlib
/mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.genomic_coordinates \
  --ucsc-sequence results/rep3_genomic_20261019_run/inputs/ucsc_hg38_chr19_11089230_11089548.json \
  --output results/coordinate_audit_reproduction
/mmfs1/home/bmallo/miniconda3/bin/python -m fsmpra.landscape_atlas \
  --config config/ldlr_landscape.json \
  --coordinate-audit results/coordinate_audit_reproduction/audit.json \
  --output results/genomic_atlas_reproduction
/mmfs1/home/bmallo/miniconda3/bin/python -m tools.validate_genomic_atlas \
  --atlas results/genomic_atlas_reproduction \
  --output results/genomic_validation_reproduction
/mmfs1/home/bmallo/miniconda3/bin/python -m tools.summarize_landscape_patterns \
  --atlas results/genomic_atlas_reproduction \
  --output results/genomic_patterns_reproduction
```

The atlas builder now requires a matching complete genomic audit. The
Snakemake landscape workflow accepts its path via `--config coordinate_audit=...`;
generate this audit first. Source guards fail if audited inputs change.

## Proposed next biological step

Add the reference sequence, A/T opportunities, and versioned candidate motif
annotations to the genomic atlas. Examine all substitutions around the four
WT neighborhoods, emphasizing chr19:11089389 and chr19:11089407–11089414 as
previously selected descriptive examples. Compare each region's local small
footprints, whole-element accessibility, and published activity before choosing
a small set of candidates for barcode and sequence-observability audits.
Motif matches and coherent mutation patterns nominate mechanisms; they do not
establish a bound protein or an activator/repressor identity.

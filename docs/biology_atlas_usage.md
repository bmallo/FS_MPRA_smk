# Annotated biological atlas

Open [the completed atlas](../results/rep3_landscape_20261018/atlas_biology_20261020/index.html).
Read the [candidate report](rep3_candidate_biology_20261020.md) for findings,
limitations and proposed discriminating experiments.

The existing arrows, variant dropdown, functional lollipop, mutation map and
tables retain their behavior. The sequence panel follows the selected variant;
its selector can instead hold a WT neighborhood, the adjacent-base candidate,
or the whole element. Green marks identify A/T opportunities and orange marks
identify the selected substitution. Motif bars display up to eight candidates;
expand the table to inspect every overlapping match. Search by name or matrix
ID and filter by collection. Motif orientations and genomic coordinates are
explicit; alleles always remain on the genomic forward strand.

The candidate panel below the comparison shows every region's summary and
all local substitutions. The audit measurement selector switches between
local small footprints, whole-element small footprints, and whole-element
MSP accessibility. When the sequence follows a variant outside the named
regions, the audit endpoint is explicitly the whole element. Missing values
are displayed as missing, including the SNV with no eligible molecules.

The viewer is a self-contained local HTML file of approximately 17 MB.
Keep its directory together for table, report and figure download links.
No server, installation, deployment, or external runtime dependency is needed.
The motif databases are existing fixed historical versions, not updated
downloads. The implementation does not send the experimental data externally.

## Completed outputs

- `results/rep3_landscape_20261018/biology_20261020_complete/`: annotations,
  region tables, 114 motif hits, 5,022 fixed-site motif/variant comparisons,
  and 4,431 barcode sensitivity rows (633 SNVs × seven endpoints).
- `results/rep3_landscape_20261018/atlas_biology_20261020/`: annotated viewer,
  original genomic exports, biological exports, report, and two new
  PNG/SVG figures: `annotated_wt` and `position_125_sensitivity`.
- `results/rep3_biology_20261020_run/`: logs, source snapshots, manifest, and
  completed artifact validation.

The earlier `biology_20261020/` directory is a development-stage result;
the `_complete` directory adds motif/phenotype correlations and final checks.
The original atlas remains preserved. Its `atlas.json` and all original
non-HTML exports are byte-identical in the annotated output.

## Validation

Seventeen Python tests passed. Motif score tails were checked against exhaustive
sequence enumeration; reverse-strand ALT scoring, shared-site exclusion and
standardization support were tested. The actual-data DOM/canvas controller test
passed all 633 navigation steps and checked synchronized annotation updates,
motif filters, region selection, accessibility audit selection, and missingness.
The two new scientific figures were visually inspected. No full browser visual
QA was available.

The artifact validator confirmed all 114 motif sequences and genomic
placements, all 5,022 alternative motif scores, the 4,431 barcode audit rows,
original-data preservation, embedded JSON consistency, and working export file
targets. Barcode effects reconcile with the rounded original profiles within
0.00001 percentage points.

These are numerical and software checks, not statistical calibration or
independent biological validation.

## Reproduction and checks

Run inside a Slurm allocation with one thread. CLI reproduction commands are
in the candidate report. After producing a new annotated atlas:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /mmfs1/home/bmallo/miniconda3/bin/python -m tools.validate_biology_artifacts \
  --atlas results/new_biology_atlas --output results/new_biology_validation
/mmfs1/home/bmallo/miniconda3/bin/node tests/test_landscape_viewer.js \
  results/new_biology_atlas/atlas.json results/new_biology_atlas/biology.json
```

The frozen source checks deliberately fail if upstream measurements, references,
or audited implementations change. Choose new output directories for re-runs.

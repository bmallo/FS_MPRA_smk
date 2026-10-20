# Portable evidence snapshot — 20 October 2026

Start with the [project handoff](../../PROJECT_HANDOFF_20261020.md). This directory contains **34 selected input/result/configuration files (~3.3 MB)** copied from local inputs and completed analyses. It allows a GitHub reader to inspect the numerical evidence without access to the large Klone datasets. It does not contain raw reads or enough intermediate data to rerun every analysis offline.

| Files | Inspect for |
|---|---|
| `design.fa.gz`, `oligo_aliases.tsv.gz` | Exact unique mutable sequences, stable IDs and original design aliases |
| `published_GRCh38_LDLR.tsv`, `published_GRCh38_LDLR.2.csv` | Existing published functional tables, preserved separately without claiming independence |
| `ledger_config.json`, `ledger_summary.json`, `ledger_validation.json` | Actual ledger run settings, molecule/barcode totals, conservation |
| `barcode_benchmark.json` | Targeted full-depth exact/quality/Pacybara benchmark |
| `footprint_join_summary.json` | Annotation join and eligible-cohort counts |
| `statistics_release_status.json`, `statistics_calibration.json` | Withheld significance decision and first calibration failures |
| `heterogeneity_summary.json` | WT/variant cohort definitions and heterogeneity summary |
| `quality_benchmark.tsv` | Initial Welch/permutation comparisons |
| `heldout_calibration.tsv`, `heldout_effect_intervals.tsv` | Expanded error rates, small-effect power and intervals |
| `caller_provenance.json` | Retained ML distributions, historical processing and score limitations |
| `genomic_audit.json`, `genomic_validation.json`, `ucsc_hg38_chr19_11089230_11089548.json` | Sequence orientation, independent genomic response and measurement preservation |
| `landscape_categories.tsv` | Per-allele annotation-category coverage and presence summaries |
| `candidate_region_summary.tsv`, `candidate_region_variants.tsv` | WT neighborhoods, individual variants and functional comparisons |
| `motif_hits.tsv`, `motif_concordance.tsv` | Fixed-version sequence matches and descriptive phenotype correlations |
| `barcode_sensitivity.tsv.gz` | All 633 SNVs × seven barcode/observability audit endpoints |
| `raw_extract_audit.json`, `raw_summary_audit.json` | Latest candidate cohort, retained ML histogram, exact analysis definitions |
| `raw_comparisons.tsv.gz`, `raw_retention.tsv` | 15,552 comparisons and 504 cohort-support records |
| `raw_resampling.tsv.gz`, `raw_joint_states.tsv.gz` | Descriptive repeated selections and methylation/protection mixtures |
| `validation_panel.tsv` | Proposed ten-construct experimental panel with exact mutable sequences |
| `raw_validation.json`, `raw_viewer_validation.json` | Numerical conservation/reconciliation and DOM/canvas navigation validation |

`manifest.json` records each source path, source SHA256, packaged SHA256, sizes and any compression. Plain copies are byte-identical. Newly compressed TSVs preserve the original bytes after decompression; gzip timestamp is fixed to zero. Use `tools/package_handoff_evidence.py` to reproduce this explicit selection on the original workspace. It refuses changed evidence and enforces a package size limit.

Scoped `.gitattributes` preserve original line endings and trailing empty TSV fields. Do not trim these files: blank fields encode missing measurements and the manifest verifies exact source bytes.

Inspect compressed TSVs with Python `gzip.open(..., 'rt')`/`csv.DictReader(..., delimiter='\t')` or `gzip -cd`. Blank numeric fields mean unavailable measurements, not zero. Candidate statistical values and sequence-match probabilities have the limitations stated in the reports; this package does not authorize discovery labels.

Verify the snapshot from the repository root without scientific dependencies:

```bash
python - <<'PY'
import gzip, hashlib, json
from pathlib import Path
p = Path('docs/evidence/20261020')
m = json.loads((p / 'manifest.json').read_text())
for r in m['files']:
    b = (p / r['file']).read_bytes()
    assert len(b) == r['bytes']
    assert hashlib.sha256(b).hexdigest() == r['sha256']
    original = gzip.decompress(b) if r['transform'].startswith('gzip') else b
    assert len(original) == r['source_bytes']
    assert hashlib.sha256(original).hexdigest() == r['source_sha256']
print('Verified', len(m['files']), 'packaged files and source-byte identities')
PY
```

Audit records retain historical absolute paths/timestamps. They are evidence, not portable substitute inputs for source-guarded stage execution. The new-checkout procedure is in the handoff. Main configs, tests, reports, dashboard assets, analysis/validation scripts and motif collections are elsewhere in Git. BAMs, full molecule ledgers, barcode arrays, per-iteration calibration files and generated viewers remain in the local dated `results/` directories; see the reports for locations and reconstruction commands.

The handoff packaging regression run is saved separately in [validation_20261020.txt](validation_20261020.txt); it is not one of the 34 frozen input/result files.

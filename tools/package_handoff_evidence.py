"""Copy a bounded, explicit set of frozen result summaries into a Git handoff.

No BAMs, per-molecule ledgers, barcode matrices or generated dashboards are
included. The manifest hashes both source bytes and packaged bytes.
"""
import gzip
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'docs/evidence/20261020'
FILES={
 'design.fa.gz':'results/rep3_ledger_20261016/catalog/design.fa',
 'oligo_aliases.tsv.gz':'results/rep3_ledger_20261016/catalog/oligo_aliases.tsv',
 'published_GRCh38_LDLR.tsv':'/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/mpra_published/GRCh38_LDLR.tsv',
 'published_GRCh38_LDLR.2.csv':'/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/mpra_published/GRCh38_LDLR.2.csv',
 'ledger_config.json':'results/rep3_ledger_20261016/config.json',
 'ledger_summary.json':'results/rep3_ledger_20261016/summary.json',
 'ledger_validation.json':'results/rep3_ledger_20261016/ledger/validation.json',
 'barcode_benchmark.json':'results/barcode_benchmark_20261016/family_comparison_final/report.json',
 'footprint_join_summary.json':'results/rep3_footprints_20261017/summary.json',
 'statistics_release_status.json':'results/rep3_statistics_20261017/inference/release_status.json',
 'statistics_calibration.json':'results/rep3_statistics_20261017/calibration/report.json',
 'heterogeneity_summary.json':'results/rep3_heterogeneity_20261017/summary/summary.json',
 'quality_benchmark.tsv':'results/rep3_footprint_quality_20261018/benchmark/summary.tsv',
 'heldout_calibration.tsv':'results/rep3_heldout_calibration_20261018/calibration/summary.tsv',
 'heldout_effect_intervals.tsv':'results/rep3_heldout_calibration_20261018/plots/effect_and_interval_summary.tsv',
 'caller_provenance.json':'results/rep3_caller_provenance_20261018/audit_v2.json',
 'genomic_audit.json':'results/rep3_landscape_20261018/genomic_coordinates_20261019/audit.json',
 'ucsc_hg38_chr19_11089230_11089548.json':'results/rep3_genomic_20261019_run/inputs/ucsc_hg38_chr19_11089230_11089548.json',
 'genomic_validation.json':'results/rep3_genomic_20261019_run/validation/audit.json',
 'landscape_categories.tsv':'results/rep3_landscape_20261018/landscape/categories.tsv',
 'candidate_region_summary.tsv':'results/rep3_landscape_20261018/biology_20261020_complete/candidate_region_summary.tsv',
 'candidate_region_variants.tsv':'results/rep3_landscape_20261018/biology_20261020_complete/candidate_region_variants.tsv',
 'motif_hits.tsv':'results/rep3_landscape_20261018/biology_20261020_complete/motif_hits.tsv',
 'motif_concordance.tsv':'results/rep3_landscape_20261018/biology_20261020_complete/motif_concordance.tsv',
 'barcode_sensitivity.tsv.gz':'results/rep3_landscape_20261018/biology_20261020_complete/barcode_sensitivity.tsv',
 'raw_extract_audit.json':'results/rep3_raw_candidates_20261020/extracted/audit.json',
 'raw_summary_audit.json':'results/rep3_raw_candidates_20261020/summary/audit.json',
 'raw_comparisons.tsv.gz':'results/rep3_raw_candidates_20261020/summary/comparisons.tsv',
 'raw_retention.tsv':'results/rep3_raw_candidates_20261020/summary/retention.tsv',
 'raw_resampling.tsv.gz':'results/rep3_raw_candidates_20261020/summary/resampling.tsv',
 'raw_joint_states.tsv.gz':'results/rep3_raw_candidates_20261020/summary/joint_molecule_states.tsv',
 'validation_panel.tsv':'results/rep3_raw_candidates_20261020/viewer/validation_panel.tsv',
 'raw_validation.json':'results/rep3_raw_candidates_20261020_run/validation.json',
 'raw_viewer_validation.json':'results/rep3_raw_candidates_20261020_run/viewer_test.json',
}


def main():
    entries=[];payloads={}
    for name,source in FILES.items():
        path=ROOT/source
        if path.stat().st_size>25_000_000:raise ValueError('Unexpectedly large source: '+source)
        original=path.read_bytes();compressed=name.endswith('.gz') and not source.endswith('.gz')
        payload=gzip.compress(original,mtime=0) if compressed else original
        entries.append(dict(file=name,source=source,source_bytes=len(original),source_sha256=hashlib.sha256(original).hexdigest(),
                            bytes=len(payload),sha256=hashlib.sha256(payload).hexdigest(),transform='gzip, mtime=0' if compressed else 'byte-identical copy'))
        payloads[name]=payload
    if sum(map(len,payloads.values()))>10_000_000:raise ValueError('Evidence package exceeds 10 MB budget')
    DEST.mkdir(parents=True,exist_ok=True)
    for name,payload in payloads.items():
        target=DEST/name
        if target.exists() and target.read_bytes()!=payload:raise ValueError('Refusing to replace changed evidence: '+name)
    for name,payload in payloads.items():(DEST/name).write_bytes(payload)
    (DEST/'manifest.json').write_text(json.dumps(dict(schema_version='fsmpra.handoff_evidence.v1',as_of='2026-10-20',
        purpose='Selected summaries for inspection; not full raw-data reproduction',files=entries),indent=2)+'\n')
    print(f'Packaged {len(entries)} files; {sum(map(len,payloads.values())):,} bytes')


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Validate completed statistics tables and report coverage/calibration counts."""
import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control


def table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def summarize(run):
    run = Path(run)
    agg, inf = run / 'aggregates', run / 'inference'
    audit = json.loads((agg / 'audit.json').read_text())
    calibration = json.loads((run / 'calibration/report.json').read_text())
    release = json.loads((inf / 'release_status.json').read_text())
    metadata = table(agg / 'barcode_metadata.tsv.gz')
    variants = table(agg / 'variant_metadata.tsv')
    effects = table(inf / 'effects.tsv.gz')
    profiles = table(inf / 'profiles.tsv.gz')
    sensitivities = table(inf / 'sensitivities.tsv.gz')
    n = np.load(agg / 'barcode_molecules.npy', mmap_mode='r')
    counts = np.load(agg / 'barcode_counts.npy', mmap_mode='r')
    assert len(metadata) == len(n) == audit['counts']['barcodes']
    assert sum(int(r['molecules']) for r in metadata) == int(n.sum()) == audit['counts']['footprint_analysis_eligible_with_barcode']
    assert int(np.load(agg / 'cohort_molecules.npy').sum()) == audit['counts']['footprint_analysis_eligible']
    for a in range(0, len(n), 10000):
        assert np.all(counts[a:a+10000] <= n[a:a+10000, None, None])
    wt_id = next(v['variant_id'] for v in variants if v['genotype_class'] == 'WT')
    wt_idx = np.array([i for i, r in enumerate(metadata) if r['variant_id'] == wt_id])
    assert len(wt_idx) == calibration['wt_barcodes']
    assert int(n[wt_idx].sum()) == calibration['wt_molecules']
    wt_profile = np.zeros((2, 211))
    for a in range(0, len(wt_idx), 10000):
        idx = wt_idx[a:a+10000]
        wt_profile += (counts[idx] / n[idx, None, None]).sum(axis=0)
    wt_profile /= len(wt_idx)
    for row in profiles:
        if row['variant_id'] == wt_id:
            track = 0 if row['track'] == 'tf_small' else 1
            np.testing.assert_allclose(float(row['mean_equal_barcode']), wt_profile[track, int(row['relative_pos'])], atol=1e-12)
    assert len(profiles) == len({(r['variant_id'], r['track'], r['relative_pos']) for r in profiles}) == 634*2*211
    assert len(effects) == len({(r['variant_id'], r['endpoint_id']) for r in effects}) == 633*22
    assert len(sensitivities) == len({(r['variant_id'], r['endpoint_id'], r['sensitivity']) for r in sensitivities}) == 633*22*3
    p = np.array([float(r['candidate_p_value']) if r['candidate_p_value'] else 1. for r in effects])
    independent_q = false_discovery_control(p, method='by')
    for row, q in zip(effects, independent_q):
        if row['status'] == 'ok':
            np.testing.assert_allclose(float(row['candidate_q_value']), q, atol=1e-14)
            assert float(row['ci_low']) <= float(row['delta']) <= float(row['ci_high'])
        else:
            assert all(row[k] == '' for k in ('candidate_p_value','candidate_q_value','ci_low','ci_high'))
        if int(row['caller_effect_supported']):
            assert calibration['status'] == 'pass' and row['endpoint_id'] in calibration['validated_endpoint_ids']
            assert q <= .05 and row['status'] == 'ok'
    candidate = [r for r in effects if r['candidate_q_value'] and float(r['candidate_q_value']) <= .05]
    supported = [r for r in effects if int(r['caller_effect_supported'])]
    assert len(supported) == release['supported_effects']
    coverage = {r['variant_id']: (int(r['n_variant_barcodes']), int(r['n_variant_molecules'])) for r in effects}
    return dict(schema_version='fsmpra.statistics_summary.v1', validation='passed', counts=audit['counts'],
        wt_barcodes=len(wt_idx), wt_molecules=int(n[wt_idx].sum()),
        variant_coverage_gate_passed=sum(b >= 50 and m >= 100 for b, m in coverage.values()),
        hypotheses=len(effects), statuses=dict(Counter(r['status'] for r in effects)),
        candidate_q_le_05=len(candidate), candidate_q_variants=len({r['variant_id'] for r in candidate}),
        candidate_by_track=dict(Counter(r['track'] for r in candidate)), supported_effects=len(supported),
        calibration_status=calibration['status'],
        calibration_by_kind={kind: dict(Counter(s['status'] for s in calibration['strata'] if s['kind'] == kind))
                             for kind in ('wt_aa', 'simulation')},
        failed_strata=[s for s in calibration['strata'] if s['status'] == 'fail'],
        unresolved_null_simulations=calibration['unresolved_null_simulations'],
        profiles=len(profiles), sensitivity_rows=len(sensitivities),
        interpretation='Candidate q-value counts are diagnostics, not released discoveries when calibration fails or remains unresolved.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = summarize(args.run)
    with open(args.output, 'x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'failed_strata'}))

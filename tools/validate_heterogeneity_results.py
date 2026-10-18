#!/usr/bin/env python3
"""Independent reconciliation to the earlier per-base barcode aggregates."""
import argparse
import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from itertools import zip_longest
from pathlib import Path

import numpy as np


def rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def validate(run, previous):
    extraction = json.loads((run / 'features/audit.json').read_text())
    summary = json.loads((run / 'summary/summary.json').read_text())
    assert extraction['complete'] and summary['complete']
    counts = np.load(previous / 'barcode_counts.npy', mmap_mode='r')
    depths = np.load(previous / 'barcode_molecules.npy', mmap_mode='r')
    assert counts.shape == (568354, 2, 211)
    total = 0
    largest_error = 0.0
    weighted = Counter()
    wt_spread = defaultdict(Counter)
    for i, (old, new) in enumerate(zip_longest(
            rows(previous / 'barcode_metadata.tsv.gz'),
            rows(run / 'summary/barcode_heterogeneity.tsv.gz'))):
        assert old is not None and new is not None
        assert int(old['barcode_index']) == i
        assert old['barcode'] == new['corrected_barcode']
        assert old['variant_id'] == new['variant_id']
        depth = int(new['molecules'])
        assert depth == int(old['molecules']) == int(depths[i])
        for track, name in enumerate(('element_tf_small_fraction', 'element_msp_fraction')):
            expected = counts[i, track].sum(dtype=np.uint64) / (depth * 211)
            observed = float(new[name + '_mean'])
            error = abs(observed - expected)
            largest_error = max(largest_error, error)
            assert error < 1e-12
            assert int(new[name + '_valid_molecules']) == depth
            assert (new[name + '_within_sd'] == '') == (depth < 2)
        for region in ('whole', 'element', 'distal'):
            lo = new[region + '_m6a_fraction_204_mean']
            hi = new[region + '_m6a_fraction_mean']
            if lo and hi:
                assert 0 <= float(lo) <= float(hi) + 1e-12 <= 1 + 1e-12
        if new['source_variant'] == 'WT':
            for feature in ('whole_nucleosome_segments', 'distal_m6a_fraction',
                            'element_tf_small_fraction', 'element_msp_fraction'):
                if int(new[feature + '_valid_molecules']) >= 2:
                    wt_spread[feature]['repeated_barcodes'] += 1
                    wt_spread[feature]['observed_sd_above_1e_minus_6'] += int(
                        float(new[feature + '_within_sd']) > 1e-6)
        total += depth
        weighted[new['variant_id']] += 1
    assert i + 1 == summary['primary_barcodes'] == 568354
    assert total == summary['primary_molecules_from_barcode_depths'] == 2102906

    totals = Counter()
    seen = set()
    for row in rows(run / 'summary/variant_heterogeneity.tsv.gz'):
        key = (row['cohort'], row['variant_id'])
        if key in seen:
            continue
        seen.add(key)
        totals[(row['cohort'], row['genotype_group'], 'molecule')] += int(row['molecules'])
        if row['barcodes']:
            assert int(row['barcodes']) == weighted[row['variant_id']]
            totals[(row['cohort'], row['genotype_group'], 'equal_barcode')] += int(row['barcodes'])
    checked = {}
    for filename, dimensions in (
            ('histograms.tsv.gz', ('feature',)),
            ('joint_histograms.tsv.gz', ('x_feature', 'y_feature')),
            ('count_histograms.tsv.gz', ('feature',))):
        sums = defaultdict(float)
        for row in rows(run / 'summary' / filename):
            prefix = (row['cohort'], row['genotype_group'], row['weighting'])
            key = prefix + tuple(row[field] for field in dimensions)
            sums[key] += float(row['count_weight']) + float(row['missing_weight'])
        for key, value in sums.items():
            assert math.isclose(value, totals[key[:3]], abs_tol=1e-5, rel_tol=1e-10), (filename, key, value)
        checked[filename] = len(sums)
    retained = Counter()
    for row in rows(run / 'features/retention.tsv'):
        for name in ('total_primary', 'sequence_qualified', 'footprint_analysis_eligible',
                     'footprint_analysis_eligible_with_barcode'):
            retained[name] += int(row[name])
    assert retained['total_primary'] == 8715827
    assert retained['sequence_qualified'] == 2660013
    assert retained['footprint_analysis_eligible'] == 2594368
    assert retained['footprint_analysis_eligible_with_barcode'] == 2102906
    return {'status': 'passed', 'barcodes_compared': i + 1, 'molecules_reconciled': total,
            'max_barcode_mean_error_against_previous_per_base_arrays': largest_error,
            'histogram_groups_conserved': checked, 'retention': dict(retained),
            'WT_repeated_barcode_variation': {key: dict(value) for key, value in wt_spread.items()}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--previous', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    report = validate(args.run, args.previous)
    with args.output.open('x') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    print(json.dumps(report))

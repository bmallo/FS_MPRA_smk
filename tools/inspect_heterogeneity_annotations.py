#!/usr/bin/env python3
"""Audit long protected annotations and retained m6A probability bytes.

These descriptions do not classify compartments or create eligibility filters.
The probability-byte audit is a bounded, nonrandom BAM sample.
"""
import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pysam


def inspect(features, previous, bam_path, wt_id, limit):
    with gzip.open(previous / 'barcode_metadata.tsv.gz', 'rt') as handle:
        depths = {r['barcode']: int(r['molecules']) for r in csv.DictReader(handle, delimiter='\t')
                  if r['variant_id'] == wt_id}
    patterns = defaultdict(Counter)
    tails = Counter()
    coverage = defaultdict(Counter)
    n = 0
    with gzip.open(features, 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            cohorts = ['sequence_qualified']
            cohorts += [flag for flag in ('footprint_analysis_eligible',
                                          'footprint_analysis_eligible_with_barcode') if row[flag] == '1']
            for cohort in cohorts:
                for region in ('whole', 'element', 'distal'):
                    coverage[(cohort, region)][float(row[region + '_coverage_fraction'])] += 1
            if row['genotype_class'] != 'WT' or row['footprint_analysis_eligible_with_barcode'] != '1':
                continue
            n += 1
            weight = 1 / depths[row['corrected_barcode']]
            calls = int(row['whole_nucleosome_segments'])
            length = float(row['whole_nucleosome_length_median']) if row['whole_nucleosome_length_median'] else 0
            assert calls > 0 or length == 0
            pattern = ('no_nucleosome_calls' if calls == 0 else
                       'one_call_longer_than_1000bp' if calls == 1 and length > 1000 else
                       'other_called_patterns')
            bucket = patterns[pattern]
            bucket['molecules'] += 1
            bucket['weight'] += weight
            for feature in ('distal_m6a_fraction', 'distal_nucleosome_fraction',
                            'element_tf_small_fraction', 'element_msp_fraction'):
                bucket[feature] += weight * float(row[feature])
            bucket['any_element_small_footprint'] += weight * int(float(row['element_tf_small_fraction']) > 0)
            for cutoff in (300, 500, 1000):
                if length > cutoff:
                    tails[str(cutoff)] += weight
    assert n == sum(depths.values())
    assert abs(sum(v['weight'] for v in patterns.values()) - len(depths)) < 1e-5
    coverage_report = []
    for (cohort, region), distribution in sorted(coverage.items()):
        total = sum(distribution.values())
        quantiles = {}
        for label, probability in (('q01', .01), ('q10', .1), ('median', .5), ('q90', .9)):
            cumulative = 0
            for value, count in sorted(distribution.items()):
                cumulative += count
                if cumulative >= probability * total:
                    quantiles[label] = value
                    break
        coverage_report.append(dict(cohort=cohort, region=region, molecules=total,
            fraction_at_least_95pct=sum(n for v, n in distribution.items() if v >= .95) / total,
            fraction_at_least_99pct=sum(n for v, n in distribution.items() if v >= .99) / total,
            **quantiles))
    values = Counter()
    examined = 0
    with pysam.AlignmentFile(str(bam_path), 'rb', threads=1) as bam:
        for read in bam.fetch(until_eof=True):
            if examined == limit:
                break
            examined += 1
            if not read.has_tag('MM') or not read.has_tag('ML'):
                continue
            ml = read.get_tag('ML')
            cursor = 0
            for group in read.get_tag('MM').split(';'):
                if not group:
                    continue
                header, *deltas = group.split(',')
                match = re.fullmatch(r'([ACGTUN])([+-])([a-z]+|[0-9]+)([.?]?)', header)
                assert match is not None
                base, strand, codes, omitted = match.groups()
                codes = [codes] if codes.isdigit() else list(codes)
                end = cursor + len(deltas) * len(codes)
                assert end <= len(ml)
                for j, code in enumerate(codes):
                    if (base, strand, code) in (('A', '+', 'a'), ('T', '-', 'a')):
                        values.update(ml[cursor + j:end:len(codes)])
                cursor = end
            assert cursor == len(ml)
    return {
        'WT_molecules': n, 'WT_barcodes': len(depths),
        'regional_alignment_coverage': coverage_report,
        'weighting': 'Each eligible supported WT barcode contributes total weight one. Length cutoffs are descriptive only.',
        'median_nucleosome_call_length_above_threshold_proportion_all_WT':
            {k: v / len(depths) for k, v in tails.items()},
        'patterns': {key: {'molecules': v['molecules'], 'weight': v['weight'],
                          'proportion': v['weight'] / len(depths),
                          'means': {f: value / v['weight'] for f, value in v.items()
                                    if f not in ('molecules', 'weight')}} for key, v in patterns.items()},
        'm6a_probability_pilot': {'records_examined': examined, 'selection': 'First BAM records, nonrandom',
                                 'byte_counts': dict(sorted(values.items())),
                                 'minimum': min(values) if values else None,
                                 'maximum': max(values) if values else None}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--features', required=True, type=Path)
    parser.add_argument('--previous', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--bam-records', type=int, default=1000)
    args = parser.parse_args()
    if args.bam_records < 1:
        raise ValueError('BAM record limit must be positive')
    cfg = json.loads(args.config.read_text())
    catalog = json.loads((Path(cfg['ledger_run_dir']) / 'catalog/catalog.json').read_text())
    wt_ids = [v['variant_id'] for v in catalog['library']['variants'] if v['class'] == 'WT']
    assert len(wt_ids) == 1
    bam_path = Path(json.loads(Path(cfg['footprint_contract']).read_text())['bam'])
    result = inspect(args.features, args.previous, bam_path, wt_ids[0], args.bam_records)
    result['sources'] = [str(path.resolve()) for path in (args.config, args.features, args.previous, bam_path)]
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result))

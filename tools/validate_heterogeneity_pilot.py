#!/usr/bin/env python3
"""Independently decode regional methylation for a bounded real-read pilot."""
import argparse
import csv
import gzip
import json
import re
from collections import Counter

import pysam


def reference_calls(read):
    forward = read.get_forward_sequence()
    ml = list(read.get_tag('ML'))
    cursor = 0
    observed = {128: set(), 204: set()}
    groups = []
    for section in read.get_tag('MM').strip(';').split(';'):
        header, *deltas = section.split(',')
        match = re.fullmatch(r'([ACGTUN])([+-])([a-z]+|[0-9]+)([.?]?)', header)
        assert match is not None
        base, strand, codes, omission = match.groups()
        codes = [codes] if codes.isdigit() else list(codes)
        positions = [i for i, b in enumerate(forward) if b == base or base == 'N']
        at_index = -1
        for delta in deltas:
            at_index += int(delta) + 1
            assert 0 <= at_index < len(positions)
            original_q = positions[at_index]
            stored_q = len(forward) - 1 - original_q if read.is_reverse else original_q
            for j, code in enumerate(codes):
                if (base, strand, code) in (('A', '+', 'a'), ('T', '-', 'a')):
                    for threshold in observed:
                        if ml[cursor+j] >= threshold:
                            observed[threshold].add(stored_q)
            cursor += len(codes)
        groups.append(header)
    assert cursor == len(ml)
    return observed, groups


def validate(config, features, limit):
    if limit < 1:
        raise ValueError('Pilot validation limit must be positive')
    spec = json.load(open(config))
    with gzip.open(features, 'rt') as handle:
        rows = {}
        for row in csv.DictReader(handle, delimiter='\t'):
            rows[row['read_id']] = row
            if len(rows) >= limit:
                break
    contract = json.load(open(spec['footprint_contract']))
    catalog = json.load(open(spec['ledger_run_dir'] + '/catalog/catalog.json'))['library']
    tally = Counter()
    groups = Counter()
    pending = set(rows)
    if not pending:
        raise ValueError('Pilot feature table is empty')
    with pysam.AlignmentFile(contract['bam'], 'rb', threads=1) as bam:
        for read in bam.fetch(catalog['contig']):
            if read.is_secondary or read.is_supplementary or read.query_name not in pending:
                continue
            pending.remove(read.query_name)
            row = rows[read.query_name]
            if row['m6a_status'] != 'ok':
                tally['skipped_unknown'] += 1
                if not pending:
                    break
                continue
            observed, headers = reference_calls(read)
            groups.update(headers)
            qref = read.get_reference_positions(full_length=True)
            assert len(observed[128]) == int(row['m6a_count'])
            for name, intervals in spec['regions'].items():
                selected = {q for q, r in enumerate(qref) if r is not None
                            and any(a <= r < b for a, b in intervals)}
                at = sum(read.query_sequence[q] in 'AT' for q in selected)
                assert len(selected) == int(row[name+'_aligned_bases'])
                assert at == int(row[name+'_at_count'])
                for threshold in (128, 204):
                    suffix = '' if threshold == 128 else '_204'
                    numerator = len(selected & observed[threshold])
                    assert numerator == int(row[name+'_m6a_count'+suffix])
                    value = row[name+'_m6a_fraction'+suffix]
                    assert (not value and not at) or abs(float(value) - numerator/at) < 1e-12
            tally['reverse' if read.is_reverse else 'forward'] += 1
            tally['checked'] += 1
            if not pending:
                break
    assert not pending
    return dict(status='passed', counts=dict(tally), mm_headers=dict(groups),
                interpretation='Independent explicit MM delta decoding and pysam CIGAR reference-map check; bounded real-read pilot, not sampling calibration.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--features', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--limit', type=int, default=1000)
    args = parser.parse_args()
    result = validate(args.config, args.features, args.limit)
    with open(args.output, 'x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result))

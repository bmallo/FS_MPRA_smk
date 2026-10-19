"""Descriptive full-cohort call landscape; no significance inference."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np

from .io import atomic_output, open_text, read_json, source_info, write_json, write_table
from .processing import check_sources

TRACKS = ('tf_small', 'tf_other', 'nuc_short', 'nuc_85_250',
          'nuc_251_500', 'nuc_501_1000', 'nuc_gt1000', 'msp')
START, END = 1234, 1445
LENGTH = END - START


def category(kind, length):
    if type(length) is not int or length < 1:
        raise ValueError('Invalid original query length')
    if kind == 'tf':
        return 0 if 10 <= length <= 60 else 1
    if kind == 'nuc':
        return 2 if length < 85 else 3 if length <= 250 else 4 if length <= 500 else 5 if length <= 1000 else 6
    if kind == 'msp':
        return 7
    raise ValueError('Unknown annotation kind')


def union(intervals):
    result = []
    for a, b in sorted(intervals):
        if type(a) is not int or type(b) is not int or not START <= a < b <= END:
            raise ValueError('Call interval outside element')
        if result and a <= result[-1][1]:
            result[-1][1] = max(b, result[-1][1])
        else:
            result.append([a, b])
    return result


def decode_calls(row):
    groups = [[] for _ in TRACKS]
    original = []
    for kind in ('tf', 'nuc', 'msp'):
        raw = row[kind + '_ref_calls0']
        if not raw:
            raise ValueError('Missing annotations for evaluated molecule')
        calls = json.loads(raw)
        if not isinstance(calls, list):
            raise ValueError('Annotations must be a list')
        for call in calls:
            length = call['query_length']
            track = category(kind, length)
            intervals = union(call['ref_intervals0'])
            if not intervals:
                raise ValueError('Empty original call projection')
            groups[track].extend(intervals)
            original.append((kind, length, intervals))
    return [union(g) for g in groups], original


def state_code(groups):
    return int(bool(groups[0])) + 2 * int(bool(groups[3])) + 4 * int(any(groups[4:7]))


def reconcile(counts, depths, variant_indices, expected):
    sums = np.zeros_like(expected)
    for first in range(0, len(depths), 4096):
        last = min(first + 4096, len(depths))
        block = np.asarray(counts[first:last], dtype=float) / depths[first:last, None, None]
        np.add.at(sums, variant_indices[first:last], block)
    error = float(np.max(np.abs(sums - expected)))
    # Accumulation order differs: raw sums can reach 268,283 barcode weights.
    if not np.allclose(sums, expected, rtol=1e-10, atol=1e-7):
        raise ValueError('Landscape differs from original TF/MSP barcode aggregates')
    return error


def run(config, output=None):
    config = Path(config)
    cfg_source = source_info(config, True)
    cfg = read_json(config)
    if cfg.get('schema_version') != 'fsmpra.landscape_config.v1' or cfg.get('sample_id') != 'HepG2_rep3':
        raise ValueError('Unexpected landscape configuration')
    agg, joined = Path(cfg['aggregates']), Path(cfg['joined_dir'])
    catalog = Path(cfg['catalog'])
    output = Path(output or Path(cfg['output_dir']) / 'landscape')
    paths = [(config, True), (catalog, True), (agg/'audit.json', True),
             (agg/'barcode_metadata.tsv.gz', True), (agg/'variant_metadata.tsv', True),
             (agg/'barcode_counts.npy', False), (agg/'barcode_molecules.npy', False),
             (joined/'audit.json', True), (joined/'molecule_footprints.tsv.gz', False),
             (Path(__file__), True), (Path(__file__).with_name('io.py'), True),
             (Path(__file__).with_name('processing.py'), True)]
    sources = [source_info(p, checksum) for p, checksum in paths]
    if sources[0] != cfg_source:
        raise ValueError('Configuration changed while loading')
    aa, ja = read_json(agg/'audit.json'), read_json(joined/'audit.json')
    if (aa.get('schema_version') != 'fsmpra.statistics_aggregate.v1' or
            aa.get('interval') != [START, END] or ja.get('interval') != [START, END] or
            aa['counts'].get('barcodes') != 568354 or
            aa['counts'].get('footprint_analysis_eligible_with_barcode') != 2102906 or
            ja['counts'].get('output_molecules') != 8715827):
        raise ValueError('Unexpected aggregate or joined cohort')
    recorded = {s['path']: s for s in aa['sources']}
    for p, checksum in [(catalog, True), (joined/'audit.json', True), (joined/'molecule_footprints.tsv.gz', False)]:
        current = source_info(p, checksum)
        if recorded.get(current['path']) != current:
            raise ValueError('Aggregate and joined provenance disagree')
    library = read_json(catalog)['library']
    if (library['mutable_ref_start0'], library['mutable_ref_end0']) != (START, END):
        raise ValueError('Catalog interval differs')
    with open(agg/'variant_metadata.tsv') as h:
        variants = [v for v in csv.DictReader(h, delimiter='\t') if v['genotype_class'] in ('WT', 'SNV')]
    if len(variants) != 634:
        raise ValueError('Expected WT and 633 SNVs')
    vidx = {v['variant_id']: i for i, v in enumerate(variants)}
    wt_indices = [i for i, v in enumerate(variants) if v['genotype_class'] == 'WT']
    if len(vidx) != len(variants) or len(wt_indices) != 1:
        raise ValueError('Variant identities invalid')
    wi = wt_indices[0]
    depths = np.load(agg/'barcode_molecules.npy', mmap_mode='r', allow_pickle=False)
    counts = np.load(agg/'barcode_counts.npy', mmap_mode='r', allow_pickle=False)
    if depths.shape != (568354,) or counts.shape != (568354, 2, LENGTH) or np.any(depths <= 0):
        raise ValueError('Unexpected barcode arrays')
    bmap, bvariant = {}, np.empty(len(depths), dtype=np.int32)
    with open_text(agg/'barcode_metadata.tsv.gz') as h:
        for j, row in enumerate(csv.DictReader(h, delimiter='\t')):
            if int(row['barcode_index']) != j or int(row['molecules']) != depths[j] or row['barcode'] in bmap:
                raise ValueError('Barcode metadata mismatch')
            bmap[row['barcode']] = j
            bvariant[j] = vidx[row['variant_id']]
    if len(bmap) != len(depths):
        raise ValueError('Incomplete barcode metadata')
    nbarcodes = np.bincount(bvariant, minlength=len(variants))
    nmolecules = np.bincount(bvariant, weights=depths, minlength=len(variants)).astype(np.int64)
    if nbarcodes[wi] != 268283 or nmolecules[wi] != 985264 or np.any(nbarcodes == 0):
        raise ValueError('WT or variant cohort mismatch')
    diff = np.zeros((len(variants), len(TRACKS), LENGTH+1), dtype=float)
    mdiff = np.zeros_like(diff, dtype=np.int64)
    presence, mpresence = np.zeros((len(variants), len(TRACKS))), np.zeros((len(variants), len(TRACKS)), dtype=np.int64)
    states, mstates = np.zeros((len(variants), 8)), np.zeros((len(variants), 8), dtype=np.int64)
    seen = np.zeros(len(depths), dtype=np.int64)
    read_ids = set()
    lengths, weighted_lengths = Counter(), Counter()
    boundaries, weighted_boundaries = Counter(), Counter()
    wt_multi_interval_tf = 0
    total = eligible = 0
    with open_text(joined/'molecule_footprints.tsv.gz') as h:
        for row in csv.DictReader(h, delimiter='\t'):
            total += 1
            if total % 1000000 == 0:
                print(f'joined rows {total:,}; eligible {eligible:,}', flush=True)
            if row['footprint_analysis_eligible_with_barcode'] != '1':
                continue
            if (row['sample_id'] != cfg['sample_id'] or row['footprint_status'] != 'evaluated' or
                    row['region_complete'] != '1' or int(row['aligned_region_bases']) != LENGTH or
                    union(json.loads(row['eligible_ref_intervals0'])) != [[START, END]]):
                raise ValueError('Eligible molecule not fully evaluated/aligned')
            bi = bmap[row['corrected_barcode']]
            vi = bvariant[bi]
            if row['assigned_variant_id'] != variants[vi]['variant_id'] or row['observed_variant_id'] != variants[vi]['variant_id']:
                raise ValueError('Molecule/barcode genotype mismatch')
            if row['read_id'] in read_ids:
                raise ValueError('Duplicate eligible molecule')
            read_ids.add(row['read_id'])
            seen[bi] += 1
            if seen[bi] > depths[bi]:
                raise ValueError('Barcode depth exceeded')
            eligible += 1
            weight = 1.0 / depths[bi]
            groups, original = decode_calls(row)
            for ti, intervals in enumerate(groups):
                if intervals:
                    presence[vi, ti] += weight
                    mpresence[vi, ti] += 1
                for a, b in intervals:
                    a, b = a-START, b-START
                    diff[vi, ti, a] += weight
                    diff[vi, ti, b] -= weight
                    mdiff[vi, ti, a] += 1
                    mdiff[vi, ti, b] -= 1
            state = state_code(groups)
            states[vi, state] += weight
            mstates[vi, state] += 1
            if vi == wi:
                for kind, length, intervals in original:
                    lengths[kind, length] += 1
                    weighted_lengths[kind, length] += weight
                    if kind == 'tf' and 10 <= length <= 60:
                        if len(intervals) == 1:
                            a, b = intervals[0]
                            boundaries[a-START, b-START] += 1
                            weighted_boundaries[a-START, b-START] += weight
                        else:
                            wt_multi_interval_tf += 1
    if total != 8715827 or eligible != 2102906 or not np.array_equal(seen, depths):
        raise ValueError('Incomplete joined stream or barcode conservation failure')
    coverage_sum = np.cumsum(diff, axis=2)[:, :, :LENGTH]
    mcoverage_sum = np.cumsum(mdiff, axis=2)[:, :, :LENGTH]
    raw_error = reconcile(counts, depths, bvariant, coverage_sum[:, [0, 7], :])
    coverage = coverage_sum / nbarcodes[:, None, None]
    mcoverage = mcoverage_sum / nmolecules[:, None, None]
    if np.min(coverage) < -1e-10 or np.max(coverage) > 1+1e-10 or np.any(mcoverage_sum < 0):
        raise ValueError('Invalid coverage')
    coverage = np.clip(coverage, 0, 1)
    p = presence / nbarcodes[:, None]
    if np.any(coverage.mean(axis=2) > p + 1e-10) or np.any(p > 1+1e-10):
        raise ValueError('Coverage/presence conservation failure')
    state_fractions = states / nbarcodes[:, None]
    if not np.allclose(state_fractions.sum(axis=1), 1, atol=1e-10, rtol=0) or not np.array_equal(mstates.sum(axis=1), nmolecules):
        raise ValueError('State conservation failure')
    audit = dict(schema_version='fsmpra.landscape.v1', complete=True, descriptive_only=True,
                 sample_id=cfg['sample_id'], interval=[START, END], tracks=list(TRACKS),
                 counts=dict(joined_rows=total, eligible_molecules=eligible, barcodes=len(depths),
                             wt_barcodes=int(nbarcodes[wi]), wt_molecules=int(nmolecules[wi]), variants=len(variants)),
                 reconciliation_max_raw_weight_error=raw_error,
                 reconciliation_fraction_error_upper_bound=raw_error / int(nbarcodes.min()),
                 wt_multi_interval_small_tf_calls=wt_multi_interval_tf, sources=sources,
                 limitations=['Observed annotation/length categories, not protein identities.',
                              'Primary interval is the 211-bp element; original lengths can extend outside it.',
                              'Categories overlap. States describe co-occurrence anywhere within element.',
                              'No sequence-bias correction or variant confidence audit in this stage.'])
    check_sources(sources)
    with atomic_output(output) as out:
        np.savez_compressed(out/'landscape.npz', variant_ids=np.array([v['variant_id'] for v in variants]),
                            source_variants=np.array([v['source_variant'] for v in variants]),
                            tracks=np.array(TRACKS), coverage=coverage, molecule_coverage=mcoverage,
                            presence=p, molecule_presence=mpresence/nmolecules[:, None],
                            states=state_fractions, molecule_states=mstates/nmolecules[:, None],
                            n_barcodes=nbarcodes, n_molecules=nmolecules)
        write_table(out/'variants.tsv', variants, list(variants[0]))
        def profile_rows():
            for vi, v in enumerate(variants):
                for ti, track in enumerate(TRACKS):
                    for pos in range(LENGTH):
                        yield dict(variant_id=v['variant_id'], source_variant=v['source_variant'], track=track,
                                   relative_pos0=pos, reference_pos0=START+pos,
                                   mean_equal_barcode=float(coverage[vi, ti, pos]),
                                   mean_molecule_weighted=float(mcoverage[vi, ti, pos]),
                                   n_barcodes=int(nbarcodes[vi]), n_molecules=int(nmolecules[vi]))
        fields=['variant_id','source_variant','track','relative_pos0','reference_pos0',
                'mean_equal_barcode','mean_molecule_weighted','n_barcodes','n_molecules']
        write_table(out/'profiles.tsv.gz', profile_rows(), fields)
        categories, state_rows = [], []
        for vi, v in enumerate(variants):
            base=dict(variant_id=v['variant_id'], source_variant=v['source_variant'],
                      n_barcodes=int(nbarcodes[vi]), n_molecules=int(nmolecules[vi]))
            for ti, track in enumerate(TRACKS):
                categories.append(dict(base, track=track, mean_coverage=float(coverage[vi,ti].mean()),
                                       any_call_fraction=float(p[vi,ti]),
                                       molecule_mean_coverage=float(mcoverage[vi,ti].mean()),
                                       molecule_any_call_fraction=float(mpresence[vi,ti]/nmolecules[vi])))
            for state in range(8):
                state_rows.append(dict(base, state_code=state, tf_small=bool(state&1), nuc_85_250=bool(state&2),
                                       nuc_gt250=bool(state&4), fraction=float(state_fractions[vi,state]),
                                       molecule_fraction=float(mstates[vi,state]/nmolecules[vi]),
                                       state_molecules=int(mstates[vi,state])))
        write_table(out/'categories.tsv', categories, list(categories[0]))
        write_table(out/'states.tsv', state_rows, list(state_rows[0]))
        rows=[dict(caller_track=k[0], original_query_length=k[1], calls=lengths[k],
                   weighted_calls_per_barcode=weighted_lengths[k]/nbarcodes[wi]) for k in sorted(lengths)]
        write_table(out/'wt_call_lengths.tsv', rows, ['caller_track','original_query_length','calls','weighted_calls_per_barcode'])
        rows=[dict(relative_start0=k[0],relative_end0=k[1],calls=boundaries[k],
                   weighted_calls_per_barcode=weighted_boundaries[k]/nbarcodes[wi]) for k in sorted(boundaries)]
        write_table(out/'wt_tf_boundaries.tsv', rows, ['relative_start0','relative_end0','calls','weighted_calls_per_barcode'])
        write_json(out/'audit.json', audit)
        check_sources(sources)
    return audit


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--output')
    args=parser.parse_args()
    print(json.dumps(run(args.config,args.output)['counts']))


if __name__ == '__main__':
    main()

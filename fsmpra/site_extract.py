"""Library-wide per-site retained m6A and annotation extraction for replicate 3.

Every sequence-qualified WT/SNV molecule that spans the element with known m6A
status is projected onto the 211-bp element. Arrays are memory-mappable .npy
files so analysis never rereads the BAM. Proxy-control classes are assigned
later, from distal methylation only. Scores are retained input calls (stored
MM/ML), not original kinetics. Never calls pysam.modified_bases.
"""
import argparse
import csv
import time
from pathlib import Path
import numpy as np
from .io import atomic_output, open_text, read_json, source_info, write_json
from .processing import check_sources
from .raw_candidate_extract import project_read

ACCEPTED = ('supported', 'corrected_low_quality_substitution')
META_FIELDS = ['read_id', 'observed_variant_id', 'corrected_barcode', 'barcode_assignment_status', 'passes',
               'read_quality', 'mapq', 'read_length', 'strand', 'circular', 'footprint_status',
               'footprint_analysis_eligible', 'footprint_analysis_eligible_with_barcode', 'whole_m6a_fraction',
               'distal_m6a_fraction', 'distal_nucleosome_fraction', 'whole_nucleosome_segments',
               'element_msp_fraction', 'element_tf_small_fraction']
ARRAYS = [('scores', np.uint8), ('known', np.bool_), ('qualities', np.uint8), ('coverage', np.uint8)]


def number(value, default=np.nan):
    return float(value) if value not in ('', None) else default


def selected(row, variants):
    return row['observed_variant_id'] in variants and row['region_complete'] == '1' and row['m6a_status'] == 'ok'


def scan_features(path, variants, limit=None):
    """First pass: select molecules and collect compact numeric metadata in file order."""
    cols = {k: [] for k in ['variant_index', 'barcode', 'supported', 'eligible', 'eligible_with_barcode', 'circular',
                            'minus', 'evaluated', 'passes', 'read_quality', 'mapq', 'whole_m6a', 'distal_m6a',
                            'distal_nucleosome', 'nucleosome_segments', 'element_msp', 'element_tf_small']}
    barcodes = {}
    with open_text(path) as h:
        for r in csv.DictReader(h, delimiter='\t'):
            if not selected(r, variants):
                continue
            bc = r['corrected_barcode']
            cols['variant_index'].append(variants[r['observed_variant_id']])
            cols['barcode'].append(barcodes.setdefault(bc, len(barcodes)) if bc else -1)
            cols['supported'].append(r['barcode_assignment_status'] in ACCEPTED)
            for key, field in [('eligible', 'footprint_analysis_eligible'), ('eligible_with_barcode', 'footprint_analysis_eligible_with_barcode'),
                               ('circular', 'circular')]:
                cols[key].append(r[field] == '1')
            cols['minus'].append(r['strand'] == '-')
            cols['evaluated'].append(r['footprint_status'] == 'evaluated')
            for key, field in [('passes', 'passes'), ('read_quality', 'read_quality'), ('mapq', 'mapq'),
                               ('whole_m6a', 'whole_m6a_fraction'), ('distal_m6a', 'distal_m6a_fraction'),
                               ('distal_nucleosome', 'distal_nucleosome_fraction'), ('nucleosome_segments', 'whole_nucleosome_segments'),
                               ('element_msp', 'element_msp_fraction'), ('element_tf_small', 'element_tf_small_fraction')]:
                cols[key].append(number(r[field]))
            if limit and len(cols['variant_index']) >= limit:
                break
    dtypes = dict(variant_index=np.int16, barcode=np.int32, supported=bool, eligible=bool, eligible_with_barcode=bool,
                  circular=bool, minus=bool, evaluated=bool, passes=np.float32, read_quality=np.float32, mapq=np.float32,
                  whole_m6a=np.float32, distal_m6a=np.float32, distal_nucleosome=np.float32, nucleosome_segments=np.float32,
                  element_msp=np.float32, element_tf_small=np.float32)
    meta = {k: np.asarray(v, dtype=dtypes[k]) for k, v in cols.items()}
    if not len(meta['variant_index']):
        raise ValueError('No molecules selected')
    # Molecules without a barcode are their own resampling cluster.
    cluster = meta['barcode'].astype(np.int64)
    lone = cluster < 0
    cluster[lone] = len(barcodes) + np.arange(lone.sum())
    meta['cluster'] = cluster
    return meta, len(barcodes)


def run(config, output, limit=None, progress_every=500000):
    import pysam
    started = time.monotonic()
    cfg = read_json(config)
    if cfg.get('schema_version') != 'fsmpra.site_methylation_config.v1':
        raise ValueError('Unexpected site-methylation configuration')
    feature_path, feature_audit = Path(cfg['features']), read_json(cfg['feature_audit'])
    check_sources(feature_audit['sources'])
    lib = read_json(cfg['catalog'])['library']
    variants = [v for v in lib['variants'] if v['class'] in ('WT', 'SNV')]
    if [v['class'] for v in variants].count('WT') != 1 or len(variants) != 634:
        raise ValueError('Expected WT plus 633 SNVs')
    byid = {v['variant_id']: i for i, v in enumerate(variants)}
    bam_path = Path(read_json(cfg['footprint_contract'])['bam'])
    frozen = next(s for s in feature_audit['sources'] if s['path'] == str(bam_path.resolve()))
    if source_info(bam_path) != frozen:
        raise ValueError('Annotation BAM differs from frozen feature extraction source')
    sources = [source_info(p, p.suffix not in ('.bam', '.bai', '.gz')) for p in
               [Path(config), Path(cfg['catalog']), Path(cfg['feature_audit']), feature_path, Path(cfg['footprint_contract']),
                bam_path, Path(str(bam_path) + '.bai')]]
    sources += [source_info(Path(__file__).with_name(n), True) for n in
                ['site_extract.py', 'raw_candidate_extract.py', 'heterogeneity.py', 'footprint_join.py', 'io.py', 'processing.py']]
    meta, n_barcodes = scan_features(feature_path, byid, limit)
    n = len(meta['variant_index'])
    print(f'Selected {n:,} molecules ({n_barcodes:,} barcodes); scanning BAM', flush=True)
    with atomic_output(output) as out:
        arrays = {name: np.lib.format.open_memmap(out / f'{name}.npy', mode='w+', dtype=dt, shape=(n, 211))
                  for name, dt in ARRAYS}
        hist = np.zeros(256, dtype=np.int64)
        k = total = 0
        with open_text(feature_path) as h, open_text(out / 'molecules.tsv.gz', 'wt') as dest, \
                pysam.AlignmentFile(str(bam_path), 'rb', threads=1) as bam:
            rows = (r for r in csv.DictReader(h, delimiter='\t') if selected(r, byid))
            writer = csv.DictWriter(dest, fieldnames=['row'] + META_FIELDS, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            pending = next(rows)
            for read in bam.fetch(lib['contig']):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                total += 1
                if read.query_name == pending['read_id']:
                    r = pending
                    if read.query_length != int(r['read_length']) or ('-' if read.is_reverse else '+') != r['strand']:
                        raise ValueError('BAM/features identity differs: ' + r['read_id'])
                    vi = byid[r['observed_variant_id']]
                    if vi != meta['variant_index'][k]:
                        raise ValueError('Feature passes disagree')
                    values = project_read(read, variants[vi]['sequence'])
                    for (name, _), value in zip(ARRAYS, values):
                        arrays[name][k] = value
                    hist += np.bincount(values[0][values[0] > 0], minlength=256)
                    if r['footprint_status'] == 'evaluated' and r['element_msp_fraction'] != '':
                        if abs(np.mean((values[3] & 2) > 0) - float(r['element_msp_fraction'])) > 1e-10:
                            raise ValueError('Footprint re-projection mismatch: ' + r['read_id'])
                    writer.writerow(dict(r, row=k))
                    k += 1
                    if k == n:
                        break
                    pending = next(rows)
                if total % progress_every == 0:
                    print(f'BAM {total:,}; extracted {k:,}/{n:,}', flush=True)
        if k != n:
            raise ValueError(f'BAM scan incomplete: {k:,}/{n:,} selected molecules found in BAM order')
        for a in arrays.values():
            a.flush()
        del arrays
        np.savez(out / 'meta.npz', **meta)
        write_json(out / 'variants.json', [{key: v.get(key) for key in
                   ('variant_id', 'source_variant', 'class', 'oligo_pos1', 'ref', 'alt', 'sequence')} for v in variants])
        check_sources(sources)
        write_json(out / 'audit.json', dict(
            schema_version='fsmpra.site_extract.v1', complete=limit is None, molecule_limit=limit, sources=sources,
            molecules=n, barcodes=n_barcodes, bam_primary_rows_scanned=total, retained_positive_ml_histogram=hist.tolist(),
            seconds=time.monotonic() - started, element_reference_interval0=[1234, 1445],
            selection='Sequence-qualified exact WT/SNV molecules from the heterogeneity feature table that span the element '
                      'with m6a_status ok. No barcode requirement; barcode and eligibility flags are kept for cohort choice.',
            arrays='scores: retained ML byte (0 = no retained call); known: A/T opportunity on the decoded strand; '
                   'qualities: base quality; coverage: bit flags 1 TF 10-60 bp, 2 MSP, 4 nucleosome 85-250 bp, 8 nucleosome >250 bp.',
            limitations='Retained MM/ML calls, not kinetics or original probabilities. No genotype, barcode or footprint re-calling.'))
    print(f'Complete: {n:,} molecules in {time.monotonic() - started:.1f}s', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--config', default='config/ldlr_site_methylation.json')
    p.add_argument('--output', required=True)
    p.add_argument('--limit-molecules', type=int, help='Smoke test: first N selected molecules only (audit marks incomplete)')
    a = p.parse_args()
    run(a.config, a.output, a.limit_molecules)


if __name__ == '__main__':
    main()

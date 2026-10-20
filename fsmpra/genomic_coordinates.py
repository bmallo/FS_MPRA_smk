"""Audit the LDLR atlas against indexed GRCh38 and a saved UCSC sequence.

This adapter deliberately accepts only the validated forward-strand, contiguous
mutable element. Plasmid backbone coordinates and indel VCF normalization are
outside its scope. All displayed genomic bases are 1-based; BED is half-open.
"""
import argparse
import csv
import hashlib
from pathlib import Path

from .io import atomic_output, read_fasta, read_json, source_info, write_json, write_table
from .processing import check_sources


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def verify_sequence(sequence, genomic, strand):
    if strand != '+':
        raise ValueError('Only forward-strand atlas coordinates are supported')
    if not sequence or sequence != genomic:
        raise ValueError('Construct does not exactly match the genomic forward strand')
    if sequence == reverse_complement(genomic):
        raise ValueError('Sequence orientation is ambiguous')


def interval_fields(mapping, start0, end0):
    """Convert an element-relative half-open interval without dropping its end."""
    if not 0 <= start0 < end0 <= mapping['length']:
        raise ValueError('Interval outside verified element')
    origin = mapping['genomic_start1'] - 1
    return dict(genomic_start0=origin + start0, genomic_end0=origin + end0,
                genomic_start1=origin + start0 + 1, genomic_end1=origin + end0)


def mapping_for(library):
    if library['genomic_strand'] != '+':
        raise ValueError('Only forward-strand atlas coordinates are supported')
    return dict(assembly=library['genomic_assembly'], chrom=library['genomic_chrom'],
                strand='+', genomic_start1=library['genomic_mutable_start1'],
                genomic_end1=library['genomic_mutable_start1'] + len(library['mutable_sequence']) - 1,
                length=len(library['mutable_sequence']),
                reference_start0=library['mutable_ref_start0'],
                oligo_start1=library['oligo_mutable_start1'])


def run(config, catalog, mpra_catalog, ucsc_sequence, legacy, output):
    import pysam
    cfg, cat, mpra = map(read_json, (config, catalog, mpra_catalog))
    lib = cat['library']
    mapping = mapping_for(lib)
    genome = Path(cfg['alignment']['reference_fasta'])
    paths = [config, catalog, mpra_catalog, ucsc_sequence, legacy, Path(__file__),
             Path(__file__).with_name('io.py'), Path(__file__).with_name('processing.py')]
    sources = [source_info(p, True) for p in paths]
    sources += [source_info(genome), source_info(str(genome) + '.fai', True)]
    sources += lib['sources'] + [a['source'] for a in mpra['mpra_audit']]
    check_sources(sources)
    for key in ('genomic_assembly', 'genomic_chrom', 'genomic_strand',
                'genomic_mutable_start1', 'mutable_ref_start0', 'mutable_ref_end0'):
        if cfg['library'][key] != lib[key]:
            raise ValueError('Processing configuration and catalog differ: ' + key)
    construct = read_fasta(lib['sources'][0]['path'], lib['contig'])
    wt = lib['mutable_sequence']
    if construct[lib['mutable_ref_start0']:lib['mutable_ref_end0']] != wt:
        raise ValueError('Physical construct reference differs from catalog')
    if hashlib.sha256(construct.encode()).hexdigest() != lib['reference_sequence_sha256']:
        raise ValueError('Physical construct checksum differs')
    ucsc = read_json(ucsc_sequence)
    if mapping['assembly'] != 'GRCh38' or ucsc['genome'] != 'hg38' or ucsc['chrom'] != mapping['chrom']:
        raise ValueError('Independent reference assembly/chromosome differs')
    if len(ucsc['dna']) != ucsc['end'] - ucsc['start']:
        raise ValueError('Incomplete independent reference response')
    start = mapping['genomic_start1'] - 1
    end = start + len(wt)
    if not ucsc['start'] <= start < end <= ucsc['end']:
        raise ValueError('Independent reference does not cover element')
    snvs = {v['variant_id']: v for v in lib['variants'] if v['class'] == 'SNV'}
    expected = {(mapping['chrom'], start+i+1, ref, alt)
                for i, ref in enumerate(wt) for alt in 'ACGT' if alt != ref}
    rows_by_source = {s['id']: [] for s in mpra['mpra_sources']}
    for row in mpra['mpra']:
        rows_by_source[row['mpra_source_id']].append(row)
    reports, raw_by_id = [], {}
    with pysam.FastaFile(str(genome)) as fa:
        sequence = fa.fetch(mapping['chrom'], start, end).upper()
        verify_sequence(wt, sequence, lib['genomic_strand'])
        if fa.fetch(ucsc['chrom'], ucsc['start'], ucsc['end']).upper() != ucsc['dna'].upper():
            raise ValueError('Local reference differs from independent UCSC reference')
        lo, hi = max(0, start-10000), end+10000
        neighborhood = fa.fetch(mapping['chrom'], lo, hi).upper()
        forward_hits = [lo+i+1 for i in range(len(neighborhood)-len(wt)+1) if neighborhood.startswith(wt, i)]
        reverse_hits = [lo+i+1 for i in range(len(neighborhood)-len(wt)+1) if neighborhood.startswith(reverse_complement(wt), i)]
        if forward_hits != [start+1] or reverse_hits:
            raise ValueError('Unexpected exact sequence matches in genomic neighborhood')
        for src in mpra['mpra_sources']:
            path = next(a['source']['path'] for a in mpra['mpra_audit'] if a['id'] == src['id'])
            with open(path) as handle:
                raw = list(csv.DictReader(handle, delimiter=src['delimiter']))
            imported = rows_by_source[src['id']]
            if len(imported) != len(raw):
                raise ValueError('MPRA import row count differs')
            observed, keyed = set(), {}
            for index, (r, joined) in enumerate(zip(raw, imported), 1):
                chrom = 'chr' + r['Chromosome'].removeprefix('chr')
                pos, ref, alt = int(r['Position']), r['Ref'].upper(), r['Alt'].upper()
                key = chrom, pos, ref, alt
                if key in keyed:
                    raise ValueError('Duplicate raw MPRA allele')
                keyed[key] = r
                if fa.fetch(chrom, pos-1, pos-1+len(ref)).upper() != ref:
                    raise ValueError('Raw MPRA reference disagrees with genome: ' + str(key))
                if (joined['chrom'], joined['genomic_pos1'], joined['ref'], joined['alt']) != key or joined['source_row'] != index:
                    raise ValueError('Imported MPRA coordinates/alleles differ from raw source')
                if joined['effect_log2'] != float(r['Value']) or joined['p_value'] != float(r['P-Value']):
                    raise ValueError('Imported MPRA measurements differ from raw source')
                if key in expected:
                    v = snvs.get(joined['variant_id'])
                    i = pos-start-1
                    if not v or v['sequence'] != wt[:i]+alt+wt[i+1:]:
                        raise ValueError('MPRA allele assigned to wrong construct sequence')
                    if (v['genomic_pos1'], v['oligo_pos1'], v['ref_pos0'], v['ref'], v['alt']) != (pos, i+mapping['oligo_start1'], i+mapping['reference_start0'], ref, alt):
                        raise ValueError('Design SNV coordinate/allele differs')
                    observed.add(key)
            if observed != expected or len(snvs) != len(expected):
                raise ValueError('Incomplete genomic SNV coverage')
            reports.append(dict(source=src['id'], raw_rows=len(raw), genomic_ref_matches=len(raw),
                                atlas_snvs=len(observed), allele_sequence_and_effect_matches=len(observed)))
            raw_by_id[src['id']] = keyed
    with open(legacy) as handle:
        old = list(csv.DictReader(handle, delimiter='\t'))
    for r in old:
        key = ('chr19', 11092732-int(r['position']), reverse_complement(r['ref']), reverse_complement(r['alt']))
        target = raw_by_id['published_LDLR_2'].get(key)
        if not target or float(r['score']) != float(target['Value']) or float(r['pvalue']) != float(target['P-Value']):
            raise ValueError('Legacy reverse-orientation mapping failed: ' + str(key))
    bases = [dict(assembly=mapping['assembly'], chrom=mapping['chrom'], strand='+',
                  relative_pos0=i, reference_pos0=mapping['reference_start0']+i,
                  oligo_pos1=mapping['oligo_start1']+i, genomic_pos1=start+i+1, ref=ref)
             for i, ref in enumerate(wt)]
    with atomic_output(output) as out:
        write_table(out/'base_coordinates.tsv', bases, list(bases[0]))
        (out/'element.bed').write_text(f"{mapping['chrom']}\t{start}\t{end}\t{lib['id']}\t0\t+\n")
        (out/'genomic_element.fa').write_text(f">{mapping['assembly']}:{mapping['chrom']}:{start+1}-{end}:+\n{sequence}\n")
        write_json(out/'ucsc_reference.json', ucsc)
        write_json(out/'audit.json', dict(schema_version='fsmpra.genomic_coordinates.v1', complete=True,
                   library_id=lib['id'], mapping=mapping, mutable_sequence=wt,
                   sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                   forward_exact_bases=len(wt), reverse_complement_exact=False,
                   local_search_interval0=[lo, hi], forward_matches1=forward_hits, reverse_matches1=reverse_hits,
                   independent_reference=dict(url=f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={ucsc['chrom']};start={ucsc['start']};end={ucsc['end']}",
                                              matched_bases=len(ucsc['dna']), download_time=ucsc.get('downloadTime')),
                   mpra=reports, legacy=dict(matched_rows=len(old), source='published_LDLR_2',
                                            position_formula='genomic_pos1 = 11092732 - legacy_position',
                                            alleles='complement both REF and ALT', effects='unchanged; do not negate'),
                   scope='Mutable element only; genomic forward alleles. Original plasmid coordinates retained. No indel VCF normalization.', sources=sources))
        check_sources(sources)
    return reports


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', default='config/ldlr_processing.json')
    p.add_argument('--catalog', default='results/rep3_ledger_20261016/catalog/catalog.json')
    p.add_argument('--mpra-catalog', default='results/pilot_20261014/catalog/catalog.json')
    p.add_argument('--ucsc-sequence', required=True)
    p.add_argument('--legacy', default='extras/data/LDLR_mpra_scores.tsv')
    p.add_argument('--output', required=True)
    args = p.parse_args()
    print(run(args.config, args.catalog, args.mpra_catalog, args.ucsc_sequence, args.legacy, args.output))

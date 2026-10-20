"""Check genomic exports against the preserved measurements and coordinate map."""
import argparse
import csv
from itertools import zip_longest
from pathlib import Path

from fsmpra.io import atomic_output, open_text, read_json, source_info, write_json
from fsmpra.processing import check_sources


def run(atlas, previous, landscape, output):
    atlas, previous, landscape = map(Path, (atlas, previous, landscape))
    paths = [atlas/'atlas.json', previous/'atlas.json', atlas/'audit.json',
             atlas/'genomic_profiles.tsv.gz', landscape/'profiles.tsv.gz',
             atlas/'region_effects.tsv.gz', previous/'region_effects.tsv.gz',
             atlas/'variant_summary.tsv', previous/'variant_summary.tsv',
             atlas/'wt_tf_boundaries_genomic.tsv', landscape/'wt_tf_boundaries.tsv',
             atlas/'functional_correlations.tsv',previous/'functional_correlations.tsv',
             atlas/'base_coordinates.tsv',atlas/'wt_regions.bed', Path(__file__)]
    sources = [source_info(p, True) for p in paths]
    check_sources(read_json(atlas/'audit.json')['sources'])
    current, old = read_json(atlas/'atlas.json'), read_json(previous/'atlas.json')
    mapping = current['coordinates']
    if len(current['variants']) != len(old['variants']):
        raise ValueError('Variant count changed')
    for before, after in zip(old['variants'], current['variants']):
        if any(after[key] != value for key, value in before.items()):
            raise ValueError('Variant data changed: ' + before['source_variant'])
    for before, after in zip(old['regions'], current['regions']):
        if any(after[key] != value for key, value in before.items()):
            raise ValueError('Region data changed')
        if after['genomic_end0']-after['genomic_start0'] != after['end0']-after['start0']:
            raise ValueError('Region interval length changed')
    def table_equal(before, after, profile=False):
        count = 0
        with open_text(before) as a, open_text(after) as b:
            for row, mapped in zip_longest(csv.DictReader(a, delimiter='\t'), csv.DictReader(b, delimiter='\t')):
                if row is None or mapped is None or any(mapped[k] != v for k,v in row.items()):
                    raise ValueError(f'Original table measurement changed: {before}, row {count}')
                if profile:
                    p = int(row['relative_pos0'])
                    if (int(mapped['genomic_pos1']) != mapping['genomic_start1']+p
                            or mapped['chrom'] != mapping['chrom'] or mapped['strand'] != '+'
                            or mapped['ref'] != current['sequence'][p]):
                        raise ValueError('Profile genomic position/reference mapping failed')
                count += 1
        return count
    counts = dict(profiles=table_equal(landscape/'profiles.tsv.gz', atlas/'genomic_profiles.tsv.gz',True),
                  region_effects=table_equal(previous/'region_effects.tsv.gz',atlas/'region_effects.tsv.gz'),
                  variants=table_equal(previous/'variant_summary.tsv',atlas/'variant_summary.tsv'),
                  wt_boundaries=table_equal(landscape/'wt_tf_boundaries.tsv',atlas/'wt_tf_boundaries_genomic.tsv'))
    if counts['profiles'] != len(current['variants'])*len(current['tracks'])*mapping['length']:
        raise ValueError('Profile grid incomplete')
    if (atlas/'functional_correlations.tsv').read_bytes() != (previous/'functional_correlations.tsv').read_bytes():
        raise ValueError('Functional correlations changed')
    with open(atlas/'base_coordinates.tsv') as h:
        bases=list(csv.DictReader(h,delimiter='\t'))
    if [int(r['genomic_pos1']) for r in bases] != list(range(mapping['genomic_start1'],mapping['genomic_end1']+1)):
        raise ValueError('Base mapping incomplete')
    bed=(atlas/'wt_regions.bed').read_text().splitlines()
    if len(bed)!=len(current['regions']):
        raise ValueError('BED region count differs')
    for line,r in zip(bed,current['regions']):
        chrom,start,end,name,score,strand=line.split('\t')
        if (chrom,int(start),int(end),name,strand)!=(mapping['chrom'],r['genomic_start1']-1,r['genomic_end1'],r['id'],'+'):
            raise ValueError('BED and inclusive region coordinates differ')
    check_sources(sources)
    with atomic_output(output) as out:
        write_json(out/'audit.json',dict(complete=True,original_measurements_unchanged=True,
                   counts=counts,genomic_bases=len(bases),sources=sources))
    print(counts)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--atlas',required=True)
    p.add_argument('--previous',default='results/rep3_landscape_20261018/atlas_navigation_20261019')
    p.add_argument('--landscape',default='results/rep3_landscape_20261018/landscape')
    p.add_argument('--output',required=True)
    a=p.parse_args();run(a.atlas,a.previous,a.landscape,a.output)

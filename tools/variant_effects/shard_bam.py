"""Write the reads of a site extraction (molecules.tsv.gz read_id column) into N shard BAMs.

Shard = crc32(read name) % N, so shards are stable across runs. Input order (coordinate-sorted)
is preserved, so each shard is indexable. Used to feed plasmid-LLR calling of every analysed
molecule.

Usage: python tools/variant_effects/shard_bam.py --molecules <extract>/molecules.tsv.gz --bam <in.bam> \
           --out-dir <workdir>/data --prefix shard --shards 4
"""
import argparse
import csv
import gzip
import zlib
from pathlib import Path

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--molecules', required=True)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--prefix', default='shard')
    ap.add_argument('--shards', type=int, default=4)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with gzip.open(a.molecules, 'rt') as h:
        wanted = {r['read_id'] for r in csv.DictReader(h, delimiter='\t')}
    src = pysam.AlignmentFile(a.bam, 'rb')
    names = [out / f'{a.prefix}{k}.bam' for k in range(a.shards)]
    tmp = [p.with_suffix('.tmp.bam') for p in names]
    files = [pysam.AlignmentFile(str(p), 'wb', template=src) for p in tmp]
    found = 0
    for rec in src.fetch(until_eof=True):
        if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or rec.query_name not in wanted:
            continue
        files[zlib.crc32(rec.query_name.encode()) % a.shards].write(rec)
        found += 1
    for f in files:
        f.close()
    for t, p in zip(tmp, names):
        t.rename(p)
        pysam.index(str(p))
    print(f'wrote {found:,} of {len(wanted):,} reads into {a.shards} shards', flush=True)
    if found != len(wanted):
        raise SystemExit('read count mismatch')


if __name__ == '__main__':
    main()

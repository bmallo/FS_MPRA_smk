"""One pass over the replicate-3 FiberHMM BAM: write per-set BAMs for the selected reads.

<set>.bam        MM/ML kept, footprint tags stripped (input for re-calling)
<set>.orig.bam   untouched, keeps the original 2.12.0 circular/linear calls

Usage: python tools/llr_compare/subset_bams.py --workdir results/llr_vs_fiberhmm_20261026 --bam <fiberhmm.bam>
"""
import argparse
from pathlib import Path

import pysam

CONTIG = 'LDLR_MPRA_NotI'
FOOTPRINT_TAGS = ('ns', 'nl', 'as', 'al', 'nq', 'aq', 'MA', 'AQ', 'AN')
SETS = ('wt_unmeth', 'wt_meth', 'wt_chrom')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workdir', required=True)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--contig', default=CONTIG)
    a = ap.parse_args()
    wd = Path(a.workdir)
    (wd / 'data').mkdir(exist_ok=True)
    owner = {}
    for s in SETS:
        for line in open(wd / 'config' / f'{s}.read_ids.txt'):
            if line.strip():
                owner[line.strip()] = s
    src = pysam.AlignmentFile(a.bam, 'rb')
    clean = {s: pysam.AlignmentFile(str(wd / 'data' / f'{s}.tmp.bam'), 'wb', template=src) for s in SETS}
    orig = {s: pysam.AlignmentFile(str(wd / 'data' / f'{s}.orig.tmp.bam'), 'wb', template=src) for s in SETS}
    found = dict.fromkeys(SETS, 0)
    for i, rec in enumerate(src.fetch(a.contig), 1):
        if rec.is_secondary or rec.is_supplementary:
            continue
        s = owner.get(rec.query_name)
        if s is None:
            continue
        orig[s].write(rec)
        for t in FOOTPRINT_TAGS:
            if rec.has_tag(t):
                rec.set_tag(t, None)
        clean[s].write(rec)
        found[s] += 1
        if i % 1_000_000 == 0:
            print(f'{i:,} records scanned; found {found}', flush=True)
    for f in list(clean.values()) + list(orig.values()):
        f.close()
    for s in SETS:
        for tmp, final in ((f'{s}.tmp.bam', f'{s}.bam'), (f'{s}.orig.tmp.bam', f'{s}.orig.bam')):
            (wd / 'data' / tmp).rename(wd / 'data' / final)
            pysam.index(str(wd / 'data' / final))
    expected = {s: sum(1 for v in owner.values() if v == s) for s in SETS}
    print('found', found, 'expected', expected, flush=True)
    if found != expected:
        raise SystemExit('read count mismatch')


if __name__ == '__main__':
    main()

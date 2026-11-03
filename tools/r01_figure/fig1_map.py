"""Figure 1: locate each DPAGT1 construct's genomic insert in hg38 (exact 25-mer anchoring, both strands).

For each construct FASTA, every 25-mer is looked up in the hg38 window; hits are collinear (genome offset - construct
offset constant, per strand) inside the insert. The largest collinear block gives the construct span, genome span
and orientation. Writes inserts.tsv: construct, contig length, construct start/end (1-based), hg38 start/end,
strand (+ = construct reads in genome direction), matched 25-mers, identity check.
Usage: python tools/r01_figure/fig1_map.py --fasta <hg38.fa> --region chr11:119080000-119120000 --constructs <fa ...>
           --output <dir>
"""
import argparse
from collections import Counter
from pathlib import Path

import pandas as pd
import pysam

K = 25
COMP = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(COMP)[::-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--region', required=True)
    ap.add_argument('--constructs', nargs='+', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    chrom, rng = a.region.split(':')
    g0, g1 = map(int, rng.replace(',', '').split('-'))
    gseq = pysam.FastaFile(a.fasta).fetch(chrom, g0 - 1, g1).upper()
    index = {}
    for i in range(len(gseq) - K + 1):
        index.setdefault(gseq[i:i + K], []).append(i)
    rows = []
    for path in a.constructs:
        fa = pysam.FastaFile(path)
        for name in fa.references:
            s = fa.fetch(name).upper()
            for strand, q in (('+', s), ('-', rc(s))):
                diag = Counter()
                hits = []
                for j in range(len(q) - K + 1):
                    for i in index.get(q[j:j + K], []):
                        diag[i - j] += 1
                        hits.append((i - j, j))
                if not diag:
                    continue
                d, n = diag.most_common(1)[0]
                js = sorted(j for dd, j in hits if dd == d)
                qs, qe = js[0], js[-1] + K - 1          # 0-based inclusive on q
                gs, ge = qs + d, qe + d
                ident = sum(x == y for x, y in zip(q[qs:qe + 1], gseq[gs:ge + 1])) / (qe - qs + 1)
                cs, ce = (qs + 1, qe + 1) if strand == '+' else (len(s) - qe, len(s) - qs)
                rows.append(dict(file=Path(path).name, construct=name, length=len(s), strand=strand, kmers=n,
                                 construct_start=cs, construct_end=ce, insert_bp=qe - qs + 1,
                                 hg38_chrom=chrom, hg38_start=g0 + gs, hg38_end=g0 + ge, identity=round(ident, 4)))
    t = pd.DataFrame(rows).sort_values(['construct', 'kmers'], ascending=[True, False]).drop_duplicates('construct')
    t.to_csv(out / 'inserts.tsv', sep='\t', index=False)
    print(t.to_string(index=False))


if __name__ == '__main__':
    main()

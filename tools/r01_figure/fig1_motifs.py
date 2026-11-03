"""Figure 1: motif matches (JASPAR 2026 + HOCOMOCO v14, FIMO-style exact P) in an hg38 window, both strands.

Background = the window's own strand-symmetric base composition. Writes motifs.tsv (P < --p) with hg38 coordinates.
Usage: python tools/r01_figure/fig1_motifs.py --fasta <hg38.fa> --region chr11:119100225-119104342
           --db results/r01_figure_20261102/motif_db --p 1e-5 --output <dir>
"""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

HERE = Path(__file__).parent


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--region', required=True)
    ap.add_argument('--db', required=True)
    ap.add_argument('--p', type=float, default=1e-5)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    s = importlib.util.spec_from_file_location('motif_scan', HERE / 'motif_scan.py')
    ms = importlib.util.module_from_spec(s)
    sys.modules['motif_scan'] = ms
    s.loader.exec_module(ms)
    out, db = Path(a.output), Path(a.db)
    out.mkdir(parents=True, exist_ok=True)
    chrom, rng = a.region.split(':')
    g0, g1 = map(int, rng.replace(',', '').split('-'))
    seq = pysam.FastaFile(a.fasta).fetch(chrom, g0 - 1, g1).upper()
    bg = np.array([seq.count(b) for b in ms.BASES], float)
    bg = (bg + bg[::-1]) / (2 * bg.sum())
    fam = {json.loads(x)['name']: json.loads(x)['masterlist_info']['tf'] for x in open(db / ms.HOCO_ANN)}
    rows = []
    for m in ms.read_motifs(db / ms.JASPAR, 'JASPAR2026') + ms.read_motifs(db / ms.HOCO, 'HOCOMOCOv14'):
        mn, sf, I = ms.pvalue_table(ms.log_odds(m['counts'], bg), bg)
        L = I.shape[0]
        pf, pr = ms.scan(seq, I, mn, sf)
        for strand, pv in (('+', pf), ('-', pr)):
            for i in np.flatnonzero(pv < a.p):
                rows.append(dict(db=m['db'], motif=m['id'], tf=fam.get(m['id'], m['name']).upper(), chrom=chrom,
                                 start=g0 + i, end=g0 + i + L - 1, strand=strand, p=pv[i]))
    t = pd.DataFrame(rows).sort_values('p')
    t.to_csv(out / 'motifs.tsv', sep='\t', index=False, float_format='%.3g')
    print(t.drop_duplicates('tf').head(25).to_string(index=False))


if __name__ == '__main__':
    main()

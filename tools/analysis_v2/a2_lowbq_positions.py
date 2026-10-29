"""Step A2: where in the element do bases fall below Q20? (the ledger stores only the barcode's qualities)

Streams primary plasmid alignments over the element, keeps a deterministic 1% (crc32(read_id) % 100 == 0),
and for reads whose element sequence exactly matches a design records each element base's quality.
Same element extraction as the ledger (fsmpra.processing.extract_interval), same Q20 threshold.

Outputs: per_position.tsv (element index, WT base, homopolymer run length, reads, reads with that base < Q20),
per_read.tsv.gz (np, rq, design class, number of bases < Q20, min quality) and summary.json.

Usage: python tools/analysis_v2/a2_lowbq_positions.py --bam <aligned.bam> --catalog <ledger>/catalog/catalog.json \
           --output results/analysis_v2_20261028/a2_accounting/HepG2_run3_lowbq
"""
import argparse
import json
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from fsmpra.io import read_fasta
from fsmpra.processing import extract_interval

Q = 20


def runs(seq):
    out = np.zeros(len(seq), int)
    i = 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        out[i:j] = j - i
        i = j
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--every', type=int, default=100)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    lib = json.load(open(a.catalog))['library']
    ref = read_fasta(Path(a.catalog).parent / 'plasmid.fa', lib['contig'])
    s0, e0 = lib['mutable_ref_start0'], lib['mutable_ref_end0']
    wt = lib['mutable_sequence']
    design = {v['sequence']: v['class'] for v in lib['variants']}
    L = len(wt)
    low = np.zeros(L, int)
    seen = np.zeros(L, int)
    rows = []
    with pysam.AlignmentFile(a.bam, 'rb') as bam:
        for r in bam.fetch(lib['contig'], s0, e0):
            if r.is_secondary or r.is_supplementary or zlib.crc32(r.query_name.encode()) % a.every:
                continue
            seq, qual, gs = extract_interval(r, s0, e0, ref)
            cls = design.get(seq) if gs == 'ok' else None
            nlow = int(sum(q < Q for q in qual)) if qual else None
            rows.append(dict(np=r.get_tag('np') if r.has_tag('np') else None, rq=r.get_tag('rq') if r.has_tag('rq') else None,
                             status=gs, design_class=cls or ('off_design' if gs == 'ok' else ''), n_low=nlow,
                             min_q=min(qual) if qual else None))
            if cls in ('WT', 'SNV') and len(qual) == L:  # same length as WT: index-aligned
                qa = np.asarray(qual)
                seen += 1
                low += qa < Q
    pr = pd.DataFrame(rows)
    pr.to_csv(out / 'per_read.tsv.gz', sep='\t', index=False)
    pp = pd.DataFrame(dict(element_index=np.arange(L), oligo_pos1=np.arange(L) + lib['oligo_mutable_start1'], wt_base=list(wt),
                           homopolymer_run=runs(wt), reads=seen, reads_below_q20=low))
    pp['fraction_below_q20'] = pp.reads_below_q20 / pp.reads
    pp.to_csv(out / 'per_position.tsv', sep='\t', index=False)
    ex = pr[pr.design_class.isin(['WT', 'SNV'])]
    by_run = pp.groupby('homopolymer_run').agg(positions=('element_index', 'size'), low=('reads_below_q20', 'sum'), reads=('reads', 'sum'))
    by_run['fraction_below_q20'] = by_run.low / by_run.reads
    summ = dict(sampled_reads=len(pr), exact_wt_snv=len(ex), exact_with_any_low=int((ex.n_low > 0).sum()),
                exact_low_with_one_base=int((ex.n_low == 1).sum()), exact_low_with_2plus=int((ex.n_low >= 2).sum()),
                low_fraction_by_homopolymer_run=by_run.reset_index().to_dict('records'),
                top_positions=pp.sort_values('fraction_below_q20', ascending=False).head(10).to_dict('records'))
    json.dump(summ, open(out / 'summary.json', 'w'), indent=1, default=float)
    print(json.dumps(summ, indent=1, default=float))


if __name__ == '__main__':
    main()

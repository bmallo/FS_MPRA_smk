"""Figure 2 QC: are FiberHMM MA spans in the SLC39A4 fp BAMs in stored-SEQ or original-read orientation?

For each primary read: m6A calls (ML >= 0.9) in stored-SEQ coordinates, and nucleosome calls (nuc track, >= 85 bp)
taken (i) as stored (SEQ orientation) and (ii) flipped for reverse reads (qlen - start - len). Correct orientation
gives few m6A inside nucleosomes on both strands; the wrong one gives ~the genome-wide m6A rate inside them on
reverse reads. Same check as fig1_profiles.py qc.tsv. Writes orientation_qc.tsv.
Usage: python tools/r01_figure/fig2_orientation_qc.py --bams <fp.bam ...> --output <dir> [--max-reads N]
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


def nuc_spans(tag):
    out = []
    for part in tag.split(';')[1:]:
        name, _, ents = part.partition(':')
        if re.sub(r'[+\-.].*$', '', name) != 'nuc':
            continue
        out += [tuple(int(x) for x in e.split('-')) for e in ents.split(',') if e]
    return [(s - 1, ln) for s, ln in out if ln >= 85]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bams', nargs='+', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--max-reads', type=int, default=3000)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for bam in a.bams:
        n = 0
        with pysam.AlignmentFile(bam) as f:
            for r in f.fetch(until_eof=True):
                if r.is_unmapped or r.is_secondary or r.is_supplementary or not r.has_tag('MA'):
                    continue
                n += 1
                if n > a.max_reads:
                    break
                ql = r.query_length
                m6a = np.array(sorted({q for (b, st, code), lst in (r.modified_bases or {}).items() if code == 'a'
                                       for q, p in lst if p >= 0.9 * 255}), int)
                if not m6a.size:
                    continue
                sp = nuc_spans(r.get_tag('MA'))
                for mode in ('as_stored', 'flipped_if_reverse'):
                    mask = np.zeros(ql, bool)
                    for s, ln in sp:
                        s_ = ql - s - ln if (mode == 'flipped_if_reverse' and r.is_reverse) else s
                        mask[max(s_, 0):min(s_ + ln, ql)] = True
                    rows.append(dict(bam=Path(bam).name, reverse=r.is_reverse, mode=mode, m6a=m6a.size,
                                     m6a_in_nuc=int(mask[m6a].sum()), nuc_bp=int(mask.sum()), qlen=ql))
    t = pd.DataFrame(rows).groupby(['bam', 'mode', 'reverse'])[['m6a', 'm6a_in_nuc', 'nuc_bp', 'qlen']].sum()
    t['share_m6a_in_nuc'] = t.m6a_in_nuc / t.m6a
    t['share_read_in_nuc'] = t.nuc_bp / t.qlen
    t.to_csv(out / 'orientation_qc.tsv', sep='\t', float_format='%.4f')
    print(t[['m6a', 'share_m6a_in_nuc', 'share_read_in_nuc']].to_string())


if __name__ == '__main__':
    main()

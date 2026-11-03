"""Figure 1: per-molecule footprints and per-position profiles on hg38 for the DPAGT1 datasets (fig1_rules.json).

For every primary read in each fiberhmm BAM: m6A calls (MM/ML, prob >= 0.9), fiberhmm MA annotations (nuc / tf
tracks, 1-based query start-length in original read orientation, flipped for reverse reads; footprint = any nuc or tf call), lifted to the reference through
the alignment and then to hg38 through inserts.tsv. Filters (rules): plasmid reads must span the construct end to end
(+/-50 bp), m6A/AT >= 0.05, largest footprint < 30% of the read; endogenous reads: m6A/AT >= 0.05.
Size classes: TF 10-40 bp, PIC 41-80 bp, nucleosome >= 85 bp.
Writes per dataset: profile.tsv (hg38 position: molecules, A/T observed, m6A fraction, coverage share per class),
molecules.parquet (kept molecules: footprints and m6A positions in hg38), filters.tsv (counts removed per filter),
plus qc.tsv (m6A inside nucleosome calls by read strand: orientation check).
Usage: python tools/r01_figure/fig1_profiles.py --calls <calls dir> --inserts <inserts.tsv> --output <dir>
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

DATASETS = {'DPAGT11': ('DPAGT11', 'plasmid'), 'DPAGT13': ('DPAGT13', 'plasmid'),
            'DPAGT1_full': ('DPAGT1', 'plasmid'), 'endogenous': (None, 'endogenous')}
WIN = ('chr11', 119100225, 119104342)
CLASSES = {'tf': (10, 40), 'pic': (41, 80), 'nuc': (85, 10_000)}
M6A_P, MIN_MF, MAX_FP_FRAC, END_SLACK = 0.9, 0.05, 0.3, 50


def parse_ma(ma):
    out = []
    for part in ma.split(';')[1:]:
        if ':' not in part:
            continue
        name, spans = part.split(':', 1)
        kind = name.split('.')[0]
        if kind not in ('nuc', 'tf') or not spans:
            continue
        for sp in spans.split(','):
            s, ln = sp.split('-')
            out.append((int(s) - 1, int(ln)))          # 0-based query start, length
    return out


def lift(qr, s, e):
    """Reference span of query [s, e) using the query->ref map (nearest aligned bases inside the span)."""
    seg = qr[s:e]
    seg = seg[seg >= 0]
    return (int(seg.min()), int(seg.max())) if seg.size else None


def hg38_map(ins):
    """{contig: (a, b)} so that hg38 = a + b * (0-based contig position); endogenous contigs map as hg38 = pos + 1."""
    return {c: ((r.hg38_start - r.construct_start + 1, 1) if r.strand == '+' else (r.hg38_end + r.construct_start - 1, -1))
            for c, r in ins.iterrows()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--calls', required=True)
    ap.add_argument('--inserts', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    ins = pd.read_csv(a.inserts, sep='\t').set_index('construct')
    hmap = hg38_map(ins)
    qc = []
    for ds, (contig, kind) in DATASETS.items():
        bam = pysam.AlignmentFile(Path(a.calls) / 'fiberhmm' / f'{ds}.fp.bam')
        clen = dict(zip(bam.references, bam.lengths))
        n0 = {k: 0 for k in ('reads', 'not_primary', 'not_spanning', 'under_methylated', 'unresolved', 'kept')}
        mols = []
        for r in bam.fetch(until_eof=True):
            n0['reads'] += 1
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                n0['not_primary'] += 1
                continue
            if kind == 'plasmid':
                if r.reference_name != contig or r.reference_start > END_SLACK or r.reference_end < clen[contig] - END_SLACK:
                    n0['not_spanning'] += 1
                    continue
            seq = r.query_sequence
            at = seq.count('A') + seq.count('T')
            mods = r.modified_bases or {}
            m6a = sorted({q for (b, st, code), lst in mods.items() if code == 'a' for q, p in lst if p >= M6A_P * 255})
            if at == 0 or len(m6a) / at < MIN_MF:
                n0['under_methylated'] += 1
                continue
            fps = parse_ma(r.get_tag('MA')) if r.has_tag('MA') else []
            if r.is_reverse:            # MA spans are in original read orientation; flip to SEQ orientation (QC-checked)
                fps = [(r.query_length - s0 - ln, ln) for s0, ln in fps]
            if fps and max(l for _, l in fps) > MAX_FP_FRAC * r.query_length:
                n0['unresolved'] += 1
                continue
            pairs = np.array(r.get_aligned_pairs(matches_only=True), int)
            qr = np.full(r.query_length, -1, int)
            qr[pairs[:, 0]] = pairs[:, 1]
            ha, hb = hmap.get(r.reference_name, (1, 1))
            to_h = lambda x: ha + hb * np.asarray(x, int)                       # noqa: E731
            fp_ref = []
            for s, ln in fps:
                sp = lift(qr, s, s + ln)
                if sp:
                    h = sorted(to_h(sp).tolist())
                    fp_ref.append((h[0], h[1], ln))
            mq = np.array(m6a, int)
            m_ref = to_h(qr[mq[qr[mq] >= 0]]).tolist() if mq.size else []
            sq = np.frombuffer(seq.encode(), 'S1')
            isat = (sq[pairs[:, 0]] == b'A') | (sq[pairs[:, 0]] == b'T')
            at_ref = to_h(pairs[isat, 1]).tolist()
            # orientation QC: share of m6A query positions that fall inside nucleosome calls (>= 85 bp)
            innuc = np.zeros(r.query_length, bool)
            for s0, ln in fps:
                if ln >= 85:
                    innuc[s0:s0 + ln] = True
            inside = int(innuc[mq].sum()) if mq.size else 0
            qc.append(dict(dataset=ds, reverse=r.is_reverse, m6a=len(m6a), m6a_in_nuc=inside))
            n0['kept'] += 1
            mols.append(dict(read=r.query_name, reverse=r.is_reverse, mf=len(m6a) / at,
                             fp=json.dumps(fp_ref), m6a=json.dumps(sorted(m_ref)), at=json.dumps(sorted(at_ref))))
        M = pd.DataFrame(mols)
        M.to_parquet(out / f'{ds}.molecules.parquet')
        pd.DataFrame([n0]).assign(dataset=ds).to_csv(out / f'{ds}.filters.tsv', sep='\t', index=False)
        # profiles on the hg38 window
        L = WIN[2] - WIN[1] + 1
        cover = np.zeros(L)
        atobs = np.zeros(L)
        meth = np.zeros(L)
        cls = {k: np.zeros(L) for k in CLASSES}
        for row in M.itertuples():
            atp = np.array(json.loads(row.at), int) - WIN[1]
            atp = atp[(atp >= 0) & (atp < L)]
            atobs[atp] += 1
            mp = np.array(json.loads(row.m6a), int) - WIN[1]
            mp = mp[(mp >= 0) & (mp < L)]
            meth[mp] += 1
            if atp.size:
                cover[atp.min():atp.max() + 1] += 1
            for s, e, ln in json.loads(row.fp):
                for k, (lo, hi) in CLASSES.items():
                    if lo <= ln <= hi:
                        a0, a1 = max(s - WIN[1], 0), min(e - WIN[1], L - 1)
                        if a1 >= a0:
                            cls[k][a0:a1 + 1] += 1
        if contig in ins.index:                         # plasmid: keep only positions inside this construct's insert
            r_ = ins.loc[contig]
            pos_ = np.arange(WIN[1], WIN[2] + 1)
            out_ = (pos_ < r_.hg38_start) | (pos_ > r_.hg38_end)
            for arr in [cover, atobs, meth] + list(cls.values()):
                arr[out_] = 0
        with np.errstate(invalid='ignore', divide='ignore'):
            prof = pd.DataFrame({'pos': np.arange(WIN[1], WIN[2] + 1), 'molecules': cover, 'at_observed': atobs,
                                 'm6a_frac': np.where(atobs > 0, meth / atobs, np.nan),
                                 **{f'{k}_share': np.where(cover > 0, v / cover, np.nan) for k, v in cls.items()}})
        prof.to_csv(out / f'{ds}.profile.tsv', sep='\t', index=False, float_format='%.5g')
        print(ds, n0, flush=True)
    q = pd.DataFrame(qc).groupby(['dataset', 'reverse'])[['m6a', 'm6a_in_nuc']].sum()
    q['share_in_nuc'] = q.m6a_in_nuc / q.m6a
    q.to_csv(out / 'qc.tsv', sep='\t')
    print(q)


if __name__ == '__main__':
    main()

"""Step A2: where every read goes, for one run's ledger (fsmpra ledger, workflow/ledger.smk).

1. Funnel: each primary plasmid molecule is assigned to the FIRST filter it fails, in this order:
     not full-length (partial or concatemer) -> low sequencing quality (np < 10, rq < 0.99 or MAPQ < 20)
     -> element unreadable (flanks incomplete or disagreeing) -> off-design element sequence
     -> exact design but a base below Q20 in the element -> designed 1-bp indel (not analysed)
     -> no usable m6A data -> analysis-eligible WT/single SNV.
   Reads before that: HiFi reads that are not primary plasmid alignments (genomic, unmapped).
2. Off-design classes: full-length, sequence-QC-passing molecules whose element matches no design are
   globally aligned to WT (mismatch 1, gap open 2, extend 0.5; Biopython). Each contiguous edit is tested on its own:
   "designed" if applying only that edit to WT gives a designed element.
3. Off-design and low-base-quality fractions by pass count (sequencing error falls with passes; real
   molecules do not).
4. Partial reads: molecules failing only the full-length test (exact WT/SNV element, Q20, sequencing QC, m6A),
   by read length; these would be usable if partial molecules were allowed.
5. Barcode outcome for eligible molecules (read-level genotype is kept either way; molecules without an
   accepted barcode are their own resampling cluster in the current effect code).

Usage: python tools/analysis_v2/a2_read_accounting.py --ledger-dir results/rep3_ledger_20261016 \
           --run HepG2_run3 --raw-reads 10513725 --output results/analysis_v2_20261028/a2_accounting/HepG2_run3
"""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align

COLS = ['observed_sequence', 'genotype_class', 'genotype_status', 'barcode_status', 'passes', 'read_quality', 'mapq',
        'plasmid_class', 'circular', 'm6a_status', 'sequence_qc', 'length_qc', 'methylation_qc', 'analysis_eligible',
        'barcode_assignment_status', 'read_length']
NP_BINS = [10, 15, 20, 30, 1000]
NP_LABELS = ['10-14', '15-19', '20-29', '30+']
STEPS = ['not full-length', 'low sequencing quality', 'element unreadable', 'off-design sequence',
         'exact design, base < Q20', 'designed 1-bp indel', 'no m6A data', 'eligible WT/SNV']


def step(d):
    s = np.full(len(d), 'eligible WT/SNV', dtype=object)
    rules = [(d.length_qc != 1, 'not full-length'), (d.sequence_qc != 1, 'low sequencing quality'),
             (d.genotype_status.isin(['incomplete_flanks', 'flank_disagreement']), 'element unreadable'),
             (d.genotype_status == 'off_design', 'off-design sequence'),
             (d.genotype_status == 'low_base_quality', 'exact design, base < Q20'),
             (d.genotype_class.isin(['INS', 'DEL']), 'designed 1-bp indel'), (d.methylation_qc != 1, 'no m6A data')]
    done = np.zeros(len(d), bool)
    for mask, lab in rules:
        m = mask.to_numpy() & ~done
        s[m] = lab
        done |= m
    return s


def aligner():
    # affine gaps so a long deletion is one event rather than scattered mismatches and gaps
    return Align.PairwiseAligner(mode='global', match_score=0, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-0.5)


def edits(al, wt, seq):
    """Contiguous edits of seq against wt: list of (wt_start, wt_end, replacement)."""
    if not seq:  # whole element deleted between intact flanks
        return [(0, len(wt), '')]
    aln = al.align(wt, seq)[0]
    ops, (tb, qb) = [], aln.aligned
    pi, pj = 0, 0
    for (i0, i1), (j0, j1) in zip(tb, qb):
        if i0 > pi or j0 > pj:
            ops.append((pi, i0, seq[pj:j0]))
        for k in range(i1 - i0):  # mismatches inside aligned blocks
            if wt[i0 + k] != seq[j0 + k]:
                ops.append((i0 + k, i0 + k + 1, seq[j0 + k]))
        pi, pj = i1, j1
    if pi < len(wt) or pj < len(seq):
        ops.append((pi, len(wt), seq[pj:]))
    merged = []  # merge adjacent ops into one contiguous edit
    for o in sorted(ops):
        if merged and o[0] <= merged[-1][1]:
            a0, a1, r = merged[-1]
            merged[-1] = (a0, max(a1, o[1]), r + o[2])
        else:
            merged.append(o)
    return merged


def homopolymer(wt, i0, i1, rep):
    """True for a 1-bp indel inside or extending a homopolymer run of length >= 2."""
    base = rep if rep else wt[i0]
    if len(base) != 1 or (i1 - i0) > 1 or (i1 - i0 == 1 and rep):
        return False
    lo, hi = i0, i0
    while lo > 0 and wt[lo - 1] == base:
        lo -= 1
    while hi < len(wt) and wt[hi] == base:
        hi += 1
    return hi - lo >= 2 or (i1 - i0 == 0 and hi - lo >= 1)


def classify(wt, designed, ed):
    kinds = []
    for i0, i1, rep in ed:
        single = wt[:i0] + rep + wt[i1:]
        size = max(i1 - i0, len(rep))
        if single in designed:
            kinds.append('designed')
        elif size >= 10:
            kinds.append('large ' + ('deletion' if i1 - i0 > len(rep) else 'insertion' if len(rep) > i1 - i0 else 'replacement'))
        elif i1 - i0 == 1 and len(rep) == 1:
            kinds.append('substitution')
        else:
            kinds.append('indel (homopolymer)' if homopolymer(wt, i0, i1, rep) else 'indel (other)')
    n_des = kinds.count('designed')
    large = [k for k in kinds if k.startswith('large')]
    if large:
        return 'large deletion/insertion (>= 10 bp)', kinds
    if len(kinds) == 2 and n_des == 2:
        return 'two designed edits', kinds
    if len(kinds) == 2 and n_des == 1:
        return 'designed edit + one other edit', kinds
    if len(kinds) == 1:
        return 'one non-designed ' + kinds[0], kinds
    if len(kinds) == 2:
        return 'two non-designed edits', kinds
    return 'three or more edits', kinds


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ledger-dir', required=True)
    ap.add_argument('--run', required=True)
    ap.add_argument('--raw-reads', type=int, required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--min-new-snv', type=int, default=50, help='molecules needed to count a non-designed SNV as measurable')
    a = ap.parse_args()
    L = Path(a.ledger_dir)
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    lib = json.load(open(L / 'catalog/catalog.json'))['library']
    wt = lib['mutable_sequence']
    designed = {v['sequence'] for v in lib['variants']}

    funnel, offseq = Counter(), Counter()
    by_np = defaultdict(Counter)
    bc = Counter()
    partial = Counter()
    n = 0
    for d in pd.read_csv(L / 'ledger/molecule_ledger.tsv.gz', sep='\t', usecols=COLS, chunksize=500_000,
                         dtype={'observed_sequence': str, 'barcode_assignment_status': str}, keep_default_na=False):
        n += len(d)
        s = step(d)
        funnel.update(s)
        npb = pd.cut(pd.to_numeric(d.passes, errors='coerce'), NP_BINS, right=False, labels=NP_LABELS).astype(str)
        # error-rate panels: full-length, readable element, MAPQ/rq passing, np >= 10
        base = (d.length_qc == 1) & (d.sequence_qc == 1) & ~d.genotype_status.isin(['incomplete_flanks', 'flank_disagreement'])
        for lab, g in d[base].groupby(npb[base]):
            by_np[lab]['molecules'] += len(g)
            by_np[lab]['off_design'] += int((g.genotype_status == 'off_design').sum())
            by_np[lab]['exact_low_bq'] += int((g.genotype_status == 'low_base_quality').sum())
            by_np[lab]['exact_ok'] += int((g.genotype_status == 'exact_design').sum())
        off = d[(s == 'off-design sequence')]
        offseq.update(off.observed_sequence)
        pm = ((d.length_qc != 1) & (d.sequence_qc == 1) & (d.genotype_status == 'exact_design')
              & d.genotype_class.isin(['WT', 'SNV']) & (d.methylation_qc == 1))
        lb = pd.cut(pd.to_numeric(d.read_length[pm], errors='coerce'), [0, 1000, 2000, 3000, 4000, 1e9], right=False,
                    labels=['<1 kb', '1-2 kb', '2-3 kb', '3-4 kb', '>=4 kb']).astype(str)
        partial.update(zip(d.genotype_class[pm], lb))
        el = d[s == 'eligible WT/SNV']
        bc.update(zip(el.genotype_class, el.barcode_status, el.barcode_assignment_status.replace('', 'none')))
        print(f'{n:,} rows', flush=True)

    plasmid = n
    fun = [dict(step='not a primary plasmid alignment (genomic or unmapped)', molecules=a.raw_reads - plasmid)]
    fun += [dict(step=k, molecules=funnel.get(k, 0)) for k in STEPS]
    fun = pd.DataFrame(fun)
    fun['fraction_of_hifi'] = fun.molecules / a.raw_reads
    fun['fraction_of_plasmid'] = np.where(fun.index == 0, np.nan, fun.molecules / plasmid)
    fun.insert(0, 'run', a.run)
    fun.to_csv(out / 'funnel.tsv', sep='\t', index=False)

    al = aligner()
    rows, new_snv = [], Counter()
    for seq, cnt in offseq.items():
        ed = edits(al, wt, seq)
        cls, kinds = classify(wt, designed, ed)
        rows.append(dict(cls=cls, molecules=cnt, n_edits=len(ed)))
        if cls == 'one non-designed substitution':
            i0, _, rep = ed[0]
            new_snv[f'{i0 + lib["oligo_mutable_start1"]}{wt[i0]}>{rep}'] += cnt
    oc = pd.DataFrame(rows)
    oc = oc.groupby('cls').agg(molecules=('molecules', 'sum'), unique_sequences=('molecules', 'size')).reset_index()
    oc['fraction_of_off_design'] = oc.molecules / oc.molecules.sum()
    oc.insert(0, 'run', a.run)
    oc.sort_values('molecules', ascending=False).to_csv(out / 'offdesign_classes.tsv', sep='\t', index=False)
    ns = pd.Series(new_snv, name='molecules', dtype=int).rename_axis('snv').reset_index().sort_values('molecules', ascending=False)
    ns.to_csv(out / 'nondesigned_snvs.tsv', sep='\t', index=False)

    npt = pd.DataFrame([dict(np_bin=k, **v) for k, v in by_np.items() if k in NP_LABELS])
    npt['off_design_fraction'] = npt.off_design / npt.molecules
    npt['low_bq_fraction_of_exact'] = npt.exact_low_bq / (npt.exact_low_bq + npt.exact_ok)
    npt.insert(0, 'run', a.run)
    npt.sort_values('np_bin', key=lambda x: x.map({k: i for i, k in enumerate(NP_LABELS)})).to_csv(out / 'error_by_passes.tsv', sep='\t', index=False)

    bct = pd.DataFrame([dict(genotype_class=g, barcode_read=b, barcode_assignment=s, molecules=c) for (g, b, s), c in bc.items()])
    bct.insert(0, 'run', a.run)
    bct.sort_values('molecules', ascending=False).to_csv(out / 'barcode_outcomes.tsv', sep='\t', index=False)

    pt = pd.DataFrame([dict(genotype_class=g, read_length=l, molecules=c) for (g, l), c in partial.items()])
    pt.insert(0, 'run', a.run)
    pt.to_csv(out / 'partial_recoverable.tsv', sep='\t', index=False)
    summ = dict(run=a.run, partial_recoverable=int(sum(partial.values())), hifi_reads=a.raw_reads, plasmid_primary=plasmid, eligible=int(funnel['eligible WT/SNV']),
                off_design_molecules=int(sum(offseq.values())), off_design_unique=len(offseq),
                nondesigned_snvs_with_min_molecules=int((ns.molecules >= a.min_new_snv).sum()) if len(ns) else 0,
                min_new_snv=a.min_new_snv,
                eligible_with_accepted_barcode=int(sum(c for (g, b, s), c in bc.items() if s in ('supported',))))
    json.dump(summ, open(out / 'summary.json', 'w'), indent=1)
    pd.set_option('display.width', 200)
    print(fun.to_string(index=False))
    print(oc.to_string(index=False))
    print(npt.to_string(index=False))
    print(bct.head(20).to_string(index=False))
    print(json.dumps(summ, indent=1))


if __name__ == '__main__':
    main()

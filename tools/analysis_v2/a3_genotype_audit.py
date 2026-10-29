"""Step A3: are genotype and barcode calls consistent and correct? (WT, SNVs and designed 1-bp indels)

Read-level genotype = exact match of the read's 211-bp element to a design (fsmpra ledger). A barcode
identifies one plasmid clone, so every read with the same barcode should carry the same genotype.

Per run (fsmpra ledger):
 1. Within-barcode concordance. "Qualified" reads: exact design, sequencing QC passed, readable N15 barcode.
    In barcodes with >= --deep qualified reads and a clear majority genotype (> 50%), a read that differs
    from the majority is discordant. Type of difference:
      one edit apart: majority or read is WT (the other a single designed edit), or both edits at the
          same position -> consistent with a sequencing, synthesis or plasmid-propagation error
      different positions: two unrelated designs -> barcode shared by two clones or a template switch
    by passes (sequencing error falls with passes) and by genotype class; also for full-length reads only
    (barcode majorities still use all qualified reads).
 2. Barcode assignment rules for analysis-ready molecules (exact design, QC, full-length, m6A; any class):
      exact (current): no competing qualified genotype, >= 2 reads
      majority 0.9 / 0.75: majority genotype holds >= that fraction of >= 2 qualified reads, read agrees
    plus molecules whose own genotype contradicts a confident majority (>= 0.9 of >= 3 reads).
Across runs (same physical library): barcodes with a confident majority (>= 0.9 of >= 2 reads) in two runs
should name the same genotype.
Across implementations (optional): per-read calls of the October genotyper (genotype_reads.py read
ledger) vs the fsmpra ledger, compared as element sequences.

Usage: python tools/analysis_v2/a3_genotype_audit.py --catalog <ledger>/catalog/catalog.json \
           --ledger HepG2_run1=<dir> --ledger HepG2_run2=<dir> --ledger HepG2_run3=<dir> \
           --second HepG2_run1=<read_ledger.tsv.gz> --second HepG2_run2=<...> --output <dir>
"""
import argparse
import itertools
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

COLS = ['read_id', 'observed_variant_id', 'genotype_status', 'barcode', 'barcode_status', 'passes', 'sequence_qc', 'length_qc',
        'methylation_qc']
NP_BINS, NP_LABELS = [0, 10, 15, 20, 30, 10000], ['<10', '10-14', '15-19', '20-29', '30+']


def catalog(path):
    lib = json.load(open(path))['library']
    v = pd.DataFrame(lib['variants'])
    v['eq'] = v.equivalent_positions.map(lambda x: frozenset(x or []))
    return lib, v.set_index('variant_id')


def relation(a, b, cls, eq):
    """Edit relation between two different designs a, b (variant ids)."""
    if cls[a] == 'WT' or cls[b] == 'WT':
        return 'one edit apart'
    return 'one edit apart' if eq[a] & eq[b] else 'different positions'


def load(ledger_dir):
    d = pd.read_csv(Path(ledger_dir) / 'ledger/molecule_ledger.tsv.gz', sep='\t', usecols=COLS, keep_default_na=False,
                    dtype={'barcode': str, 'observed_variant_id': str})
    d['qualified'] = (d.genotype_status == 'exact_design') & (d.sequence_qc == 1) & (d.barcode_status == 'ok')
    d['ready'] = (d.genotype_status == 'exact_design') & (d.sequence_qc == 1) & (d.length_qc == 1) & (d.methylation_qc == 1)
    d['np_bin'] = pd.cut(pd.to_numeric(d.passes, errors='coerce'), NP_BINS, right=False, labels=NP_LABELS).astype(str)
    return d


def barcode_table(q):
    """Per barcode: qualified reads, majority genotype and its share."""
    c = q.groupby(['barcode', 'observed_variant_id'], observed=True).size().rename('n').reset_index()
    c = c.sort_values(['barcode', 'n'], ascending=[True, False])
    tot = c.groupby('barcode').n.sum()
    top = c.drop_duplicates('barcode').set_index('barcode')
    return pd.DataFrame(dict(reads=tot, majority=top.observed_variant_id, majority_n=top.n,
                             genotypes=c.groupby('barcode').size())).assign(share=lambda t: t.majority_n / t.reads)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--ledger', action='append', required=True, help='RUN=ledger_dir')
    ap.add_argument('--second', action='append', default=[], help='RUN=genotype_reads read_ledger.tsv.gz')
    ap.add_argument('--deep', type=int, default=5)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    lib, cat = catalog(a.catalog)
    cls, eq = cat['class'].to_dict(), cat['eq'].to_dict()
    disc_rows, cls_rows, rule_rows, maj, second_rows = [], [], [], {}, []
    for spec in a.ledger:
        run, path = spec.split('=', 1)
        d = load(path)
        q = d[d.qualified]
        bt = barcode_table(q)
        maj[run] = bt
        x = q.join(bt, on='barcode')
        deep = x[(x.reads >= a.deep) & (x.share > 0.5)].copy()
        deep['discordant'] = deep.observed_variant_id != deep.majority
        dd = deep[deep.discordant]
        rel = pd.Series([relation(r, m, cls, eq) for r, m in zip(dd.observed_variant_id, dd.majority)], index=dd.index, dtype=object)
        deep['relation'] = 'agrees'
        deep.loc[dd.index, 'relation'] = rel
        for reads_set, sub in (('all qualified', deep), ('full-length only', deep[deep.length_qc == 1])):
            for npb, g in sub.groupby('np_bin'):
                vc = g.relation.value_counts()
                disc_rows.append(dict(run=run, reads_set=reads_set, np_bin=npb, reads=len(g), one_edit_apart=int(vc.get('one edit apart', 0)),
                                      different_positions=int(vc.get('different positions', 0))))
        deep['read_class'] = deep.observed_variant_id.map(cls)
        deep['majority_class'] = deep.majority.map(cls)
        for mc, g in deep.groupby('majority_class'):
            vc = g.relation.value_counts()
            row = dict(run=run, majority_class=mc, reads=len(g), one_edit_apart=int(vc.get('one edit apart', 0)),
                       different_positions=int(vc.get('different positions', 0)))
            for rc in ('WT', 'SNV', 'DEL', 'INS'):
                row['discordant_read_is_' + rc] = int(((g.relation != 'agrees') & (g.read_class == rc)).sum())
            cls_rows.append(row)
        # assignment rules on analysis-ready molecules
        r = d[d.ready].join(bt, on='barcode')
        r['cls'] = r.observed_variant_id.map(cls)
        bcok = r.barcode_status == 'ok'
        agree = r.observed_variant_id == r.majority
        exact = bcok & (r.genotypes == 1) & (r.reads >= 2)
        contra = bcok & (r.reads >= 3) & (r.share >= 0.9) & ~agree
        for c, g in r.groupby('cls'):
            m = r.cls == c
            rule_rows.append(dict(run=run, genotype_class=c, ready_molecules=int(m.sum()), barcode_unreadable=int((m & ~bcok).sum()),
                                  singleton_barcode=int((m & bcok & (r.reads == 1)).sum()),
                                  exact_rule=int((m & exact).sum()),
                                  majority_90=int((m & bcok & (r.reads >= 2) & (r.share >= 0.9) & agree).sum()),
                                  majority_75=int((m & bcok & (r.reads >= 2) & (r.share >= 0.75) & agree).sum()),
                                  contradicts_confident_majority=int((m & contra).sum())))
        n_bc = int((bt.reads >= 1).sum())
        print(run, len(d), 'qualified', len(q), 'barcodes', n_bc, flush=True)
        # second genotyper
        sp = dict(s.split('=', 1) for s in a.second)
        if run in sp:
            wt = lib['mutable_sequence']
            o1 = lib['oligo_mutable_start1']
            seq2id = dict(zip(cat.sequence, cat.index))

            def to_seq(name):
                if name == 'WT':
                    return wt
                m = re.fullmatch(r'(\d+)([ACGT])>([ACGT])', name)
                if m:
                    i = int(m.group(1)) - o1
                    return wt[:i] + m.group(3) + wt[i + 1:] if wt[i] == m.group(2) else None
                m = re.fullmatch(r'(\d+)del([ACGT])', name)
                if m:
                    i = int(m.group(1)) - o1
                    return wt[:i] + wt[i + 1:] if wt[i] == m.group(2) else None
                m = re.fullmatch(r'(\d+)_(\d+)ins([ACGT]+)', name)
                if m:
                    i = int(m.group(2)) - o1
                    return wt[:i] + m.group(3) + wt[i:]
                return None
            s2 = pd.read_csv(sp[run], sep='\t', usecols=['read', 'class', 'variant'], keep_default_na=False)
            names = s2.variant.unique()
            name2id = {n: seq2id.get(to_seq(n)) for n in names}
            s2['id2'] = s2.variant.map(name2id)
            j = d[(d.sequence_qc == 1)][['read_id', 'observed_variant_id', 'genotype_status']].merge(
                s2[['read', 'class', 'id2']], left_on='read_id', right_on='read', how='inner')
            j['fsmpra_design'] = j.genotype_status.isin(['exact_design', 'low_base_quality'])
            j['second_design'] = j.id2.notna()
            both = j[j.fsmpra_design & j.second_design]
            second_rows.append(dict(run=run, reads_compared=len(j), both_call_design=len(both),
                                    same_design=int((both.observed_variant_id == both.id2).sum()),
                                    different_design=int((both.observed_variant_id != both.id2).sum()),
                                    only_fsmpra=int((j.fsmpra_design & ~j.second_design).sum()),
                                    only_second=int((~j.fsmpra_design & j.second_design).sum()),
                                    unmapped_second_names=int(sum(v is None for v in name2id.values()))))
            ex = both[both.observed_variant_id != both.id2].head(20)
            ex.assign(fsmpra=ex.observed_variant_id.map(cat.source_variant), second=ex.id2.map(cat.source_variant))[
                ['read_id', 'fsmpra', 'second']].to_csv(out / f'second_genotyper_disagreements_{run}.tsv', sep='\t', index=False)
        del d, q, x, deep, r

    disc = pd.DataFrame(disc_rows)
    disc['one_edit_rate'] = disc.one_edit_apart / disc.reads
    disc['different_positions_rate'] = disc.different_positions / disc.reads
    disc.to_csv(out / 'discordance_by_passes.tsv', sep='\t', index=False)
    ct = pd.DataFrame(cls_rows)
    ct['discordant_rate'] = (ct.one_edit_apart + ct.different_positions) / ct.reads
    ct.to_csv(out / 'discordance_by_class.tsv', sep='\t', index=False)
    rules = pd.DataFrame(rule_rows)
    rules.to_csv(out / 'assignment_rules.tsv', sep='\t', index=False)
    cross = []
    for r1, r2 in itertools.combinations(maj, 2):
        b1 = maj[r1][(maj[r1].reads >= 2) & (maj[r1].share >= 0.9)]
        b2 = maj[r2][(maj[r2].reads >= 2) & (maj[r2].share >= 0.9)]
        s = b1.join(b2, how='inner', lsuffix='_a', rsuffix='_b')
        same = s.majority_a == s.majority_b
        rel = pd.Series([relation(x, y, cls, eq) for x, y in zip(s.majority_a[~same], s.majority_b[~same])], dtype=object)
        cross.append(dict(run_a=r1, run_b=r2, confident_barcodes_a=len(b1), confident_barcodes_b=len(b2), shared=len(s),
                          same_genotype=int(same.sum()), differ_one_edit=int((rel == 'one edit apart').sum()),
                          differ_positions=int((rel == 'different positions').sum())))
    cr = pd.DataFrame(cross)
    cr['agreement'] = cr.same_genotype / cr.shared
    cr.to_csv(out / 'cross_run_barcodes.tsv', sep='\t', index=False)
    space = 4 ** 15
    coll = pd.DataFrame([dict(run=k, barcodes=len(v), barcodes_ge2_reads=int((v.reads >= 2).sum()),
                              multi_genotype_ge2=int(((v.reads >= 2) & (v.genotypes > 1)).sum()),
                              two_strong_genotypes=int(((v.reads - v.majority_n) >= 2).sum()),
                              expected_random_shared_pairs=len(v) ** 2 / (2 * space)) for k, v in maj.items()])
    coll.to_csv(out / 'barcode_collisions.tsv', sep='\t', index=False)
    if second_rows:
        pd.DataFrame(second_rows).to_csv(out / 'second_genotyper.tsv', sep='\t', index=False)
    pd.set_option('display.width', 250)
    for t in (disc, ct, rules, cr, coll, pd.DataFrame(second_rows)):
        print(t.round(5).to_string(index=False))


if __name__ == '__main__':
    main()

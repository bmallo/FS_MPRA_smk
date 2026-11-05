"""Candidate pathogenic and benign truth sets per element (definitions: truth_set_rules.json).

Inputs under --output (from fetch.py, ag_score.py, assemble.py): clinvar/<window>.json, variants.tsv, all_snvs.tsv,
alphagenome/tracks.tsv. Pathogenic = ClinVar P/LP small variants in the window. Benign candidates = B1 (common in any
ancestry group) and B2 (not a GTEx V10 eQTL for any gene) and B3 (no large AlphaGenome RNA and DNase change, top-20%
rule within the element, and PromoterAI |score| < 0.1).
Writes to <output>/truth_sets/: clinvar_variants.tsv, b3_snvs.tsv, candidates.tsv, truth_sets.tsv, summary.json.
Usage: python tools/promoter_variant_census/truth_sets.py --output <census dir>
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import census_lib as cl  # noqa: E402

ELEMENTS = {'LDLR': ['LDLR'], 'HBB': ['HBB'], 'GATA2_IS': ['GATA2_IS'], 'GATA2_IG': ['GATA2_IG'],
            'SLC39A4': ['SLC39A4_W1', 'SLC39A4_W2'], 'GATA2_E95': ['GATA2_E95_W1', 'GATA2_E95_W2']}
NO_PAI = {'GATA2_E95'}                           # enhancer: PromoterAI not applicable (not counted as failing)
AG_TRACKS = {   # (RNA_SEQ curies, DNASE curies); LDLR = HepG2 as Figure 3 panel b
    'LDLR': (['EFO:0001187'], ['EFO:0001187']),
    'HBB': (['EFO:0002067', 'UBERON:0013756'], ['EFO:0002067']),
    'GATA2': (['EFO:0002067', 'CL:0001059', 'CL:0000837'], ['EFO:0002067', 'CL:0001059', 'CL:0000837']),
    'SLC39A4': (['UBERON:0002108', 'EFO:0001099'], ['UBERON:0002108', 'EFO:0001099'])}
TOP, PAI = 0.8, 0.1
MAX_SMALL = 50                                   # ClinVar records longer than this are structural, not counted


def cv_class(desc):
    d = desc.lower()
    if 'conflicting' in d:
        return 'conflicting'
    if 'pathogenic' in d and 'benign' not in d and 'uncertain' not in d:
        return 'plp'
    if 'benign' in d and 'pathogenic' not in d and 'uncertain' not in d:
        return 'blb'
    if 'uncertain' in d:
        return 'vus'
    return 'other'


def clinvar_table(out):
    w = cl.windows().set_index('window')
    rows = []
    for win in w.index:
        x = w.loc[win]
        res = json.load(open(out / 'clinvar' / f'{win}.json'))['result']
        for uid in res['uids']:
            r = res[uid]
            vs = r['variation_set'][0]
            spdi = vs.get('canonical_spdi', '')
            gc = r.get('germline_classification', {})
            row = dict(window=win, gene=x.gene, uid=uid, accession=r['accession'], title=r['title'],
                       clinvar_desc=gc.get('description', ''), review=gc.get('review_status', ''),
                       cls=cv_class(gc.get('description', '')), obj_type=r.get('obj_type', ''), spdi=spdi,
                       traits=';'.join(sorted({t['trait_name'] for t in gc.get('trait_set', [])})))
            if not spdi:
                rows.append(dict(row, kind='structural/no SPDI', in_window=False))
                continue
            _, p0, d, i = spdi.split(':')
            p = int(p0) + 1
            while d and i and d[0] == i[0]:                       # trim shared bases
                d, i, p = d[1:], i[1:], p + 1
            while d and i and d[-1] == i[-1]:
                d, i = d[:-1], i[:-1]
            if max(len(d), len(i)) > MAX_SMALL:
                rows.append(dict(row, kind='structural/no SPDI', in_window=False))
                continue
            inw = (x.start <= p + len(d) - 1 and p <= x.end) if d else (x.start <= p <= x.end)
            kind = 'SNV' if len(d) == 1 and len(i) == 1 else 'indel'
            key = f'{x.chrom}:{p}:{d}:{i}' if kind == 'SNV' else ''
            rows.append(dict(row, kind=kind, in_window=inw, pos=p, ref=d, alt=i, key=key))
    return pd.DataFrame(rows).astype({'pos': 'Int64'})


def ag_b3(out):
    """Per element and possible SNV: mean raw AlphaGenome RNA (target gene) and DNase over the element's tracks,
    large-change flags (top 20% of |mean| among the element's possible SNVs)."""
    t = pd.read_csv(out / 'alphagenome/tracks.tsv', sep='\t',
                    usecols=['key', 'context', 'output_type', 'ontology_curie', 'raw_score'])
    s = pd.read_csv(out / 'all_snvs.tsv', sep='\t')
    s['key'] = s.chrom + ':' + s.pos.astype(str) + ':' + s.ref + ':' + s.alt
    s['element'] = s.window.map({wn: e for e, ws in ELEMENTS.items() for wn in ws})
    s = s.drop_duplicates(['element', 'key'])
    rows = []
    for el, g in s.groupby('element'):
        ctx = cl.CONTEXT[ELEMENTS[el][0]]
        gene = g.gene.iloc[0]
        rna, dn = AG_TRACKS[gene]
        tt = t[t.context == ctx]
        r = tt[(tt.output_type == 'RNA_SEQ') & tt.ontology_curie.isin(rna)].groupby('key').raw_score.mean()
        dd = tt[(tt.output_type == 'DNASE') & tt.ontology_curie.isin(dn)].groupby('key').raw_score.mean()
        e = g[['element', 'gene', 'key', 'pai_score']].copy()
        e['ag_rna'] = e.key.map(r)
        e['ag_dnase'] = e.key.map(dd)
        for c in ('ag_rna', 'ag_dnase'):
            cut = e[c].abs().quantile(TOP)
            e[f'{c}_cut'] = cut
            e[f'{c}_large'] = e[c].abs() >= cut
        e['pai_applicable'] = el not in NO_PAI
        e['pai_large'] = e.pai_score.abs() >= PAI if el not in NO_PAI else False
        e['b3'] = ~e.ag_rna_large & ~e.ag_dnase_large & ~e.pai_large
        rows.append(e)
    return pd.concat(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    d = out / 'truth_sets'
    d.mkdir(exist_ok=True)
    cv = clinvar_table(out)
    b3 = ag_b3(out)
    checks = dict(b3_missing_ag=int(b3[['ag_rna', 'ag_dnase']].isna().any(axis=1).sum()),
                  b3_missing_pai=int((b3.pai_score.isna() & b3.pai_applicable).sum()))
    v = pd.read_csv(out / 'variants.tsv', sep='\t')
    v = v[v.in_gnomad].copy()
    v['key'] = v.chrom + ':' + v.pos.astype(str) + ':' + v.ref + ':' + v.alt
    v['element'] = v.window.map({wn: e for e, ws in ELEMENTS.items() for wn in ws})
    v = v.drop_duplicates(['element', 'key'])
    sig_any = v.gtex_v10_sig_any_genes.fillna('') != ''
    v['b2_class'] = np.where(sig_any, 'significant', np.where(v.gtex_v10_maf01, 'tested_not_significant',
                                                              'not_tested'))
    v['b1'] = v.af_gt_1pct_any_group
    v['b1_overall'] = v.af_gt_1pct
    v['b2'] = v.b2_class != 'significant'
    v['b2_strict'] = v.b2_class == 'tested_not_significant'
    v = v.merge(b3[['element', 'key', 'ag_rna', 'ag_dnase', 'ag_rna_large', 'ag_dnase_large', 'pai_large', 'b3']],
                on=['element', 'key'], how='left')
    v['b3'] = v.b3.fillna(False).astype(bool)                   # indels: no PromoterAI score -> cannot pass B3
    cvs = cv[cv.in_window & (cv.kind != 'structural/no SPDI')].copy()
    cvs['element'] = cvs.window.map({wn: e for e, ws in ELEMENTS.items() for wn in ws})
    cvs = cvs.drop_duplicates(['element', 'accession'])
    v = v.merge(cvs[cvs.kind == 'SNV'][['element', 'key', 'cls', 'accession']].rename(
        columns={'cls': 'clinvar_cls', 'accession': 'clinvar'}), on=['element', 'key'], how='left')
    v['clinvar_cls'] = v.clinvar_cls.fillna('not in ClinVar')
    # ClinVar variants with census evidence (gnomAD, B1-B3) for the SNVs
    cvs = cvs.merge(v[['element', 'key', 'g_af', 'g_grpmax_af', 'b1', 'b2_class', 'b3']], on=['element', 'key'],
                    how='left').merge(b3[['element', 'key', 'ag_rna_large', 'ag_dnase_large', 'pai_large',
                                          'pai_score']].rename(columns={'pai_large': 'pai_large_any'}),
                                      on=['element', 'key'], how='left')
    rows = []
    for el in ELEMENTS:
        V, C, B = v[v.element == el], cvs[cvs.element == el], b3[b3.element == el]
        s, i = V[V.kind == 'SNV'], V[V.kind == 'indel']
        r = dict(element=el, possible_snvs=len(B),
                 clinvar_small=len(C), clinvar_structural=int(((cv.window.isin(ELEMENTS[el])) &
                                                             (cv.kind == 'structural/no SPDI')).sum()),
                 plp_snv=int(((C.cls == 'plp') & (C.kind == 'SNV')).sum()),
                 plp_indel=int(((C.cls == 'plp') & (C.kind == 'indel')).sum()),
                 blb_snv=int(((C.cls == 'blb') & (C.kind == 'SNV')).sum()),
                 blb_indel=int(((C.cls == 'blb') & (C.kind == 'indel')).sum()),
                 conflicting=int((C.cls == 'conflicting').sum()), vus=int((C.cls == 'vus').sum()),
                 other=int((C.cls == 'other').sum()),
                 plp_in_gnomad=int((C.cls == 'plp').sum() and C[(C.cls == 'plp')].g_af.notna().sum()),
                 gnomad_snvs=len(s), gnomad_indels=len(i),
                 b1_snvs=int(s.b1.sum()), b1_indels=int(i.b1.sum()), b1_overall_snvs=int(s.b1_overall.sum()),
                 b1_b2_snvs=int((s.b1 & s.b2).sum()), b1_b2_tested_snvs=int((s.b1 & s.b2_strict).sum()),
                 b1_b2_b3_snvs=int((s.b1 & s.b2 & s.b3).sum()),
                 b1_b2_b3_strict_snvs=int((s.b1 & s.b2_strict & s.b3).sum()),
                 b1_b2_b3_overall_snvs=int((s.b1_overall & s.b2 & s.b3).sum()),
                 b3_possible_snvs=int(B.b3.sum()), b3_gnomad_snvs=int(s.b3.sum()),
                 b1_eqtl_any_gene=int((s.b1 & (s.b2_class == 'significant')).sum()),
                 b1_untested=int((s.b1 & (s.b2_class == 'not_tested')).sum()),
                 candidates_also_clinvar_blb=int((s.b1 & s.b2 & s.b3 & (s.clinvar_cls == 'blb')).sum()),
                 b1_also_clinvar_plp=int((s.b1 & (s.clinvar_cls == 'plp')).sum()))
        rows.append(r)
    ts = pd.DataFrame(rows)
    ts.to_csv(d / 'truth_sets.tsv', sep='\t', index=False)
    cvs.drop(columns=['uid']).sort_values(['element', 'pos']).to_csv(d / 'clinvar_variants.tsv', sep='\t',
                                                                      index=False, float_format='%.4g')
    b3.to_csv(d / 'b3_snvs.tsv', sep='\t', index=False, float_format='%.4g')
    cols = ['element', 'window', 'gnomad_id', 'rsid', 'kind', 'g_af', 'g_grpmax_af', 'g_grpmax_group', 'b1',
            'b1_overall', 'b2_class', 'gtex_v10_sig_any_genes', 'ag_rna', 'ag_rna_large', 'ag_dnase', 'ag_dnase_large',
            'pai_score', 'pai_large', 'b3', 'clinvar_cls', 'clinvar']
    v[v.b1][cols].sort_values(['element', 'gnomad_id']).to_csv(d / 'candidates.tsv', sep='\t', index=False,
                                                               float_format='%.4g')
    cut = b3.groupby('element')[['ag_rna_cut', 'ag_dnase_cut']].first().to_dict('index')
    json.dump(dict(checks=checks, ag_cutoffs=cut, clinvar_query_date='2026-11-04',
                   clinvar_terms={wn: json.load(open(out / 'clinvar' / f'{wn}.json'))['term']
                                  for ws in ELEMENTS.values() for wn in ws}),
              open(d / 'summary.json', 'w'), indent=1)
    print(json.dumps(checks))
    print(ts.set_index('element').T.to_string())


if __name__ == '__main__':
    main()

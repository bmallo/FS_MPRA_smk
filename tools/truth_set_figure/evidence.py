"""Per-SNV evidence table and truth-set counts for the truth-set figure (definitions: rules.json).

Elements: LDLR, HBB, SLC39A4 (promoter variant census), CARD11, JAK1, JAK3, STAT1, STAT3 (immune promoter scan) and
TERT (TERT promoter scan). Universe = every possible SNV. Arm A = confident (>= 2-star) ClinVar P/LP or B/LB. Arm B =
everything else, split into likely benign candidates (common, not an eQTL for any gene, no PromoterAI or AlphaGenome
change), likely pathogenic candidates (rare, PromoterAI |score| >= 0.2 and AlphaGenome large change) and no call.
Writes to --output: snv_evidence.tsv, clinvar_confident.tsv, counts.tsv, funnel.tsv, summary.json
--no-alphagenome: rules.json draft1_no_alphagenome (candidates without the AlphaGenome condition).
--min-stars 1 --benign in_gnomad: rules.json draft3 (>= 1-star ClinVar; benign = in gnomAD at any frequency, no GTEx).
Usage: python tools/truth_set_figure/evidence.py --output results/truth_set_figure_<date> [--elements ...]
       [--no-alphagenome] [--min-stars {1,2}] [--benign {common_not_eqtl,in_gnomad}]
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CENSUS = ROOT / 'results/promoter_variant_census_20261104'
IMMUNE = ROOT / 'results/immune_promoter_scan_20261105'
TERT = ROOT / 'results/tert_promoter_scan_20261105'
IMMUNE_SET = ['CL:0000084', 'CL:0000624', 'CL:0000625', 'CL:0000236', 'CL:0001054', 'EFO:0002784', 'CLO:0007045']
TERT_SET = ['EFO:0003042', 'EFO:0003045', 'EFO:0002067', 'EFO:0001187', 'EFO:0002791']
# element -> (result dir, windows, AlphaGenome context, RNA_SEQ curies, DNASE curies)
ELEMENTS = {
    'LDLR': (CENSUS, ['LDLR'], 'LDLR', ['EFO:0001187'], ['EFO:0001187']),
    'HBB': (CENSUS, ['HBB'], 'HBB', ['EFO:0002067', 'UBERON:0013756'], ['EFO:0002067']),
    'SLC39A4': (CENSUS, ['SLC39A4_W1', 'SLC39A4_W2'], 'SLC39A4', ['UBERON:0002108', 'EFO:0001099'],
                ['UBERON:0002108', 'EFO:0001099']),
    **{g: (IMMUNE, [g], g, IMMUNE_SET, IMMUNE_SET) for g in ('CARD11', 'JAK1', 'JAK3', 'STAT1', 'STAT3')},
    'TERT': (TERT, ['TERT'], 'TERT', TERT_SET, TERT_SET)}
TOP, PAI_BENIGN, PAI_PATH, RARE = 0.8, 0.1, 0.2, 1e-4


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


def stars(review):
    r = str(review)
    if r == 'practice guideline':
        return 4
    if r == 'reviewed by expert panel':
        return 3
    if r == 'criteria provided, multiple submitters, no conflicts':
        return 2
    if r.startswith('criteria provided'):
        return 1
    return 0


def load_clinvar(d, el, wins, min_stars):
    if d == CENSUS:
        c = pd.read_csv(d / 'truth_sets/clinvar_variants.tsv', sep='\t')
        c = c[c.element == el]
    else:
        c = pd.read_csv(d / 'clinvar_variants.tsv', sep='\t')
        c = c[c.window.isin(wins)]
        chrom = pd.read_csv(d / 'variants.tsv', sep='\t', usecols=['window', 'chrom']).groupby('window').chrom.first()
        chrom = c.window.map(chrom)
        c = c.assign(key=np.where(c.kind == 'SNV', chrom + ':' + c.pos.astype('Int64').astype(str) + ':' + c.ref.fillna('')
                                  + ':' + c.alt.fillna(''), ''))
    c = c[c.in_window.astype(bool)].drop_duplicates('accession').copy()
    c['element'] = el
    c['stars'] = c.review.map(stars)
    c['confident'] = c.cls.isin(['plp', 'blb']) & (c.stars >= min_stars)
    return c


def ag_large(d, ctx, rna, dn, keys):
    t = pd.read_csv(d / 'alphagenome/tracks.tsv', sep='\t',
                    usecols=['key', 'context', 'output_type', 'ontology_curie', 'raw_score'])
    t = t[t.context == ctx]
    r = t[(t.output_type == 'RNA_SEQ') & t.ontology_curie.isin(rna)].groupby('key').raw_score.mean()
    n = t[(t.output_type == 'DNASE') & t.ontology_curie.isin(dn)].groupby('key').raw_score.mean()
    e = pd.DataFrame({'key': keys})
    e['ag_rna'], e['ag_dnase'] = e.key.map(r), e.key.map(n)
    for c in ('ag_rna', 'ag_dnase'):
        cut = e[c].abs().quantile(TOP)
        e[f'{c}_cut'] = cut
        e[f'{c}_large'] = e[c].abs() >= cut
    e['ag_large'] = e.ag_rna_large | e.ag_dnase_large
    e['ag_ntracks_rna'] = t[(t.output_type == 'RNA_SEQ') & t.ontology_curie.isin(rna)].ontology_curie.nunique()
    e['ag_ntracks_dnase'] = t[(t.output_type == 'DNASE') & t.ontology_curie.isin(dn)].ontology_curie.nunique()
    return e


def element_table(el, use_ag=True, min_stars=2, benign='common_not_eqtl'):
    d, wins, ctx, rna, dn = ELEMENTS[el]
    s = pd.read_csv(d / 'all_snvs.tsv', sep='\t')
    s = s[s.window.isin(wins)].copy()
    s['key'] = s.chrom + ':' + s.pos.astype(str) + ':' + s.ref + ':' + s.alt
    s = s.drop_duplicates('key')[['key', 'chrom', 'pos', 'ref', 'alt', 'pai_score']]
    g = pd.read_csv(d / 'variants.tsv', sep='\t', low_memory=False)
    g = g[g.window.isin(wins) & (g.kind == 'SNV')].copy()
    g['key'] = g.chrom + ':' + g.pos.astype(str) + ':' + g.ref + ':' + g.alt
    g = g.drop_duplicates('key').set_index('key')
    ga = g.g_ac.where(g.g_filters.fillna('') == '', 0)            # FILTER failures count as absent
    ea = g.e_ac.where(g.e_filters.fillna('') == '', 0)
    s['joint_ac'] = s.key.map(ga + ea).fillna(0).astype(int)
    s['joint_af'] = s.key.map((ga + ea) / (g.g_an + g.e_an).replace(0, np.nan)).fillna(0.0)
    s['in_gnomad_genomes'] = s.key.map(g.in_gnomad).fillna(False).astype(bool)
    s['g_af'] = s.key.map(g.g_af)
    s['g_grpmax_af'] = s.key.map(g.g_grpmax_af)
    s['g_grpmax_group'] = s.key.map(g.g_grpmax_group)
    s['b1_common'] = s.key.map(g.af_gt_1pct_any_group).fillna(False).astype(bool)
    s['gtex_v10_class'] = s.key.map(g.gtex_v10_class)
    s['gtex_any_gene'] = s.key.map(g.gtex_v10_sig_any_genes).fillna('')
    s['b2_not_eqtl'] = s.gtex_any_gene == ''
    s['rare'] = (s.joint_ac == 0) | (s.joint_af < RARE)
    s['in_gnomad_any'] = s.joint_ac > 0                          # genomes or exomes, FILTER PASS
    if use_ag:
        s = s.merge(ag_large(d, ctx, rna, dn, s.key), on='key', how='left')
    else:
        s = s.assign(ag_rna=np.nan, ag_dnase=np.nan, ag_large=False, ag_ntracks_rna=0, ag_ntracks_dnase=0)
    s['pai_scored'] = s.pai_score.notna()
    s['pai_noeffect'] = s.pai_scored & (s.pai_score.abs() < PAI_BENIGN)
    s['pai_large'] = s.pai_scored & (s.pai_score.abs() >= PAI_PATH)
    c = load_clinvar(d, el, wins, min_stars)
    cs = c[c.kind == 'SNV'].sort_values('stars', ascending=False)
    first = cs.drop_duplicates('key').set_index('key')
    conf = cs[cs.confident].drop_duplicates('key').set_index('key')
    s['clinvar_cls'] = s.key.map(first.cls).fillna('')
    s['clinvar_stars'] = s.key.map(first.stars)
    s['clinvar'] = s.key.map(first.accession).fillna('')
    s['arm_a'] = s.key.map(conf.cls).fillna('')
    s['b3_no_effect'] = s.pai_noeffect & ~s.ag_large
    pop = s.b1_common & s.b2_not_eqtl if benign == 'common_not_eqtl' else s.in_gnomad_any
    s['lb_candidate'] = (s.arm_a == '') & pop & s.b3_no_effect
    s['lp_candidate'] = (s.arm_a == '') & s.rare & s.pai_large & (s.ag_large | (not use_ag))
    s['call'] = np.select([s.arm_a == 'plp', s.arm_a == 'blb', s.lp_candidate, s.lb_candidate],
                          ['ClinVar P/LP', 'ClinVar B/LB', 'likely pathogenic candidate', 'likely benign candidate'],
                          'no call')
    s.insert(0, 'element', el)
    return s, c


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output', required=True)
    ap.add_argument('--elements', nargs='*', default=list(ELEMENTS))
    ap.add_argument('--no-alphagenome', action='store_true')
    ap.add_argument('--min-stars', type=int, default=2, choices=[1, 2], help='ClinVar review stars for arm A')
    ap.add_argument('--benign', default='common_not_eqtl', choices=['common_not_eqtl', 'in_gnomad'],
                    help='population part of the likely benign rule')
    a = ap.parse_args()
    use_ag = not a.no_alphagenome
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    tabs, cvs = zip(*[element_table(el, use_ag, a.min_stars, a.benign) for el in a.elements])
    S, C = pd.concat(tabs), pd.concat(cvs)
    checks = dict(alphagenome_used=use_ag,
                  snvs_missing_ag=int(S[['ag_rna', 'ag_dnase']].isna().any(axis=1).sum()) if use_ag else None,
                  snvs_missing_pai=int((~S.pai_scored).sum()),
                  b1_but_rare=int((S.b1_common & S.rare).sum()))
    rows, fun = [], []
    for el, s in S.groupby('element', sort=False):
        c = C[C.element == el]
        B = s[s.arm_a == '']
        r = dict(element=el, possible_snvs=len(s), ag_tracks_rna=int(s.ag_ntracks_rna.iloc[0]),
                 ag_tracks_dnase=int(s.ag_ntracks_dnase.iloc[0]), pai_unscored=int((~s.pai_scored).sum()))
        for cls, lab in (('plp', 'plp'), ('blb', 'blb')):
            for k in ('SNV', 'indel'):
                x = c[(c.cls == cls) & (c.kind == k)]
                r[f'clinvar_{lab}_{k.lower()}_2star'] = int((x.stars >= 2).sum())
                r[f'clinvar_{lab}_{k.lower()}_confident'] = int(x.confident.sum())
                r[f'clinvar_{lab}_{k.lower()}_1star'] = int((x.stars == 1).sum())
                r[f'clinvar_{lab}_{k.lower()}_0star'] = int((x.stars == 0).sum())
        r['clinvar_vus_snv'] = int(((c.cls == 'vus') & (c.kind == 'SNV')).sum())
        r['clinvar_conflicting_snv'] = int(((c.cls == 'conflicting') & (c.kind == 'SNV')).sum())
        r['arm_a_plp_snvs'] = int((s.arm_a == 'plp').sum())
        r['arm_a_blb_snvs'] = int((s.arm_a == 'blb').sum())
        r['arm_b_snvs'] = len(B)
        r['arm_b_with_other_clinvar'] = int((B.clinvar_cls != '').sum())
        r['lb_candidates'] = int(B.lb_candidate.sum())
        r['lb_candidates_not_tested_gtex'] = int((B.lb_candidate & (B.gtex_v10_class == 'not_tested')).sum())
        r['lp_candidates'] = int(B.lp_candidate.sum())
        r['lp_candidates_absent_gnomad'] = int((B.lp_candidate & (B.joint_ac == 0)).sum())
        r['lp_candidates_pai_negative'] = int((B.lp_candidate & (B.pai_score < 0)).sum())
        r['no_call'] = int((B.call == 'no call').sum())
        rows.append(r)
        if a.benign == 'common_not_eqtl':
            steps = [('benign', 'not confident ClinVar', B), ('benign', 'common (B1)', B[B.b1_common])]
            pop = [('not an eQTL, any gene (B2)', 'b2_not_eqtl')]
        else:
            steps = [('benign', 'not confident ClinVar', B), ('benign', 'in gnomAD (any AF)', B[B.in_gnomad_any])]
            pop = []
        x = steps[-1][2]
        for lab, col in pop + [('PromoterAI |score| < 0.1', 'pai_noeffect')]:
            x = x[x[col]]
            steps.append(('benign', lab, x))
        if use_ag:
            steps.append(('benign', 'no AlphaGenome large change', x[~x.ag_large]))
        x = B[B.rare]
        steps += [('pathogenic', 'not confident ClinVar', B), ('pathogenic', 'rare (absent or AF < 0.01%)', x)]
        x = x[x.pai_large]
        steps.append(('pathogenic', 'PromoterAI |score| >= 0.2', x))
        if use_ag:
            steps.append(('pathogenic', 'AlphaGenome large change', x[x.ag_large]))
        fun += [dict(element=el, side=sd, step=i, filter=lab, n=len(v)) for i, (sd, lab, v) in enumerate(steps)]
    counts, funnel = pd.DataFrame(rows), pd.DataFrame(fun)
    S.to_csv(out / 'snv_evidence.tsv', sep='\t', index=False, float_format='%.6g')
    C[C.confident].to_csv(out / 'clinvar_confident.tsv', sep='\t', index=False)
    counts.to_csv(out / 'counts.tsv', sep='\t', index=False)
    funnel.to_csv(out / 'funnel.tsv', sep='\t', index=False)
    json.dump(dict(elements=a.elements, alphagenome=use_ag, min_stars=a.min_stars, benign=a.benign, checks=checks,
                   md5={p: md5(out / p) for p in ('snv_evidence.tsv', 'clinvar_confident.tsv', 'counts.tsv',
                                                  'funnel.tsv')}),
              open(out / 'summary.json', 'w'), indent=1)
    print(json.dumps(checks))
    print(counts.T.to_string())


if __name__ == '__main__':
    main()

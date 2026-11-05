"""Assemble the promoter variant census from the raw pulls (fetch.py, ag_score.py); definitions in rules.json.

Writes to --output:
  variants.tsv  one row per gnomAD variant in a window (genome or exome AC >= 1), all flags and scores
  all_snvs.tsv  every possible SNV per window with gnomAD flags, AlphaGenome and PromoterAI calls
  counts.tsv    main table per window, per element (SLC39A4 W1+W2 counted once) and per gene (GATA2 IS+IG)
  funnel.tsv    count left after each filter step (SNVs; indels through GTEx), both MAF versions
  summary.json  resource versions, md5s, sanity checks
Usage: python tools/promoter_variant_census/assemble.py --output <dir>
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import census_lib as cl  # noqa: E402

AG_SETS = {
    'LDLR': ['EFO:0001187', 'UBERON:0002107', 'UBERON:0001114', 'UBERON:0001115'],
    'HBB': ['EFO:0002067', 'UBERON:0013756', 'UBERON:0000178', 'UBERON:0002371', 'UBERON:0012168', 'CL:0002092'],
    'GATA2': ['EFO:0002067', 'CL:0001059', 'CL:0000837', 'UBERON:0002371', 'CL:0002092', 'UBERON:0012168'],
    'SLC39A4': ['UBERON:0002108', 'UBERON:0002114', 'EFO:0001099', 'CL:0002563']}
EXPT = ['EFO:0001187', 'EFO:0001182']                    # HepG2, HEK293 (current experiment's cell lines)
Q1, Q2, P1, P2 = 0.99, 0.95, 0.1, 0.2
UNITS = {'LDLR': ['LDLR'], 'HBB': ['HBB'], 'GATA2_IS': ['GATA2_IS'], 'GATA2_IG': ['GATA2_IG'],
         'SLC39A4_W1': ['SLC39A4_W1'], 'SLC39A4_W2': ['SLC39A4_W2'], 'SLC39A4 (W1+W2)': ['SLC39A4_W1', 'SLC39A4_W2'],
         'GATA2 (IS+IG)': ['GATA2_IS', 'GATA2_IG'], 'GATA2_E95_W1': ['GATA2_E95_W1'], 'GATA2_E95_W2': ['GATA2_E95_W2'],
         'GATA2 +9.5 (W1+W2)': ['GATA2_E95_W1', 'GATA2_E95_W2']}
NO_PAI = {'GATA2_E95_W1', 'GATA2_E95_W2'}                  # enhancer: PromoterAI not applicable (truth_set_rules.json)


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


def ag_calls(out):
    """Per (context, key): max |quantile| over the gene's biosample set (RNA_SEQ target gene + CAGE), over all
    tracks (RNA_SEQ + CAGE, from maxabs), over HepG2/HEK293, and the DNase max over the gene set."""
    t = pd.read_csv(out / 'alphagenome/tracks.tsv', sep='\t')
    m = pd.read_csv(out / 'alphagenome/maxabs.tsv', sep='\t')
    t['aq'] = t.quantile_score.abs()
    t['gene'] = t.context.str.split('_').str[0]
    t['in_set'] = [c in AG_SETS[g] for c, g in zip(t.ontology_curie, t.gene)]
    expr = t[t.output_type.isin(['RNA_SEQ', 'CAGE'])]
    k = ['context', 'key']
    a = expr[expr.in_set].groupby(k).agg(ag_set_maxq=('aq', 'max'), ag_set_ntracks=('aq', 'size'))
    j = expr[expr.in_set].loc[expr[expr.in_set].groupby(k).aq.idxmax()].set_index(k)
    a['ag_set_track_at_max'] = j.output_type + ' ' + j.biosample_name
    a['ag_set_q_at_max'] = j.quantile_score
    a['ag_set_rna_maxq'] = expr[expr.in_set & (expr.output_type == 'RNA_SEQ')].groupby(k).aq.max()
    a['ag_set_cage_maxq'] = expr[expr.in_set & (expr.output_type == 'CAGE')].groupby(k).aq.max()
    a['ag_expt_maxq'] = expr[expr.ontology_curie.isin(EXPT)].groupby(k).aq.max()
    a['ag_set_dnase_maxq'] = t[t.in_set & (t.output_type == 'DNASE')].groupby(k).aq.max()
    mm = m[m.output_type.isin(['RNA_SEQ', 'CAGE'])].groupby(k).max_abs_q.max()
    a['ag_all_maxq'] = mm
    a = a.reset_index()
    for q, s in ((Q1, '099'), (Q2, '095')):
        a[f'ag_noeffect_{s}'] = a.ag_set_maxq < q
    a['ag_all_noeffect_099'] = a.ag_all_maxq < Q1
    a['ag_expt_noeffect_099'] = a.ag_expt_maxq < Q1
    return a


def pai_calls(out, w):
    r = pd.read_csv(out / 'promoterai/rows.tsv', sep='\t')
    r['tx'] = r.transcript_id.str.split('.').str[0]
    r['key'] = r.chrom + ':' + r.pos.astype(str) + ':' + r.ref + ':' + r.alt
    rows = []
    for x in w.itertuples():
        s = r[(r.chrom == x.chrom) & r.pos.between(x.start, x.end) & (r.gene == x.gene)]
        own = s[s.tx == x.transcript].set_index('key').promoterAI
        anyt = s.assign(a=s.promoterAI.abs()).groupby('key').a.max()
        rows.append(pd.DataFrame({'window': x.window, 'pai_score': own, 'pai_maxabs_any_tx': anyt}).rename_axis('key'))
    p = pd.concat(rows).reset_index()
    p['pai_noeffect_01'] = p.pai_score.abs() < P1
    p['pai_noeffect_02'] = p.pai_score.abs() < P2
    return p, r


def gtex_calls(out, g):
    d = out / 'gtex'
    sig = pd.read_csv(d / 'significant_target_gene.tsv', sep='\t')
    fm = pd.read_csv(d / 'finemapping_v8_target_gene.tsv', sep='\t')
    lk = pd.read_csv(d / 'variant_lookup.tsv', sep='\t').set_index('variant_id')
    s10 = sig[sig.datasetId == 'gtex_v10'].groupby(['variantId', 'target_gene']).tissueSiteDetailId.apply(
        lambda s: ';'.join(sorted(set(s))))
    s8 = sig[sig.datasetId == 'gtex_v8'].groupby(['variantId', 'target_gene']).size()
    f8 = fm.groupby(['variantId', 'target_gene']).method.apply(lambda s: ';'.join(sorted(set(s))))
    ix = list(zip(g.gtex_id, g.gene))
    g['gtex_v10_sig_tissues'] = [s10.get(i, '') for i in ix]
    g['gtex_v10_sig_target'] = g.gtex_v10_sig_tissues != ''
    g['gtex_v8_sig_target'] = [i in s8.index for i in ix]
    g['gtex_v8_finemap_target'] = [f8.get(i, '') for i in ix]
    for c in ('gtex_v10_present', 'gtex_v10_maf01', 'gtex_v8_present', 'gtex_v8_maf01'):
        g[c] = g.gtex_id.map(lk[c]).fillna(False).astype(bool)
    g['gtex_v10_sig_any_genes'] = g.gtex_id.map(lk['v10_sig_genes']).fillna('') if 'v10_sig_genes' in lk else ''
    g['gtex_v10_class'] = np.where(g.gtex_v10_sig_target, 'significant',
                                   np.where(g.gtex_v10_maf01, 'tested_not_significant', 'not_tested'))
    g['gtex_looked_up'] = g.gtex_id.isin(lk.index)
    return g


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    w = cl.windows()
    checks = {}

    # gnomAD
    g = cl.gnomad_table(out)
    g = g[(g.g_ac > 0) | (g.e_ac > 0)].copy()
    g['g_pass'] = g.g_filters == ''
    g['in_gnomad'] = (g.g_ac > 0) & g.g_pass
    g['g_fail_filter'] = (g.g_ac > 0) & ~g.g_pass
    g['exome_only'] = (g.g_ac == 0) & (g.e_ac > 0) & (g.e_filters == '')
    g['af_gt_1pct'] = g.in_gnomad & (g.g_af > 0.01)
    g['af_gt_1pct_any_group'] = g.in_gnomad & ((g.g_af > 0.01) | (g.g_grpmax_af > 0.01))
    g['ref_is_minor'] = g.in_gnomad & (g.g_af > 0.99)
    g['singleton'] = g.in_gnomad & (g.g_ac == 1)
    g = gtex_calls(out, g)
    g['not_eqtl'] = ~g.gtex_v10_sig_target
    g['not_eqtl_strict'] = g.gtex_v10_class == 'tested_not_significant'
    checks['gtex_lookup_missing_for_genome_ac1'] = int((~g.gtex_looked_up & (g.g_ac > 0)).sum())

    # reference check for gnomAD SNVs
    snv = cl.all_snvs(out)
    refmap = dict(zip(snv.chrom + ':' + snv.pos.astype(str), snv.ref))
    gs = g[g.kind == 'SNV']
    checks['gnomad_snv_ref_mismatch'] = int(sum(refmap[f'{c}:{p}'] != r for c, p, r in zip(gs.chrom, gs.pos, gs.ref)))

    # predictors
    ag = ag_calls(out)
    snv['context'] = snv.window.map(cl.CONTEXT)
    g['context'] = g.window.map(cl.CONTEXT)
    pai, praw = pai_calls(out, w)
    checks['promoterai_ref_mismatch'] = int(sum(refmap.get(f'{c}:{p}', r) != r for c, p, r in
                                                zip(praw.chrom, praw.pos, praw.ref)))
    snv = snv.merge(ag, on=['context', 'key'], how='left').merge(pai, on=['window', 'key'], how='left')
    g = g.merge(ag, on=['context', 'key'], how='left').merge(pai, on=['window', 'key'], how='left')
    checks['snvs_missing_alphagenome'] = int(snv.ag_set_maxq.isna().sum())
    checks['snvs_missing_promoterai'] = int((snv.pai_score.isna() & ~snv.window.isin(NO_PAI)).sum())
    checks['gnomad_indels_missing_alphagenome'] = int((g.in_gnomad & (g.kind == 'indel') & g.ag_set_maxq.isna()).sum())
    gk = g[g.in_gnomad & (g.kind == 'SNV')].set_index(['window', 'key'])
    for c in ('in_gnomad', 'af_gt_1pct', 'af_gt_1pct_any_group', 'not_eqtl', 'not_eqtl_strict', 'gtex_v10_class',
              'g_af', 'g_grpmax_af', 'rsid'):
        snv[c] = [gk[c].get((wn, k), False if gk[c].dtype == bool else '') for wn, k in zip(snv.window, snv.key)]
    for c in ('ag_noeffect_099', 'ag_noeffect_095', 'ag_all_noeffect_099', 'ag_expt_noeffect_099',
              'pai_noeffect_01', 'pai_noeffect_02'):
        na = c.startswith('pai_')                                 # not applicable counts as no effect
        snv[c] = np.where(snv.window.isin(NO_PAI) & na, True, snv[c].fillna(False)).astype(bool)
        g[c] = np.where(g.window.isin(NO_PAI) & na & (g.kind == 'SNV'), True, g[c].fillna(False)).astype(bool)
    snv['both_noeffect'] = snv.ag_noeffect_099 & snv.pai_noeffect_01
    g['both_noeffect'] = g.ag_noeffect_099 & g.pai_noeffect_01

    # counts and funnel per unit (dedupe by key across W1/W2 of SLC39A4)
    rows, fun = [], []
    for unit, wins in UNITS.items():
        S = snv[snv.window.isin(wins)].drop_duplicates('key')
        G = g[g.window.isin(wins)].drop_duplicates('key')
        Gs, Gi = G[G.in_gnomad & (G.kind == 'SNV')], G[G.in_gnomad & (G.kind == 'indel')]
        bp = len(set(p for x in w[w.window.isin(wins)].itertuples() for p in range(x.start, x.end + 1)))
        c = Gs.af_gt_1pct
        r = dict(unit=unit, bp=bp, all_possible_snvs=len(S),
                 gnomad_snvs=len(Gs), gnomad_indels=len(Gi),
                 excluded_filter_snvs=int((G.g_fail_filter & (G.kind == 'SNV')).sum()),
                 excluded_filter_indels=int((G.g_fail_filter & (G.kind == 'indel')).sum()),
                 exome_only_snvs=int((G.exome_only & (G.kind == 'SNV')).sum()),
                 gnomad_snv_singletons=int(Gs.singleton.sum()),
                 af1_overall_snvs=int(c.sum()), af1_overall_indels=int(Gi.af_gt_1pct.sum()),
                 af1_anygroup_snvs=int(Gs.af_gt_1pct_any_group.sum()),
                 af1_anygroup_indels=int(Gi.af_gt_1pct_any_group.sum()),
                 ref_is_minor_snvs=int(Gs.ref_is_minor.sum()))
        for k_, V in (('snvs', Gs), ('indels', Gi)):
            for cls in ('significant', 'tested_not_significant', 'not_tested'):
                r[f'gtex_{cls}_{k_}'] = int((V.gtex_v10_class == cls).sum())
        r['gtex_v8_sig_target_snvs'] = int(Gs.gtex_v8_sig_target.sum())
        r['gtex_v8_finemap_target_snvs'] = int((Gs.gtex_v8_finemap_target != '').sum())
        r['gtex_v10_sig_any_gene_snvs'] = int((Gs.gtex_v10_sig_any_genes != '').sum())
        for nm, col in (('ag_noeffect', 'ag_noeffect_099'), ('ag_noeffect_q095', 'ag_noeffect_095'),
                        ('ag_noeffect_alltracks', 'ag_all_noeffect_099'), ('pai_noeffect', 'pai_noeffect_01'),
                        ('pai_noeffect_02', 'pai_noeffect_02'), ('both_noeffect', 'both_noeffect')):
            r[f'{nm}_allsnvs'] = int(S[col].sum())
            r[f'{nm}_gnomad_snvs'] = int(Gs[col].sum())
            r[f'{nm}_af1_snvs'] = int(Gs[col][c].sum())
        r['ag_noeffect_gnomad_indels'] = int(Gi.ag_noeffect_099.sum())
        for mv, mcol in (('overall', 'af_gt_1pct'), ('any_group', 'af_gt_1pct_any_group')):
            for ev, ecol in (('', 'not_eqtl'), ('_strict', 'not_eqtl_strict')):
                X = Gs[Gs[mcol] & Gs[ecol] & Gs.ag_noeffect_099 & Gs.pai_noeffect_01]
                r[f'intersection_{mv}{ev}'] = len(X)
        rows.append(r)
        for mv, mcol in (('overall', 'af_gt_1pct'), ('any_group', 'af_gt_1pct_any_group')):
            steps = [('all possible SNVs', S.key), ('in gnomAD genomes (PASS)', Gs.key)]
            X = Gs[Gs[mcol]]
            steps.append((f'AF > 1% ({mv})', X.key))
            X = X[X.not_eqtl]
            steps.append(('not GTEx V10 eQTL for target gene', X.key))
            X = X[X.ag_noeffect_099]
            steps.append(('AlphaGenome no effect (|q| < 0.99, gene set)', X.key))
            X = X[X.pai_noeffect_01]
            steps.append(('PromoterAI no effect (|score| < 0.1)', X.key))
            fun += [dict(unit=unit, maf_version=mv, kind='SNV', step=i + 1, filter=s, n=len(k_))
                    for i, (s, k_) in enumerate(steps)]
            X = Gi[Gi[mcol]]
            fun += [dict(unit=unit, maf_version=mv, kind='indel', step=i + 1, filter=s, n=n_) for i, (s, n_) in
                    enumerate([('in gnomAD genomes (PASS)', len(Gi)), (f'AF > 1% ({mv})', len(X)),
                               ('not GTEx V10 eQTL for target gene', int(X.not_eqtl.sum()))])]
    counts, funnel = pd.DataFrame(rows), pd.DataFrame(fun)
    counts.to_csv(out / 'counts.tsv', sep='\t', index=False)
    funnel.to_csv(out / 'funnel.tsv', sep='\t', index=False)
    vcols = ['window', 'gene', 'gnomad_id', 'rsid', 'kind', 'chrom', 'pos', 'ref', 'alt', 'g_ac', 'g_an', 'g_af',
             'g_hom', 'g_filters', 'g_grpmax_af', 'g_grpmax_group', 'g_groups_af_gt_1pct', 'e_ac', 'e_an', 'e_filters',
             'in_gnomad', 'g_fail_filter', 'exome_only', 'af_gt_1pct', 'af_gt_1pct_any_group', 'ref_is_minor',
             'gtex_v10_class', 'gtex_v10_sig_tissues', 'gtex_v10_present', 'gtex_v10_maf01', 'gtex_v10_sig_any_genes',
             'gtex_v8_sig_target', 'gtex_v8_finemap_target', 'not_eqtl', 'not_eqtl_strict',
             'ag_set_maxq', 'ag_set_q_at_max', 'ag_set_track_at_max', 'ag_set_rna_maxq', 'ag_set_cage_maxq',
             'ag_set_dnase_maxq', 'ag_expt_maxq', 'ag_all_maxq', 'ag_noeffect_099', 'ag_noeffect_095',
             'ag_all_noeffect_099', 'pai_score', 'pai_maxabs_any_tx', 'pai_noeffect_01', 'pai_noeffect_02',
             'both_noeffect']
    g.sort_values(['window', 'pos', 'ref', 'alt'])[vcols].to_csv(out / 'variants.tsv', sep='\t', index=False,
                                                               float_format='%.6g')
    scols = ['window', 'gene', 'chrom', 'pos', 'ref', 'alt', 'in_gnomad', 'rsid', 'g_af', 'g_grpmax_af', 'af_gt_1pct',
             'af_gt_1pct_any_group', 'gtex_v10_class', 'ag_set_maxq', 'ag_set_track_at_max', 'ag_all_maxq',
             'ag_noeffect_099', 'ag_noeffect_095', 'ag_all_noeffect_099', 'pai_score', 'pai_noeffect_01',
             'pai_noeffect_02', 'both_noeffect']
    snv[scols].to_csv(out / 'all_snvs.tsv', sep='\t', index=False, float_format='%.6g')
    pai_src = Path('/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/promoterai/promoterAI_tss500.tsv.gz')
    summary = dict(
        resources={'gnomAD': 'v4.1 genomes (GraphQL gnomad_r4), queried 2026-11-04',
                   'GTEx': 'Analysis V10 (gtex_v10, GENCODE v39) + V8 (gtex_v8) via Portal API v2, queried 2026-11-04',
                   'AlphaGenome': 'API, client 0.10.0, 1 Mb intervals, recommended RNA_SEQ/CAGE/DNASE scorers',
                   'PromoterAI': f'{pai_src.name} md5 f6ac61c2ea66313222c4ba2d1c94033e (Jaganathan et al. Science 2025)',
                   'TSS': 'Ensembl 116 MANE', 'reference': 'UCSC hg38'},
        gencode_ids=json.load(open(out / 'gtex/gencode_ids.json')), checks=checks,
        md5={p.name: md5(p) for p in [out / 'counts.tsv', out / 'funnel.tsv', out / 'variants.tsv',
                                       out / 'all_snvs.tsv']})
    json.dump(summary, open(out / 'summary.json', 'w'), indent=1)
    print(json.dumps(checks, indent=1))


if __name__ == '__main__':
    main()

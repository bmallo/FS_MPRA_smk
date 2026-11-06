"""gnomAD, GTEx and PromoterAI scan of five immune-gene promoters (settings: rules.json, windows: windows.tsv).

Reuses the promoter variant census code (gnomAD parsing, GTEx and PromoterAI calls); AlphaGenome not applied yet.
Raw pulls first, with tools/promoter_variant_census/fetch.py {gnomad,promoterai,gtex --batched}
--windows-tsv tools/immune_promoter_scan/windows.tsv --output <dir>.
ClinVar (fetch.py clinvar) parsed with the census truth-set code; context only, not a filter.
Writes to --output: variants.tsv, all_snvs.tsv, counts.tsv, common.tsv, eqtl_target.tsv, clinvar_variants.tsv,
summary.json
Usage: python tools/immune_promoter_scan/scan.py --output results/immune_promoter_scan_20261105 [--windows-tsv <file>]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / 'promoter_variant_census'))
import census_lib as cl  # noqa: E402
import assemble as asm  # noqa: E402
import truth_sets as ts  # noqa: E402

GROUP_COLS = ['g_af', 'g_grpmax_af', 'g_grpmax_group']


def offset(pos, w):
    return (pos - w.tss) * (1 if w.strand == '+' else -1)       # upstream negative


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output', required=True)
    ap.add_argument('--windows-tsv', default=str(HERE / 'windows.tsv'), help='windows_tert.tsv for TERT')
    a = ap.parse_args()
    cl.WINDOWS_TSV = Path(a.windows_tsv)
    out = Path(a.output)
    w = cl.windows()
    W = w.set_index('window')
    checks = {}

    # gnomAD (as assemble.py) + joint genomes/exomes AF
    g = cl.gnomad_table(out)
    g = g[(g.g_ac > 0) | (g.e_ac > 0)].copy()
    g['in_gnomad'] = (g.g_ac > 0) & (g.g_filters == '')
    g['g_fail_filter'] = (g.g_ac > 0) & (g.g_filters != '')
    g['exome_only'] = (g.g_ac == 0) & (g.e_ac > 0) & (g.e_filters == '')
    g['joint_af'] = (g.g_ac + g.e_ac) / (g.g_an + g.e_an)
    g['af_gt_1pct'] = g.in_gnomad & (g.g_af > 0.01)
    g['af_gt_1pct_any_group'] = g.in_gnomad & ((g.g_af > 0.01) | (g.g_grpmax_af > 0.01))
    g['singleton'] = g.in_gnomad & (g.g_ac == 1)
    g = asm.gtex_calls(out, g)
    g['not_eqtl'] = ~g.gtex_v10_sig_target
    g['not_eqtl_strict'] = g.gtex_v10_class == 'tested_not_significant'
    checks['gtex_lookup_missing_for_genome_ac1'] = int((~g.gtex_looked_up & (g.g_ac > 0)).sum())
    g['tss_offset'] = [offset(p, W.loc[x]) for p, x in zip(g.pos, g.window)]

    # every possible SNV + PromoterAI (window transcript; max over transcripts as sensitivity)
    snv = cl.all_snvs(out)
    refmap = dict(zip(snv.chrom + ':' + snv.pos.astype(str), snv.ref))
    gs = g[g.kind == 'SNV']
    checks['gnomad_snv_ref_mismatch'] = int(sum(refmap[f'{c}:{p}'] != r for c, p, r in zip(gs.chrom, gs.pos, gs.ref)))
    pai, praw = asm.pai_calls(out, w)
    checks['promoterai_ref_mismatch'] = int(sum(refmap.get(f'{c}:{p}', r) != r for c, p, r in
                                                zip(praw.chrom, praw.pos, praw.ref)))
    snv = snv.merge(pai, on=['window', 'key'], how='left')
    g = g.merge(pai, on=['window', 'key'], how='left')
    for d in (snv, g):
        d['pai_scored'] = d.pai_score.notna()
        d['pai_noeffect_01'] = d.pai_scored & (d.pai_score.abs() < 0.1)
    # ClinVar small variants in the window; class attached to gnomAD SNVs and possible SNVs by key
    cv = ts.clinvar_table(out)
    cv['tss_offset'] = [offset(p, W.loc[x]) if pd.notna(p) else pd.NA for p, x in zip(cv.pos, cv.window)]
    cs = cv[cv.in_window & (cv.kind == 'SNV')]
    checks['clinvar_snv_ref_mismatch'] = int(sum(refmap[f'{c}:{p}'] != r for c, p, r in
                                                 zip(cs.key.str.split(':').str[0], cs.pos, cs.ref)))
    ckey = cs.set_index(['window', 'key'])
    for d in (snv, g):
        ix = list(zip(d.window, d.key))
        d['clinvar_cls'] = [ckey.cls.get(i, '') for i in ix]
        d['clinvar_review'] = [ckey.review.get(i, '') for i in ix]
        d['clinvar'] = [ckey.accession.get(i, '') for i in ix]
    gk = g[g.in_gnomad & (g.kind == 'SNV')].set_index(['window', 'key'])
    snv['in_gnomad'] = [(x, k) in gk.index for x, k in zip(snv.window, snv.key)]
    snv['tss_offset'] = [offset(p, W.loc[x]) for p, x in zip(snv.pos, snv.window)]

    # GTEx V10 target-gene eQTLs (any variant, gnomAD or not)
    sig = pd.read_csv(out / 'gtex/significant_target_gene.tsv', sep='\t')
    s10 = sig[sig.datasetId == 'gtex_v10']
    ev = []
    for x in w.itertuples():
        e = s10[(s10.target_gene == x.gene) & (s10.chromosome == x.chrom)]
        e = e.groupby(['variantId', 'snpId', 'pos']).agg(
            tissues=('tissueSiteDetailId', lambda s: ';'.join(sorted(set(s)))), n_tissues=('tissueSiteDetailId', 'nunique'),
            min_p=('pValue', 'min')).reset_index()
        e['window'] = x.window
        e['tss_offset'] = [offset(p, x) for p in e.pos]
        e['in_window'] = e.pos.between(x.start, x.end)
        ev.append(e)
    ev = pd.concat(ev)

    rows = []
    for x in w.itertuples():
        S, G = snv[snv.window == x.window], g[g.window == x.window]
        C = cv[(cv.window == x.window)]
        Cs = C[C.in_window & (C.kind != 'structural/no SPDI')]
        Gs, Gi = G[G.in_gnomad & (G.kind == 'SNV')], G[G.in_gnomad & (G.kind == 'indel')]
        E = ev[ev.window == x.window]
        n = E.loc[E.tss_offset.abs().idxmin()] if len(E) else None
        r = dict(window=x.window, bp=x.end - x.start + 1, possible_snvs=len(S), gnomad_snvs=len(Gs),
                 gnomad_indels=len(Gi), gnomad_snv_singletons=int(Gs.singleton.sum()),
                 exome_only_snvs=int((G.exome_only & (G.kind == 'SNV')).sum()),
                 excluded_filter=int(G.g_fail_filter.sum()),
                 af1_overall_snvs=int(Gs.af_gt_1pct.sum()), af1_anygroup_snvs=int(Gs.af_gt_1pct_any_group.sum()),
                 af1_anygroup_indels=int(Gi.af_gt_1pct_any_group.sum()),
                 v10_target_eqtl_variants=len(E), v10_target_eqtl_tissues=s10[s10.target_gene == x.gene]
                 .tissueSiteDetailId.nunique(), v10_target_eqtl_in_window=int(E.in_window.sum()),
                 nearest_target_eqtl=None if n is None else f'{n.snpId} ({int(n.tss_offset):+d} bp)')
        for cls in ('significant', 'tested_not_significant', 'not_tested'):
            r[f'gtex_{cls}_gnomad'] = int((G[G.in_gnomad].gtex_v10_class == cls).sum())
        r['gtex_v10_table_gnomad'] = int(G[G.in_gnomad].gtex_v10_present.sum())
        r['gtex_v10_any_gene_eqtl_gnomad'] = int((G[G.in_gnomad].gtex_v10_sig_any_genes != '').sum())
        r['gtex_v8_sig_target_gnomad'] = int(G[G.in_gnomad].gtex_v8_sig_target.sum())
        r['gtex_v8_finemap_target_gnomad'] = int((G[G.in_gnomad].gtex_v8_finemap_target != '').sum())
        sc = S[S.pai_scored]
        r['pai_scored_bp'] = sc.pos.nunique()
        r['pai_unscored_bp'] = r['bp'] - r['pai_scored_bp']
        for t in (0.1, 0.2, 0.5):
            r[f'pai_le_-{t}'] = int((sc.pai_score <= -t).sum())
            r[f'pai_ge_+{t}'] = int((sc.pai_score >= t).sum())
        gsc = Gs[Gs.pai_scored]
        r['gnomad_snvs_pai_scored'] = len(gsc)
        r['gnomad_snvs_pai_abs_ge_0.1'] = int((gsc.pai_score.abs() >= 0.1).sum())
        b1 = Gs[Gs.af_gt_1pct_any_group]
        r['cand_b1'] = len(b1)
        r['cand_b1_b2'] = int(b1.not_eqtl.sum())
        r['cand_b1_b2_pai'] = int((b1.not_eqtl & b1.pai_noeffect_01).sum())
        r['cand_b1_b2strict_pai'] = int((b1.not_eqtl_strict & b1.pai_noeffect_01).sum())
        r['clinvar_records'] = len(C)
        r['clinvar_structural'] = int((C.kind == 'structural/no SPDI').sum())
        for cls in ('plp', 'blb', 'vus', 'conflicting', 'other'):
            for k_ in ('SNV', 'indel'):
                r[f'clinvar_{cls}_{k_.lower()}'] = int(((Cs.cls == cls) & (Cs.kind == k_)).sum())
        r['clinvar_plp_snv_in_gnomad'] = int((Gs.clinvar_cls == 'plp').sum())
        r['clinvar_blb_snv_in_gnomad'] = int((Gs.clinvar_cls == 'blb').sum())
        r['clinvar_plp_snv_b1'] = int((b1.clinvar_cls == 'plp').sum())
        cand = b1[b1.not_eqtl & b1.pai_noeffect_01]
        r['cand_b1_b2_pai_clinvar_blb'] = int((cand.clinvar_cls == 'blb').sum())
        rows.append(r)
    counts = pd.DataFrame(rows)

    vcols = ['window', 'gene', 'gnomad_id', 'rsid', 'kind', 'chrom', 'pos', 'tss_offset', 'ref', 'alt', 'g_ac', 'g_an',
             'g_af', 'g_hom', 'g_filters', 'g_grpmax_af', 'g_grpmax_group', 'g_groups_af_gt_1pct', 'e_ac', 'e_an',
             'e_filters', 'joint_af', 'in_gnomad', 'g_fail_filter', 'exome_only', 'af_gt_1pct', 'af_gt_1pct_any_group',
             'gtex_v10_class', 'gtex_v10_sig_tissues', 'gtex_v10_present', 'gtex_v10_maf01', 'gtex_v10_sig_any_genes',
             'gtex_v8_sig_target', 'gtex_v8_finemap_target', 'not_eqtl', 'not_eqtl_strict', 'pai_score',
             'pai_maxabs_any_tx', 'pai_scored', 'pai_noeffect_01', 'clinvar_cls', 'clinvar_review', 'clinvar']
    g = g.sort_values(['window', 'pos', 'ref', 'alt'])[vcols]
    g.to_csv(out / 'variants.tsv', sep='\t', index=False, float_format='%.6g')
    common = g[g.af_gt_1pct_any_group]
    common.to_csv(out / 'common.tsv', sep='\t', index=False, float_format='%.6g')
    snv[['window', 'gene', 'chrom', 'pos', 'tss_offset', 'ref', 'alt', 'in_gnomad', 'pai_score', 'pai_maxabs_any_tx',
         'pai_scored', 'clinvar_cls', 'clinvar']].to_csv(out / 'all_snvs.tsv', sep='\t', index=False, float_format='%.6g')
    ev.sort_values(['window', 'pos']).to_csv(out / 'eqtl_target.tsv', sep='\t', index=False, float_format='%.4g')
    counts.to_csv(out / 'counts.tsv', sep='\t', index=False)
    cv.drop(columns='key').to_csv(out / 'clinvar_variants.tsv', sep='\t', index=False)
    json.dump(dict(resources={'gnomAD': 'v4.1 (gnomad_r4) GraphQL, queried 2026-11-05',
                              'GTEx': 'V10 + V8 Portal API v2 (batched variant-table lookup), queried 2026-11-05',
                              'ClinVar': 'NCBI E-utilities esearch + esummary, queried 2026-11-05',
                              'PromoterAI': 'promoterAI_tss500.tsv.gz md5 f6ac61c2ea66313222c4ba2d1c94033e',
                              'reference': 'UCSC hg38'},
                   gencode_ids=json.load(open(out / 'gtex/gencode_ids.json')), checks=checks,
                   md5={p: asm.md5(out / p) for p in ['variants.tsv', 'all_snvs.tsv', 'counts.tsv', 'common.tsv',
                                                       'eqtl_target.tsv', 'clinvar_variants.tsv']}),
              open(out / 'summary.json', 'w'), indent=1)
    print(json.dumps(checks, indent=1))
    print(counts.T.to_string())


if __name__ == '__main__':
    main()

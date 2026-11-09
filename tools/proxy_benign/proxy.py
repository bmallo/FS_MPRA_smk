"""Proxy benign evidence level per possible SNV for the 16 elements (definitions: rules.json).

Classes: in gnomAD (any frequency); tested in GTEx V10 and not an eQTL for any gene in any tissue; no large AlphaGenome
change (element's track set, top-20% rule). SNVs with a >= 1-star ClinVar P/LP or B/LB record are the benchmark and are
left out. Inputs under --root (fetch.py gnomad/gtex --batched/clinvar, ag_array.sbatch): gnomad/, gtex/, clinvar/,
reference/, ag/<element>/alphagenome/tracks.tsv.
Writes to --root: snv_proxy.tsv (one row per possible SNV), counts.tsv (per element), summary.json
--benchmark <clinvar_annotated.tsv> (clinvar_noncoding.py): benchmark = its 'benchmark' rows (v3: protein-changing
records dropped); --output: where to write (default --root).
Usage: python tools/proxy_benign/proxy.py --root results/proxy_benign_20261106 [--benchmark <tsv> --output <dir>]
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / 'promoter_variant_census'))
import census_lib as cl  # noqa: E402
import truth_sets as ts  # noqa: E402

cl.WINDOWS_TSV = HERE / 'elements.tsv'
TOP = 0.8


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


def stars(review):
    r = str(review)
    if r in ('practice guideline', 'reviewed by expert panel', 'criteria provided, multiple submitters, no conflicts'):
        return 2
    return 1 if r.startswith('criteria provided') else 0


def ag_class(root, w, keys):
    t = pd.read_csv(root / 'ag' / w.window / 'alphagenome/tracks.tsv', sep='\t',
                    usecols=['key', 'output_type', 'gene_name', 'ontology_curie', 'raw_score'])
    rna, dn = str(w.ag_rna).split(';'), str(w.ag_dnase).split(';')
    r = t[(t.output_type == 'RNA_SEQ') & (t.gene_name == w.gene) & t.ontology_curie.isin(rna)]
    n = t[(t.output_type == 'DNASE') & t.ontology_curie.isin(dn)]
    e = pd.DataFrame({'key': keys})
    e['ag_rna'] = e.key.map(r.groupby('key').raw_score.mean())
    e['ag_dnase'] = e.key.map(n.groupby('key').raw_score.mean())
    for c in ('ag_rna', 'ag_dnase'):
        e[f'{c}_large'] = e[c].abs() >= e[c].abs().quantile(TOP)
    e['ag_no_change'] = ~e.ag_rna_large & ~e.ag_dnase_large
    return e, dict(rna_tracks=int(r.ontology_curie.nunique()), dnase_tracks=int(n.ontology_curie.nunique()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--benchmark', help='clinvar_annotated.tsv from clinvar_noncoding.py')
    ap.add_argument('--output', help='output directory (default --root)')
    a = ap.parse_args()
    root = Path(a.root)
    out_dir = Path(a.output) if a.output else root
    out_dir.mkdir(parents=True, exist_ok=True)
    W = cl.windows()
    snv = cl.all_snvs(root)
    g = cl.gnomad_table(root)
    g = g[g.kind == 'SNV'].drop_duplicates(['window', 'key']).set_index(['window', 'key'])
    ac = g.g_ac.where(g.g_filters == '', 0) + g.e_ac.where(g.e_filters == '', 0)   # FILTER PASS alleles only
    lk = pd.read_csv(root / 'gtex/variant_lookup.tsv', sep='\t').set_index('variant_id')
    cv = ts.clinvar_table(root)
    cv['stars'] = cv.review.map(stars)
    cv = cv[cv.in_window.astype(bool)]
    bench = cv[cv.cls.isin(['plp', 'blb']) & (cv.stars >= 1)]
    if a.benchmark:
        ann = pd.read_csv(a.benchmark, sep='\t')
        cv = ann
        bench = ann[ann.benchmark]
    out, rows, ntracks = [], [], {}
    for w in W.itertuples():
        s = snv[snv.window == w.window].copy()
        ix = list(zip(s.window, s.key))
        s['gnomad_ac'] = [int(ac.get(i, 0)) for i in ix]
        s['in_gnomad'] = s.gnomad_ac > 0
        gid = s.chrom + '_' + s.pos.astype(str) + '_' + s.ref + '_' + s.alt + '_b38'
        s['gtex_tested'] = gid.map(lk.gtex_v10_maf01).fillna(False).astype(bool)
        s['gtex_eqtl_any'] = gid.map(lk.v10_sig_pairs).fillna(0) > 0
        s['gtex_no_eqtl'] = s.gtex_tested & ~s.gtex_eqtl_any
        e, ntracks[w.window] = ag_class(root, w, s.key)
        s = s.merge(e, on='key', how='left')
        b = bench[(bench.window == w.window) & (bench.kind == 'SNV')]
        s['clinvar_benchmark'] = s.key.map(b.drop_duplicates('key').set_index('key').cls).fillna('')
        s['level'] = s[['in_gnomad', 'gtex_no_eqtl', 'ag_no_change']].sum(axis=1)
        s.loc[s.clinvar_benchmark != '', 'level'] = np.nan
        out.append(s)
        c = cv[cv.window == w.window]
        r = dict(element=w.window, possible_snvs=len(s), benchmark_snvs=int((s.clinvar_benchmark != '').sum()),
                 rna_tracks=ntracks[w.window]['rna_tracks'], dnase_tracks=ntracks[w.window]['dnase_tracks'])
        for cls in ('plp', 'blb'):
            x = c[(c.cls == cls) & (c.stars >= 1)]
            if 'protein_changing' in x:                     # v3: protein-changing records are not benchmark
                x = x[~x.protein_changing]
            r[f'clinvar_{cls}_1star_snv'] = int((x.kind == 'SNV').sum())
            r[f'clinvar_{cls}_1star_all'] = len(x)
        P = s[s.clinvar_benchmark == '']
        for k in ('in_gnomad', 'gtex_no_eqtl', 'ag_no_change'):
            r[f'class_{k}'] = int(P[k].sum())
        for lv in range(4):
            r[f'level{lv}'] = int((P.level == lv).sum())
            r[f'level{lv}_pct'] = round(100 * r[f'level{lv}'] / len(s), 2)
        rows.append(r)
    S, counts = pd.concat(out), pd.DataFrame(rows)
    checks = dict(snvs_missing_ag=int(S[['ag_rna', 'ag_dnase']].isna().any(axis=1).sum()),
                  benchmark_snvs_unmatched=int((~bench[bench.kind == 'SNV'].key.isin(S.key)).sum()),
                  gtex_lookup_missing=int((~(S.chrom + '_' + S.pos.astype(str) + '_' + S.ref + '_' + S.alt + '_b38')
                                           .isin(lk.index) & S.in_gnomad & S.key.isin(
                                               g[g.g_ac > 0].index.get_level_values(1))).sum()))
    S.to_csv(out_dir / 'snv_proxy.tsv', sep='\t', index=False, float_format='%.6g')
    counts.to_csv(out_dir / 'counts.tsv', sep='\t', index=False)
    json.dump(dict(benchmark=a.benchmark or 'built in (>= 1-star P/LP and B/LB)', checks=checks, tracks=ntracks,
                   md5={p: md5(out_dir / p) for p in ('snv_proxy.tsv', 'counts.tsv')}),
              open(out_dir / 'summary.json', 'w'), indent=1)
    print(json.dumps(checks))
    print(counts.to_string(index=False))


if __name__ == '__main__':
    main()

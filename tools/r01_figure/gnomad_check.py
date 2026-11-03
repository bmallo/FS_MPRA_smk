"""R01 Aim 1B preview: do gnomAD-observed LDLR promoter SNVs behave as benign in the existing data?
Settings in tools/r01_figure/gnomad_rules.json (written before the download).

Inputs: gnomAD v4.1 region JSON (GraphQL download), the step-3 clinical table, Kircher 2019 activity, fp_features
effects, and the nucleosome-size endpoint computed here for all 633 SNVs with repeat3.hepg2 (A7 and +/-5 masks).
Writes snvs.tsv (one row per designed SNV with gnomAD fields, classes and measures), classes.tsv (abnormal counts and
Wilson 95% CI per class and measure), nucleosome_all.tsv (repeat3 endpoint rows for all SNVs), summary.json.
Usage: python tools/r01_figure/gnomad_check.py --root results/analysis_v2_20261028 --gnomad <json>
           --clinical <step3 dir> --procs 4 --output <dir>
"""
import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent


def _load(name):
    s = importlib.util.spec_from_file_location(name, HERE / f'{name}.py')
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m                    # so worker processes can pickle its functions
    s.loader.exec_module(m)
    return m


RULES = json.load(open(HERE / 'gnomad_rules.json'))
OLIGO_OFFSET = 11089264
P_ACT, Z_NUC = 1e-5, 3.29


def gnomad_table(path):
    v = json.load(open(path))['data']['region']['variants']
    rows = []
    for x in v:
        if len(x['ref']) != 1 or len(x['alt']) != 1:
            continue
        ds = [d for d in (x['genome'], x['exome']) if d]
        ac = sum(d['ac'] for d in ds)
        an = sum(d['an'] for d in ds)
        rows.append(dict(name=f"{x['pos'] - OLIGO_OFFSET}{x['ref']}>{x['alt']}", oligo=x['pos'] - OLIGO_OFFSET,
                         gnomad_id=x['variant_id'], rsid=';'.join(x['rsids']), ac=ac, an=an,
                         af=ac / an if an else np.nan, hom=sum(d['homozygote_count'] for d in ds),
                         filtered=any(d['filters'] for d in ds if d['ac'] > 0),
                         in_genomes=bool(x['genome'] and x['genome']['ac'] > 0),
                         in_exomes=bool(x['exome'] and x['exome']['ac'] > 0)))
    return pd.DataFrame(rows)


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / d
    return c - h, c + h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--gnomad', required=True)
    ap.add_argument('--clinical', required=True)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    t0 = time.time()
    ro, r3 = _load('reorient'), _load('repeat3')

    fn = pd.read_csv(root / 'a11_motifs_function/function_snvs.tsv', sep='\t').set_index('snv')
    snvs = fn.index.tolist()
    g = gnomad_table(a.gnomad)
    n_raw = len(g)
    g0 = g[g.ac == 0]
    g = g[g.ac > 0]
    in_el = g[g.oligo.between(21, 231)]
    undesigned = sorted(set(in_el.name) - set(snvs))
    g = in_el.set_index('name')

    t = fn[['pos', 'act', 'act2', 'C1']].join(ro.activity()[['act_p', 'act2_p']])
    t = t.join(g[['gnomad_id', 'rsid', 'ac', 'an', 'af', 'hom', 'filtered', 'in_genomes', 'in_exomes']])
    t['ac'] = t.ac.fillna(0).astype(int)
    c = pd.read_csv(Path(a.clinical) / 'clinical_variants.tsv', sep='\t')
    c = c[c.type == 'SNV'].set_index('designed').clinvar_class.fillna('')
    t['clinvar'] = t.index.map(c).fillna('')

    nuc, _ = r3.hepg2(root, {'ALL': snvs}, a.procs, t0)
    nuc.to_csv(out / 'nucleosome_all.tsv', sep='\t', index=False, float_format='%.5g')
    n7 = nuc[(nuc.endpoint == 'nucleosome') & (nuc.calls == 'confident') & (nuc['mask'] == 'A7')].set_index('snv')
    t['nuc_D_pts'] = 100 * n7.D
    t['nuc_se_pts'] = 100 * n7.se
    t['nuc_z'] = n7.D / n7.se
    fe = pd.read_parquet(root / 'fp_features/effects.parquet').query('mask == "A7" and calls == "confident" and sig')
    t['fp_sig'] = t.index.isin(fe.snv)
    t['fp_sig_list'] = t.index.map(fe.assign(x=fe.footprint + ':' + (100 * fe.D).round(2).astype(str))
                                   .groupby('snv').x.agg(';'.join)).fillna('')

    t['ab_activity'] = (t.act_p < P_ACT) & (np.sign(t.act) == np.sign(t.act2))
    t['ab_nucleosome'] = t.nuc_z.abs() > Z_NUC
    t['ab_footprint'] = t.fp_sig
    t['ab_either'] = t.ab_activity | t.ab_nucleosome
    t['ab_chromatin_any'] = t.ab_nucleosome | t.ab_footprint

    af = t.af.fillna(0)
    classes = {
        'clinvar_plp': t.clinvar.isin(['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic']),
        'clinvar_blb': t.clinvar.isin(['Benign', 'Likely benign']),
        'gnomad_any': t.ac >= 1, 'gnomad_singleton': t.ac == 1, 'gnomad_ac2plus': t.ac >= 2,
        'gnomad_af_1e-4': af >= 1e-4, 'gnomad_af_1e-3': af >= 1e-3, 'gnomad_common_1pct': af >= 0.01,
        'gnomad_homozygote': t.hom.fillna(0) >= 1, 'not_in_gnomad': t.ac == 0,
        'all': pd.Series(True, index=t.index)}
    for k, m in classes.items():
        t[f'class_{k}'] = m.values
    t.to_csv(out / 'snvs.tsv', sep='\t', float_format='%.5g', index_label='snv')

    rows = []
    for k, m in classes.items():
        s = t[m.values]
        r = dict(cls=k, n=len(s), act_median=s.act.median(), nuc_median_pts=s.nuc_D_pts.median(),
                 members=';'.join(s.index) if len(s) <= 12 else '')
        for meas in ('activity', 'nucleosome', 'footprint', 'either', 'chromatin_any'):
            kk = int(s[f'ab_{meas}'].sum())
            lo, hi = wilson(kk, len(s))
            r.update({f'{meas}_abnormal': kk, f'{meas}_frac': kk / len(s) if len(s) else np.nan,
                      f'{meas}_lo': lo, f'{meas}_hi': hi})
        rows.append(r)
    pd.DataFrame(rows).to_csv(out / 'classes.tsv', sep='\t', index=False, float_format='%.4g')
    json.dump(dict(gnomad_snv_records=n_raw, ac0_dropped=len(g0), in_element=len(in_el),
                   in_element_not_designed=undesigned, filtered_in_element=int(in_el.filtered.sum()),
                   seconds=round(time.time() - t0)), open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

"""One table of every evidence type for the 633 designed LDLR promoter SNVs (oligo 21-231). No new tests.

Database / predictor evidence: ClinVar + LOVD (step 3), gnomAD v4.1 (gnomad_check), GTEx eQTLs (v8 and v10 API
downloads), and optional AlphaGenome / PromoterAI score files (TSV with columns snv and one or more score columns;
merged with the prefix ag_ / pai_ when given).
Functional evidence:
  HepG2 MPRA        Kircher 2019 log2 effect and P (LDLR, LDLR.2)
  HepG2 Fiber-seq   nucleosome-size change (gnomad_check, A7), fp_features significant footprints and F04 D,
                    A11 element-wide shift C1, local accessibility = mean A10 E1 D (inaccessible fraction, A7,
                    testable adenines) within 10 bp of the SNV, A10 global selection and hit count
  HEK293T Fiber-seq per SNV, molecules and clones (labeling > 0.01), and the mean labeling-standardised
                    inaccessible-fraction residual vs HEK WT over molecules (local: adenines within 10 bp; far:
                    adenines > 50 bp away; edited A/T unobserved). Descriptive only: most SNVs are one clone.
  Visibility        A4: caller drops the edited base; A/T count within the local window
Writes evidence.tsv and summary.json to --output.
Usage: python tools/r01_figure/merge_evidence.py --root results/analysis_v2_20261028 --work results/r01_figure_20261102
           [--alphagenome f.tsv] [--promoterai f.tsv] --output <dir>
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
OLIGO_OFFSET, C_OFFSET = 11089264, 285


def hek_per_snv(snvs):
    from tools.analysis_v2 import access
    a9f = access.a9f
    d = access.load_element('HEK293T')
    ref = access.E1Reference(d)
    use = d.distal > a9f.E1_MIN
    pm = ref.p[a9f.strata(d.distal)]
    resid = np.where((d.X >= 0) & use[:, None] & np.isfinite(pm), (d.X == 0) - np.nan_to_num(pm), np.nan)
    rows = []
    for v in snvs:
        idx = np.flatnonzero((d.label == v) & use)
        if idx.size == 0:
            rows.append(dict(snv=v, hek_molecules=0, hek_clones=0))
            continue
        k, rb = int(v[:-3]), v[-3]
        R = resid[idx].copy()
        if rb in 'AT':
            R[:, d.oligo == k] = np.nan
        loc = np.abs(d.oligo - k) <= 10
        far = np.abs(d.oligo - k) > 50
        rows.append(dict(snv=v, hek_molecules=int(idx.size), hek_clones=int(np.unique(d.clone[idx]).size),
                         hek_local_pts=100 * np.nanmean(R[:, loc]) if loc.any() else np.nan,
                         hek_far_pts=100 * np.nanmean(R[:, far]) if far.any() else np.nan))
    return pd.DataFrame(rows).set_index('snv')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--work', required=True, help='results/r01_figure_<date> (step3_clinical, gnomad_check, predictors_raw)')
    ap.add_argument('--alphagenome')
    ap.add_argument('--promoterai')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, w, out = Path(a.root), Path(a.work), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)

    g = pd.read_csv(w / 'gnomad_check/snvs.tsv', sep='\t').set_index('snv')
    t = g[['pos', 'act', 'act_p', 'act2', 'act2_p', 'ac', 'af', 'hom', 'clinvar', 'nuc_D_pts', 'nuc_z', 'fp_sig_list',
           'C1']].rename(columns={'pos': 'oligo', 'ac': 'gnomad_ac', 'af': 'gnomad_af', 'hom': 'gnomad_hom'})
    t['c_dot'] = t.oligo - C_OFFSET
    t['hg38'] = t.oligo + OLIGO_OFFSET
    c = pd.read_csv(w / 'step3_clinical/clinical_variants.tsv', sep='\t')
    c = c[c.type == 'SNV'].set_index('designed')
    t['lovd_effect'] = c.lovd_effect.reindex(t.index)

    eq = []
    for f in sorted((w / 'predictors_raw/gtex').glob('eqtl_*.json')):
        eq += [x for x in json.load(open(f)).get('data', [])]
    hit = {}
    for x in eq:
        k = x['pos'] - OLIGO_OFFSET
        if 21 <= k <= 231:
            hit.setdefault(f"{k}{x['variantId'].split('_')[2]}>{x['variantId'].split('_')[3]}", []).append(
                x['tissueSiteDetailId'])
    t['gtex_eqtl_tissues'] = t.index.map(lambda v: ';'.join(sorted(set(hit.get(v, [])))))
    gtex = dict(eqtl_rows=len(eq), in_element=sum(len(v) for v in hit.values()),
                nearest_bp=int(min(abs(x['pos'] - (OLIGO_OFFSET + 126)) for x in eq)) - 105 if eq else None,
                tissues=sorted({x['tissueSiteDetailId'] for x in eq}))

    fe = pd.read_parquet(root / 'fp_features/effects.parquet').query('mask == "A7" and calls == "confident"')
    f4 = fe[fe.footprint == 'F04'].set_index('snv')
    t['F04_D_pts'], t['F04_p'] = 100 * f4.D, f4.p
    e = pd.read_parquet(root / 'a10_effects/HepG2_run3/effects.parquet', columns=['snv', 'endpoint', 'oligo', 'snv_oligo',
                        'D', 'p'], filters=[('endpoint', '==', 'E1')])
    e = e[e.p.notna() & ((e.oligo - e.snv_oligo).abs() <= 10)]
    t['local_inacc_pts'] = 100 * e.groupby('snv').D.mean()
    s7 = pd.read_csv(root / 'a10_effects/HepG2_run3/snvs.tsv', sep='\t').query('endpoint == "E1"').set_index('snv')
    t['a10_selected'], t['a10_hits'] = s7.selected, s7.hits
    vis = pd.read_csv(root / 'a4_informative/visibility/variants_visibility.tsv', sep='\t').query('cls == "SNV"')
    vis = vis.set_index('source_variant')
    t['drops_its_base'] = vis.caller_drops_at.reindex(t.index).astype(bool)
    t['at_local_wt'] = vis.at_local_wt.reindex(t.index)
    d = pd.read_csv(root / 'a11_motifs_function/dose_snvs.tsv', sep='\t')
    t['sp1_r3_dpwm'] = d[d.motif == 'Repeat 3 (Sp1)'].set_index('snv').d_pwm.reindex(t.index)

    t = t.join(hek_per_snv(t.index.tolist()))
    for path, pre in ((a.alphagenome, 'ag_'), (a.promoterai, 'pai_')):
        if path:
            p = pd.read_csv(path, sep='\t').set_index('snv')
            t = t.join(p.add_prefix(pre))
    t.to_csv(out / 'evidence.tsv', sep='\t', float_format='%.5g', index_label='snv')
    json.dump(dict(snvs=len(t), gtex=gtex, hek_snvs=int((t.hek_molecules > 0).sum()),
                   hek_ge3_clones=int((t.hek_clones >= 3).sum()), alphagenome=bool(a.alphagenome),
                   promoterai=bool(a.promoterai)), open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

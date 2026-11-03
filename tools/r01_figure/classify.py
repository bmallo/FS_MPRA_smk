"""Mechanism categories for LDLR promoter SNVs from MPRA + plasmid Fiber-seq + motif evidence (classify_rules.json).

Writes variants.tsv (every designed SNV: evidence columns, category, evidence level), summary.tsv (counts per ClinVar
class x category, primary and |log2| >= 1 sensitivity) and summary.json to --output.
Usage: python tools/r01_figure/classify.py --root results/analysis_v2_20261028 --work results/r01_figure_20261102
           --db results/r01_figure_20261102/motif_db --output <dir>
"""
import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).parent
RULES = json.load(open(HERE / 'classify_rules.json'))
WIN, P_MOTIF, FOLD = 10, 1e-4, 10
ACTIVATOR_FAMILY = re.compile(r'^(SP[1-4]|KLF\d*|SREBF[12]|SRBP[12])$')
PLP = ['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic']
BLB = ['Benign', 'Likely benign']


def _load(name):
    s = importlib.util.spec_from_file_location(name, HERE / f'{name}.py')
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m
    s.loader.exec_module(m)
    return m


def motif_alleles(snvs, ref, db):
    ms, cl = _load('motif_scan'), _load('clinical')
    bg = np.array([ref.count(b) for b in ms.BASES], float)
    bg = (bg + bg[::-1]) / (2 * len(ref))
    fam = {}
    for line in open(db / ms.HOCO_ANN):
        r = json.loads(line)
        fam[r['name']] = r['masterlist_info']['tf']
    motifs = ms.read_motifs(db / ms.JASPAR, 'JASPAR2026') + ms.read_motifs(db / ms.HOCO, 'HOCOMOCOv14')
    for m in motifs:
        m['mn'], m['sf'], m['int'] = ms.pvalue_table(ms.log_odds(m['counts'], bg), bg)
        m['tf'] = fam.get(m['id'], m['name']).upper()
    wt_best = {}
    rows = []
    for v in snvs:
        k, alt = int(v[:-3]), v[-1]
        seq = cl.edit(ref, k, k, alt)
        brk, crt = [], []
        for m in motifs:
            L = m['int'].shape[0]
            res = []
            for which, sq in (('wt', ref), ('var', seq)):
                if which == 'wt' and (m['id'], k) in wt_best:
                    res.append(wt_best[(m['id'], k)])
                    continue
                pf, pr = ms.scan(sq, m['int'], m['mn'], m['sf'])
                st = np.arange(len(pf)) + cl.LO
                cov = (st <= k) & (st + L - 1 >= k)
                b = min(pf[cov].min(), pr[cov].min()) if cov.any() else 1.0
                if which == 'wt':
                    wt_best[(m['id'], k)] = b
                res.append(b)
            pw, pv = res
            if pw < P_MOTIF and pv >= FOLD * pw and ACTIVATOR_FAMILY.match(m['tf']):
                brk.append((pw, m['tf']))
            if pv < P_MOTIF and pw >= FOLD * pv:
                crt.append((pv, m['tf']))
        rows.append(dict(snv=v, breaks_activator_motif=bool(brk), broken=';'.join(sorted({t for _, t in brk})),
                         creates_motif=bool(crt), created=';'.join(t for _, t in sorted(crt)[:3])))
    return pd.DataFrame(rows).set_index('snv')


def footprint_signal(root, snvs):
    e = pd.read_parquet(root / 'a10_effects/HepG2_run3/effects.parquet',
                        columns=['snv', 'endpoint', 'oligo', 'snv_oligo', 'bin', 'D', 'se', 'masked', 'hit'],
                        filters=[('endpoint', '==', 'E2'), ('bin', 'in', ['sub_TF', 'TF'])])
    e = e[e.snv.isin(snvs) & ((e.oligo - e.snv_oligo).abs() <= WIN) & ~e.masked & e.se.gt(0)]
    e['z'] = e.D / e.se
    fp = pd.read_csv(root / 'fp_features/footprints.tsv', sep='\t')
    tf = fp[(fp.len_lo <= 45) & (fp.len_hi >= 17)].set_index('feature')
    fe = pd.read_parquet(root / 'fp_features/effects.parquet').query('mask == "A7" and calls == "confident" and sig')
    fe = fe[fe.footprint.isin(tf.index)].join(tf[['span_lo', 'span_hi']], on='footprint')
    fe['pos'] = fe.snv.str[:-3].astype(int)
    fe = fe[(fe.pos >= fe.span_lo - WIN) & (fe.pos <= fe.span_hi + WIN)]
    rows = []
    for v in snvs:
        g = e[e.snv == v]
        h = g[g.hit]
        f = fe[fe.snv == v]
        level, sign = 'none', 0
        if len(h) or len(f):
            s = np.sign(h.D).sum()
            sign = int(np.sign(s)) if s != 0 else int(np.sign(f.D.sum()))
            level = 'significant'
        elif len(g) and abs(g.z.mean()) >= 1:
            level, sign = 'supporting', int(np.sign(g.z.mean()))
        rows.append(dict(snv=v, fp_level=level, fp_sign=sign, fp_mean_z=g.z.mean() if len(g) else np.nan,
                         fp_rows=len(g), fp_hits=len(h), fp_features_sig=';'.join(f.footprint)))
    return pd.DataFrame(rows).set_index('snv')


def classify(r, min_abs):
    expr = 0
    if abs(r.act) >= min_abs and r.act_p < 1e-5 and np.sign(r.act) == np.sign(r.act2):
        expr = int(np.sign(r.act))
    sig, sup = r.fp_level == 'significant', r.fp_level == 'supporting'
    if expr < 0:
        if sig and r.fp_sign < 0:
            return 'Activator disrupted', 'footprint (significant)'
        if sig and r.fp_sign > 0 and not r.breaks_activator_motif:
            return 'Repressor gained', 'footprint (significant)'
        if r.breaks_activator_motif:
            return 'Activator disrupted', 'motif broken'
        if sup and r.fp_sign < 0:
            return 'Activator disrupted', 'footprint (supporting)'
        if sup and r.fp_sign > 0:
            return 'Repressor gained', 'footprint (supporting)'
        return 'Expression down, mechanism unresolved', 'nucleosome closes' if r.nuc_call == 'closes' else '-'
    if expr > 0:
        if sig and r.fp_sign > 0:
            return 'Activator gained', 'footprint (significant)'
        if sig and r.fp_sign < 0:
            return 'Repressor lost', 'footprint (significant)'
        if r.creates_motif:
            return 'Activator gained', 'motif created'
        if sup and r.fp_sign > 0:
            return 'Activator gained', 'footprint (supporting)'
        if sup and r.fp_sign < 0:
            return 'Repressor lost', 'footprint (supporting)'
        return 'Expression up, mechanism unresolved', 'nucleosome opens' if r.nuc_call == 'opens' else '-'
    if r.nuc_call != 'none' or sig:
        return 'Chromatin change only', ('nucleosome ' + r.nuc_call) if r.nuc_call != 'none' else 'footprint (significant)'
    if r.poorly_visible:
        return 'Not assessable by Fiber-seq', 'few A/T within 5 bp'
    return 'No detectable effect', '-'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--work', required=True)
    ap.add_argument('--db', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, w, out = Path(a.root), Path(a.work), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    cl = _load('clinical')
    ref, _ = cl.reference(root)
    t = pd.read_csv(w / 'evidence_v2/evidence.tsv', sep='\t').set_index('snv')
    t['clin'] = np.select([t.clinvar.isin(PLP), t.clinvar.isin(BLB), t.clinvar.fillna('').str.startswith('Conflicting'),
                           t.clinvar == 'Uncertain significance'], ['P/LP', 'B/LB', 'Conflicting', 'VUS'], 'not in ClinVar')
    snvs = t.index.tolist()

    cv = t.clin != 'not in ClinVar'
    p = pd.Series(2 * stats.norm.sf(t.nuc_z.abs()), index=t.index)
    t['nuc_q'] = np.nan
    t.loc[cv, 'nuc_q'] = stats.false_discovery_control(p[cv].fillna(1))
    t['nuc_call'] = np.where(t.nuc_q < 0.05, np.where(t.nuc_z > 0, 'closes', 'opens'), 'none')
    t['at_pm5'] = [sum(ref[j - cl.LO] in 'AT' for j in range(max(cl.LO, k - 5), min(cl.HI, k + 5) + 1) if j != k)
                   for k in t.oligo]
    t['poorly_visible'] = t.at_pm5 <= 2
    t = t.join(footprint_signal(root, snvs)).join(motif_alleles(snvs, ref, Path(a.db)))
    for tag, m in (('', 0.5), ('_abs1', 1.0)):
        cats = t.apply(lambda r: classify(r, m), axis=1)
        t[f'category{tag}'] = [c[0] for c in cats]
        t[f'evidence{tag}'] = [c[1] for c in cats]
    t.to_csv(out / 'variants.tsv', sep='\t', float_format='%.5g', index_label='snv')
    rows = []
    for tag in ('', '_abs1'):
        x = t[cv].groupby(['clin', f'category{tag}']).size().rename('n').reset_index()
        x = x.rename(columns={f'category{tag}': 'category'}).assign(threshold='|log2|>=0.5' if not tag else '|log2|>=1')
        rows.append(x)
    pd.concat(rows).to_csv(out / 'summary.tsv', sep='\t', index=False)
    json.dump(dict(clinvar_snvs=int(cv.sum()), poorly_visible=int(t[cv].poorly_visible.sum()),
                   fp_levels=t[cv].fp_level.value_counts().to_dict(),
                   breaks_activator=int(t[cv].breaks_activator_motif.sum()), creates=int(t[cv].creates_motif.sum()),
                   nuc_calls=t[cv].nuc_call.value_counts().to_dict()), open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

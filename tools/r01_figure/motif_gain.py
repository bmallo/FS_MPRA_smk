"""Motifs created by one SNV, with the element-wide base rate (exploratory; motif_scan.py scoring).

For every designed SNV: each motif's best match P covering the SNV, WT vs variant (JASPAR 2026 + HOCOMOCO v14,
FIMO-style exact P). 'Created' = variant P < 1e-4 and >= 10-fold better than WT. Writes:
  all_snvs.tsv   per SNV: number of created motifs, best created P and TF
  target.tsv     every motif the target SNV improves to P < 1e-4 (any fold), with the variant match span (c.),
                 fold change, HOCOMOCO TF class/family, and overlap with --protected (oligo lo-hi)
Usage: python tools/r01_figure/motif_gain.py --snv 129C>T --protected 120 142 --root results/analysis_v2_20261028
           --db results/r01_figure_20261102/motif_db --output <dir>
"""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
P_HIT, FOLD, C0 = 1e-4, 10, 285


def _load(name):
    s = importlib.util.spec_from_file_location(name, HERE / f'{name}.py')
    m = importlib.util.module_from_spec(s)
    sys.modules[name] = m
    s.loader.exec_module(m)
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--snv', required=True)
    ap.add_argument('--protected', nargs=2, type=int, required=True)
    ap.add_argument('--root', required=True)
    ap.add_argument('--db', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out, db = Path(a.output), Path(a.db)
    out.mkdir(parents=True, exist_ok=False)
    ms, cl = _load('motif_scan'), _load('clinical')
    ref, vis = cl.reference(Path(a.root))
    snvs = vis[vis.cls == 'SNV'].source_variant.tolist()
    bg = np.array([ref.count(b) for b in ms.BASES], float)
    bg = (bg + bg[::-1]) / (2 * len(ref))
    ann = {}
    for line in open(db / ms.HOCO_ANN):
        r = json.loads(line)
        mi = r['masterlist_info']
        ann[r['name']] = (mi['tf'], mi.get('tfclass_class', ''), mi.get('tfclass_family', ''))
    motifs = ms.read_motifs(db / ms.JASPAR, 'JASPAR2026') + ms.read_motifs(db / ms.HOCO, 'HOCOMOCOv14')
    for m in motifs:
        m['mn'], m['sf'], m['int'] = ms.pvalue_table(ms.log_odds(m['counts'], bg), bg)
        m['tf'], m['class'], m['family'] = ann.get(m['id'], (m['name'], '', ''))
        pf, pr = ms.scan(ref, m['int'], m['mn'], m['sf'])
        m['wt'] = np.minimum(pf, pr)

    def best(m, seq, k):
        L = m['int'].shape[0]
        pf, pr = ms.scan(seq, m['int'], m['mn'], m['sf'])
        p = np.minimum(pf, pr)
        st = np.arange(len(p)) + cl.LO
        cov = (st <= k) & (st + L - 1 >= k)
        i = np.flatnonzero(cov)[np.argmin(p[cov])]
        return p[i], st[i], st[i] + L - 1, i, cov

    rows, target = [], []
    for v in snvs:
        k, alt = int(v[:-3]), v[-1]
        seq = cl.edit(ref, k, k, alt)
        made = []
        for m in motifs:
            pv, s0, s1, i, cov = best(m, seq, k)
            pw = m['wt'][cov].min()
            if pv < P_HIT and pw >= FOLD * pv:
                made.append((pv, m['tf']))
            if v == a.snv and pv < P_HIT and pv < pw:
                target.append(dict(db=m['db'], motif=m['id'], tf=m['tf'], tf_class=m['class'], tf_family=m['family'],
                                   p_wt=pw, p_variant=pv, fold=pw / pv, oligo_lo=s0, oligo_hi=s1,
                                   c=f'c.{s0 - C0}..c.{s1 - C0}',
                                   overlap_protected=int(max(0, min(s1, a.protected[1]) - max(s0, a.protected[0]) + 1))))
        made.sort()
        rows.append(dict(snv=v, created=len(made), best_created_p=made[0][0] if made else np.nan,
                         best_created_tf=made[0][1] if made else ''))
    al = pd.DataFrame(rows)
    al.to_csv(out / 'all_snvs.tsv', sep='\t', index=False, float_format='%.3g')
    tg = pd.DataFrame(target).sort_values('p_variant')
    tg.to_csv(out / 'target.tsv', sep='\t', index=False, float_format='%.3g')
    me = al.set_index('snv').loc[a.snv]
    json.dump(dict(snv=a.snv, created=int(me.created), best_created_p=float(me.best_created_p),
                   element_created_median=float(al.created.median()),
                   element_created_pct_rank=float((al.created <= me.created).mean()),
                   best_p_pct_rank=float((al.best_created_p.fillna(1) >= me.best_created_p).mean()),
                   snvs_with_any_created=int((al.created > 0).sum()), snvs=len(al)),
              open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

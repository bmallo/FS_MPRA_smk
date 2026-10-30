"""Step A10 sensitivity (post hoc): empirical null (Efron 2004) for the HepG2 per-SNV tests.

Far from each SNV (> 50 bp) the z values of real SNVs spread wider than on A9's pseudo-variants. Here the null
SD sigma0 per endpoint is the robust SD (1.4826 * MAD) of those far z values, every z is divided by sigma0, and
the Simes / Benjamini-Bogomolov procedure of a10_effects.py is rerun. Writes snvs.tsv, hits.tsv and summary.json.

Usage: python tools/analysis_v2/a10_sensitivity.py --run <a10 HepG2 run dir> --output <new dir>
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

_spec = importlib.util.spec_from_file_location('a10', Path(__file__).with_name('a10_effects.py'))
a10 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a10)
FAR = 50


def robust_sd(x):
    x = np.asarray(x, float)
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    e = pd.read_parquet(Path(a.run) / 'effects.parquet')
    e['dist'] = (e.oligo - e.snv_oligo).abs()
    sig = {ep: robust_sd(g[g.testable & (g.dist > FAR)].z) for ep, g in e.groupby('endpoint')}
    e['z0'] = e.z / e.endpoint.map(sig)
    e['p0'] = np.where(e.testable, 2 * stats.norm.sf(np.abs(e.z0)), np.nan)
    e['hit0'] = False
    rows, summ = [], {}
    for (ep, grp), g in e.groupby(['endpoint', 'group']):
        fam, feat, idx = {}, {}, {}
        for v, h in g.groupby('snv'):
            fam[v], feat[v], idx[v] = a10.simes(h.p0.to_numpy()), h.p0.to_numpy(), h.index
        chosen, rej = a10.bb(fam, feat)
        for v in fam:
            e.loc[idx[v][rej[v]], 'hit0'] = True
            rows.append(dict(snv=v, endpoint=ep, group=grp, simes_p=fam[v], selected=v in chosen, hits=int(rej[v].sum())))
        h = e.loc[g.index][lambda x: x.hit0]
        summ[f'{ep} {grp}'] = dict(snvs=len(fam), selected=len(chosen), hit_features=int(len(h)),
                                   hits_within_10bp=int((h.dist <= 10).sum()), hits_beyond_50bp=int((h.dist > FAR).sum()))
    pd.DataFrame(rows).to_csv(out / 'snvs.tsv', sep='\t', index=False)
    e[e.hit0].drop(columns=['hit']).to_csv(out / 'hits.tsv', sep='\t', index=False)
    info = dict(run=a.run, far_bp=FAR, sigma0=sig, results=summ)
    json.dump(info, open(out / 'summary.json', 'w'), indent=1, default=float)
    print(json.dumps(info, indent=1, default=float))


if __name__ == '__main__':
    main()

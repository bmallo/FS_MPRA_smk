"""Figure 1 (WASF1, EF1a): per-base FIRE share on hg38 from fig1_fire.sh output (fig1_rules.json: wasf1.draft1, ef1a).

Share = kept molecules (fig1_profiles parquet read IDs) with a FIRE element over the base / kept molecules covering it
(ft extract st-en), as the LRSS 07_profiles.py FIRE column. Plasmid reads are lifted to hg38 with inserts.tsv and limited
to the construct insert. Writes fire_profiles.tsv.gz (track, pos, coverage_fire, FIRE_share) in the LRSS layout.
Usage: python tools/r01_figure/fig1_fire_profiles.py --work results/r01_figure_20261102/fig1_wasf1 [--gene EF1a]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

GENES = {   # window (hg38) and tracks {track: (dataset, construct contig or None for endogenous)}
    'WASF1': ((110176895, 110181442), {'genome': ('endogenous', None), '300bp': ('WASF12', 'WASF12')}),
    'EF1a': ((73520048, 73521236), {'genome': ('endogenous', None), 'plasmid': ('pEFGFP', '#1')}),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--work', required=True)
    ap.add_argument('--gene', default='WASF1', choices=sorted(GENES))
    a = ap.parse_args()
    (X0, X1), TRACKS = GENES[a.gene]
    w = Path(a.work)
    ins = pd.read_csv(w / 'inserts.tsv', sep='\t').set_index('construct')
    L = X1 - X0 + 1
    out = []
    for track, (ds, contig) in TRACKS.items():
        kept = set(pd.read_parquet(w / f'profiles/{ds}.molecules.parquet', columns=['read']).read)
        x = pd.read_csv(w / f'fire/{ds}.extract.tsv.gz', sep='\t',
                        usecols=['#ct', 'st', 'en', 'fiber', 'ref_fire_starts', 'ref_fire_lengths'])
        n_all = len(x)
        x = x[x.fiber.isin(kept)].drop_duplicates('fiber')
        if contig:
            r = ins.loc[contig]
            off, lo, hi = int(r.hg38_start - r.construct_start), int(r.hg38_start), int(r.hg38_end)
            x = x[x['#ct'] == contig]
        else:
            off, lo, hi = 0, X0, X1
        cov, cnt = np.zeros(L), np.zeros(L)
        n_fire = 0
        for row in x.itertuples():
            a0, a1 = max(row.st + off + 1, lo, X0), min(row.en + off, hi, X1)   # ft extract: 0-based half-open
            if a1 < a0:
                continue
            cov[a0 - X0:a1 - X0 + 1] += 1
            st = [int(v) for v in str(row.ref_fire_starts).split(',') if v not in ('', '.', 'nan')]
            ln = [int(v) for v in str(row.ref_fire_lengths).split(',') if v not in ('', '.', 'nan')]
            hit = np.zeros(L, bool)
            for s, l_ in zip(st, ln):
                if s < 0 or l_ <= 0:
                    continue
                s_, e_ = max(s + off + 1, a0), min(s + off + l_, a1)
                if e_ >= s_:
                    hit[s_ - X0:e_ - X0 + 1] = True
            n_fire += hit.any()
            cnt += hit
        print(f'{track}: extract rows {n_all}, kept molecules {len(kept)}, matched {len(x)}, '
              f'with a FIRE element in range {n_fire}', flush=True)
        out.append(pd.DataFrame({'track': track, 'pos': np.arange(X0, X1 + 1), 'coverage_fire': cov,
                                 'FIRE_share': np.where(cov > 0, cnt / np.maximum(cov, 1), np.nan)}))
    pd.concat(out).to_csv(w / 'fire/fire_profiles.tsv.gz', sep='\t', index=False, float_format='%.5g')


if __name__ == '__main__':
    main()

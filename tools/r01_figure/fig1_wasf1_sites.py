"""Figure 1 (WASF1): how well supported are the footprint peaks in the peak-4 300-bp element? (exploratory;
fig1_rules.json: wasf1.site_check)

Sites: the GC box (reference), the 31-50 bp peak on the left of the element and the 20-30 bp peak on the right (centres =
maxima of the smoothed call-midpoint density of that bin in WASF12 molecules). Per site:
- motif matches overlapping it (P < 1e-5 from motifs/, P < 1e-4 from motifs_peak4_1e4/);
- sequence, A/T count and the longest run without A/T (Hia5 can only mark A/T, so an A/T-poor stretch can look protected);
- per dataset (WASF12 replicate 1 and 2, 4.5-kb construct, endogenous): molecules spanning the site, share with a call of
  the site's bin whose midpoint is within +/- WIN bp of the centre, and the element-wide median share for that bin
  (same rule at every centre in the element) as a background;
- per call at the site (WASF12): observed A/T inside the call and how many of them carry m6A (evidence the HMM had).
Usage: python tools/r01_figure/fig1_wasf1_sites.py --work results/r01_figure_20261102/fig1_wasf1
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy.ndimage import gaussian_filter1d

HG38 = '/mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.fa'
PEAK = (110179911, 110180211)
WIN, SPAN_SLACK = 10, 20
BIN = {'20-30': (20, 30), '31-50': (31, 50), '20-80': (20, 80)}


def density_peak(mol, lo_ln, hi_ln, x0, x1):
    mid = np.zeros(PEAK[1] - PEAK[0] + 1)
    for r in mol.itertuples():
        for s, e, ln in json.loads(r.fp):
            if lo_ln <= ln <= hi_ln and PEAK[0] <= (s + e) / 2 <= PEAK[1]:
                mid[int((s + e) / 2) - PEAK[0]] += 1
    sm = gaussian_filter1d(mid, 3)
    i0, i1 = x0 - PEAK[0], x1 - PEAK[0]
    return int(np.argmax(sm[i0:i1 + 1])) + x0


def shares(mol, centres, lo_ln, hi_ln):
    """Per centre: molecules spanning centre +/- (WIN + SPAN_SLACK) and those with a call of the bin whose midpoint is
    within +/- WIN of it."""
    n, k = np.zeros(len(centres)), np.zeros(len(centres))
    c = np.asarray(centres)
    for r in mol.itertuples():
        at = json.loads(r.at)
        if not at:
            continue
        span = (c - WIN - SPAN_SLACK >= min(at)) & (c + WIN + SPAN_SLACK <= max(at))
        mids = np.array([(s + e) / 2 for s, e, ln in json.loads(r.fp) if lo_ln <= ln <= hi_ln])
        hit = (np.abs(mids[:, None] - c[None, :]) <= WIN).any(0) if mids.size else np.zeros(len(c), bool)
        n += span
        k += span & hit
    return n, k


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--work', required=True)
    a = ap.parse_args()
    w = Path(a.work)
    out = w / 'site_check'
    out.mkdir(exist_ok=True)
    mol = {d: pd.read_parquet(w / f'profiles/{d}.molecules.parquet') for d in ('WASF12', 'WASF1_full', 'endogenous')}
    for rep in (1, 2):
        names = set((w / f'bams/WASF12_rep{rep}.names').read_text().split())
        mol[f'WASF12_rep{rep}'] = mol['WASF12'][mol['WASF12'].read.isin(names)]
    sites = {'GC box': ((110180054 + 110180097) // 2, '20-80'),
             '31-50 bp peak (left)': (density_peak(mol['WASF12'], 31, 50, PEAK[0] + 20, 110180040), '31-50'),
             '20-30 bp peak (right)': (density_peak(mol['WASF12'], 20, 30, 110180110, PEAK[1] - 20), '20-30')}
    fa = pysam.FastaFile(HG38)
    m5 = pd.read_csv(w / 'motifs/motifs.tsv', sep='\t')
    m4 = pd.read_csv(w / 'motifs_peak4_1e4/motifs.tsv', sep='\t')
    rows, mot_rows, ev_rows = [], [], []
    for name, (c, b) in sites.items():
        lo, hi = c - 12, c + 12                                  # ~25 bp around the centre for sequence and motifs
        seq = fa.fetch('chr6', lo - 1, hi).upper()
        at = [ch in 'AT' for ch in seq]
        runs = [len(x) for x in ''.join('x' if v else '.' for v in at).split('x')]
        for tag, m in (('1e-5', m5), ('1e-4', m4)):
            ov = m[(m.end >= c - WIN) & (m.start <= c + WIN)].sort_values('p')
            best = ov.groupby('tf', as_index=False).first().sort_values('p')
            mot_rows += [dict(site=name, centre=c, threshold=tag, tf=r.tf, db=r.db, start=r.start, end=r.end,
                              strand=r.strand, p=r.p) for r in best.itertuples()]
        lo_ln, hi_ln = BIN[b]
        for d in ('WASF12_rep1', 'WASF12_rep2', 'WASF1_full', 'endogenous'):
            n, k = shares(mol[d], [c], lo_ln, hi_ln)
            cs = np.arange(PEAK[0] + WIN + SPAN_SLACK, PEAK[1] - WIN - SPAN_SLACK + 1, 5)
            nb, kb = shares(mol[d], cs, lo_ln, hi_ln)
            bg = np.median(kb[nb > 0] / nb[nb > 0]) if (nb > 0).any() else np.nan
            rows.append(dict(site=name, centre=c, bin=b, dataset=d, molecules=int(n[0]), with_call=int(k[0]),
                             share=round(k[0] / n[0], 3) if n[0] else np.nan, element_median_share=round(bg, 3),
                             seq=seq, at_count=sum(at), longest_no_at=max(runs)))
        for r in mol['WASF12'].itertuples():                    # evidence inside each call at the site
            atp, m6 = set(json.loads(r.at)), set(json.loads(r.m6a))
            for s, e, ln in json.loads(r.fp):
                if lo_ln <= ln <= hi_ln and abs((s + e) / 2 - c) <= WIN:
                    inside = [p_ for p_ in range(s, e + 1) if p_ in atp]
                    ev_rows.append(dict(site=name, length=ln, at_observed=len(inside),
                                        m6a_inside=sum(p_ in m6 for p_ in inside)))
    R, M, E = pd.DataFrame(rows), pd.DataFrame(mot_rows), pd.DataFrame(ev_rows)
    R.to_csv(out / 'site_shares.tsv', sep='\t', index=False)
    M.to_csv(out / 'site_motifs.tsv', sep='\t', index=False)
    E.to_csv(out / 'site_call_evidence.tsv', sep='\t', index=False)
    print(R.drop(columns='seq').to_string(index=False))
    for name, (c, b) in sites.items():
        print(f'\n{name} centre {c}: sequence {c - 12}-{c + 12}', R[R.site == name].seq.iloc[0])
        mm = M[M.site == name]
        for tag in ('1e-5', '1e-4'):
            x = mm[mm.threshold == tag]
            print(f'  motifs P < {tag}: {len(x)} TFs;', ', '.join(f'{r.tf} {r.start}-{r.end} ({r.p:.1e})'
                                                             for r in x.head(8).itertuples()))
        e = E[E.site == name]
        print(f'  WASF12 calls at site: {len(e)}; observed A/T inside a call: median {e.at_observed.median():.0f} '
              f'(IQR {e.at_observed.quantile(.25):.0f}-{e.at_observed.quantile(.75):.0f}); calls with <= 2 A/T: '
              f'{(e.at_observed <= 2).mean():.0%}; m6A inside: {(e.m6a_inside > 0).mean():.0%} of calls')


if __name__ == '__main__':
    main()

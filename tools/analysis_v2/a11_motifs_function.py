"""Step A11: motif dose-response and the link to published MPRA activity (HepG2; rules.json motifs_function).

Dose-response: per motif, the WT best PWM placement near the window; each SNV in the window or site gets
delta PWM (variant - WT log2-odds at that placement) and two chromatin changes from the A10 +/-5 bp mask run:
  primary   mean change in 10-40 bp call coverage over the site positions (10-19 + 20-40 bin D)
  secondary mean change in unmarked fraction over testable adenines within the site +/-10 bp
  model     endpoint ~ delta PWM + drops-its-base, OLS with HC3 errors; BH across motifs per set and endpoint.
Function: Kircher 2019 log2 effects (LDLR primary, LDLR.2 replication) against
  C1 mean E1 change > 50 bp from the SNV (element-wide shift), C2 mean E1 change 11-50 bp away (local);
  F1 activity ~ C1 + position, F2 activity ~ C2 + position, F3 |activity| ~ RMS(C2) + position,
  OLS with errors clustered by position; also without the position term (description).

Writes dose_snvs.tsv, dose.tsv, function_snvs.tsv, function.tsv, summary.json.

Usage: python tools/analysis_v2/a11_motifs_function.py --a10 <a10 HepG2_run3_pm5 dir> --catalog <catalog.json> \
           --jaspar extras/data/JASPAR2024_CORE_vertebrates_nr.jaspar --ldlr <GRCh38_LDLR.tsv> \
           --ldlr2 <GRCh38_LDLR.2.csv> --output <new dir>
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location('tfa', ROOT / 'extras' / 'fs_mpra_tf_annotate.py')
tfa = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tfa)

OLIGO0 = 21
EXTEND = 5
E1_FLANK = 10
FAR, LOCAL = 50, (11, 50)
PSEUDOCOUNT = 0.8
MOTIFS = [('Repeat 1 (Sp1)', 89, 104, 'MA0079.5', 'primary'), ('Repeat 2 (SRE-1)', 124, 139, 'MA0596.1', 'primary'),
          ('Repeat 3 (Sp1)', 140, 155, 'MA0079.5', 'primary'), ('ZNF652-like', 52, 62, 'MA1657.2', 'secondary'),
          ('SP4', 82, 102, 'MA0685.2', 'secondary'), ('CREB1/CREM', 182, 192, 'MA0018.5', 'secondary')]
REGION = (143, 151)
REGION_MIN_FAR_HITS, REGION_FAR = 20, 30
COMP = str.maketrans('ACGT', 'TGCA')


def score(lod, seq):
    idx = tfa._seq_idx(seq)
    return float(lod[np.arange(lod.shape[0]), idx].sum())


def best_site(lod, wt, lo, hi):
    """WT best placement (score, start oligo, strand) with its start in [lo - EXTEND, hi + EXTEND - L + 1]."""
    L = lod.shape[0]
    best = (-1e18, None, None)
    for start in range(max(OLIGO0, lo - EXTEND), min(OLIGO0 + len(wt) - L, hi + EXTEND - L + 1) + 1):
        w = wt[start - OLIGO0:start - OLIGO0 + L]
        for strand, s in (('+', w), ('-', w.translate(COMP)[::-1])):
            sc = score(lod, s)
            if sc > best[0]:
                best = (sc, start, strand)
    return best


def site_score(lod, seq, start, strand):
    L = lod.shape[0]
    w = seq[start - OLIGO0:start - OLIGO0 + L]
    return score(lod, w if strand == '+' else w.translate(COMP)[::-1])


def ols(formula, data, cluster=None):
    m = smf.ols(formula, data=data)
    return m.fit(cov_type='cluster', cov_kwds=dict(groups=data[cluster])) if cluster else m.fit(cov_type='HC3')


def bh(p):
    p = np.asarray(p, float)
    out = np.full(p.size, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        q = p[ok]
        o = np.argsort(q)
        adj = np.minimum.accumulate((q[o] * q.size / np.arange(1, q.size + 1))[::-1])[::-1]
        r = np.empty(q.size)
        r[o] = np.minimum(adj, 1)
        out[ok] = r
    return out


def activity(path, sep, catalog_snvs):
    t = pd.read_csv(path, sep=sep)
    t = t.rename(columns={'Position': 'genomic_pos1', 'Ref': 'ref', 'Alt': 'alt', 'Value': 'act'})
    m = catalog_snvs.merge(t[['genomic_pos1', 'ref', 'alt', 'act']], on=['genomic_pos1', 'ref', 'alt'], how='left')
    return m.set_index('source_variant').act


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--a10', required=True)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--jaspar', required=True)
    ap.add_argument('--ldlr', required=True)
    ap.add_argument('--ldlr2', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    lib = json.load(open(a.catalog))['library']
    wt = lib['mutable_sequence']
    cat = pd.DataFrame(lib['variants'])
    snvs = cat[cat['class'] == 'SNV'][['source_variant', 'oligo_pos1', 'ref', 'alt', 'genomic_pos1']].copy()
    snvs['oligo_pos1'] = snvs.oligo_pos1.astype(int)
    pwms = {p['id']: p for p in tfa.load_pwms(a.jaspar, None, PSEUDOCOUNT, [0.25] * 4)}
    e = pd.read_parquet(Path(a.a10) / 'effects.parquet')
    e['dist'] = (e.oligo - e.snv_oligo).abs()
    group = e.drop_duplicates('snv').set_index('snv').group

    # dose-response
    drows, mrows = [], []
    for name, lo, hi, mid, mset in MOTIFS:
        lod = pwms[mid]['lod']
        L = lod.shape[0]
        wt_sc, start, strand = best_site(lod, wt, lo, hi)
        s_lo, s_hi = start, start + L - 1
        inside = snvs[((snvs.oligo_pos1 >= lo) & (snvs.oligo_pos1 <= hi)) |
                      ((snvs.oligo_pos1 >= s_lo) & (snvs.oligo_pos1 <= s_hi))]
        for v in inside.itertuples():
            k = v.oligo_pos1 - OLIGO0
            var = wt[:k] + v.alt + wt[k + 1:]
            d_pwm = site_score(lod, var, start, strand) - wt_sc if s_lo <= v.oligo_pos1 <= s_hi else 0.0
            g = e[e.snv == v.source_variant]
            e2 = g[(g.endpoint == 'E2') & g['bin'].isin(['sub_TF', 'TF']) & g.oligo.between(s_lo, s_hi) & g.testable]
            per = e2.groupby('oligo').agg(D=('D', 'sum'), k=('D', 'size'))
            per = per[per.k == 2]
            e1 = g[(g.endpoint == 'E1') & g.testable & g.oligo.between(s_lo - E1_FLANK, s_hi + E1_FLANK)]
            drows.append(dict(motif=name, set=mset, snv=v.source_variant, oligo=v.oligo_pos1,
                              in_site=bool(s_lo <= v.oligo_pos1 <= s_hi), d_pwm=d_pwm,
                              drops=int(group.get(v.source_variant, '') == 'drops its base'),
                              calls_D=per.D.mean() if len(per) else np.nan, calls_positions=len(per),
                              unmarked_D=e1.D.mean() if len(e1) else np.nan, unmarked_adenines=len(e1)))
        d = pd.DataFrame([r for r in drows if r['motif'] == name])
        for ep in ('calls_D', 'unmarked_D'):
            x = d.dropna(subset=[ep])
            row = dict(motif=name, set=mset, matrix=mid, tf=pwms[mid]['tf'], site=f'{s_lo}-{s_hi}', strand=strand,
                       wt_score=round(wt_sc, 2), endpoint=ep, snvs=len(x), snvs_in_site=int(x.in_site.sum()))
            if len(x) >= 4 and x.d_pwm.std() > 0:
                f = ols(f'{ep} ~ d_pwm' + (' + drops' if x.drops.nunique() > 1 else ''), x)
                ci = f.conf_int().loc['d_pwm']
                row.update(slope=f.params['d_pwm'], se=f.bse['d_pwm'], ci_lo=ci[0], ci_hi=ci[1], p=f.pvalues['d_pwm'],
                           r_pwm=float(np.corrcoef(x.d_pwm, x[ep])[0, 1]))
            mrows.append(row)
    dose_snvs = pd.DataFrame(drows)
    dose = pd.DataFrame(mrows)
    dose['q'] = np.nan
    for (_, _), g in dose.groupby(['set', 'endpoint']):
        dose.loc[g.index, 'q'] = bh(g.get('p', pd.Series(np.nan, index=g.index)).to_numpy())
    dose_snvs.to_csv(out / 'dose_snvs.tsv', sep='\t', index=False)
    dose.to_csv(out / 'dose.tsv', sep='\t', index=False)

    # function
    e1 = e[(e.endpoint == 'E1') & e.testable]
    c1 = e1[e1.dist > FAR].groupby('snv').D.mean()
    loc = e1[e1.dist.between(*LOCAL)]
    c2 = loc.groupby('snv').D.mean()
    rms = loc.groupby('snv').D.apply(lambda x: float(np.sqrt(np.mean(x ** 2))))
    far_hits = e[(e.endpoint == 'E1') & e.hit & (e.dist > REGION_FAR)].groupby('snv').size()
    f = snvs.set_index('source_variant')[['oligo_pos1']].rename(columns={'oligo_pos1': 'pos'})
    f['act'] = activity(a.ldlr, '\t', snvs)
    f['act2'] = activity(a.ldlr2, ',', snvs)
    f['C1'], f['C2'], f['rms'] = c1, c2, rms
    f['far_hits'] = far_hits.reindex(f.index).fillna(0).astype(int)
    f['region'] = f.pos.between(*REGION) & (f.far_hits >= REGION_MIN_FAR_HITS)
    for col in ('C1', 'C2', 'rms'):
        f[col] = 100 * f[col]                                 # percentage points
    f.index.name = 'snv'
    f.reset_index().to_csv(out / 'function_snvs.tsv', sep='\t', index=False)
    frows = []
    for src in ('act', 'act2'):
        for model, y, x in (('F1', src, 'C1'), ('F2', src, 'C2'), ('F3', f'abs_{src}', 'rms')):
            d = f.dropna(subset=[src, x]).copy()
            d[f'abs_{src}'] = d[src].abs()
            for fe in (True, False):
                fit = ols(f'{y} ~ {x}' + (' + C(pos)' if fe else ''), d, cluster='pos')
                ci = fit.conf_int().loc[x]
                frows.append(dict(source='LDLR' if src == 'act' else 'LDLR.2', model=model, outcome=y.replace(src, 'activity'),
                                  chromatin=x, position_fe=fe, n=len(d), positions=d.pos.nunique(), slope=fit.params[x],
                                  se=fit.bse[x], ci_lo=ci[0], ci_hi=ci[1], p=fit.pvalues[x],
                                  spearman=float(d[[y, x]].corr('spearman').iloc[0, 1])))
    func = pd.DataFrame(frows)
    func.to_csv(out / 'function.tsv', sep='\t', index=False)
    reg = f[f.region]
    summ = dict(settings=vars(a), motifs=dose.to_dict('records'), function=func.to_dict('records'),
                region=dict(snvs=reg.index.tolist(), act_median=float(reg.act.median()),
                            act_median_others=float(f[~f.region].act.median()),
                            act2_median=float(reg.act2.median()), act2_median_others=float(f[~f.region].act2.median())))
    json.dump(summ, open(out / 'summary.json', 'w'), indent=1, default=float)
    pd.set_option('display.width', 250)
    print(dose.round(4).to_string())
    print(func.round(4).to_string())
    print(json.dumps(summ['region'], default=float))


if __name__ == '__main__':
    main()

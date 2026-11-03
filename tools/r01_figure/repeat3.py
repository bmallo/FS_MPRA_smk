"""R01 figure: Repeat 3 Sp1 site (oligo 143-151), clinical classes, chromatin and activity (settings in
tools/r01_figure/repeat3_rules.json, written before this was run).

HepG2 run 3
  1. Per SNV in the groups: occupancy of two call boxes (site: confident 10-40 bp calls with midpoint 143-151;
     nucleosome-size: confident >= 81 bp calls with midpoint in the element), variant vs labeling-matched WT with
     the A9 random-effects test. Masked adenines (A7 primary; +/-5 bp sensitivity) are set unobserved in the SNV's
     molecules and in WT before calling, exactly as fp_features does; unmasked SNVs use the stored calls.
  2. Per group: DerSimonian-Laird pool across SNVs of each endpoint (SNVs with n >= 10, C >= 3), and of the stored
     A10 E1 D per adenine (A7 mask, testable rows).
HEK293T (display only)
  3. Per clone (>= 20 molecules with labeling > 0.01): mean labeling-standardised inaccessible-fraction residual vs
     HEK WT at adenines 140-156 and at adenines > 50 bp from 147 (edited A/T unobserved). WT clones for reference.

Writes snv_endpoints.tsv, groups.tsv, profile.tsv, hek_clones.tsv, hek_check.json, summary.json to --output.
Usage: python tools/r01_figure/repeat3.py --root results/analysis_v2_20261028 --clinical <step3 dir> --procs 4 --output <dir>
"""
import argparse
import importlib.util
import json
import multiprocessing as mp
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
_spec = importlib.util.spec_from_file_location('fp_features', HERE.parent / 'analysis_v2/fp_features.py')
fpf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fpf)
a10, a9f, a7, a6x = fpf.a10, fpf.a9f, fpf.a7, fpf.a6x
from tools.analysis_v2 import access  # noqa: E402

RULES = json.load(open(HERE / 'repeat3_rules.json'))
SITE = tuple(RULES['site']['oligo'])
CATALOG = 'results/rep3_ledger_20261016/catalog/catalog.json'
VARIANTS = 'results/analysis_v2_20261028/a4_informative/visibility/variants_visibility.tsv'
FEATS = [dict(name='site', mid_lo=SITE[0] - fpf.OLIGO0, mid_hi=SITE[1] - fpf.OLIGO0, len_lo=10, len_hi=40),
         dict(name='nucleosome', mid_lo=0, mid_hi=fpf.WIDTH - 1, len_lo=81, len_hi=10_000)]
MIN_N, MIN_C, HEK_MIN_MOL = 10, 3, 20


def dersimonian_laird(est, se):
    est, se = np.asarray(est, float), np.asarray(se, float)
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    est, se = est[ok], se[ok]
    k = est.size
    if k == 0:
        return dict(k=0, est=np.nan, se=np.nan, lo=np.nan, hi=np.nan, tau2=np.nan, I2=np.nan, p=np.nan)
    w = 1 / se ** 2
    fixed = (w * est).sum() / w.sum()
    Q = (w * (est - fixed) ** 2).sum()
    tau2 = max(0.0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0.0
    ws = 1 / (se ** 2 + tau2)
    m = (ws * est).sum() / ws.sum()
    s = np.sqrt(1 / ws.sum())
    from scipy import stats
    return dict(k=k, est=m, se=s, lo=m - 1.96 * s, hi=m + 1.96 * s, tau2=tau2,
                I2=max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0, p=2 * stats.norm.sf(abs(m / s)))


def groups(root, clinical):
    fn = pd.read_csv(root / 'a11_motifs_function/function_snvs.tsv', sep='\t').set_index('snv')
    d = pd.read_csv(root / 'a11_motifs_function/dose_snvs.tsv', sep='\t').query('motif == "Repeat 3 (Sp1)"')
    d = d.set_index('snv')
    g = RULES['groups']
    out = {}
    g1 = g['G1_clinvar_plp_conflicting_in_site']['snvs']
    site = [v for v in fn.index if SITE[0] <= fn.pos[v] <= SITE[1]]
    out['G1'] = g1
    out['G2'] = [v for v in site if v not in g1 and d.d_pwm.get(v, 0) <= -2]
    out['G3'] = [v for v in site if d.d_pwm.get(v, 0) > -2]
    assert sorted(out['G3']) == sorted(g['G3_site_preserving']['snvs']), out['G3']
    out['G4'] = [v for v in fn.index if 138 <= fn.pos[v] <= 142 or 152 <= fn.pos[v] <= 157]
    out['G5'] = g['G5_clinvar_benign_lb']['snvs']
    out['G6'] = g['G6_clinvar_plp_outside_R3']['snvs']
    # check the clinical lists against the step-3 table
    c = pd.read_csv(Path(clinical) / 'clinical_variants.tsv', sep='\t')
    c = c[c.type == 'SNV'].set_index('designed').clinvar_class.fillna('')
    plp = {v for v, k in c.items() if k in ('Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic')}
    conf = {v for v, k in c.items() if k.startswith('Conflicting')}
    ben = {v for v, k in c.items() if k in ('Benign', 'Likely benign')}
    in_site = lambda v: SITE[0] <= int(v[:-3]) <= SITE[1]
    assert set(out['G1']) == {v for v in plp | conf if in_site(v)}
    assert set(out['G5']) == ben and set(out['G6']) == {v for v in plp if not in_site(v)}
    return out, fn, c


def snv_task(task):
    v, mask = task
    g = fpf._G
    rv = g['snv_rows'][v]
    cols = g['masks'][mask][v]['mask']
    if cols.size:
        tw = fpf.call_table(fpf._masked(g['X'][g['wt2']], cols), g['hit'], g['miss'], g['lp'], g['pos'], *fpf.EL,
                            fpf.KMAX, a7.MIN_LLR, a7.MIN_OPPS, a7.DRAWDOWN)
        yw = fpf.occupancy(tw, g['distal'][g['wt2']], FEATS)
        tv = fpf.call_table(fpf._masked(g['X'][rv], cols), g['hit'], g['miss'], g['lp'], g['pos'], *fpf.EL,
                            fpf.KMAX, a7.MIN_LLR, a7.MIN_OPPS, a7.DRAWDOWN)
        yv = fpf.occupancy(tv, g['distal'][rv], FEATS)
    else:
        yw, yv = g['occ_wt'], g['occ_full'][rv]
    out = []
    for kind, bit in (('all', 1), ('confident', 2)):
        Yw = ((yw & bit) > 0).astype(float)
        Sy, So = a10.stratum_sums(Yw, g['st'][g['wt2']])
        ref = g['ref'][kind]
        p = np.divide(Sy, So, out=np.full_like(Sy, np.nan), where=So > 0)
        Y = (yv & bit) > 0
        r = a9f.contrast(Y, np.ones_like(Y), g['st'][rv], g['clone'][rv], p, So,
                         SimpleNamespace(V=ref.V, W2=ref.W2, tau2=ref.tau2))
        for f, b in enumerate(FEATS):
            out.append(dict(snv=v, mask=mask, calls=kind, endpoint=b['name'], masked_adenines=int(cols.size),
                            occ_variant=float(Y[:, f].mean()), wt_occ=float(np.nanmean(Yw[:, f])),
                            n=r['n'][f], C=r['C'][f], D=r['D'][f], se=r['se'][f], p=r['p'][f],
                            expected=r['expected'][f]))
    return out


def hepg2(root, gl, procs, t0):
    ext = root / 'a9_statistics/extract/HepG2_run3'
    pos, X, distal, _, _ = a6x.load(ext)
    meta = np.load(ext / 'meta.npz')
    label, cls, clone = meta['label'], meta['cls'], meta['clone']
    hit, miss = a7.emission_steps(fpf.EMISSIONS, pos)
    t = pd.read_csv(fpf.EMISSIONS, sep='\t')
    t = t[t.fold == 'A']
    pa = dict(zip(t.ref_pos - 1, t.p_a))
    lp = np.log(1 - np.clip([pa.get(int(p), 0.0) for p in pos], 0, 0.9999))
    k2 = np.flatnonzero(distal >= a9f.E2_MIN)
    z = np.load(root / 'fp_features/calls.npz')
    assert np.array_equal(z['rows'], k2)
    tab = (z['mid'], z['length'], z['start'], z['end'], z['S'])
    occ_full = np.zeros((distal.size, len(FEATS)), np.uint8)
    occ_full[k2] = fpf.occupancy(tab, distal[k2], FEATS)
    wt2 = k2[cls[k2] == 'WT']
    st = a9f.strata(distal)
    ow = occ_full[wt2]
    refs = {kind: a9f.Reference((ow & bit) > 0, np.ones(ow.shape, bool), st[wt2], clone[wt2])
            for kind, bit in (('all', 1), ('confident', 2))}
    masks = {m: a10.snv_masks(CATALOG, VARIANTS, pos, *fpf.EL, False, pm5=(m == 'pm5')) for m in ('A7', 'pm5')}
    snvs = sorted({v for g in gl.values() for v in g})
    snv2 = k2[cls[k2] == 'SNV']
    fpf._G.update(X=X, hit=hit, miss=miss, lp=lp, pos=pos, distal=distal, st=st, clone=clone, wt2=wt2,
                  occ_wt=ow, occ_full=occ_full, ref=refs, masks=masks,
                  snv_rows={v: snv2[label[snv2] == v] for v in snvs})
    print(f'loaded, {len(snvs)} SNVs, WT {wt2.size} molecules, {time.time() - t0:.0f} s', flush=True)
    tasks = [(v, m) for v in snvs for m in ('A7', 'pm5')]
    rows = []
    with mp.get_context('fork').Pool(procs) as pool:
        for i, r in enumerate(pool.imap_unordered(snv_task, tasks)):
            rows += r
            if i % 20 == 0:
                print(f'{i + 1}/{len(tasks)}, {time.time() - t0:.0f} s', flush=True)
    wt_occ = {f'{kind}_{b["name"]}': float((((ow & bit) > 0)[:, f]).mean())
              for kind, bit in (('all', 1), ('confident', 2)) for f, b in enumerate(FEATS)}
    return pd.DataFrame(rows), wt_occ


def hek(root, gl):
    d = access.load_element('HEK293T')
    ref = access.E1Reference(d)
    base = json.load(open(root / 'explorer_v6/base.json'))['HEK293T']
    wt_ours = ref.wt_fraction()
    check = dict(adenines_match=bool(np.array_equal(np.asarray(base['adenines']), d.oligo)),
                 max_abs_diff_unmarked=float(np.nanmax(np.abs(np.asarray(base['unmarked']) / 1e4 - wt_ours))))
    use = d.distal > a9f.E1_MIN
    st = a9f.strata(d.distal)
    pm = ref.p[st]
    obs = (d.X >= 0) & use[:, None] & np.isfinite(pm)
    resid = np.where(obs, (d.X == 0) - np.nan_to_num(pm), np.nan)
    site_cols = (d.oligo >= 140) & (d.oligo <= 156)
    far_cols = np.abs(d.oligo - 147) > 50
    lab = pd.Series(d.label)
    grp = {v: k for k, vs in gl.items() for v in vs}
    rows = []
    for (lb, cl), idx in lab.groupby([lab, pd.Series(d.clone)]).groups.items():
        idx = np.asarray(idx)
        idx = idx[use[idx]]
        if idx.size < HEK_MIN_MOL or (lb != 'WT' and lb not in grp):
            continue
        R = resid[idx].copy()
        if lb != 'WT':
            k, ref_b = int(lb[:-3]), lb[-3]
            if ref_b in 'AT':
                R[:, d.oligo == k] = np.nan
        rows.append(dict(label=lb, group=grp.get(lb, 'WT'), clone=int(cl), molecules=int(idx.size),
                         site_resid_pts=100 * np.nanmean(R[:, site_cols]),
                         far_resid_pts=100 * np.nanmean(R[:, far_cols])))
    return pd.DataFrame(rows), check


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--clinical', required=True)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    t0 = time.time()
    gl, fn, clin = groups(root, a.clinical)
    json.dump(gl, open(out / 'groups_membership.json', 'w'), indent=1)

    se, wt_occ = hepg2(root, gl, a.procs, t0)
    se = se.join(fn[['pos', 'act', 'act2', 'C1']], on='snv')
    se['group'] = se.snv.map({v: k for k, vs in gl.items() for v in vs})
    se['clinvar'] = se.snv.map(clin).fillna('')
    se.to_csv(out / 'snv_endpoints.tsv', sep='\t', index=False, float_format='%.5g')

    grows = []
    for gname, vs in gl.items():
        for (mask, calls, ep), t in se[se.snv.isin(vs)].groupby(['mask', 'calls', 'endpoint']):
            t = t.drop_duplicates('snv')
            t = t[(t.n >= MIN_N) & (t.C >= MIN_C)]
            r = dersimonian_laird(100 * t.D, 100 * t.se)
            grows.append(dict(group=gname, measure=f'{ep}_{calls}_{mask}', **r))
        f = fn.reindex(vs)
        for col in ('act', 'act2', 'C1'):
            x = f[col].dropna()
            grows.append(dict(group=gname, measure=col, k=len(x), est=x.mean(), lo=x.min(), hi=x.max(),
                              median=x.median()))
    pd.DataFrame(grows).to_csv(out / 'groups.tsv', sep='\t', index=False, float_format='%.5g')

    snvs = sorted({v for g in gl.values() for v in g})
    e = pd.read_parquet(root / 'a10_effects/HepG2_run3/effects.parquet', columns=['snv', 'endpoint', 'oligo', 'D', 'se', 'p'],
                        filters=[('endpoint', '==', 'E1'), ('snv', 'in', snvs)])
    e = e[e.p.notna()]
    prow = []
    for gname, vs in gl.items():
        for o, t in e[e.snv.isin(vs)].groupby('oligo'):
            prow.append(dict(group=gname, oligo=o, **dersimonian_laird(100 * t.D, 100 * t.se)))
    pd.DataFrame(prow).to_csv(out / 'profile.tsv', sep='\t', index=False, float_format='%.5g')

    hc, check = hek(root, gl)
    hc.to_csv(out / 'hek_clones.tsv', sep='\t', index=False, float_format='%.5g')
    json.dump(check, open(out / 'hek_check.json', 'w'), indent=1)
    json.dump(dict(groups={k: len(v) for k, v in gl.items()}, wt_occupancy=wt_occ, hek_check=check,
                   seconds=round(time.time() - t0)), open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

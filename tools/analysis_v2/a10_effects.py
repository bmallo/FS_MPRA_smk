"""Step A10: variant effects per cell (settings: rules.json variant_effects; tests from A9, masks from A7).

HepG2 (--mode snv): every SNV against the WT reference with the A9 random-effects test, per feature.
  E1 (primary): element adenine unmarked, labeling > 0.01; the SNV's A7 primary mask is not tested.
  E2 (secondary): coverage by a call of each size bin, labeling >= 0.2; masked adenines set unobserved before
  calling, in the SNV's molecules and the WT reference alike (tau2 and clone weights from the unmasked reference).
  Global test per SNV = Simes over its testable features; Benjamini-Bogomolov hierarchical FDR (q = 0.05) per
  A7 caller group ('drops its base' apart from the rest) and endpoint. Secondary: max |z| against the A9
  pseudo-variants of the same clone tier.
HEK293T (--mode motif): per motif window and caller group, the pooled SNV molecules' mean motif residual against
  1,000 pseudo-pools (prespecified; invalid here, see rules.json HEK293T amendment_1) and the clone-mean test:
  mean of SNV clone means against WT clones of >= 100 scoring molecules (10,000 draws); BH across tests.

Writes effects.parquet + snvs.tsv + hits.tsv (snv mode) or motifs.tsv (motif mode), and summary.json.

Usage: python tools/analysis_v2/a10_effects.py --mode snv|motif --extract <a9 extract dir> \
           --emissions <fold-A emissions.tsv.gz> --element 1234:1445 [--reverse] --catalog <catalog.json> \
           --variants <A4 variants_visibility.tsv> [--a9 <a9 run dir>] --procs 4 --output <new dir>
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

_spec = importlib.util.spec_from_file_location('a9f', Path(__file__).with_name('a9_framework.py'))
a9f = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a9f)
a8r, a7, a6x = a9f.a8r, a9f.a7, a9f.a6x

Q = 0.05
WIDTH = 211
OLIGO0 = 21
BIN_NAMES = [b[0] for b in a7.BINS]
B_LO = np.array([b[1] for b in a7.BINS], np.int64)
B_HI = np.array([b[2] for b in a7.BINS], np.int64)
MOTIFS = [('ZNF652-like', 52, 62), ('SP4', 82, 102), ('SREBF1/2', 127, 136), ('SP1/SP3/KLF', 140, 156),
          ('CREB1/CREM', 182, 192)]
MOTIF_BINS = (0, 1, 2)                          # 10-80 bp calls for the HEK293T E2 motif score
N_NULL = 1000
N_CLONE_NULL = 10000
MIN_CLONE_MOLECULES = 100     # HEK293T clone-mean test (amendment 1)
SEED = 20261029


def caller_group(row):
    if str(row.caller_drops_at) == 'True':
        return 'drops its base'
    return 'gains an A/T' if row.at_count_change > 0 else 'neither'


def snv_masks(catalog, variants, pos, el_s, el_e, reverse, pm5=False):
    """Per SNV label: caller group, oligo position and the reference A/T columns of its A7 primary mask (with pm5:
    plus every A/T within +/-5 bp of the edited base, the prespecified sensitivity mask)."""
    lib = json.load(open(catalog))['library']
    wt = lib['mutable_sequence']
    seqs = {v['variant_id']: v['sequence'] for v in lib['variants']}
    var = pd.read_csv(variants, sep='\t')
    var = var[var.cls == 'SNV'].reset_index(drop=True)
    col_of_ref = {int(r): c for c, r in enumerate(pos)}
    el_ref = [(el_e - 1 - j) if reverse else (el_s + j) for j in range(WIDTH)]
    out = {}
    for v, m in zip(var.itertuples(), a7.mask_sets(var, wt, seqs)):
        idx = set(m['primary']) | (set(a7.pm_set(wt, int(v.index0))) if pm5 else set())
        cols = sorted(col_of_ref[el_ref[k]] for k in idx)
        out[v.source_variant] = dict(group=caller_group(v), oligo=int(v.oligo_pos1), mask=np.array(cols, np.int64))
    return out


def simes(p):
    p = np.sort(p[np.isfinite(p)])
    if p.size == 0:
        return np.nan
    return float(np.min(p.size * p / np.arange(1, p.size + 1)))


def bh(p, q):
    """Boolean rejections (NaN never rejected)."""
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    out = np.zeros(p.size, bool)
    pv = p[ok]
    m = pv.size
    if m == 0:
        return out
    o = np.argsort(pv)
    k = np.flatnonzero(pv[o] <= q * np.arange(1, m + 1) / m)
    if k.size:
        rej = np.zeros(m, bool)
        rej[o[:k.max() + 1]] = True
        out[np.flatnonzero(ok)] = rej
    return out


def bb(snv_p, feat_p, q=Q):
    """Benjamini-Bogomolov: BH on family p-values selects R of N; within each selected family BH at q R / N.
    snv_p: dict family -> p; feat_p: dict family -> array of feature p. Returns selected set, dict of rejections."""
    fams = list(snv_p)
    sel = bh(np.array([snv_p[f] for f in fams]), q)
    N, R = len(fams), int(sel.sum())
    chosen = {f for f, s in zip(fams, sel) if s}
    rej = {f: (bh(feat_p[f], q * R / N) if f in chosen else np.zeros(len(feat_p[f]), bool)) for f in fams}
    return chosen, rej


def masked_calls(X, cols, hit, miss, pos, el_s, el_e, reverse):
    Xm = X.copy()
    if cols.size:
        Xm[:, cols] = -1
    return a8r.call_bits(Xm, hit, miss, pos, el_s, el_e, reverse, B_LO, B_HI).reshape(X.shape[0], B_LO.size * WIDTH)


def _calls(args):
    return a9f._calls(args)


def stratum_sums(Y, st):
    Sy = np.zeros((a9f.N_STRATA, Y.shape[1]))
    So = np.zeros((a9f.N_STRATA, Y.shape[1]))
    for s in np.unique(st):
        m = st == s
        Sy[s] = Y[m].sum(0)
        So[s] = m.sum()
    return Sy, So


_G = {}


def _e2_masked(label):
    """E2 for one masked SNV: WT reference and the SNV's molecules called with the mask unobserved."""
    g = _G
    cols = g['masks'][label]['mask']
    bw = masked_calls(g['Xw2'], cols, g['hit'], g['miss'], g['pos'], g['el_s'], g['el_e'], g['reverse'])
    Sy, So = stratum_sums(bw, g['stw2'])
    p = np.divide(Sy, So, out=np.full_like(Sy, np.nan), where=So > 0)
    rows = g['snv_rows2'][label]
    bv = masked_calls(g['Xs'][rows], cols, g['hit'], g['miss'], g['pos'], g['el_s'], g['el_e'], g['reverse'])
    r = a9f.contrast(bv, np.ones_like(bv), g['sts'][rows], g['clone_s'][rows], p, So, g['ref2'])
    return label, r


def run_snv(a, out, t0):
    el_s, el_e = (int(x) for x in a.element.split(':'))
    pos, X, distal, passes, strand = a6x.load(a.extract)
    meta = np.load(Path(a.extract) / 'meta.npz')
    label, cls, clone = meta['label'], meta['cls'], meta['clone']
    wt, snv = np.flatnonzero(cls == 'WT'), np.flatnonzero(cls == 'SNV')
    Xw, Xs = X[wt], X[snv]
    del X
    dw, ds = distal[wt], distal[snv]
    stw, sts = a9f.strata(dw), a9f.strata(ds)
    cw, cs = clone[wt], clone[snv]
    lab_s = label[snv]
    masks = snv_masks(a.catalog, a.variants, pos, el_s, el_e, a.reverse, a.pm5)
    hit, miss = a7.emission_steps(a.emissions, pos)
    ecol = np.flatnonzero((pos >= el_s) & (pos < el_e))
    e1_oligo = ((el_e - 1 - pos[ecol]) if a.reverse else (pos[ecol] - el_s)) + OLIGO0
    ref1 = a9f.Reference(Xw[:, ecol] == 0, (Xw[:, ecol] >= 0) & (dw > a9f.E1_MIN)[:, None], stw, cw)
    p1, So1 = ref1.rates()
    k2 = np.flatnonzero(dw >= a9f.E2_MIN)
    chunks = np.array_split(k2, max(1, a.procs * 8))
    with mp.get_context('fork').Pool(a.procs) as pool:
        bw = np.concatenate(pool.map(_calls, [(Xw[c], hit, miss, pos, el_s, el_e, a.reverse) for c in chunks]))
    bw = bw.reshape(k2.size, -1)
    ref2 = a9f.Reference(bw, np.ones_like(bw), stw[k2], cw[k2])
    p2, So2 = ref2.rates()
    s2 = np.flatnonzero(ds >= a9f.E2_MIN)
    chunks = np.array_split(s2, max(1, a.procs * 8))
    with mp.get_context('fork').Pool(a.procs) as pool:
        bs = np.concatenate(pool.map(_calls, [(Xs[c], hit, miss, pos, el_s, el_e, a.reverse) for c in chunks]))
    bs = bs.reshape(s2.size, -1)
    pos_in_s2 = np.full(snv.size, -1)
    pos_in_s2[s2] = np.arange(s2.size)
    print(f'calls: {k2.size} WT, {s2.size} SNV molecules, {time.time() - t0:.0f} s', flush=True)

    snvs = sorted(set(lab_s) & set(masks))
    rows_of = {v: np.flatnonzero(lab_s == v) for v in snvs}
    res = {}
    for v in snvs:
        rr = rows_of[v]
        Y = Xs[rr][:, ecol] == 0
        O = (Xs[rr][:, ecol] >= 0) & (ds[rr] > a9f.E1_MIN)[:, None]
        r = a9f.contrast(Y, O, sts[rr], cs[rr], p1, So1, ref1)
        masked = np.isin(ecol, masks[v]['mask'])
        r['p'] = np.where(masked, np.nan, r['p'])
        r['masked'] = masked
        res[(v, 'E1')] = r
        if masks[v]['mask'].size == 0:
            r2 = rr[pos_in_s2[rr] >= 0]
            Yb = bs[pos_in_s2[r2]]
            res[(v, 'E2')] = a9f.contrast(Yb, np.ones_like(Yb), sts[r2], cs[r2], p2, So2, ref2)
    print(f'E1 and unmasked E2 done, {time.time() - t0:.0f} s', flush=True)
    todo = [v for v in snvs if masks[v]['mask'].size]
    _G.update(Xw2=Xw[k2], stw2=stw[k2], Xs=Xs, sts=sts, clone_s=cs, hit=hit, miss=miss, pos=pos, el_s=el_s,
              el_e=el_e, reverse=a.reverse, masks=masks, ref2=SimpleNamespace(V=ref2.V, W2=ref2.W2, tau2=ref2.tau2),
              snv_rows2={v: rows_of[v][ds[rows_of[v]] >= a9f.E2_MIN] for v in todo})
    with mp.get_context('fork').Pool(a.procs) as pool:
        for i, (v, r) in enumerate(pool.imap_unordered(_e2_masked, todo)):
            res[(v, 'E2')] = r
            if i % 50 == 0:
                print(f'masked E2 {i + 1}/{len(todo)}, {time.time() - t0:.0f} s', flush=True)

    frames = []
    e2_bin = np.repeat(BIN_NAMES, WIDTH)
    e2_oligo = np.tile(np.arange(WIDTH) + OLIGO0, len(BIN_NAMES))
    for (v, ep), r in res.items():
        f = pd.DataFrame({k: r[k] for k in ('n', 'C', 'D', 'se', 'p', 'expected')})
        f['feature'] = np.arange(f.shape[0])
        f['snv'], f['endpoint'], f['group'] = v, ep, 'drops its base' if masks[v]['group'] == 'drops its base' else 'rest'
        f['caller_class'], f['snv_oligo'] = masks[v]['group'], masks[v]['oligo']
        f['oligo'] = e1_oligo if ep == 'E1' else e2_oligo
        f['bin'] = 'adenine' if ep == 'E1' else e2_bin
        f['masked'] = r.get('masked', np.zeros(f.shape[0], bool))
        frames.append(f)
    eff = pd.concat(frames, ignore_index=True)
    eff['z'] = eff.D / eff.se
    eff['testable'] = eff.p.notna()
    mde = pd.read_parquet(Path(a.a9) / 'mde.parquet', columns=['snv', 'endpoint', 'feature', 'mde'])
    eff = eff.merge(mde, on=['snv', 'endpoint', 'feature'], how='left')

    # global tests and hierarchical FDR
    rows, eff['hit'] = [], False
    tiers = max_null(Path(a.a9))
    for (ep, grp), g in eff.groupby(['endpoint', 'group']):
        fam_p, feat_p, idx = {}, {}, {}
        for v, h in g.groupby('snv'):
            fam_p[v], feat_p[v], idx[v] = simes(h.p.to_numpy()), h.p.to_numpy(), h.index
        chosen, rej = bb(fam_p, feat_p)
        for v in fam_p:
            eff.loc[idx[v][rej[v]], 'hit'] = True
            h = g.loc[idx[v]]
            zmax = np.nanmax(np.abs(h.z.where(h.testable))) if h.testable.any() else np.nan
            C_v = int(np.unique(cs[rows_of[v]]).size)
            null = tiers.get((ep, a9f.tier(C_v)), np.array([]))
            rows.append(dict(snv=v, endpoint=ep, group=grp, caller_class=masks[v]['group'],
                             snv_oligo=masks[v]['oligo'], molecules=rows_of[v].size, clones=C_v,
                             testable=int(h.testable.sum()), masked=int(h.masked.sum()), simes_p=fam_p[v],
                             selected=v in chosen, hits=int(rej[v].sum()), max_abs_z=zmax,
                             maxz_p=(1 + (null >= zmax).sum()) / (1 + null.size) if null.size else np.nan))
    snv_t = pd.DataFrame(rows).sort_values(['endpoint', 'group', 'simes_p'])
    snv_t.to_csv(out / 'snvs.tsv', sep='\t', index=False)
    eff.to_parquet(out / 'effects.parquet', index=False)
    hits = eff[eff.hit].copy()
    hits['distance'] = hits.oligo - hits.snv_oligo
    hits.to_csv(out / 'hits.tsv', sep='\t', index=False)
    summ = {}
    for (ep, grp), g in snv_t.groupby(['endpoint', 'group']):
        h = hits[(hits.endpoint == ep) & (hits.group == grp)]
        summ[f'{ep} {grp}'] = dict(snvs=int(len(g)), selected=int(g.selected.sum()), hit_features=int(len(h)),
                                   maxz_p_lt_0_05=int((g.maxz_p < 0.05).sum()),
                                   median_abs_distance=float(h.distance.abs().median()) if len(h) else None)
    return summ, dict(wt=int(wt.size), snv_molecules=int(snv.size), snvs=len(snvs), masked=len(todo))


def max_null(a9dir):
    """Per (endpoint, clone tier): max |z| over each A9 pseudo-variant's testable features."""
    c = pd.read_parquet(a9dir / 'calibration.parquet', columns=['endpoint', 'snv', 'draw', 'D', 'se', 'p', 'tier'])
    c = c[c.p.notna()]
    c['az'] = (c.D / c.se).abs()
    m = c.groupby(['endpoint', 'tier', 'snv', 'draw']).az.max().reset_index()
    return {k: g.az.to_numpy() for k, g in m.groupby(['endpoint', 'tier'])}


def molecule_scores(Y, O, p, st):
    """Per molecule: mean of (y - p_WT) over observed columns (NaN if none)."""
    r = np.where(O, Y - p[st], 0.0).sum(1)
    n = O.sum(1)
    return np.divide(r, n, out=np.full(r.size, np.nan), where=n > 0)


def run_motif(a, out, t0):
    el_s, el_e = (int(x) for x in a.element.split(':'))
    pos, X, distal, passes, strand = a6x.load(a.extract)
    meta = np.load(Path(a.extract) / 'meta.npz')
    label, cls, clone = meta['label'], meta['cls'], meta['clone']
    wt, snv = np.flatnonzero(cls == 'WT'), np.flatnonzero(cls == 'SNV')
    Xw, Xs = X[wt], X[snv]
    del X
    dw, ds = distal[wt], distal[snv]
    stw, sts = a9f.strata(dw), a9f.strata(ds)
    cw, cs, lab_s = clone[wt], clone[snv], label[snv]
    masks = snv_masks(a.catalog, a.variants, pos, el_s, el_e, a.reverse)
    hit, miss = a7.emission_steps(a.emissions, pos)
    ecol = np.flatnonzero((pos >= el_s) & (pos < el_e))
    e1_oligo = ((el_e - 1 - pos[ecol]) if a.reverse else (pos[ecol] - el_s)) + OLIGO0
    k1w, k1s = dw > a9f.E1_MIN, ds > a9f.E1_MIN
    Sy, So = stratum_sums(((Xw[:, ecol] == 0) & (Xw[:, ecol] >= 0) & k1w[:, None]).astype(float), stw)
    So = np.zeros_like(Sy)
    for s in np.unique(stw):
        m = stw == s
        So[s] = ((Xw[m][:, ecol] >= 0) & k1w[m][:, None]).sum(0)
    p1 = np.divide(Sy, So, out=np.full_like(Sy, np.nan), where=So > 0)
    k2w, k2s = dw >= a9f.E2_MIN, ds >= a9f.E2_MIN
    # WT clone pools by size (A9 rule)
    u, inv = np.unique(cw, return_inverse=True)
    order = np.argsort(inv, kind='stable')
    pools = {}
    for g in np.split(order, np.cumsum(np.bincount(inv))[:-1]):
        pools.setdefault(g.size, []).append(g)
    rng = np.random.default_rng(SEED)
    rows = []
    cache = {}
    for name, lo, hi in MOTIFS:
        mcol = (e1_oligo >= lo) & (e1_oligo <= hi)
        win = np.arange(lo - OLIGO0, hi - OLIGO0 + 1)
        in_motif = [v for v in sorted(set(lab_s)) if v in masks and lo <= masks[v]['oligo'] <= hi]
        for grp in ('rest', 'drops its base'):
            pool_v = [v for v in in_motif if (masks[v]['group'] == 'drops its base') == (grp == 'drops its base')]
            if not pool_v:
                continue
            for ep in ('E1', 'E2'):
                obs, oclones, wt_scores, sizes = [], [], {}, {}
                for v in pool_v:
                    rr = np.flatnonzero(lab_s == v)
                    sizes[v] = np.bincount(np.unique(cs[rr], return_inverse=True)[1])
                    cols = mcol & ~np.isin(ecol, masks[v]['mask'])
                    if ep == 'E1':
                        Yv = Xs[rr][:, ecol[cols]] == 0
                        Ov = (Xs[rr][:, ecol[cols]] >= 0) & k1s[rr][:, None]
                        obs.append(molecule_scores(Yv, Ov, p1[:, cols], sts[rr]))
                        oclones.append(cs[rr])
                        Yw = Xw[:, ecol[cols]] == 0
                        Ow = (Xw[:, ecol[cols]] >= 0) & k1w[:, None]
                        wt_scores[v] = molecule_scores(Yw, Ow, p1[:, cols], stw)
                    else:
                        key = tuple(masks[v]['mask'])
                        if key not in cache:
                            bw = masked_calls(Xw[k2w], masks[v]['mask'], hit, miss, pos, el_s, el_e, a.reverse)
                            bw = bw.reshape(-1, len(BIN_NAMES), WIDTH)[:, list(MOTIF_BINS)].any(1)
                            Sy2, So2 = stratum_sums(bw.astype(float), stw[k2w])
                            cache[key] = (bw, np.divide(Sy2, So2, out=np.full_like(Sy2, np.nan), where=So2 > 0))
                        bw, p2 = cache[key]
                        r2 = rr[k2s[rr]]
                        bv = masked_calls(Xs[r2], masks[v]['mask'], hit, miss, pos, el_s, el_e, a.reverse)
                        bv = bv.reshape(-1, len(BIN_NAMES), WIDTH)[:, list(MOTIF_BINS)].any(1)
                        obs.append(molecule_scores(bv[:, win], np.ones((r2.size, win.size), bool), p2[:, win],
                                                   sts[r2]))
                        oclones.append(cs[r2])
                        sc = np.full(wt.size, np.nan)
                        sc[k2w] = molecule_scores(bw[:, win], np.ones((bw.shape[0], win.size), bool), p2[:, win],
                                                  stw[k2w])
                        wt_scores[v] = sc
                o = np.concatenate(obs)
                D = float(np.nanmean(o)) if np.isfinite(o).any() else np.nan
                null = np.empty(N_NULL)
                for b in range(N_NULL):
                    vals = [wt_scores[v][a9f.draw_clones(sizes[v], pools, rng)[0]] for v in pool_v]
                    null[b] = np.nanmean(np.concatenate(vals))
                pe = (1 + (np.abs(null) >= abs(D)).sum()) / (1 + N_NULL) if np.isfinite(D) else np.nan
                cm = clone_mean_test(pool_v, obs, oclones, wt_scores, cw, rng)
                rows.append(dict(motif=name, start=lo, end=hi, group=grp, endpoint=ep, snvs=len(pool_v),
                                 snv_list=','.join(pool_v), molecules=int(np.isfinite(o).sum()),
                                 clones=int(sum(s.size for s in sizes.values())), D=D, null_sd=float(np.std(null)),
                                 null_lo=float(np.quantile(null, 0.025)), null_hi=float(np.quantile(null, 0.975)),
                                 p=pe, **cm))
                print(f'{name} {grp} {ep}: D={D:.4f} p={pe:.4f}; clone-mean effect={cm["cm_effect"]:.4f} '
                      f'p={cm["cm_p"]:.4f} ({cm["cm_snv_clones"]} clones), {time.time() - t0:.0f} s', flush=True)
    t = pd.DataFrame(rows)
    t['bh'], t['cm_bh'] = False, False
    for ep, g in t.groupby('endpoint'):
        t.loc[g.index, 'bh'] = bh(g.p.to_numpy(), Q)
        t.loc[g.index, 'cm_bh'] = bh(g.cm_p.to_numpy(), Q)
    t.to_csv(out / 'motifs.tsv', sep='\t', index=False)
    return {ep: dict(tests=int(g.cm_p.notna().sum()), clone_mean_bh=int(g.cm_bh.sum()),
                     prespecified_bh_invalid=int(g.bh.sum())) for ep, g in t.groupby('endpoint')}, \
        dict(wt=int(wt.size), snv_molecules=int(snv.size))


def clone_means(scores, clones, min_n=MIN_CLONE_MOLECULES):
    """Mean score per clone with >= min_n finite scores: dict clone -> mean."""
    ok = np.isfinite(scores)
    f = pd.DataFrame(dict(c=clones[ok], s=scores[ok])).groupby('c').s.agg(['size', 'mean'])
    f = f[f['size'] >= min_n]
    return dict(zip(f.index, f['mean']))


def clone_mean_test(pool_v, obs, oclones, wt_scores, cw, rng):
    """HEK293T amendment 1: mean of the pool's SNV clone means against the same number of WT clone means (WT clones
    with >= MIN_CLONE_MOLECULES scoring molecules, without replacement, scored with each SNV's mask)."""
    slots, vals = [], []
    wm = [clone_means(wt_scores[v], cw) for v in pool_v]
    for k, (o, c) in enumerate(zip(obs, oclones)):
        for m in clone_means(o, c).values():
            slots.append(k)
            vals.append(m)
    elig = sorted(set.intersection(*[set(w) for w in wm])) if wm else []
    nan = dict(cm_snv_clones=len(vals), cm_wt_clones=len(elig), cm_stat=np.nan, cm_null_mean=np.nan,
               cm_null_sd=np.nan, cm_effect=np.nan, cm_p=np.nan)
    if not vals or len(elig) < len(vals):
        return nan
    W = np.array([[w[c] for w in wm] for c in elig])          # WT clones x SNVs in the pool
    slots = np.array(slots)
    stat = float(np.mean(vals))
    null = np.array([W[rng.permutation(len(elig))[:slots.size], slots].mean() for _ in range(N_CLONE_NULL)])
    mu = null.mean()
    eff = stat - mu
    return dict(cm_snv_clones=len(vals), cm_wt_clones=len(elig), cm_stat=stat, cm_null_mean=float(mu),
                cm_null_sd=float(null.std()), cm_effect=float(eff),
                cm_p=float((1 + (np.abs(null - mu) >= abs(eff)).sum()) / (1 + N_CLONE_NULL)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mode', choices=['snv', 'motif'], required=True)
    ap.add_argument('--extract', required=True)
    ap.add_argument('--emissions', required=True)
    ap.add_argument('--element', required=True)
    ap.add_argument('--reverse', action='store_true')
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--a9', help='A9 run dir of the same cell (snv mode: MDE and pseudo-variant max |z|)')
    ap.add_argument('--pm5', action='store_true', help='snv mode: +/-5 bp sensitivity mask (rules.json)')
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    t0 = time.time()
    summ, counts = (run_snv if a.mode == 'snv' else run_motif)(a, out, t0)
    info = dict(settings=vars(a), q=Q, counts=counts, results=summ, seconds=round(time.time() - t0))
    json.dump(info, open(out / 'summary.json', 'w'), indent=1, default=float)
    print(json.dumps(info, default=float))


if __name__ == '__main__':
    main()

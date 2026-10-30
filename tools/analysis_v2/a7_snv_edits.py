"""Step A7: what a designed variant does to plasmid-LLR footprint calls by itself (settings: rules.json snv_edits).

WT molecules (A6 extraction, classes C and M) are re-called with plasmid-LLR's own Kadane scan (global mode,
min_llr 5, drawdown 25, min_opps 3) after an in-silico edit, and compared with the same molecules unedited:
  drop           remove one element A/T observation: what the caller does to any SNV at that base (to G/C or the
                 A<->T swap) and to its deletion
  primary        remove a variant's symmetric-mask set (edited base when dropped + re-contexted adenines)
  pm5            remove every A/T within +/-5 bp of an element position (sensitivity mask)
  context        SNVs, model-based: drop, plus rescale methylation at re-contexted adenines by the naked pEF-GFP
                 relative-rate ratio within 0.05 distal strata (seeded)
Artifact = edited - unedited coverage per element position (oligo order), class and call-size bin.

Writes to --output: baseline.tsv, edits.tsv, artifacts.npz (edit x class x bin x element position), variants.tsv
(which edit each designed variant maps to), and, with --features, features.tsv (v1 feature occupancy artifacts) and,
with --v1-tests, v1_check.tsv. summary.json holds counts and settings.

Usage: python tools/analysis_v2/a7_snv_edits.py --extract <a6 extract dir> --element 1234:1445 [--reverse] \
           --emissions <fold-A emissions.tsv.gz> --catalog <catalog.json> --variants <a4 variants_visibility.tsv> \
           --naked-lookup <a4 naked_lookup.tsv> [--features <v1 features.tsv> --v1-tests <v1 tests.tsv.gz>] \
           --procs 4 --output <new dir>
"""
import argparse
import importlib.util
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from numba import njit

sys.path.insert(0, '/mmfs1/gscratch/stergachislab/bmallo/large_home/code/mine/plasmid-llr/llr')
import kadane  # noqa: E402


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(file))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


a6x = _load('a6x', 'a6_extract.py')
a4v = _load('a4v', 'a4_variant_visibility.py')

MIN_LLR, DRAWDOWN, MIN_OPPS = 5.0, 25.0, 3
U_MAX, M_MIN, STRATUM_WIDTH = 0.01, 0.60, 0.05
BINS = [('sub_TF', 10, 19), ('TF', 20, 40), ('PIC', 41, 80), ('NUC', 81, 100000)]
CHECK_BIN = ('20-80', 20, 80)               # for the comparison with the existing calls (A6 m_calls.tsv)
PM = 5
WIDTH = 211
OLIGO0 = 21
CLASSES = ['C', 'M']
SEED = 20261029

_scan = kadane._scan


@njit
def scan_cov(X, hit, miss, pos, ov_cols, ov_vals, el_s, el_e, reverse, b_lo, b_hi, f_mlo, f_mhi, f_llo, f_lhi):
    """Call every row of X (molecule x A/T column, -1/0/1) with columns ov_cols replaced by ov_vals (row x k).
    Returns coverage counts (bin x element position, oligo order) and feature counts (molecules with a call whose
    midpoint element index is in [f_mlo, f_mhi] and length in [f_llo, f_lhi])."""
    n, P = X.shape
    nb = b_lo.size
    nf = f_mlo.size
    cov = np.zeros((nb, el_e - el_s), np.int64)
    feat = np.zeros(nf, np.int64)
    fl = np.zeros(nf, np.bool_)
    step = np.zeros(P)
    ih = np.zeros(P, np.bool_)
    ii = np.zeros(P, np.bool_)
    for r in range(n):
        for j in range(P):
            x = X[r, j]
            ii[j] = x >= 0
            ih[j] = x == 1
            step[j] = hit[j] if x == 1 else (miss[j] if x == 0 else 0.0)
        for k in range(ov_cols.size):
            j = ov_cols[k]
            x = ov_vals[r, k]
            ii[j] = x >= 0
            ih[j] = x == 1
            step[j] = hit[j] if x == 1 else (miss[j] if x == 0 else 0.0)
        lo = -1
        hi = 0
        for j in range(P):
            if ii[j]:
                if lo < 0:
                    lo = j
                hi = j + 1
        if lo < 0:
            continue
        s, e, _l, _o, _la, _ra = _scan(step, ih, ii, lo, hi, MIN_LLR, MIN_OPPS, DRAWDOWN)
        for f in range(nf):
            fl[f] = False
        for c in range(s.size):
            a = pos[s[c]]
            b = pos[e[c] - 1] + 1
            length = b - a
            lo2 = max(a, el_s)
            hi2 = min(b, el_e)
            for q in range(nb):
                if length >= b_lo[q] and length <= b_hi[q]:
                    for rr in range(lo2, hi2):
                        jj = (el_e - 1 - rr) if reverse else (rr - el_s)
                        cov[q, jj] += 1
            mid = np.rint((a + b - 1) / 2.0)
            mi = int((el_e - 1) - mid) if reverse else int(mid - el_s)
            for f in range(nf):
                if f_mlo[f] <= mi <= f_mhi[f] and f_llo[f] <= length <= f_lhi[f]:
                    fl[f] = True
        for f in range(nf):
            if fl[f]:
                feat[f] += 1
    return cov, feat


def classes(distal):
    return np.where(distal <= U_MAX, 'U', np.where(distal >= M_MIN, 'M', 'C'))


def emission_steps(path, pos):
    t = pd.read_csv(path, sep='\t')
    t = t[t.fold == 'A']
    m = dict(zip(t.ref_pos.to_numpy() - 1, zip(t.p_a.to_numpy(), t.p_p.to_numpy())))
    missing = [p for p in pos if p not in m]
    if missing:
        raise ValueError(f'{len(missing)} A/T positions have no emission')
    pa = np.array([m[p][0] for p in pos])
    pp = np.array([m[p][1] for p in pos])
    return np.log(pp / pa), np.log((1 - pp) / (1 - pa))


def recontexted(wt, var_seq, cls, i):
    """[(wt index, wt class, variant class)] for adenines other than the edited base whose (run, position) class
    changes (a4_variant_visibility.compare's rule, returning the positions)."""
    cw, cv = a4v.adenine_classes(wt), a4v.adenine_classes(var_seq)

    def to_wt(k):
        if cls == 'DEL':
            return k + 1 if k >= i else k
        if cls == 'INS':
            return None if k == i else (k - 1 if k > i else k)
        return k
    out = []
    for (s, k), c in cv.items():
        kw = to_wt(k)
        if kw is None or (cls == 'SNV' and kw == i):
            continue
        if (s, kw) in cw and cw[(s, kw)] != c:
            out.append((kw, cw[(s, kw)], c))
    return sorted(out)


def mask_sets(var, wt, seqs):
    """Per variant: dropped index (or None), re-contexted list, primary mask (sorted element indices)."""
    rows = []
    for v in var.itertuples():
        drop = int(v.index0) if str(v.caller_drops_at) == 'True' else None
        rc = recontexted(wt, seqs[v.variant_id], v.cls, int(v.index0))
        prim = sorted(set([k for k, _, _ in rc] + ([drop] if drop is not None else [])))
        rows.append(dict(source_variant=v.source_variant, cls=v.cls, index0=int(v.index0), drop=drop,
                         recontext=rc, primary=prim))
    return rows


def pm_set(wt, j, pm=PM):
    return [k for k in range(max(0, j - pm), min(len(wt), j + pm + 1)) if wt[k] in 'AT']


_G = {}


def _work(task):
    """task = (edit index, class, element indices to override, mode, payload); returns counts."""
    eid, c, idx, mode, payload = task
    g = _G[c]
    cols = np.array([g['col_of_el'][k] for k in idx], np.int64)
    X = g['X']
    if mode == 'context':
        vals = np.full((X.shape[0], cols.size), -1, np.int8)
        rng = np.random.default_rng([SEED, payload['seed_key'], CLASSES.index(c)])
        for q, (k, ratio) in enumerate(zip(idx, payload['ratios'])):
            if ratio is None:                    # the dropped base
                continue
            col = cols[q]
            x = X[:, col].copy()
            if ratio < 1:
                flip = (x == 1) & (rng.random(x.size) < 1 - ratio)
                x[flip] = 0
            elif ratio > 1:
                p = g['rate'][col][g['stratum']]
                pn = np.minimum(0.999, ratio * p)
                q_ = np.where(p < 1, (pn - p) / np.maximum(1 - p, 1e-12), 0.0)
                flip = (x == 0) & (rng.random(x.size) < q_)
                x[flip] = 1
            vals[:, q] = x
    else:
        vals = np.full((X.shape[0], cols.size), -1, np.int8)
    cov, feat = scan_cov(X, g['hit'], g['miss'], g['pos'], cols, vals, g['el_s'], g['el_e'], g['reverse'],
                         g['b_lo'], g['b_hi'], g['f_mlo'], g['f_mhi'], g['f_llo'], g['f_lhi'])
    return eid, c, cov, feat


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--extract', required=True)
    ap.add_argument('--element', required=True, help='start0:end0 of the 211-bp element on the reference')
    ap.add_argument('--reverse', action='store_true', help='element is reverse-complemented on the reference')
    ap.add_argument('--emissions', required=True)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--naked-lookup', required=True)
    ap.add_argument('--features')
    ap.add_argument('--v1-tests')
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--max-molecules', type=int, default=0, help='0 = all (for tests only)')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    t0 = time.time()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    el_s, el_e = map(int, a.element.split(':'))
    pos, X, distal, _passes, _strand = a6x.load(a.extract)
    ref = json.load(open(Path(a.extract) / 'extract.json'))['reference']
    seq = ''.join(l.strip() for l in open(ref) if not l.startswith('>')).upper()
    lib = json.load(open(a.catalog))['library']
    wt = lib['mutable_sequence']
    el_seq = seq[el_s:el_e]
    comp = str.maketrans('ACGT', 'TGCA')
    if (el_seq.translate(comp)[::-1] if a.reverse else el_seq) != wt:
        raise ValueError('reference element differs from the catalog WT design')
    seqs = {v['variant_id']: v['sequence'] for v in lib['variants']}
    hit, miss = emission_steps(a.emissions, pos)
    col_of_ref = {int(p): c for c, p in enumerate(pos)}
    el_ref = [(el_e - 1 - j) if a.reverse else (el_s + j) for j in range(WIDTH)]
    col_of_el = {j: col_of_ref[r] for j, r in enumerate(el_ref) if r in col_of_ref}
    if sorted(col_of_el) != [j for j in range(WIDTH) if wt[j] in 'AT']:
        raise ValueError('element A/T columns do not match the WT design')

    cls = classes(distal)
    lut = {tuple([int(k.split('|')[0]), k.split('|')[1]]): float(v)
           for k, v in pd.read_csv(a.naked_lookup, sep='\t', index_col=0).p_rel.items()}
    feats = pd.read_csv(a.features, sep='\t') if a.features else None
    f_arr = (feats.mid_lo.to_numpy(np.int64), feats.mid_hi.to_numpy(np.int64), feats.len_lo.to_numpy(np.int64),
             feats.len_hi.to_numpy(np.int64)) if feats is not None else tuple(np.zeros(0, np.int64) for _ in range(4))
    bins = BINS + [CHECK_BIN]
    b_lo = np.array([b[1] for b in bins], np.int64)
    b_hi = np.array([b[2] for b in bins], np.int64)
    n_mol = {}
    for c in CLASSES:
        rows = np.flatnonzero(cls == c)
        if a.max_molecules:
            rows = rows[:a.max_molecules]
        Xc = np.ascontiguousarray(X[rows])
        st = np.minimum((distal[rows] / STRATUM_WIDTH).astype(int), int(round(1 / STRATUM_WIDTH)) - 1)
        rate = {}
        for col in col_of_el.values():
            v = Xc[:, col]
            num = np.bincount(st, weights=(v == 1), minlength=20)
            den = np.bincount(st, weights=(v >= 0), minlength=20)
            rate[col] = num / np.maximum(den, 1)
        _G[c] = dict(X=Xc, hit=hit, miss=miss, pos=pos.astype(np.int64), col_of_el=col_of_el, el_s=el_s, el_e=el_e,
                     reverse=a.reverse, b_lo=b_lo, b_hi=b_hi, f_mlo=f_arr[0], f_mhi=f_arr[1], f_llo=f_arr[2],
                     f_lhi=f_arr[3], stratum=st, rate=rate)
        n_mol[c] = int(rows.size)
    del X

    var = pd.read_csv(a.variants, sep='\t')
    ms = mask_sets(var, wt, seqs)
    edits = [dict(kind='baseline', centre=-1, idx=(), mode='mask', payload=None)]
    key_of = {}

    def add(kind, centre, idx, mode='mask', payload=None):
        k = (kind, tuple(idx)) if mode == 'mask' else (kind, centre, payload['seed_key'])
        if k not in key_of:
            key_of[k] = len(edits)
            edits.append(dict(kind=kind, centre=centre, idx=tuple(idx), mode=mode, payload=payload))
        return key_of[k]
    for j in col_of_el:
        add('drop', j, [j])
    for j in range(WIDTH):
        if pm_set(wt, j):
            add('pm5', j, pm_set(wt, j))
    vrows = []
    for n, m in enumerate(ms):
        r = dict(source_variant=m['source_variant'], cls=m['cls'], index0=m['index0'],
                 caller_drops=m['drop'] is not None, n_recontexted=len(m['recontext']), n_primary=len(m['primary']),
                 primary_oligo=','.join(str(k + OLIGO0) for k in m['primary']))
        r['edit_drop'] = key_of[('drop', (m['drop'],))] if m['drop'] is not None else 0
        r['edit_primary'] = add('primary', m['index0'], m['primary']) if m['primary'] else 0
        pm = pm_set(wt, m['index0'])
        r['edit_pm5'] = add('pm5', m['index0'], pm) if pm else 0
        r['edit_context'] = -1
        if m['cls'] == 'SNV' and m['recontext']:
            idx = ([m['drop']] if m['drop'] is not None else []) + [k for k, _, _ in m['recontext']]
            ratios = ([None] if m['drop'] is not None else []) + [lut[cv] / lut[cw] for _, cw, cv in m['recontext']]
            r['edit_context'] = add('context', m['index0'], idx, 'context', dict(ratios=ratios, seed_key=n))
        elif m['cls'] == 'SNV':
            r['edit_context'] = r['edit_drop']           # no re-contexted adenine: the context edit is the drop
        vrows.append(r)

    tasks = [(e, c, ed['idx'], ed['mode'], ed['payload']) for e, ed in enumerate(edits) for c in CLASSES]
    # compile once before forking
    g = _G['C']
    scan_cov(g['X'][:2], g['hit'], g['miss'], g['pos'], np.zeros(0, np.int64), np.zeros((2, 0), np.int8), el_s,
             el_e, a.reverse, b_lo, b_hi, *f_arr)
    print(f'{len(edits)} edits x {len(CLASSES)} classes; molecules {n_mol}', flush=True)
    cov = np.zeros((len(edits), len(CLASSES), len(bins), WIDTH), np.float64)
    fea = np.zeros((len(edits), len(CLASSES), f_arr[0].size), np.float64)
    with mp.get_context('fork').Pool(a.procs) as pool:
        for i, (e, c, cv, ft) in enumerate(pool.imap_unordered(_work, tasks, chunksize=4)):
            k = CLASSES.index(c)
            cov[e, k] = cv / n_mol[c]
            fea[e, k] = ft / n_mol[c]
            if i % 200 == 0:
                print(f'  {i}/{len(tasks)} tasks, {time.time() - t0:.0f} s', flush=True)

    base = cov[0]
    delta = (cov - base[None])[:, :, :len(BINS)].astype(np.float32)
    np.savez_compressed(out / 'artifacts.npz', delta=delta, classes=np.array(CLASSES),
                        bins=np.array([b[0] for b in BINS]))
    bl = []
    for k, c in enumerate(CLASSES):
        for q, b in enumerate(bins):
            bl.append(pd.DataFrame(dict(cls=c, bin=b[0], index0=np.arange(WIDTH), oligo_pos1=np.arange(WIDTH) + OLIGO0,
                                        coverage=base[k, q])))
    pd.concat(bl).to_csv(out / 'baseline.tsv', sep='\t', index=False)
    er = []
    for e, ed in enumerate(edits):
        r = dict(edit=e, kind=ed['kind'], centre_oligo=ed['centre'] + OLIGO0 if ed['centre'] >= 0 else -1,
                 n_removed=len(ed['idx']) if ed['mode'] == 'mask' else int(sum(x is None for x in ed['payload']['ratios'])),
                 n_rescaled=0 if ed['mode'] == 'mask' else int(sum(x is not None for x in ed['payload']['ratios'])),
                 positions_oligo=','.join(str(k + OLIGO0) for k in ed['idx']))
        far = np.abs(np.arange(WIDTH) - ed['centre']) > PM if ed['centre'] >= 0 else np.ones(WIDTH, bool)
        for k, c in enumerate(CLASSES):
            d = delta[e, k]
            r[f'{c}_max_abs'] = float(np.abs(d).max())
            r[f'{c}_max_abs_far'] = float(np.abs(d[:, far]).max()) if far.any() else 0.0
            r[f'{c}_sum_bp'] = float(d.sum())
            for q, b in enumerate(BINS):
                r[f'{c}_{b[0]}_sum_bp'] = float(d[q].sum())
        er.append(r)
    pd.DataFrame(er).to_csv(out / 'edits.tsv', sep='\t', index=False)
    vt = pd.DataFrame(vrows)
    vt.to_csv(out / 'variants.tsv', sep='\t', index=False)
    summary = dict(molecules=n_mol, edits={k: int(sum(ed['kind'] == k for ed in edits)) for k in
                                         ['drop', 'pm5', 'primary', 'context']},
                   settings=dict(min_llr=MIN_LLR, drawdown=DRAWDOWN, min_opps=MIN_OPPS, bins=BINS, pm=PM, seed=SEED,
                                 element=[el_s, el_e], reverse=a.reverse, emissions=a.emissions,
                                 extract=a.extract))
    if feats is not None:
        fr = []
        for e, ed in enumerate(edits):
            if ed['kind'] not in ('drop', 'context'):
                continue
            for f, name in enumerate(feats.feature):
                fr.append(dict(edit=e, kind=ed['kind'], centre_oligo=ed['centre'] + OLIGO0, feature=name,
                               dC=fea[e, 0, f] - fea[0, 0, f], dM=fea[e, 1, f] - fea[0, 1, f]))
        ft = pd.DataFrame(fr)
        ft.to_csv(out / 'features.tsv', sep='\t', index=False)
        pd.DataFrame(dict(feature=feats.feature, C=fea[0, 0], M=fea[0, 1])).to_csv(out / 'features_baseline.tsv',
                                                                                   sep='\t', index=False)
        if a.v1_tests:
            summary['v1_check'] = v1_check(a.v1_tests, ft, vt, var, out)
    summary['seconds'] = round(time.time() - t0)
    json.dump(summary, open(out / 'summary.json', 'w'), indent=1, default=str)
    print(json.dumps(summary, indent=1, default=str))


def caller_class(v):
    return 'drops its base' if str(v.caller_drops_at) == 'True' else ('gains an A/T' if v.at_count_change > 0
                                                                     else 'neither')


def v1_check(path, ft, vt, var, out):
    """Pooled over v1 (SNV, feature) pairs with the SNV inside the feature: A7 artifact vs v1 observed residual."""
    t = pd.read_csv(path, sep='\t')
    cc = {v.source_variant: caller_class(v) for v in var.itertuples()}
    ed = vt.set_index('source_variant')
    look = ft.set_index(['edit', 'feature'])
    rows = []
    for r in t.itertuples():
        if r.source_variant not in ed.index:
            continue
        v = ed.loc[r.source_variant]
        out_r = dict(source_variant=r.source_variant, feature=r.feature, inside=bool(r.snv_in_footprint),
                     caller_class=cc[r.source_variant], v1_resid=r.occupancy_resid_delta, slope_b=r.slope_b)
        for kind, e in (('drop', v.edit_drop), ('context', v.edit_context)):
            if e > 0:
                x = look.loc[(e, r.feature)]
                out_r[f'{kind}_dC'], out_r[f'{kind}_dM'] = x.dC, x.dM
                out_r[f'{kind}_resid'] = x.dC - r.slope_b * x.dM
            else:
                out_r[f'{kind}_dC'] = out_r[f'{kind}_dM'] = out_r[f'{kind}_resid'] = 0.0
        rows.append(out_r)
    d = pd.DataFrame(rows)
    d.to_csv(out / 'v1_pairs.tsv', sep='\t', index=False)
    s = []
    for (inside, c), g in d.groupby(['inside', 'caller_class']):
        s.append(dict(inside=inside, caller_class=c, pairs=len(g), v1_resid_mean_pp=100 * g.v1_resid.mean(),
                      drop_resid_mean_pp=100 * g.drop_resid.mean(), context_resid_mean_pp=100 * g.context_resid.mean(),
                      r_v1_vs_context=float(np.corrcoef(g.v1_resid, g.context_resid)[0, 1])
                      if g.context_resid.std() > 0 else float('nan')))
    s = pd.DataFrame(s)
    s.to_csv(out / 'v1_check.tsv', sep='\t', index=False)
    return s.round(4).to_dict(orient='records')


if __name__ == '__main__':
    main()

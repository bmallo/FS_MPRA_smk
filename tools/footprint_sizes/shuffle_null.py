"""Within-molecule co-protection test: are footprint calls more than per-base sequence bias would give?

A methyltransferase sequence preference acts base by base: it lowers the methylation rate at some
adenines on every molecule, independently. A bound protein protects a run of adjacent adenines on
the same molecule. Permuting each base's methylation state across molecules (within a labeling
class/stratum) keeps every per-base rate, so all sequence bias, labeling level and A/T spacing,
and destroys only the within-molecule co-protection. Re-running the same LLR scan on the permuted
molecules gives the calls that sequence bias + labeling alone produce in each feature box.

Observed and permuted molecules are scanned over the same window (element +- flank) with the
same emissions, min_llr, min_opps and drawdown as the production calls.

Usage: python tools/footprint_sizes/shuffle_null.py --bam wt_subset.bam --reference ref.fa \
           --emissions emissions.tsv.gz --extract <hepg2 extract_llr> --features features.tsv \
           --config config/ldlr_footprint_sizes.json --output <dir>
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from fsmpra.footprint_sizes import WIDTH, classify
from fsmpra.io import open_text, read_json

LLR_REPO = Path('/mmfs1/gscratch/stergachislab/bmallo/large_home/code/mine/plasmid-llr')
sys.path.insert(0, str(LLR_REPO / 'llr'))
import kadane  # noqa: E402
import llrlib  # noqa: E402


def emissions(path, L, fold='A'):
    t = pd.read_csv(path, sep='\t')
    t = t[t.fold == fold]
    hit, miss = np.zeros(L), np.zeros(L)
    p = t.ref_pos.to_numpy() - 1
    hit[p] = np.log(t.p_p.to_numpy() / t.p_a.to_numpy())
    miss[p] = np.log((1 - t.p_p.to_numpy()) / (1 - t.p_a.to_numpy()))
    return hit, miss


def box_hits(X, lo, hit, miss, feats, element, min_llr, min_opps, drawdown):
    """Per molecule x feature: any call in the box. X codes -1 uninformative, 0 unmethylated, 1 m6A."""
    n, W = X.shape
    out = np.zeros((n, len(feats)), bool)
    step = np.zeros(W)
    for i in range(n):
        x = X[i]
        info = x >= 0
        isHit = x == 1
        step[:] = 0.0
        step[isHit] = hit[isHit]
        step[info & ~isHit] = miss[info & ~isHit]
        s, e, *_ = kadane.scan(step, isHit, info, 0, W, min_llr, min_opps, drawdown)
        if not len(s):
            continue
        s0, e0 = s + lo, e + lo
        mid = np.rint((s0 + e0 - 1) / 2 - element[0]).astype(int)
        length = e0 - s0
        for k, f in enumerate(feats):
            out[i, k] = bool(((mid >= f['mid_lo']) & (mid <= f['mid_hi']) & (length >= f['len_lo']) & (length <= f['len_hi'])).any())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ('--bam', '--reference', '--emissions', '--extract', '--features', '--config', '--output'):
        ap.add_argument(k, required=True)
    ap.add_argument('--element', nargs=2, type=int, default=[1234, 1445])
    ap.add_argument('--flank', type=int, default=150)
    ap.add_argument('--shuffles', type=int, default=3)
    ap.add_argument('--min-at', type=int, default=500)
    ap.add_argument('--min-opps', type=int, default=3)
    ap.add_argument('--seed', type=int, default=20261030)
    a = ap.parse_args()
    cfg = read_json(a.config)
    spec = read_json(Path(cfg['base_config']))['analysis']
    min_llr, dd = cfg['calls']['min_llr'], cfg['calls']['drawdown']
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)

    seq = ''.join(l.strip() for l in open(a.reference) if not l.startswith('>')).upper()
    ref_arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    at_mask = (ref_arr == ord('A')) | (ref_arr == ord('T'))
    hit, miss = emissions(a.emissions, len(seq))
    lo, hi = a.element[0] - a.flank, a.element[1] + a.flank
    hit_w, miss_w = hit[lo:hi], miss[lo:hi]

    extract = Path(a.extract)
    _, _, _, groups, weights, _, _ = classify(extract, spec)
    with open_text(extract / 'molecules.tsv.gz') as h:
        index = {r['read_id']: int(r['row']) for r in csv.DictReader(h, delimiter='\t')}
    feats = pd.read_csv(a.features, sep='\t').to_dict('records')

    rows, mats = [], []
    for rec in pysam.AlignmentFile(a.bam, 'rb'):
        if not llrlib.usable(rec) or rec.query_name not in index:
            continue
        res = llrlib.read_observations(rec, ref_arr, at_mask)
        if res is None or res[0].size < a.min_at:
            continue
        rp, obs = res
        x = np.full(hi - lo, -1, np.int8)
        k = (rp >= lo) & (rp < hi)
        x[rp[k] - lo] = obs[k].astype(np.int8)
        g = groups[index[rec.query_name]]
        if g < 0:
            continue
        rows.append(g)
        mats.append(x)
    g = np.array(rows)
    X = np.vstack(mats)
    print(f'{len(g):,} molecules; per group {np.bincount(g, minlength=12).tolist()}', flush=True)

    scan = dict(lo=lo, hit=hit_w, miss=miss_w, feats=feats, element=a.element, min_llr=min_llr, min_opps=a.min_opps, drawdown=dd)
    obs_hits = box_hits(X, **scan)
    rng = np.random.default_rng(a.seed)
    shuf_hits = []
    for r in range(a.shuffles):
        Y = X.copy()
        for grp in np.unique(g):
            idx = np.flatnonzero(g == grp)
            perm = np.argsort(rng.random((len(idx), X.shape[1])), axis=0)
            Y[idx] = np.take_along_axis(X[idx], perm, axis=0)
        shuf_hits.append(box_hits(Y, **scan))
        print(f'shuffle {r + 1} done', flush=True)

    def occ(H):
        frac = np.array([H[g == k].mean(0) if (g == k).any() else np.zeros(H.shape[1]) for k in range(2 + len(weights))])
        return frac[2:].T @ weights, frac[1], frac[0]

    oC, oM, oU = occ(obs_hits)
    S = np.array([occ(H) for H in shuf_hits])
    res = pd.DataFrame(dict(feature=[f['feature'] for f in feats], centre_oligo=[f['mid_peak_oligo'] for f in feats],
                            length_peak=[f['len_peak'] for f in feats],
                            observed_C_std=oC, shuffled_C_std=S[:, 0].mean(0), shuffled_C_std_sd=S[:, 0].std(0),
                            observed_M=oM, shuffled_M=S[:, 1].mean(0), shuffled_M_sd=S[:, 1].std(0), observed_U=oU, shuffled_U=S[:, 2].mean(0)))
    res['C_excess_over_shuffle'] = res.observed_C_std - res.shuffled_C_std
    res['M_excess_over_shuffle'] = res.observed_M - res.shuffled_M
    res['C_fold'] = res.observed_C_std / res.shuffled_C_std.where(res.shuffled_C_std > 0)
    res['M_fold'] = res.observed_M / res.shuffled_M.where(res.shuffled_M > 0)
    res.to_csv(out / 'shuffle_null.tsv', sep='\t', index=False)
    json.dump(dict(molecules=int(len(g)), per_group=np.bincount(g, minlength=12).tolist(), window=[lo, hi], shuffles=a.shuffles,
                   min_llr=min_llr, drawdown=dd, min_opps=a.min_opps, seed=a.seed, bam=a.bam, emissions=a.emissions,
                   rule='each A/T column permuted across molecules within its U / M / C-stratum class'),
              open(out / 'summary.json', 'w'), indent=1)
    pd.set_option('display.width', 220)
    print(res.round(5).to_string(index=False))


if __name__ == '__main__':
    main()

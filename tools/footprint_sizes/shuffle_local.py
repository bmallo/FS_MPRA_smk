"""Locally conditioned co-protection test for each size-resolved footprint.

shuffle_null.py permutes every base across molecules. On chromatin molecules that turns nucleosome
protection into scattered unmethylated bases, which inflates the chance rate of short calls, so it is
not a sequence-bias null there. Here, per feature:
  * molecules = those with no call >= 60 bp overlapping the feature span +- 20 bp (locally free of
    nucleosome-sized protection), within each class (M; each C labeling stratum);
  * only the A/T columns inside span +- 20 bp are permuted across those molecules (within stratum);
    all other bases keep their real values, so flanking context and the scan are unchanged;
  * the same LLR scan is re-run and the feature box counted.
Per-base rates, sequence bias and labeling are kept; only within-molecule co-protection at the
footprint is destroyed. Observed / null occupancy is reported for C_std and M.

Usage: python tools/footprint_sizes/shuffle_local.py <same arguments as shuffle_null.py>
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from shuffle_null import emissions, kadane, llrlib  # noqa: E402
from fsmpra.footprint_sizes import classify  # noqa: E402
from fsmpra.io import open_text, read_json  # noqa: E402

PAD, LONG = 20, 60


def scan_one(x, hit, miss, min_llr, min_opps, dd):
    info = x >= 0
    is_hit = x == 1
    step = np.zeros(len(x))
    step[is_hit] = hit[is_hit]
    step[info & ~is_hit] = miss[info & ~is_hit]
    s, e, *_ = kadane.scan(step, is_hit, info, 0, len(x), min_llr, min_opps, dd)
    return s, e


def in_box(s, e, lo, el0, f):
    s0, e0 = s + lo, e + lo
    mid = np.rint((s0 + e0 - 1) / 2 - el0).astype(int)
    ln = e0 - s0
    return bool(((mid >= f['mid_lo']) & (mid <= f['mid_hi']) & (ln >= f['len_lo']) & (ln <= f['len_hi'])).any())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ('--bam', '--reference', '--emissions', '--extract', '--features', '--config', '--output'):
        ap.add_argument(k, required=True)
    ap.add_argument('--element', nargs=2, type=int, default=[1234, 1445])
    ap.add_argument('--flank', type=int, default=150)
    ap.add_argument('--shuffles', type=int, default=2)
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
    el0 = a.element[0]
    lo, hi = el0 - a.flank, a.element[1] + a.flank
    hit_w, miss_w = hit[lo:hi], miss[lo:hi]
    _, _, _, groups, weights, _, _ = classify(Path(a.extract), spec)
    with open_text(Path(a.extract) / 'molecules.tsv.gz') as h:
        index = {r['read_id']: int(r['row']) for r in csv.DictReader(h, delimiter='\t')}
    feats = pd.read_csv(a.features, sep='\t').to_dict('records')

    g, mats = [], []
    for rec in pysam.AlignmentFile(a.bam, 'rb'):
        if not llrlib.usable(rec) or rec.query_name not in index:
            continue
        res = llrlib.read_observations(rec, ref_arr, at_mask)
        if res is None or res[0].size < a.min_at:
            continue
        grp = groups[index[rec.query_name]]
        if grp < 1:  # U molecules carry no information here
            continue
        rp, obs = res
        x = np.full(hi - lo, -1, np.int8)
        k = (rp >= lo) & (rp < hi)
        x[rp[k] - lo] = obs[k].astype(np.int8)
        g.append(grp)
        mats.append(x)
    g, X = np.array(g), np.vstack(mats)
    print(f'{len(g):,} M/C molecules', flush=True)
    calls = [scan_one(X[i], hit_w, miss_w, min_llr, a.min_opps, dd) for i in range(len(g))]
    print('observed scan done', flush=True)

    rng = np.random.default_rng(a.seed)
    rows = []
    for f in feats:
        sa, sb = f['mid_peak'] - f['len_peak'] / 2, f['mid_peak'] + f['len_peak'] / 2
        wa, wb = int(np.floor(sa)) - PAD + el0 - lo, int(np.ceil(sb)) + PAD + el0 - lo + 1  # window columns [wa, wb)
        free = np.array([not (((e - s) >= LONG) & (s < wb) & (e > wa)).any() for s, e in calls])
        obs = np.array([in_box(s, e, lo, el0, f) for s, e in calls])
        null = np.zeros((a.shuffles, len(g)))
        for r in range(a.shuffles):
            for grp in np.unique(g):
                idx = np.flatnonzero((g == grp) & free)
                if len(idx) < 2:
                    continue
                Y = X[idx].copy()
                blk = Y[:, wa:wb]
                perm = np.argsort(rng.random(blk.shape), axis=0)
                Y[:, wa:wb] = np.take_along_axis(blk, perm, axis=0)
                for j, i in enumerate(idx):
                    s, e = scan_one(Y[j], hit_w, miss_w, min_llr, a.min_opps, dd)
                    null[r, i] = in_box(s, e, lo, el0, f)

        def occ(v):
            fr = np.array([v[(g == k) & free].mean() if ((g == k) & free).any() else np.nan for k in range(2 + len(weights))])
            return float(np.nansum(fr[2:] * weights)), float(fr[1])
        oc, om = occ(obs.astype(float))
        nn = np.array([occ(null[r]) for r in range(a.shuffles)])
        rows.append(dict(feature=f['feature'], centre_oligo=f['mid_peak_oligo'], length_peak=f['len_peak'],
                         free_fraction_C=float(free[g >= 2].mean()), free_fraction_M=float(free[g == 1].mean()),
                         observed_C_std_free=oc, null_C_std_free=nn[:, 0].mean(), observed_M_free=om, null_M_free=nn[:, 1].mean(),
                         C_fold=oc / nn[:, 0].mean() if nn[:, 0].mean() > 0 else np.nan,
                         M_fold=om / nn[:, 1].mean() if nn[:, 1].mean() > 0 else np.nan))
        print(rows[-1], flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(out / 'shuffle_local.tsv', sep='\t', index=False)
    json.dump(dict(molecules=int(len(g)), pad_bp=PAD, long_call_bp=LONG, shuffles=a.shuffles, seed=a.seed, window=[lo, hi],
                   rule='per feature: molecules without a >=60 bp call overlapping span +-20 bp; A/T columns inside span +-20 bp '
                        'permuted across them within class/stratum; rest of the molecule unchanged; same LLR scan'),
              open(out / 'summary.json', 'w'), indent=1)
    pd.set_option('display.width', 220)
    print(t.round(5).to_string(index=False))


if __name__ == '__main__':
    main()

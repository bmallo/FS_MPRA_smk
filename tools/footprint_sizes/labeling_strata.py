"""Footprint occupancy by per-molecule labeling level, HepG2 (5 min Hia5) vs HEK293T (10 min).

Descriptive. For each cell, WT molecules are split into U, M and the 10 C labeling strata (each
cell's own distal-m6A deciles), and the HepG2 feature boxes are counted per group. A footprint
that is stable over the reaction should stay roughly flat with labeling once above the detection
floor; one from a fast-exchanging factor would thin out as labeling (enzyme exposure) rises.
Detection sensitivity also rises with labeling, so the shape is read alongside the shuffle null.

Usage: python tools/footprint_sizes/labeling_strata.py --config config/ldlr_footprint_sizes.json \
           --features <hepg2_map>/features.tsv --output <dir>
"""
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd

from fsmpra.footprint_sizes import box_mass, classify, density, load_calls
from fsmpra.io import read_json

V = 'results/variant_effects_20261028'
CELLS = {'HepG2': (f'{V}/hepg2/extract_llr', f'{V}/hepg2/llr/results/calls/shard*.llr_global_foldA_pp0_dd25.parquet', (1234, 1445), False),
         'HEK293T': (f'{V}/hek293t/extract_llr', f'{V}/hek293t/llr/results/calls/shard*.llr_global_foldA_pp0_dd25.parquet', (3236, 3447), True)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ('--config', '--features', '--output'):
        ap.add_argument(k, required=True)
    a = ap.parse_args()
    cfg = read_json(a.config)
    spec = read_json(Path(cfg['base_config']))['analysis']
    lr = tuple(cfg['calls']['length_range'])
    feats = pd.read_csv(a.features, sep='\t').to_dict('records')
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for cell, (ext, pat, element, rev) in CELLS.items():
        ext = Path(ext)
        _, wt_i, vi, groups, weights, edges, _ = classify(ext, spec)
        meta = np.load(ext / 'meta.npz')
        dm = meta['distal_m6a'].astype(float)
        r, mid, length, _ = load_calls(ext, sorted(glob.glob(pat)), element, rev, cfg['calls']['min_llr'], lr)
        sel = (vi == wt_i) & (groups >= 0)
        dens, mol = density(r, mid, length, groups, sel, 2 + len(weights), lr)
        for g in range(dens.shape[0]):
            label = 'U' if g == 0 else 'M' if g == 1 else f'C{g - 1}'
            med = float(np.median(dm[sel & (groups == g)]))
            for f in feats:
                rows.append(dict(cell=cell, group=label, molecules=int(mol[g]), median_distal_m6a=med, feature=f['feature'],
                                 centre_oligo=f['mid_peak_oligo'], length_peak=f['len_peak'], occupancy=box_mass(dens[g], f, lr[0])))
        print(cell, 'done', flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(out / 'occupancy_by_stratum.tsv', sep='\t', index=False)
    w = t.pivot_table(index=['feature', 'centre_oligo', 'length_peak'], columns=['cell', 'group'], values='occupancy')
    pd.set_option('display.width', 300)
    print((100 * w).round(2).to_string())
    print(t.drop_duplicates(['cell', 'group'])[['cell', 'group', 'molecules', 'median_distal_m6a']].to_string(index=False))


if __name__ == '__main__':
    main()

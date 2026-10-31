"""Data files for the promoter explorer page (after Phase A; no new statistics).

base.json (loaded once)
  sequence, motif windows, published MPRA (Kircher 2019, HepG2) per SNV for both elements.
  HepG2 run 3: per-SNV summaries (A10/A11); WT unmarked fraction per adenine (labeling > 0.01) and call coverage
  per size bin (labeling >= 0.2), overall and per labeling stratum (width 0.05), so the page can show WT at
  any variant's labeling mix; a sample of WT molecules (footprints and m6A) for the pileup.
  HEK293T: the same WT tracks and pileup, A10 clone-mean motif results, v1 footprints K01-K23.
chunks/hepg2_NN.json (loaded on demand; 10 positions each)
  per SNV summaries also carry near-variant measures (+/-10 bp, +/-5 bp mask run) and the published MPRA p-values.
  per SNV: its own unmarked fraction and call coverage (as read, no mask), its labeling-stratum weights,
  a sample of its molecules, and the A10 test results per feature (E1 and E2, both masks).
Fractions and differences are stored as integers in units of 1e-4.

Usage: python tools/analysis_v2/explorer_data.py --root results/analysis_v2_20261028 --catalog <catalog.json> \
           --procs 4 --output <new dir>
"""
import argparse
import importlib.util
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('a10', Path(__file__).with_name('a10_effects.py'))
a10 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a10)
a9f, a7, a6x = a10.a9f, a10.a7, a10.a6x

WIDTH, OLIGO0 = 211, 21
SCALE = 1e4
NB = len(a10.BIN_NAMES)
PILEUP = 40
NEAR = 10                                    # bp around the variant for the near-variant lollipop measures
MPRA = {'act_p': ('LDLR', '/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/mpra_published/GRCh38_LDLR.tsv', '\t'),
        'act2_p': ('LDLR.2', '/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/mpra_published/GRCh38_LDLR.2.csv', ',')}
CHUNK = 10
SEED = 20261030
CELLS = {
    'HepG2': dict(extract='a9_statistics/extract/HepG2_run3', element=(1234, 1445), reverse=False,
                  emissions='results/llr_vs_fiberhmm_20261026/results/emissions.tsv.gz'),
    'HEK293T': dict(extract='a9_statistics/extract/HEK293T_run1', element=(3236, 3447), reverse=True,
                    emissions='results/hek293t_llr_20261028/results/emissions.tsv.gz'),
}
LIT = [('Repeat 1 (Sp1)', 89, 104), ('Repeat 2 (SRE-1)', 124, 139), ('Repeat 3 (Sp1)', 140, 155)]
V1 = [('ZNF652-like', 52, 62), ('SP4', 82, 102), ('SREBF1/2', 127, 136), ('SP1/SP3/KLF', 140, 156),
      ('CREB1/CREM', 182, 192)]


def _calls(args):
    return a10.a9f._calls(args)


def r(x, k=4):
    return None if x is None or not np.isfinite(x) else round(float(x), k)


def ints(a):
    a = np.asarray(a, float)
    return [None if not np.isfinite(v) else int(round(v * SCALE)) for v in a.ravel()]


class Cell:
    """One cell's molecules, element columns in oligo order, and footprint calls for labeling >= 0.2."""

    def __init__(self, root, cfg, procs):
        self.el_s, self.el_e = cfg['element']
        pos, X, distal, _, _ = a6x.load(root / cfg['extract'])
        meta = np.load(root / cfg['extract'] / 'meta.npz')
        self.label, self.cls = meta['label'], meta['cls']
        ecol = np.flatnonzero((pos >= self.el_s) & (pos < self.el_e))
        oligo = ((self.el_e - 1 - pos[ecol]) if cfg['reverse'] else (pos[ecol] - self.el_s)) + OLIGO0
        o = np.argsort(oligo)
        self.ecol, self.adenines = ecol[o], oligo[o]
        self.E = X[:, self.ecol]                                   # molecules x element adenines (-1/0/1)
        self.distal = distal
        self.st = a9f.strata(distal)
        self.k2 = np.flatnonzero(distal >= a9f.E2_MIN)
        hit, miss = a7.emission_steps(cfg['emissions'], pos)
        chunks = np.array_split(self.k2, max(1, procs * 16))
        with mp.get_context('fork').Pool(procs) as pool:
            bits = pool.map(_calls, [(X[c], hit, miss, pos, self.el_s, self.el_e, cfg['reverse']) for c in chunks])
        self.bits = np.concatenate(bits).reshape(self.k2.size, NB * WIDTH)
        self.row2 = np.full(distal.size, -1)
        self.row2[self.k2] = np.arange(self.k2.size)

    def unmarked(self, rows):
        x = self.E[rows]
        obs = (x >= 0).sum(0)
        return np.where(obs > 0, (x == 0).sum(0) / np.maximum(obs, 1), np.nan)

    def coverage(self, rows):
        b = self.bits[self.row2[rows]]
        return b.mean(0) if len(b) else np.full(NB * WIDTH, np.nan)

    def strata_e1(self, rows):
        out = np.full((a9f.N_STRATA, self.E.shape[1]), np.nan)
        for s in range(a9f.N_STRATA):
            m = rows[self.st[rows] == s]
            if m.size:
                out[s] = self.unmarked(m)
        return out

    def strata_e2(self, rows):
        out = np.full((a9f.N_STRATA, NB * WIDTH), np.nan)
        for s in range(a9f.N_STRATA):
            m = rows[self.st[rows] == s]
            if m.size:
                out[s] = self.bits[self.row2[m]].mean(0)
        return out

    def weights(self, rows):
        return np.bincount(self.st[rows], minlength=a9f.N_STRATA).tolist()

    def pileup(self, rows, rng):
        """Up to PILEUP molecules with labeling >= 0.2, sorted by labeling: footprints [start, end, bin] in oligo
        coordinates (runs of covered positions per size bin) and indices of marked adenines."""
        rows = rows[self.row2[rows] >= 0]
        if rows.size > PILEUP:
            rows = rng.choice(rows, PILEUP, replace=False)
        rows = rows[np.argsort(self.distal[rows])]
        out = []
        for m in rows:
            b = self.bits[self.row2[m]].reshape(NB, WIDTH)
            segs = []
            for k in range(NB):
                on = np.flatnonzero(b[k])
                if on.size:
                    br = np.flatnonzero(np.diff(on) > 1)
                    for s, e in zip(np.r_[on[0], on[br + 1]], np.r_[on[br], on[-1]]):
                        segs.append([int(s) + OLIGO0, int(e) + OLIGO0, k])
            out.append(dict(l=round(float(self.distal[m]), 2), s=segs, m=np.flatnonzero(self.E[m] == 1).tolist(),
                            u=np.flatnonzero(self.E[m] == 0).tolist()))
        return out


def wt_block(c, rng, strata):
    wt = np.flatnonzero(c.cls == 'WT')
    k1 = wt[c.distal[wt] > 0.01]
    k2 = wt[c.distal[wt] >= a9f.E2_MIN]
    out = dict(adenines=c.adenines.tolist(), unmarked=ints(c.unmarked(k1)),
               coverage={b: ints(v) for b, v in zip(a10.BIN_NAMES, c.coverage(k2).reshape(NB, WIDTH))},
               pileup=c.pileup(wt, rng), molecules=dict(wt=int(wt.size), unmarked=int(k1.size), calls=int(k2.size)))
    if strata:
        out['strata_unmarked'] = ints(c.strata_e1(k1))
        out['strata_coverage'] = ints(c.strata_e2(k2))
    return out


def near_measures(eff, v):
    """Variant - WT within NEAR bp of the variant (+/-5 bp mask run): mean change per size bin of call coverage and
    of the unmarked fraction (adenines 6-10 bp away), in points, and whether any feature there is a hit."""
    g = eff[(eff.snv == v) & ((eff.oligo - eff.snv_oligo).abs() <= NEAR) & eff.testable]
    out = {}
    for b in a10.BIN_NAMES:
        h = g[(g.endpoint == 'E2') & (g['bin'] == b)]
        out[f'near_{b}'] = r(100 * h.D.mean(), 3) if len(h) else None
        out[f'near_{b}_hit'] = bool(h.hit.any())
    h = g[g.endpoint == 'E1']
    out['near_unmarked'] = r(100 * h.D.mean(), 3) if len(h) else None
    out['near_unmarked_hit'] = bool(h.hit.any())
    return out


def mpra_p(catalog_snvs):
    out = {}
    for key, (_, path, sep) in MPRA.items():
        t = pd.read_csv(path, sep=sep).rename(columns={'Position': 'genomic_pos1', 'Ref': 'ref', 'Alt': 'alt'})
        m = catalog_snvs.merge(t[['genomic_pos1', 'ref', 'alt', 'P-Value']], on=['genomic_pos1', 'ref', 'alt'], how='left')
        out[key] = dict(zip(m.source_variant, m['P-Value']))
    return out


def summaries(root, catalog_snvs, eff_pm5):
    runs = {'A7': root / 'a10_effects/HepG2_run3', 'pm5': root / 'a10_effects/HepG2_run3_pm5'}
    tabs = {k: pd.read_csv(v / 'snvs.tsv', sep='\t') for k, v in runs.items()}
    emp = pd.read_csv(root / 'a10_effects/HepG2_run3_pm5_empirical_null/snvs.tsv', sep='\t')
    fn = pd.read_csv(root / 'a11_motifs_function/function_snvs.tsv', sep='\t').set_index('snv')
    bal = pd.read_csv(root / 'a9_statistics/HepG2_run3/balance.tsv', sep='\t').query('draw == -1').set_index('snv')
    mde = pd.read_parquet(root / 'a9_statistics/HepG2_run3/mde.parquet', columns=['snv', 'endpoint', 'mde'])
    mde1 = mde[mde.endpoint == 'E1'].groupby('snv').mde.median()
    cat = catalog_snvs.set_index('source_variant')
    pv = mpra_p(catalog_snvs)
    rows = []
    for v in sorted(cat.index, key=lambda s: (int(cat.loc[s, 'oligo_pos1']), s)):
        t = {k: tabs[k][tabs[k].snv == v].set_index('endpoint') for k in tabs}
        if 'E1' not in t['pm5'].index:
            continue
        e = emp[emp.snv == v].set_index('endpoint')
        row = dict(id=v, pos=int(cat.loc[v, 'oligo_pos1']), ref=cat.loc[v, 'ref'], alt=cat.loc[v, 'alt'],
                   caller=t['pm5'].loc['E1', 'caller_class'], clones=int(t['pm5'].loc['E1', 'clones']),
                   molecules=int(t['pm5'].loc['E1', 'molecules']),
                   act=r(fn.act.get(v), 3), act2=r(fn.act2.get(v), 3), shift=r(fn.C1.get(v), 3),
                   local=r(fn.C2.get(v), 3), local_rms=r(fn.rms.get(v), 3), region=bool(fn.region.get(v, False)),
                   mde_e1=r(100 * mde1.get(v, np.nan), 2), smd_labeling=r(bal.smd_distal.get(v), 3),
                   smd_passes=r(bal.smd_passes.get(v), 3), smd_strand=r(bal.smd_strand.get(v), 3))
        for k in ('A7', 'pm5'):
            for ep in ('E1', 'E2'):
                s = t[k].loc[ep] if ep in t[k].index else None
                tag = f'{k}_{ep}'
                row[f'sel_{tag}'] = bool(s.selected) if s is not None else None
                row[f'hits_{tag}'] = int(s.hits) if s is not None else None
                row[f'logp_{tag}'] = r(-np.log10(max(s.simes_p, 1e-300)), 2) if s is not None else None
        for ep in ('E1', 'E2'):
            row[f'sel_emp_{ep}'] = bool(e.loc[ep, 'selected']) if ep in e.index else None
        for k in pv:
            row[k] = r(pv[k].get(v), 5)
        row.update(near_measures(eff_pm5, v))
        rows.append(row)
    return rows


def tests(root, v):
    out = {}
    for k, run in (('A7', 'HepG2_run3'), ('pm5', 'HepG2_run3_pm5')):
        g = tests.eff[k]
        g = g[(g.snv == v) & (g.endpoint == 'E1')].sort_values('oligo')
        out[f'e1_{k}'] = dict(D=ints(g.D.where(g.testable)), se=ints(g.se.where(g.testable)),
                              hit=np.flatnonzero(g.hit.to_numpy()).tolist(),
                              masked=np.flatnonzero(g.masked.to_numpy()).tolist())
    for k in ('A7', 'pm5'):
        g = tests.eff[k]
        g = g[(g.snv == v) & (g.endpoint == 'E2')]
        out[f'e2_{k}'] = {}
        for b in a10.BIN_NAMES:
            h = g[g['bin'] == b].sort_values('oligo')
            out[f'e2_{k}'][b] = dict(D=ints(h.D.where(h.testable)), se=ints(h.se.where(h.testable)),
                                     hit=np.flatnonzero(h.hit.to_numpy()).tolist())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'chunks').mkdir()
    root = Path(a.root)
    rng = np.random.default_rng(SEED)
    lib = json.load(open(a.catalog))['library']
    cat = pd.DataFrame(lib['variants'])
    snvs = cat[cat['class'] == 'SNV'][['source_variant', 'oligo_pos1', 'ref', 'alt', 'genomic_pos1']].copy()
    snvs['oligo_pos1'] = snvs.oligo_pos1.astype(int)
    fn = pd.read_csv(root / 'a11_motifs_function/function_snvs.tsv', sep='\t').set_index('snv')
    base = dict(sequence=lib['mutable_sequence'], oligo0=OLIGO0, scale=SCALE, bins=a10.BIN_NAMES, chunk=CHUNK,
                n_strata=a9f.N_STRATA, motifs=dict(literature=LIT, v1=V1),
                activity=[dict(id=v, pos=int(p), ref=rf, alt=al, act=r(fn.act.get(v), 3), act2=r(fn.act2.get(v), 3))
                          for v, p, rf, al in zip(snvs.source_variant, snvs.oligo_pos1, snvs.ref, snvs.alt)])

    c = Cell(root, CELLS['HepG2'], a.procs)
    print('HepG2 calls done', flush=True)
    hep = wt_block(c, rng, strata=True)
    tests.eff = {k: pd.read_parquet(root / 'a10_effects' / run / 'effects.parquet')
                 for k, run in (('A7', 'HepG2_run3'), ('pm5', 'HepG2_run3_pm5'))}
    hep['snvs'] = summaries(root, snvs, tests.eff['pm5'])
    chunks = {}
    for s in hep['snvs']:
        v = s['id']
        rows = np.flatnonzero((c.cls == 'SNV') & (c.label == v))
        k1 = rows[c.distal[rows] > 0.01]
        k2 = rows[c.distal[rows] >= a9f.E2_MIN]
        d = dict(n1=int(k1.size), n2=int(k2.size), w1=c.weights(k1), w2=c.weights(k2),
                 unmarked=ints(c.unmarked(k1)), coverage=ints(c.coverage(k2)), pileup=c.pileup(rows, rng))
        d.update(tests(root, v))
        chunks.setdefault((s['pos'] - OLIGO0) // CHUNK, {})[v] = d
    for k, d in chunks.items():
        json.dump(d, open(out / 'chunks' / f'hepg2_{k:02d}.json', 'w'), separators=(',', ':'))
    base['HepG2'] = hep
    del c

    h = Cell(root, CELLS['HEK293T'], a.procs)
    hek = wt_block(h, rng, strata=False)
    m = pd.read_csv(root / 'a10_effects/HEK293T_run1/motifs.tsv', sep='\t')
    hek['motif_tests'] = [dict(motif=x.motif, start=int(x.start), end=int(x.end), group=x.group, endpoint=x.endpoint,
                               snvs=int(x.snvs), clones=int(x.cm_snv_clones), effect=r(100 * x.cm_effect, 3),
                               null_sd=r(100 * x.cm_null_sd, 3), p=r(x.cm_p, 5), bh=bool(x.cm_bh))
                          for x in m.itertuples()]
    fp = json.load(open('results/footprint_sizes_20261028/browser/data.json'))['cells']['HEK293T']['features']
    hek['footprints'] = [dict(id=f['id'], mid=f['mid'], len=f['len'], flag=f['flag']) for f in fp]
    base['HEK293T'] = hek
    json.dump(base, open(out / 'base.json', 'w'), separators=(',', ':'))
    sizes = {p.name: round(p.stat().st_size / 1e6, 2) for p in [out / 'base.json'] + sorted((out / 'chunks').iterdir())}
    print(json.dumps(dict(files=len(sizes), total_mb=round(sum(sizes.values()), 1), base_mb=sizes['base.json'],
                          largest_chunk_mb=max(v for k, v in sizes.items() if k != 'base.json'))))


if __name__ == '__main__':
    main()

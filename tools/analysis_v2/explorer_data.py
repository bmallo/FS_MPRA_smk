"""Data files for the promoter explorer page (after Phase A; no new statistics).

base.json (loaded once)
  sequence, motif windows, published MPRA (Kircher 2019, HepG2) per SNV for both elements.
  HepG2 run 3: per-SNV summaries (A10/A11, near-variant measures, published MPRA p-values); WT inaccessible
  (unmarked) fraction per adenine (labeling > 0.01), overall and per labeling stratum (width 0.05) so the page can
  show WT at any variant's labeling mix; WT footprints and per-variant footprint tests; WT pileup molecules.
  HEK293T: WT inaccessible fraction and pileup, A10 clone-mean motif results, v1 footprints K01-K23.
chunks/hepg2_NN.json (loaded on demand; 10 positions each)
  per SNV: its own inaccessible fraction, its labeling-stratum weights, a sample of its molecules, and the A10
  per-adenine test results (A7 and +/-5 bp masks).
Fractions and differences are stored as integers in units of 1e-4.

Confidence (display approximation, docs/footprint_confidence_design.md): each footprint call gets
P_chance, the probability that its unmarked adenines would all be unmarked with nothing bound, given the
molecule's labeling: log P = (d / d_ref) * sum log(1 - p_a) with p_a the run's fold-A exposed rates (labeling
>= 0.6 molecules, d_ref at labeling 0.70) and d = -ln(1 - labeling). Calls with P_chance < 1e-3 are 'confident';
pileup calls carry -log10 P_chance.

Footprints (HepG2): the WT footprint map and per-variant footprint occupancy tests from fp_features.py (--features);
pileup calls carry [start, end, kind, -log10 P_chance, footprint index], kind 0 = inside a WT footprint box,
1 = other 10-80 bp call, 2 = 81+ bp call (HEK293T: no footprint map, kind 1 or 2 only).

Usage: python tools/analysis_v2/explorer_data.py --root results/analysis_v2_20261028 --catalog <catalog.json> \
           --features results/analysis_v2_20261028/fp_features --procs 4 --output <new dir>
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
REF_LABELING = 0.70                          # typical labeling of the molecules the exposed rates come from
CONF_LOG10 = -3.0                            # confident call: P_chance < 1e-3
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


def dose(labeling):
    L = np.clip(labeling, 0.01, 0.99)
    return np.log(1 - L) / np.log(1 - REF_LABELING)


def segments(b):
    """Footprint calls of one molecule from its covered positions: (start, end, size bin), element index."""
    out = []
    for k in range(NB):
        on = np.flatnonzero(b[k])
        if on.size:
            br = np.flatnonzero(np.diff(on) > 1)
            out += [(int(s), int(e), k) for s, e in zip(np.r_[on[0], on[br + 1]], np.r_[on[br], on[-1]])]
    return out


def chance_log10(e_row, f, s, e, ad_idx, lp):
    """log10 P_chance of a call over element positions s..e: unmarked adenines there, scaled by the dose."""
    lo, hi = np.searchsorted(ad_idx, s, 'left'), np.searchsorted(ad_idx, e, 'right')
    seg = e_row[lo:hi]
    return float(f * np.where(seg == 0, lp[lo:hi], 0.0).sum() / np.log(10))




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
        t = pd.read_csv(cfg['emissions'], sep='\t')
        t = t[t.fold == 'A']
        pa = dict(zip(t.ref_pos - 1, t.p_a))
        self.ad_idx = self.adenines - OLIGO0
        self.lp = np.log(1 - np.clip([pa[int(p)] for p in pos[self.ecol]], 0, 0.9999))

    def unmarked(self, rows):
        x = self.E[rows]
        obs = (x >= 0).sum(0)
        return np.where(obs > 0, (x == 0).sum(0) / np.maximum(obs, 1), np.nan)


    def strata_e1(self, rows):
        out = np.full((a9f.N_STRATA, self.E.shape[1]), np.nan)
        for s in range(a9f.N_STRATA):
            m = rows[self.st[rows] == s]
            if m.size:
                out[s] = self.unmarked(m)
        return out


    def weights(self, rows):
        return np.bincount(self.st[rows], minlength=a9f.N_STRATA).tolist()

    def pileup(self, rows, rng):
        """Up to PILEUP molecules with labeling >= 0.2, sorted by labeling: footprints [start, end, bin, -log10 P_chance]
        in oligo coordinates (runs of covered positions per size bin) and indices of marked adenines."""
        rows = rows[self.row2[rows] >= 0]
        if rows.size > PILEUP:
            rows = rng.choice(rows, PILEUP, replace=False)
        rows = rows[np.argsort(self.distal[rows])]
        out = []
        for m in rows:
            b = self.bits[self.row2[m]].reshape(NB, WIDTH)
            f = dose(self.distal[m])
            segs = [[s + OLIGO0, e + OLIGO0, 2 if k == NB - 1 else 1,
                     round(-chance_log10(self.E[m], f, s, e, self.ad_idx, self.lp), 1), -1] for s, e, k in segments(b)]
            out.append(dict(l=round(float(self.distal[m]), 2), s=segs, m=np.flatnonzero(self.E[m] == 1).tolist(),
                            u=np.flatnonzero(self.E[m] == 0).tolist()))
        return out


def wt_block(c, rng, strata):
    wt = np.flatnonzero(c.cls == 'WT')
    k1 = wt[c.distal[wt] > 0.01]
    k2 = wt[c.distal[wt] >= a9f.E2_MIN]
    out = dict(adenines=c.adenines.tolist(), unmarked=ints(c.unmarked(k1)),
               pileup=c.pileup(wt, rng), molecules=dict(wt=int(wt.size), unmarked=int(k1.size), calls=int(k2.size)))
    if strata:
        out['strata_unmarked'] = ints(c.strata_e1(k1))
    return out


class FootprintCalls:
    """fp_features.py outputs: exact calls per molecule (labeling >= 0.2), WT footprints and per-variant tests."""

    def __init__(self, d):
        z = np.load(Path(d) / 'calls.npz')
        self.row = {int(r): i for i, r in enumerate(z['rows'])}
        self.mid, self.len, self.start, self.end, self.S = z['mid'], z['length'], z['start'], z['end'], z['S']
        self.fp = pd.read_csv(Path(d) / 'footprints.tsv', sep='\t')
        self.eff = pd.read_parquet(Path(d) / 'effects.parquet')

    def box(self, m, ln):
        for i, f in enumerate(self.fp.itertuples()):
            if f.mid_lo + OLIGO0 <= m <= f.mid_hi + OLIGO0 and f.len_lo <= ln <= f.len_hi:
                return i
        return -1

    def pileup(self, c, rows, rng):
        rows = rows[c.distal[rows] >= a9f.E2_MIN]
        if rows.size > PILEUP:
            rows = rng.choice(rows, PILEUP, replace=False)
        rows = rows[np.argsort(c.distal[rows])]
        out = []
        for m in rows:
            i, f, segs = self.row[int(m)], dose(c.distal[m]), []
            for k in np.flatnonzero(self.mid[i] >= 0):
                ln, fp = int(self.len[i, k]), self.box(int(self.mid[i, k]), int(self.len[i, k]))
                segs.append([max(int(self.start[i, k]), OLIGO0), min(int(self.end[i, k]), OLIGO0 + WIDTH - 1),
                             0 if fp >= 0 else (2 if ln >= 81 else 1), round(float(-f * self.S[i, k] / np.log(10)), 1), fp])
            out.append(dict(l=round(float(c.distal[m]), 2), s=segs, m=np.flatnonzero(c.E[m] == 1).tolist(),
                            u=np.flatnonzero(c.E[m] == 0).tolist()))
        return out

    def footprints(self):
        cols = ['feature', 'oligo_mid', 'len_peak', 'oligo_lo', 'oligo_hi', 'len_lo', 'len_hi', 'span_lo', 'span_hi',
                'wt_occupancy', 'wt_occupancy_confident', 'confident_share', 'c_over_m', 'flag',
                'wt_occ_0.2', 'wt_occ_0.4', 'wt_occ_0.6']
        return [{k: (r(v, 5) if isinstance(v, float) else (int(v) if isinstance(v, (np.integer,)) else v)) for k, v in row.items()}
                for row in self.fp[cols].to_dict('records')]

    def effects(self, snv_ids):
        """Per mask and call set: per SNV (in snv_ids order) x footprint, occupancy, D, SE (integers, 1e-4) and
        indices of significant footprints."""
        names = self.fp.feature.tolist()
        out = {}
        for (mask, calls), g in self.eff.groupby(['mask', 'calls']):
            g = g.set_index(['snv', 'footprint'])
            blk = dict(occ=[], D=[], se=[], sig=[])
            for v in snv_ids:
                h = g.loc[v].reindex(names) if v in g.index.get_level_values(0) else None
                if h is None:
                    for k in ('occ', 'D', 'se'):
                        blk[k].append([None] * len(names))
                    blk['sig'].append([])
                    continue
                blk['occ'].append(ints(h.occ.to_numpy()))
                blk['D'].append(ints(h.D.to_numpy()))
                blk['se'].append(ints(h.se.to_numpy()))
                blk['sig'].append(np.flatnonzero(h.sig.fillna(False).to_numpy(bool)).tolist())
            out[f'{mask}|{calls}'] = blk
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
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--features', required=True, help='fp_features.py output dir (HepG2 run 3)')
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'chunks').mkdir()
    fc = FootprintCalls(a.features)
    root = Path(a.root)
    rng = np.random.default_rng(SEED)
    lib = json.load(open(a.catalog))['library']
    cat = pd.DataFrame(lib['variants'])
    snvs = cat[cat['class'] == 'SNV'][['source_variant', 'oligo_pos1', 'ref', 'alt', 'genomic_pos1']].copy()
    snvs['oligo_pos1'] = snvs.oligo_pos1.astype(int)
    fn = pd.read_csv(root / 'a11_motifs_function/function_snvs.tsv', sep='\t').set_index('snv')
    base = dict(sequence=lib['mutable_sequence'], oligo0=OLIGO0, scale=SCALE, bins=a10.BIN_NAMES, chunk=CHUNK,
                confidence=dict(ref_labeling=REF_LABELING, log10_threshold=CONF_LOG10),
                n_strata=a9f.N_STRATA, motifs=dict(literature=LIT, v1=V1),
                activity=[dict(id=v, pos=int(p), ref=rf, alt=al, act=r(fn.act.get(v), 3), act2=r(fn.act2.get(v), 3))
                          for v, p, rf, al in zip(snvs.source_variant, snvs.oligo_pos1, snvs.ref, snvs.alt)])

    c = Cell(root, CELLS['HepG2'], a.procs)
    print('HepG2 calls done', flush=True)
    hep = wt_block(c, rng, strata=True)
    hep['pileup'] = fc.pileup(c, np.flatnonzero(c.cls == 'WT'), rng)
    tests.eff = {k: pd.read_parquet(root / 'a10_effects' / run / 'effects.parquet')
                 for k, run in (('A7', 'HepG2_run3'), ('pm5', 'HepG2_run3_pm5'))}
    hep['snvs'] = summaries(root, snvs, tests.eff['pm5'])
    hep['footprints'] = fc.footprints()
    hep['fp_effects'] = fc.effects([x['id'] for x in hep['snvs']])
    chunks = {}
    for s in hep['snvs']:
        v = s['id']
        rows = np.flatnonzero((c.cls == 'SNV') & (c.label == v))
        k1 = rows[c.distal[rows] > 0.01]
        k2 = rows[c.distal[rows] >= a9f.E2_MIN]
        d = dict(n1=int(k1.size), n2=int(k2.size), w1=c.weights(k1), unmarked=ints(c.unmarked(k1)),
                 pileup=fc.pileup(c, rows, rng))
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

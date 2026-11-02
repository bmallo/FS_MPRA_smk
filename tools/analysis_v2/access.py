"""Data access for follow-up analyses of analysis v2 (docs/analysis_v2/DATA_GUIDE.md). No new statistics.

load_element(cell)    per-molecule m6A calls at the element's 101 adenines, with labeling, clone and genotype,
                      from the A9 extract (the analysis-ready WT + SNV molecules used by A9-A11).
a7_masks(cell)        per SNV, the oligo positions of its A7 primary mask (pm5=True adds every A/T within +/-5 bp).
E1Reference(d)        WT inaccessible (unmarked) fraction per adenine and labeling stratum, with the WT clone
                      variance (A9 random-effects test).
ref.test(d, snv)      A10 E1 test of one SNV against WT: per adenine D (variant - WT, labeling-standardised),
                      random-effects se, p (NaN = not testable or masked), z.

Run on the standing node or in Slurm, with PYTHONPATH=. and BLAS threads = 1. Loading HepG2 run 3 takes ~1 min and
~0.5 GB (element columns only; the extract itself is 255 MB packed, ~3.6 GB unpacked).

Example:
    from tools.analysis_v2 import access
    d = access.load_element('HepG2')
    ref = access.E1Reference(d)
    t = ref.test(d, '79C>A', mask=access.a7_masks('HepG2')['79C>A'])
"""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('a10', Path(__file__).with_name('a10_effects.py'))
a10 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a10)
a9f = a10.a9f

ROOT = 'results/analysis_v2_20261028'
CELLS = {
    'HepG2': dict(extract='a9_statistics/extract/HepG2_run3', element=(1234, 1445), reverse=False,
                  catalog='results/rep3_ledger_20261016/catalog/catalog.json',
                  variants='results/analysis_v2_20261028/a4_informative/visibility/variants_visibility.tsv',
                  emissions='results/llr_vs_fiberhmm_20261026/results/emissions.tsv.gz'),
    'HEK293T': dict(extract='a9_statistics/extract/HEK293T_run1', element=(3236, 3447), reverse=True,
                    catalog=None, variants=None,
                    emissions='results/hek293t_llr_20261028/results/emissions.tsv.gz'),
}
CHUNK = 100_000


def oligo_of(pos, element, reverse):
    """Oligo coordinate (Kircher 2019 numbering; the element spans oligo 21-231) of plasmid positions."""
    el_s, el_e = element
    return ((el_e - 1 - pos) if reverse else (pos - el_s)) + a10.OLIGO0


def load_element(cell, root=ROOT):
    """SimpleNamespace: oligo (101 element adenines, ascending), X (molecules x 101 int8: -1 not observed,
    0 unmarked = inaccessible, 1 marked = accessible), distal (labeling: m6A fraction on the backbone), passes,
    strand, label (genotype: 'WT' or an SNV such as '79C>A'), cls ('WT'/'SNV'), clone, read_quality."""
    c = CELLS[cell]
    ext = Path(root) / c['extract']
    z = np.load(ext / 'molecules.npz')
    pos = z['pos']
    el_s, el_e = c['element']
    ecol = np.flatnonzero((pos >= el_s) & (pos < el_e))
    obs, meth = z['observed'], z['methylated']
    X = np.empty((obs.shape[0], ecol.size), np.int8)
    for i in range(0, obs.shape[0], CHUNK):
        o = np.unpackbits(obs[i:i + CHUNK], axis=1, count=pos.size)[:, ecol].astype(bool)
        m = np.unpackbits(meth[i:i + CHUNK], axis=1, count=pos.size)[:, ecol].astype(np.int8)
        X[i:i + CHUNK] = np.where(o, m, np.int8(-1))
    oligo = oligo_of(pos[ecol], c['element'], c['reverse'])
    order = np.argsort(oligo)
    meta = np.load(ext / 'meta.npz')
    return SimpleNamespace(cell=cell, oligo=oligo[order], pos=pos[ecol][order], X=X[:, order],
                           distal=z['distal_frac'], passes=z['passes'], strand=z['strand'], label=meta['label'],
                           cls=meta['cls'], clone=meta['clone'], read_quality=meta['read_quality'])


def a7_masks(cell='HepG2', pm5=False, root=ROOT):
    """{snv: array of masked oligo positions} (A7 primary mask; pm5 adds +/-5 bp). HepG2 only (HEK293T has no
    per-SNV catalog in this form; its tests are motif-level)."""
    c = CELLS[cell]
    if c['catalog'] is None:
        raise ValueError(f'no per-SNV masks for {cell}')
    pos = np.load(Path(root) / c['extract'] / 'molecules.npz')['pos']
    m = a10.snv_masks(c['catalog'], c['variants'], pos, *c['element'], c['reverse'], pm5)
    return {v: np.sort(oligo_of(pos[r['mask']], c['element'], c['reverse'])) for v, r in m.items()}


class E1Reference:
    """WT reference for E1 (molecules with labeling > 0.01), as in A10: per labeling stratum (width 0.05) and adenine,
    the WT inaccessible fraction, and the between-clone variance tau2 from WT clones."""

    def __init__(self, d):
        wt = np.flatnonzero(d.cls == 'WT')
        Y, O = self._yo(d, wt)
        self.ref = a9f.Reference(Y, O, a9f.strata(d.distal[wt]), d.clone[wt])
        self.p, self.So = self.ref.rates()

    @staticmethod
    def _yo(d, rows):
        X = d.X[rows]
        return X == 0, (X >= 0) & (d.distal[rows] > a9f.E1_MIN)[:, None]

    def wt_fraction(self, weights=None):
        """WT inaccessible fraction per adenine, standardised to the given stratum weights (default: WT's own)."""
        w = self.So.sum(1) if weights is None else np.asarray(weights, float)
        w = w / w.sum()
        return np.nansum(self.p * w[:, None], 0)

    def test(self, d, snv, mask=()):
        """Per adenine: n molecules, C clones, D (variant - WT inaccessible fraction), se, p, z, testable, masked."""
        rows = np.flatnonzero(d.label == snv)
        if rows.size == 0:
            raise KeyError(snv)
        Y, O = self._yo(d, rows)
        r = a9f.contrast(Y, O, a9f.strata(d.distal[rows]), d.clone[rows], self.p, self.So, self.ref)
        masked = np.isin(d.oligo, np.asarray(mask))
        t = pd.DataFrame({'oligo': d.oligo, **{k: r[k] for k in ('n', 'C', 'D', 'se', 'p', 'expected')}})
        t['p'] = np.where(masked, np.nan, t.p)
        t['z'] = t.D / t.se
        t['masked'], t['testable'] = masked, t.p.notna()
        return t

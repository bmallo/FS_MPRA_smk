import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location('a10', Path(__file__).parents[1] / 'tools/analysis_v2/a10_effects.py')
a10 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a10)
ROOT = Path(__file__).parents[1]
CATALOG = ROOT / 'results/rep3_ledger_20261016/catalog/catalog.json'
VARIANTS = ROOT / 'results/analysis_v2_20261028/a4_informative/visibility/variants_visibility.tsv'


def test_simes_and_bh():
    assert np.isclose(a10.simes(np.array([0.01, 0.04, 0.5, np.nan])), 0.03)
    assert a10.bh(np.array([0.001, 0.02, 0.5, np.nan]), 0.05).tolist() == [True, True, False, False]


def test_benjamini_bogomolov_scales_q_by_selected_share():
    fam = {'a': 0.001, 'b': 0.002, 'c': 0.9, 'd': 0.8}
    feat = {'a': np.array([0.001, 0.02]), 'b': np.array([0.002, 0.3]), 'c': np.array([0.9]), 'd': np.array([0.8])}
    chosen, rej = a10.bb(fam, feat, 0.05)
    assert chosen == {'a', 'b'}
    # within families BH at 0.05 * 2 / 4 = 0.025: a keeps 0.001 and 0.02 (0.02 <= 0.025), b keeps only 0.002
    assert rej['a'].tolist() == [True, True] and rej['b'].tolist() == [True, False] and not rej['c'].any()


def _layout(reverse):
    lib = json.load(open(CATALOG))['library']
    wt = lib['mutable_sequence']
    el_s, el_e = (3236, 3447) if reverse else (1234, 1445)
    comp = str.maketrans('ACGT', 'TGCA')
    ref_el = wt.translate(comp)[::-1] if reverse else wt
    pos = np.array([el_s + j for j, b in enumerate(ref_el) if b in 'AT'], np.int64)
    return pos, el_s, el_e


def test_masks_point_at_the_edited_adenine_in_both_orientations():
    var = pd.read_csv(VARIANTS, sep='\t')
    v = var[(var.cls == 'SNV') & (var.caller_drops_at.astype(str) == 'True')].iloc[0]
    for reverse in (False, True):
        pos, el_s, el_e = _layout(reverse)
        m = a10.snv_masks(CATALOG, VARIANTS, pos, el_s, el_e, reverse)[v.source_variant]
        assert m['group'] == 'drops its base' and m['oligo'] == v.oligo_pos1
        j = int(v.index0)
        ref = (el_e - 1 - j) if reverse else (el_s + j)
        assert ref in set(pos[m['mask']].tolist())


def test_masked_calls_equal_calls_on_masked_molecules():
    rng = np.random.default_rng(1)
    pos = np.arange(0, 600, 3).astype(np.int64)
    X = np.where(rng.random((30, pos.size)) < 0.6, 1, 0).astype(np.int8)
    X[:, 100:110] = 0
    hit = np.full(pos.size, np.log(0.01 / 0.7))
    miss = np.full(pos.size, np.log(0.99 / 0.3))
    cols = np.array([104, 105], np.int64)
    got = a10.masked_calls(X, cols, hit, miss, pos, 200, 411, False)
    Xm = X.copy()
    Xm[:, cols] = -1
    exp = a10.a8r.call_bits(Xm, hit, miss, pos, 200, 411, False, a10.B_LO, a10.B_HI).reshape(30, -1)
    assert np.array_equal(got, exp) and got.any()
    assert X[0, 104] == 0                                               # the input is not modified


def test_molecule_scores():
    Y = np.array([[1, 0], [0, 0]], bool)
    O = np.array([[True, True], [False, False]])
    p = np.array([[0.5, 0.5]])
    s = a10.molecule_scores(Y, O, p, np.array([0, 0]))
    assert np.isclose(s[0], 0.0) and np.isnan(s[1])


def test_clone_mean_test_centres_on_wt_clones():
    rng = np.random.default_rng(2)
    cw = np.repeat(np.arange(50), 120)
    wt = rng.normal(0, 0.01, cw.size) + np.repeat(rng.normal(0.003, 0.002, 50), 120)   # WT clones off zero
    snv_scores = [np.full(150, 0.003) + rng.normal(0, 0.01, 150)]
    r = a10.clone_mean_test(['v'], snv_scores, [np.zeros(150, int)], {'v': wt}, cw, rng)
    assert r['cm_snv_clones'] == 1 and r['cm_wt_clones'] == 50
    assert abs(r['cm_null_mean'] - wt.mean()) < 0.001 and r['cm_p'] > 0.05
    far = [np.full(150, 0.05)]
    assert a10.clone_mean_test(['v'], far, [np.zeros(150, int)], {'v': wt}, cw, rng)['cm_p'] < 0.01
    assert np.isnan(a10.clone_mean_test(['v'], [np.full(50, 0.1)], [np.zeros(50, int)], {'v': wt}, cw, rng)['cm_p'])


def test_pm5_mask_covers_the_window():
    var = pd.read_csv(VARIANTS, sep='\t')
    v = var[(var.cls == 'SNV')].iloc[40]
    pos, el_s, el_e = _layout(False)
    prim = a10.snv_masks(CATALOG, VARIANTS, pos, el_s, el_e, False)[v.source_variant]['mask']
    pm = a10.snv_masks(CATALOG, VARIANTS, pos, el_s, el_e, False, pm5=True)[v.source_variant]['mask']
    j = int(v.index0)
    near = {int(r) for r in pos if abs(int(r) - (el_s + j)) <= 5}
    assert set(prim) <= set(pm) and set(pos[pm].tolist()) >= near

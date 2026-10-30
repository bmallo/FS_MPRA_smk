import importlib.util
import sys
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location('a7', Path(__file__).parents[1] / 'tools/analysis_v2/a7_snv_edits.py')
a7 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a7)
kadane = sys.modules['kadane']

B_LO = np.array([1, 10], np.int64)
B_HI = np.array([9, 100000], np.int64)
NOF = tuple(np.zeros(0, np.int64) for _ in range(4))


def synthetic(seed=3, n=60, L=400):
    rng = np.random.default_rng(seed)
    pos = np.sort(rng.choice(L, size=150, replace=False)).astype(np.int64)
    X = np.where(rng.random((n, pos.size)) < 0.5, 1, 0).astype(np.int8)
    X[:, 60:80] = 0                                     # a protected stretch on every molecule
    X[rng.random(X.shape) < 0.05] = -1
    hit = np.log(rng.uniform(0.01, 0.05, pos.size) / rng.uniform(0.5, 0.8, pos.size))
    miss = np.log((1 - rng.uniform(0.01, 0.05, pos.size)) / (1 - rng.uniform(0.5, 0.8, pos.size)))
    return pos, X, hit, miss, L


def dense_coverage(pos, X, hit, miss, L, el_s, el_e, reverse=False):
    """Reference implementation: call.py's dense reference-coordinate scan, then per-bin coverage."""
    cov = np.zeros((2, el_e - el_s), np.int64)
    for x in X:
        step = np.zeros(L)
        ih = np.zeros(L, bool)
        ii = np.zeros(L, bool)
        o = x >= 0
        ii[pos[o]] = True
        ih[pos[x == 1]] = True
        step[pos[x == 1]] = hit[x == 1]
        step[pos[x == 0]] = miss[x == 0]
        rp = pos[o]
        s, e, *_ = kadane.scan(step, ih, ii, int(rp[0]), int(rp[-1]) + 1, a7.MIN_LLR, a7.MIN_OPPS, a7.DRAWDOWN)
        for a, b in zip(s, e):
            q = 0 if b - a <= 9 else 1
            for r in range(max(a, el_s), min(b, el_e)):
                cov[q, (el_e - 1 - r) if reverse else (r - el_s)] += 1
    return cov


def test_scan_cov_matches_dense_reference_scan():
    pos, X, hit, miss, L = synthetic()
    for rev in (False, True):
        cov, _ = a7.scan_cov(X, hit, miss, pos, np.zeros(0, np.int64), np.zeros((X.shape[0], 0), np.int8),
                             100, 300, rev, B_LO, B_HI, *NOF)
        assert np.array_equal(cov, dense_coverage(pos, X, hit, miss, L, 100, 300, rev))
        assert cov.sum() > 0


def test_override_equals_editing_the_molecules():
    pos, X, hit, miss, L = synthetic()
    cols = np.array([65, 70], np.int64)
    vals = np.full((X.shape[0], 2), -1, np.int8)
    vals[:, 1] = 1
    cov, _ = a7.scan_cov(X, hit, miss, pos, cols, vals, 100, 300, False, B_LO, B_HI, *NOF)
    Y = X.copy()
    Y[:, cols] = vals
    ref, _ = a7.scan_cov(Y, hit, miss, pos, np.zeros(0, np.int64), np.zeros((X.shape[0], 0), np.int8),
                         100, 300, False, B_LO, B_HI, *NOF)
    assert np.array_equal(cov, ref)
    base, _ = a7.scan_cov(X, hit, miss, pos, np.zeros(0, np.int64), np.zeros((X.shape[0], 0), np.int8),
                          100, 300, False, B_LO, B_HI, *NOF)
    assert not np.array_equal(cov, base)


def test_feature_counts_molecules_once():
    pos = np.arange(0, 100, 2).astype(np.int64)        # 50 A/T, every other base
    X = np.ones((3, 50), np.int8)
    X[:2, 20:30] = 0                                    # molecules 0-1: one call over ref 40-58
    hit = np.full(50, -3.0)
    miss = np.full(50, 1.0)
    f = (np.array([8, 0], np.int64), np.array([10, 2], np.int64), np.array([15, 15], np.int64),
         np.array([25, 25], np.int64))
    _, feat = a7.scan_cov(X, hit, miss, pos, np.zeros(0, np.int64), np.zeros((3, 0), np.int8), 40, 70, False,
                          B_LO, B_HI, *f)
    assert feat.tolist() == [2, 0]                      # call [40, 59): length 19, midpoint index 9


def test_recontexted_and_masks():
    wt = 'GCAAAAGCTG'
    # A->G inside the run: the run AAAA (2-5) splits into AA (2-3) and A (5)
    rc = a7.recontexted(wt, 'GCAAGAGCTG', 'SNV', 4)
    assert [k for k, _, _ in rc] == [2, 3, 5]
    assert rc[-1][1] == (4, "3' end") and rc[-1][2] == (1, 'isolated')
    # G/C -> G/C far from any run: nothing changes
    assert a7.recontexted(wt, 'GCAAAAGGTG', 'SNV', 7) == []
    # deleting one A of the run: remaining adenines shift class
    assert [k for k, _, _ in a7.recontexted(wt, 'GCAAAGCTG', 'DEL', 2)] == [3, 4, 5]
    assert a7.pm_set(wt, 0) == [2, 3, 4, 5]
    assert a7.pm_set(wt, 9) == [4, 5, 8]


def test_classes():
    assert a7.classes(np.array([0.0, 0.01, 0.3, 0.6, 0.9])).tolist() == ['U', 'U', 'C', 'M', 'M']

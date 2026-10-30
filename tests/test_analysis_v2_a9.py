import importlib.util
from pathlib import Path

from types import SimpleNamespace

import numpy as np

spec = importlib.util.spec_from_file_location('a9f', Path(__file__).parents[1] / 'tools/analysis_v2/a9_framework.py')
a9f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a9f)
spec = importlib.util.spec_from_file_location('a9x', Path(__file__).parents[1] / 'tools/analysis_v2/a9_extract.py')
a9x = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a9x)


def wt_rates(F=3, p=0.3, tau2=0.0):
    S = a9f.N_STRATA
    ref = SimpleNamespace(V=np.zeros((S, F)), W2=np.zeros((S, F)), tau2=np.full(F, tau2))
    return np.full((S, F), p), np.full((S, F), 1e6), ref


def test_clone_robust_variance_matches_direct_sandwich():
    rng = np.random.default_rng(1)
    n, F = 60, 3
    Y = rng.random((n, F)) < 0.4
    O = rng.random((n, F)) < 0.9
    st = rng.integers(0, 4, n)
    clone = rng.integers(0, 12, n)
    p, So, V = wt_rates(F)
    r = a9f.contrast(Y, O, st, clone, p, So, V)
    for f in range(F):
        use = O[:, f]
        res = Y[use, f] - 0.3
        D = res.mean()
        cs = np.array([(res[clone[use] == c] - D).sum() for c in np.unique(clone[use])])
        C = cs.size
        assert np.isclose(r['D'][f], D)
        assert np.isclose(r['se_sandwich'][f] ** 2, C / (C - 1) * (cs ** 2).sum() / use.sum() ** 2)
        assert r['C'][f] == C


def test_one_molecule_per_clone_is_residual_variance():
    rng = np.random.default_rng(2)
    n = 40
    Y = rng.random((n, 1)) < 0.5
    O = np.ones((n, 1), bool)
    p, So, V = wt_rates(1, 0.5)
    r = a9f.contrast(Y, O, np.zeros(n, int), np.arange(n), p, So, V)
    res = Y[:, 0] - 0.5
    assert np.isclose(r['se_sandwich'][0] ** 2, n / (n - 1) * ((res - res.mean()) ** 2).sum() / n ** 2)


def test_standardisation_removes_a_labeling_difference():
    """WT rate depends on labeling stratum only; a 'variant' with another labeling mix but the same per-stratum
    rates has D near 0, while the raw difference in means is large."""
    rng = np.random.default_rng(3)
    rate = np.linspace(0.9, 0.1, a9f.N_STRATA)
    st_wt = rng.integers(0, a9f.N_STRATA, 200000)
    Yw = (rng.random(st_wt.size) < rate[st_wt])[:, None]
    ref = a9f.Reference(Yw, np.ones_like(Yw), st_wt, np.arange(st_wt.size))
    st_v = rng.integers(15, 20, 5000)                                   # heavily labeled only
    Yv = (rng.random(st_v.size) < rate[st_v])[:, None]
    p, So = ref.rates()
    r = a9f.contrast(Yv, np.ones_like(Yv), st_v, np.arange(st_v.size), p, So, ref)
    assert abs(Yv.mean() - Yw.mean()) > 0.2
    assert abs(r['D'][0]) < 3 * r['se'][0]


def test_reference_rates_remove_given_molecules():
    rng = np.random.default_rng(4)
    Y = rng.random((100, 2)) < 0.5
    O = np.ones_like(Y)
    st = rng.integers(0, 3, 100)
    ref = a9f.Reference(Y, O, st, np.arange(100))
    p, So = ref.rates(Y[:10], O[:10], st[:10])
    full = a9f.Reference(Y[10:], O[10:], st[10:], np.arange(90))
    p2, So2 = full.rates()
    assert np.allclose(np.nan_to_num(p), np.nan_to_num(p2)) and np.allclose(So, So2)


def test_draw_clones_matches_sizes_without_reuse():
    rng = np.random.default_rng(5)
    rows, start = [], 0
    pools = {}
    for size, count in ((1, 50), (2, 30), (5, 10), (9, 3)):
        for _ in range(count):
            pools.setdefault(size, []).append(np.arange(start, start + size))
            start += size
    sizes = np.array([1, 1, 1, 2, 3, 5, 9])
    rr, lab, short = a9f.draw_clones(sizes, pools, rng)
    got = np.sort(np.bincount(lab))
    assert np.array_equal(got, np.sort(sizes)) and short == 0
    assert np.unique(rr).size == rr.size                                # no molecule (so no clone) twice
    rr, lab, short = a9f.draw_clones(np.array([12]), pools, rng)        # bigger than any WT clone
    assert short == 1 and rr.size == 9


def test_clone_keys_drop_contradicting_reads():
    reads = np.array([f'r{i}' for i in range(7)])
    geno = np.array(['A', 'A', 'A', 'A', 'B', 'B', 'A'])
    bc = np.array(['x'] * 5 + ['y', ''])
    q = np.ones(7, bool)
    key, drop = a9x.clone_keys(reads, geno, bc, q)
    # barcode x: 4 A, 1 B -> share 0.8 < 0.9, nothing dropped
    assert not drop.any() and key[0] == 'x|A' and key[4] == 'x|B' and key[6].startswith('#')
    geno2 = np.array(['A'] * 10 + ['B'])
    key, drop = a9x.clone_keys(np.array([f'r{i}' for i in range(11)]), geno2, np.array(['x'] * 11), np.ones(11, bool))
    assert drop.tolist() == [False] * 10 + [True]                       # 10/11 = 0.91 >= 0.9


def test_bh_any():
    assert a9f.bh_any(np.array([0.001, 0.5, 0.9]))
    assert not a9f.bh_any(np.array([0.04, 0.5, 0.9]))


def test_random_effects_variance_adds_tau2_times_clone_weights():
    rng = np.random.default_rng(6)
    n, F = 50, 2
    Y = rng.random((n, F)) < 0.3
    O = np.ones((n, F), bool)
    clone = np.repeat(np.arange(5), 10)
    p, So, ref0 = wt_rates(F, 0.3, 0.0)
    _, _, ref1 = wt_rates(F, 0.3, 0.01)
    r0 = a9f.contrast(Y, O, np.zeros(n, int), clone, p, So, ref0)
    r1 = a9f.contrast(Y, O, np.zeros(n, int), clone, p, So, ref1)
    assert np.allclose(r0['se'], r0['se_naive'])                        # tau2 = 0: random effects = binomial
    assert np.allclose(r1['se'] ** 2 - r0['se'] ** 2, 0.01 * 5 * (10 / 50) ** 2)


def test_tau2_recovers_simulated_between_clone_variance():
    rng = np.random.default_rng(7)
    C, k, tau = 3000, 8, 0.08
    mu = np.clip(0.4 + rng.normal(0, tau, C), 0.01, 0.99)
    clone = np.repeat(np.arange(C), k)
    Y = (rng.random(C * k) < mu[clone])[:, None]
    O = np.ones_like(Y)
    pm = np.full(Y.shape, Y.mean())
    assert abs(a9f.tau2(Y, O, pm, clone)[0] - tau ** 2) < 0.25 * tau ** 2
    Y0 = (rng.random(C * k) < 0.4)[:, None]
    # no clone effect: the untruncated estimator is unbiased with SD ~ sqrt(2 / C) * p(1-p) / k = 0.00075 here
    assert a9f.tau2(Y0, O, np.full(Y0.shape, Y0.mean()), clone)[0] < 3 * 0.00075


def test_features_with_few_expected_minority_events_are_not_tested():
    n = 400
    Y = np.ones((n, 2), bool)
    Y[:3, 0] = False
    O = np.ones((n, 2), bool)
    p, So, ref = wt_rates(2, 0.999)                                     # expected minority count 0.4
    r = a9f.contrast(Y, O, np.zeros(n, int), np.arange(n), p, So, ref)
    assert np.isnan(r['p']).all() and np.isclose(r['expected'][0], 0.4)
    p2, So2, ref2 = wt_rates(2, 0.9)
    assert np.isfinite(a9f.contrast(Y, O, np.zeros(n, int), np.arange(n), p2, So2, ref2)['p']).all()

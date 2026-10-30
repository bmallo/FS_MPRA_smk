import importlib.util
from pathlib import Path

import numpy as np
import pysam

spec = importlib.util.spec_from_file_location('a8s', Path(__file__).parents[1] / 'tools/analysis_v2/a8_spikes.py')
a8s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a8s)
spec = importlib.util.spec_from_file_location('a8i', Path(__file__).parents[1] / 'tools/analysis_v2/a8_ids.py')
a8i = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a8i)

B_LO = np.array([10, 20, 41, 81], np.int64)
B_HI = np.array([19, 40, 80, 100000], np.int64)


def molecules(seed=5, n=80):
    rng = np.random.default_rng(seed)
    pos = np.arange(0, 600, 3).astype(np.int64)          # an A/T every 3 bp
    X = np.where(rng.random((n, pos.size)) < 0.7, 1, 0).astype(np.int8)
    X[rng.random(X.shape) < 0.03] = -1
    hit = np.full(pos.size, np.log(0.01 / 0.7))
    miss = np.full(pos.size, np.log(0.99 / 0.3))
    return pos, X, hit, miss


def run(pos, X, hit, miss, s0, L, spiked, keep=0.0, u=None):
    c_lo, c_hi = int(np.searchsorted(pos, s0)), int(np.searchsorted(pos, s0 + L))
    if u is None:
        u = np.ones((X.shape[0], max(c_hi - c_lo, 1)))
    centre = s0 + (L - 1) // 2
    return a8s.spike_scan(X, hit, miss, pos, np.ones(pos.size, np.bool_), c_lo, c_hi, keep, u, s0, s0 + L, centre,
                          L, 200, 411, False, B_LO, B_HI, spiked)


def test_unspiked_scan_matches_stored_segments():
    pos, X, hit, miss = molecules()
    ss, ee, off = a8s.base_segments(X, hit, miss, pos, np.ones(pos.size, np.bool_))
    for s0, L in [(300, 30), (250, 147), (330, 15)]:
        d, k, _, _ = run(pos, X, hit, miss, s0, L, False)
        sb = a8s.size_bin(L)
        d2, k2 = a8s.detect(ss, ee, off, s0, s0 + L, s0 + (L - 1) // 2, L, B_LO[sb], B_HI[sb])
        assert np.array_equal(d, d2) and np.array_equal(k, k2)


def test_perfect_spike_is_found_with_the_right_size():
    pos, X, hit, miss = molecules()
    X[:, :] = np.where(X >= 0, 1, X)                     # fully marked molecules: no footprint anywhere
    d0, _, _, _ = run(pos, X, hit, miss, 300, 30, False)
    d1, k1, l1, _ = run(pos, X, hit, miss, 300, 30, True)
    assert not d0.any() and d1.all() and k1.all()
    # calls run from the first to the last observed unmarked adenine (one every 3 bp): 28 bp, or 25 if an edge
    # adenine is unobserved
    assert set(l1.tolist()) <= {25, 28} and (l1 == 28).mean() > 0.8


def test_leaky_spike_with_keep_one_changes_nothing():
    pos, X, hit, miss = molecules()
    d0, k0, _, c0 = run(pos, X, hit, miss, 300, 40, False)
    d1, k1, _, c1 = run(pos, X, hit, miss, 300, 40, True, keep=1.0, u=np.zeros((X.shape[0], 14)))
    assert np.array_equal(d0, d1) and np.array_equal(c0, c1)


def test_size_bin():
    assert [a8s.size_bin(x) for x in (10, 15, 20, 40, 41, 80, 81, 147, 9)] == [0, 0, 1, 1, 2, 2, 3, 3, -1]


def _read(cigar, seq, start=0):
    h = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'c', 'LN': 100}]})
    r = pysam.AlignedSegment(h)
    r.query_name, r.reference_id, r.reference_start = 'r', 0, start
    r.cigarstring, r.query_sequence = cigar, seq
    return r


def test_barcode_reader():
    seq = 'AAAAACCCCCGGGGGTTTTT'
    assert a8i.barcode(_read('20M', seq), 5, 10) == 'CCCCC'
    assert a8i.barcode(_read('7M1I12M', 'AAAAACCTCCCGGGGGTTTT'), 5, 10) == ''      # insertion inside
    assert a8i.barcode(_read('7M2D11M', 'AAAAACCCCCGGGGGTTT'), 5, 10) == ''        # deletion inside
    assert a8i.barcode(_read('10M1I9M', 'AAAAACCCCCAGGGGGTTTT'), 5, 10) == 'CCCCC'  # insertion right after


spec = importlib.util.spec_from_file_location('a8r', Path(__file__).parents[1] / 'tools/analysis_v2/a8_repro.py')
a8r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a8r)


def test_clone_ids_group_barcodes_and_split_unreadable():
    bc = np.array(['AAA', 'AAA', '', '', 'CCC'])
    crc = np.array([1, 2, 3, 3, 4], np.uint32)
    ids, half = a8r.clone_ids(bc, crc)
    assert ids[0] == ids[1] and len({ids[1], ids[2], ids[4]}) == 3
    assert ids[2] == ids[3]                              # same name CRC -> same key (one molecule in practice)
    assert half[0] == half[1] and half[2] == 3 % 2


def test_standardised_weights_strata():
    rng = np.random.default_rng(0)
    strata = np.array([0] * 50 + [1] * 50)
    num = np.concatenate([np.zeros((50, 1)), np.ones((50, 1))])
    w = np.zeros(a8r.N_STRATA)
    w[0], w[1] = 0.25, 0.75
    est, boot = a8r.standardised(num, np.ones(100), strata, np.arange(100), w, 20, rng)
    assert np.isclose(est[0], 0.75) and boot.shape == (20, 1)
    # per-feature denominators: a feature observed only in stratum 1 molecules
    den = np.concatenate([np.zeros((50, 1)), np.ones((50, 1))])
    w2 = np.zeros(a8r.N_STRATA)
    w2[1] = 1.0
    est2, _ = a8r.standardised(num, den, strata, np.arange(100), w2, 5, rng)
    assert np.isclose(est2[0], 1.0)


def test_emissions_use_fold_a_m_and_u_with_jeffreys():
    X = np.array([[1, 0], [1, -1], [0, 0], [1, 1]], np.int8)
    distal = np.array([0.7, 0.7, 0.0, 0.7])            # M, M, U, M
    fold = np.array([0, 0, 0, 1])                       # last M molecule is held out
    p_a, p_p, n_m, n_u = a8r.emissions(X, distal, fold)
    assert np.allclose(p_a, [(2 + 0.5) / 3, (0 + 0.5) / 2]) and np.allclose(p_p, [0.5 / 2, 0.5 / 2])
    assert (n_m, n_u) == (2, 1)

import array
import importlib.util
from pathlib import Path

import numpy as np
import pysam

spec = importlib.util.spec_from_file_location('a5x', Path(__file__).parents[1] / 'tools/analysis_v2/a5_m6a_extract.py')
a5x = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a5x)
spec = importlib.util.spec_from_file_location('a5s', Path(__file__).parents[1] / 'tools/analysis_v2/a5_m6a_summary.py')
a5s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a5s)


def test_bins_classes_and_distal_mask():
    assert [a5x.pass_bin(p) for p in (3, 9, 10, 14, 15, 29, 30, 80)] == [0, 0, 1, 1, 2, 3, 4, 4]
    assert [a5x.klass(f) for f in (0.0, 0.01, 0.0101, 0.59, 0.6, 0.9)] == [0, 0, 1, 1, 2, 2]
    m = a5x.distal_mask(1000, (400, 500), (800, 810))
    assert not m[300:600].any() and m[299] and m[600] and not m[775:835].any() and m[774] and m[835]


def test_wt_exact():
    ref = np.frombuffer(b'ACGTACGTAC', np.uint8)
    q = np.arange(10)
    r = np.arange(10)
    quals = np.full(10, 30)
    assert a5x.wt_exact(q, r, ref.copy(), quals, ref, (2, 8))
    mism = ref.copy()
    mism[4] = ord('G')
    assert not a5x.wt_exact(q, r, mism, quals, ref, (2, 8))
    assert a5x.wt_exact(q, r, mism, quals, ref, (5, 8))          # mismatch outside the region
    low = quals.copy()
    low[3] = 10
    assert not a5x.wt_exact(q, r, ref.copy(), low, ref, (2, 8))
    # deletion of reference base 4: r skips 4
    rd = np.array([0, 1, 2, 3, 5, 6, 7, 8, 9])
    assert not a5x.wt_exact(np.arange(9), rd, ref.copy(), quals, ref, (2, 8))
    # insertion after query 4: q jumps inside the region
    qi = np.array([0, 1, 2, 3, 4, 6, 7, 8, 9, 10])
    assert not a5x.wt_exact(qi, r, np.frombuffer(b'ACGTAAACGTAC', np.uint8)[:11], np.full(11, 30), ref, (2, 8))


def _record(seq, mm, ml, ref_start=0, reverse=False, np_=20, rq=0.999):
    h = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'p', 'LN': len(seq)}]})
    r = pysam.AlignedSegment(h)
    r.query_name = 'r1'
    r.query_sequence = seq
    r.reference_id = 0
    r.reference_start = ref_start
    r.cigartuples = [(0, len(seq))]
    r.mapping_quality = 60
    r.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
    r.is_reverse = reverse
    r.set_tag('MM', mm)
    r.set_tag('ML', array.array('B', ml))
    r.set_tag('np', np_)
    r.set_tag('rq', rq)
    return r


def test_ml_array_and_process():
    seq = 'GATCAAGTTC' * 100                                     # 1,000 bp, 600 A/T (>= MIN_DISTAL_AT)
    ref_arr = np.frombuffer(seq.encode(), np.uint8)
    at = (ref_arr == 65) | (ref_arr == 84)
    # every A listed (ML 250 for the first A, 10 for the rest), T channel present with no calls
    n_a = seq.count('A')
    ml = [250] + [10] * (n_a - 1)
    rec = _record(seq, 'A+a.,' + ','.join(['0'] * n_a) + ';T-a.;', ml)
    arr = a5x.ml_array(rec)
    assert arr[1] == 250 and arr[4] == 10 and arr[0] == -1 and (arr >= 0).sum() == n_a
    st, x = a5x.process(rec, ref_arr, at, len(seq), np.ones(len(seq), bool), (100, 110))
    assert st == 'ok' and x['wt'] and x['passes'] == 20 and x['strand'] == 0
    assert x['pos'].size == at.sum()
    assert abs(x['distal_frac'] - 1 / at.sum()) < 1e-12              # one call >= 128 out of all A/T
    st, _ = a5x.process(_record(seq, 'A+a.;', [], rq=0.98), ref_arr, at, len(seq), np.ones(len(seq), bool), (100, 110))
    assert st == 'read_quality'
    # m6A channel present with no calls = an unlabeled molecule, kept with distal fraction 0
    st, x = a5x.process(_record(seq, 'A+a.;C+m?;T-a.;', []), ref_arr, at, len(seq), np.ones(len(seq), bool), (100, 110))
    assert st == 'ok' and x['distal_frac'] == 0 and (x['ml'] == -1).all()
    # no m6A channel at all (only 5mC) = not callable
    st, _ = a5x.process(_record(seq, 'C+m?;', []), ref_arr, at, len(seq), np.ones(len(seq), bool), (100, 110))
    assert st == 'no_m6a_tags'


def test_gatc_sites_and_controls():
    ref = 'CCGATCCGATACCATCGG'
    top, bottom = a5s.gatc_sites(ref)
    assert top == [3] and bottom == [4] and ref[3] == 'A' and ref[4] == 'T'
    ct, cb = a5s.same_context_controls(ref)
    assert ct == [8]                                              # A of GAT in GATA (index 8), not the GATC
    assert cb == [14]                                             # T of ATC in CATC (bottom-strand adenine)
    assert a5s.revcomp('AACG') == 'CGTT'


def test_rate_pools_counts():
    c = np.zeros((3, 6), int)
    c[0, [1, 2, 4]] = [10, 30, 0]
    c[1, [1, 2]] = [1, 9]
    pooled, per, n = a5s.rate(c, [1, 2, 4])
    assert abs(pooled - 0.25) < 1e-12 and n == 40
    assert per[0] == 0.1 and per[1] == 0.3 and np.isnan(per[2])

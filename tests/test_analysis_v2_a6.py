import array
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

spec = importlib.util.spec_from_file_location('a6x', Path(__file__).parents[1] / 'tools/analysis_v2/a6_extract.py')
a6x = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a6x)
spec = importlib.util.spec_from_file_location('a6p', Path(__file__).parents[1] / 'tools/analysis_v2/a6_proxies.py')
a6p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a6p)


def test_span_cover():
    pos = np.array([0, 5, 10, 30, 35, 60, 90])
    X = np.array([[1, 0, 0, 1, 0, 0, 1],       # spans: (0,30) 29 bp; (30,90) 59 bp
                  [0, 0, 1, 0, -1, 0, 1],      # first two unbounded; (10,90) 79 bp, pos 35 not observed
                  [1, 0, 0, 0, 0, 0, 0]], np.int8)  # nothing bounded on the right
    c20 = a6p.span_cover(X, pos, 20)
    assert c20[0].tolist() == [False, True, True, False, True, True, False]
    assert c20[1].tolist() == [False, False, False, True, False, True, False]
    assert not c20[2].any()
    c60 = a6p.span_cover(X, pos, 60)
    assert not c60[0].any() and c60[1].tolist() == [False, False, False, True, False, True, False]
    assert not a6p.span_cover(X, pos, 80).any()


def test_permute_within_keeps_columns_per_group():
    rng = np.random.default_rng(1)
    X = rng.integers(-1, 2, size=(200, 7)).astype(np.int8)
    g = np.repeat(['a', 'b'], 100)
    Y = a6p.permute_within(X, g, rng)
    for k in ['a', 'b']:
        for j in range(7):
            assert sorted(X[g == k, j]) == sorted(Y[g == k, j])
    assert not np.array_equal(X, Y)


def test_clustering_null_detects_linkage():
    # independent calls: excess ~ 0; whole-molecule blocks of unmethylation: excess > 0
    rng = np.random.default_rng(2)
    pos = np.arange(0, 400, 4)
    X = (rng.random((3000, pos.size)) < 0.7).astype(np.int8)
    distal = np.full(3000, 0.8)
    ind = a6p.clustering(pos, X, distal, None, 3)
    Xb = X.copy()
    Xb[:1000, 40:50] = 0                             # a 40-bp footprint on a third of the molecules
    blk = a6p.clustering(pos, Xb, distal, None, 3)
    sel = lambda t: t[(t.cls == 'M') & (t.span == 20) & t.ref_pos0.between(160, 199)].excess.mean()
    assert abs(sel(ind)) < 0.01 and sel(blk) > 0.1
    assert set(blk.cls) == {'M'}


def test_strata_and_classes():
    d = np.array([0.0, 0.01, 0.0101, 0.049, 0.05, 0.6, 0.999, 1.0])
    assert a6p.classes(d).tolist() == ['U', 'U', 'C', 'C', 'C', 'M', 'M', 'M']
    assert a6p.strata(d).tolist() == [0, 0, 0, 0, 1, 11, 19, 19] or a6p.strata(d).tolist() == [0, 0, 0, 0, 1, 12, 19, 19]


def test_tiles_align_to_element():
    t = a6p.tiles(4647, (1234, 1445))
    assert (1234, 1445) in t and all(e - s == 211 for s, e in t) and t[0][0] == 1234 % 211 and t[-1][1] <= 4647
    assert a6p.tiles(500, None) == [(0, 211), (211, 422)]


def test_llr_weight_shrinks_when_m_is_bound(tmp_path):
    p = tmp_path / 'e.tsv'
    pd.DataFrame(dict(fold=['A', 'A', 'B'], ref_pos=[1235, 2000, 1236], strand='top', context='x',
                      p_a=[0.8, 0.8, 0.5], p_p=[0.0, 0.0, 0.0])).to_csv(p, sep='\t', index=False)
    w = a6p.llr_weight(p, (1234, 1445))
    assert len(w) == 1 and w.ref_pos0.iloc[0] == 1234
    assert abs(w.weight.iloc[0] - np.log(1 / 0.2)) < 1e-12
    assert abs(w['weight_bound_0.5'].iloc[0] - np.log(1 / 0.6)) < 1e-12


def _record(seq, mm, ml, np_=20):
    h = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'p', 'LN': len(seq)}]})
    r = pysam.AlignedSegment(h)
    r.query_name = 'r1'
    r.query_sequence = seq
    r.reference_id = 0
    r.reference_start = 0
    r.cigartuples = [(0, len(seq))]
    r.mapping_quality = 60
    r.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
    r.set_tag('MM', mm)
    r.set_tag('ML', array.array('B', ml))
    r.set_tag('np', np_)
    r.set_tag('rq', 0.999)
    return r


def test_molecule_vector():
    seq = 'GATCAAGTTC' * 100
    ref_arr = np.frombuffer(seq.encode(), np.uint8)
    at_idx = np.flatnonzero((ref_arr == 65) | (ref_arr == 84))
    n_a = seq.count('A')
    rec = _record(seq, 'A+a.,' + ','.join(['0'] * n_a) + ';T-a.;', [250] + [10] * (n_a - 1))
    st, x = a6x.molecule(rec, ref_arr, at_idx, len(seq), np.ones(len(seq), bool), (100, 110))
    assert st == 'ok' and x['vec'][0] == 1 and (x['vec'] == 1).sum() == 1 and (x['vec'] == 0).sum() == at_idx.size - 1
    assert a6x.molecule(_record(seq, 'A+a.;', [], np_=9), ref_arr, at_idx, len(seq),
                        np.ones(len(seq), bool), None)[0] == 'low_passes'
    mut = seq[:105] + ('G' if seq[105] != 'G' else 'C') + seq[106:]
    assert a6x.molecule(_record(mut, 'A+a.;', []), ref_arr, at_idx, len(seq), np.ones(len(seq), bool),
                        (100, 110))[0] == 'not_wt'
    st, x = a6x.molecule(_record(mut, 'A+a.;', []), ref_arr, at_idx, len(seq), np.ones(len(seq), bool), None)
    assert st == 'ok'

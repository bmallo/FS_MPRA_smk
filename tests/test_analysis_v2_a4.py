import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location('a4', Path(__file__).parents[1] / 'tools/analysis_v2/a4_informative_bases.py')
a4 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a4)


def brute_cover(w, need):
    n = len(w)
    out = np.full(n, np.nan)
    for s in range(n):
        for e in range(s, n):
            if sum(w[s:e + 1]) >= need - 1e-12:
                for i in range(s, e + 1):
                    if np.isnan(out[i]) or out[i] > e - s + 1:
                        out[i] = e - s + 1
    return out


def test_shortest_cover_matches_brute_force():
    rng = np.random.default_rng(0)
    for _ in range(50):
        w = rng.random(30) * (rng.random(30) < 0.4)
        need = rng.uniform(0.3, 3)
        np.testing.assert_array_equal(a4.shortest_cover(w, need), brute_cover(list(w), need))


def test_shortest_cover_counts_and_gaps():
    at = np.array([1, 0, 0, 0, 1, 1, 0, 1], bool)
    # 3 A/T needed: [0,5] covers 0-2; position 3 via [3,7] (5 bp); positions 4-7 via [4,7]
    np.testing.assert_array_equal(a4.shortest_cover(at, 3), [6, 6, 6, 5, 4, 4, 4, 4])
    np.testing.assert_array_equal(a4.gap_length(at), [0, 3, 3, 3, 0, 0, 1, 0])
    np.testing.assert_array_equal(a4.window_sum(at, 1), [1, 1, 0, 1, 2, 2, 2, 1])


def test_seq_class():
    assert a4.seq_class('SNV', 'A', 'T') == 'A<->T swap'
    assert a4.seq_class('SNV', 'A', 'G') == 'A/T loss'
    assert a4.seq_class('SNV', 'C', 'T') == 'A/T gain'
    assert a4.seq_class('SNV', 'C', 'G') == 'G<->C swap'
    assert a4.seq_class('DEL', 'T', '') == 'DEL of A/T'
    assert a4.seq_class('INS', '', 'G') == 'INS of G/C'


def test_variants_table_caller_view():
    wt = 'ACGTAC' * 10
    pos = pd.DataFrame(dict(is_at=[b in 'AT' for b in wt], llr_miss=[1.0 if b in 'AT' else 0.0 for b in wt],
                            min_fp_seq_bp=np.ones(len(wt)), min_fp_llr_bp=np.ones(len(wt))))
    o = a4.OLIGO_OFFSET
    vs = [dict(variant_id='v1', source_variant='swap', **{'class': 'SNV'}, ref='A', alt='T', oligo_pos1=o + 4,
               sequence=wt[:4] + 'T' + wt[5:]),
          dict(variant_id='v2', source_variant='gain', **{'class': 'SNV'}, ref='C', alt='A', oligo_pos1=o + 1,
               sequence=wt[:1] + 'A' + wt[2:]),
          dict(variant_id='v3', source_variant='del', **{'class': 'DEL'}, ref='T', alt='', oligo_pos1=o + 3,
               sequence=wt[:3] + wt[4:]),
          dict(variant_id='v4', source_variant='ins', **{'class': 'INS'}, ref='', alt='A', oligo_pos1=o + 2,
               sequence=wt[:2] + 'A' + wt[2:])]
    t = a4.variants_table(vs, wt, pos).set_index('source_variant')
    assert t.loc['swap', 'caller_drops_at'] and t.loc['swap', 'at_count_change'] == 0
    assert not t.loc['gain', 'caller_drops_at'] and t.loc['gain', 'at_count_change'] == 1
    assert t.loc['del', 'caller_drops_at'] and t.loc['del', 'at_count_change'] == -1
    assert not t.loc['ins', 'caller_drops_at'] and t.loc['ins', 'at_count_change'] == 1


spec_c = importlib.util.spec_from_file_location('a4c', Path(__file__).parents[1] / 'tools/analysis_v2/a4_context_bias.py')
a4c = importlib.util.module_from_spec(spec_c)
spec_c.loader.exec_module(a4c)


def test_run_context():
    assert a4c.run_context('GCAAAAAGC', 'top') == (5, 3)
    assert a4c.run_context('GCGAAAAAG', 'top') == (5, 2)      # centre is the 2nd A from the 5' (left) end
    assert a4c.run_context('GCGCTTTCC', 'bottom') == (3, 3)      # leftmost T is the bottom strand's 3' end
    assert a4c.run_context('GCGCTTCCC', 'bottom') == (2, 2)   # bottom strand reads right to left
    assert a4c.run_context('GCGCACGCG', 'top') == (1, 1)


spec_v = importlib.util.spec_from_file_location('a4v', Path(__file__).parents[1] / 'tools/analysis_v2/a4_variant_visibility.py')
a4v = importlib.util.module_from_spec(spec_v)
spec_v.loader.exec_module(a4v)


def test_adenine_classes_and_compare():
    c = a4v.adenine_classes('GAAAAAGTTC')
    assert c[('top', 1)] == (5, "5' end") and c[('top', 3)] == (5, 'interior') and c[('top', 5)] == (5, "3' end")
    assert c[('bottom', 7)] == (2, "3' end") and c[('bottom', 8)] == (2, "5' end")
    lut = {(5, "5' end"): 0.2, (5, 'interior'): 0.02, (5, "3' end"): 0.7, (2, "5' end"): 0.8, (2, "3' end"): 0.9,
           (1, 'isolated'): 1.0, (3, "5' end"): 0.2, (3, 'interior'): 0.05, (3, "3' end"): 0.75}
    wt = 'GAAAAAGTTC'
    # A>G in the middle of the run: two isolated... runs AA + AA -> classes (2, ends)
    snv = 'GAAGAAGTTC'
    vw, vv, ch = a4v.compare(wt, snv, 'SNV', 3, lut)
    assert abs(vw - (0.2 + 3 * 0.02 + 0.7 + 0.8 + 0.9)) < 1e-9
    assert abs(vv - (0.8 + 0.9 + 0.8 + 0.9 + 0.8 + 0.9)) < 1e-9
    assert ch == 4
    # deleting the G between runs is mapped back to WT coordinates without spurious changes far away
    _, _, ch_del = a4v.compare('GAAGCCCTTC', 'GAACCCTTC', 'DEL', 3, {**lut, (2, "5' end"): 0.8})
    assert ch_del == 0

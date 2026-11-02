import importlib.util
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location('fpf', Path(__file__).parents[1] / 'tools/analysis_v2/fp_features.py')
fpf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fpf)


def molecule():
    pos = np.arange(1200, 1480, 3).astype(np.int64)            # an A/T every 3 bp around the element
    x = np.ones((1, pos.size), np.int8)
    prot = (pos >= 1300) & (pos < 1330)                         # protected stretch
    x[0, prot] = 0
    hit = np.full(pos.size, np.log(0.01 / 0.7))
    miss = np.full(pos.size, np.log(0.99 / 0.3))
    lp = np.full(pos.size, np.log(1 - 0.7))
    return pos, x, hit, miss, lp, prot


def test_call_table_reports_midpoint_length_and_chance_sum():
    pos, x, hit, miss, lp, prot = molecule()
    mid, ln, st, en, S = fpf.call_table(x, hit, miss, lp, pos, 1234, 1445, 8, 5.0, 3, 25.0)
    a, b = pos[prot][0], pos[prot][-1] + 1
    assert ln[0, 0] == b - a and mid[0, 0] == (a + b - 1) // 2 - 1234 + 21
    assert st[0, 0] == a - 1234 + 21 and en[0, 0] == b - 1 - 1234 + 21
    assert np.isclose(S[0, 0], prot.sum() * np.log(0.3)) and mid[0, 1] == -1


def test_occupancy_flags_all_and_confident():
    pos, x, hit, miss, lp, prot = molecule()
    tab = fpf.call_table(x, hit, miss, lp, pos, 1234, 1445, 8, 5.0, 3, 25.0)
    m, L = tab[0][0, 0] - 21, tab[1][0, 0]
    feats = [dict(mid_lo=m - 1, mid_hi=m + 1, len_lo=L - 2, len_hi=L + 2), dict(mid_lo=0, mid_hi=5, len_lo=10, len_hi=20)]
    occ = fpf.occupancy(tab, np.array([0.7]), feats)
    assert occ[0, 0] == 3 and occ[0, 1] == 0                   # in box and confident (0.3^10 << 1e-3)
    assert fpf.occupancy(tab, np.array([0.05]), feats)[0, 0] == 1   # weakly labeled: same call, not confident

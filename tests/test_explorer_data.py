import importlib.util
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location('exd', Path(__file__).parents[1] / 'tools/analysis_v2/explorer_data.py')
exd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exd)


def test_segments_split_runs_per_size_bin():
    b = np.zeros((exd.NB, exd.WIDTH), bool)
    b[1, 10:30] = True
    b[1, 40:45] = True
    b[3, 0:100] = True
    assert exd.segments(b) == [(10, 29, 1), (40, 44, 1), (0, 99, 3)]


def test_chance_counts_only_unmarked_adenines_in_the_call():
    ad_idx = np.array([5, 10, 15, 20, 25])
    p = np.array([0.9, 0.5, 0.9, 0.9, 0.9])
    lp = np.log(1 - p)
    e_row = np.array([0, 0, 1, 0, 0], np.int8)           # adenine at 15 marked
    q = exd.chance_log10(e_row, 1.0, 8, 22, ad_idx, lp)  # adenines 10, 15, 20 inside; 10 and 20 unmarked
    assert np.isclose(q, np.log10(0.5 * 0.1))


def test_dose_scales_marking_with_labeling():
    assert np.isclose(exd.dose(exd.REF_LABELING), 1.0)
    assert exd.dose(0.3) < 1 < exd.dose(0.9)
    # at dose f, P(unmarked) = (1 - p_ref)^f: a weak dose leaves more adenines unmarked by chance
    p_ref = 0.9
    assert (1 - p_ref) ** exd.dose(0.3) > (1 - p_ref)

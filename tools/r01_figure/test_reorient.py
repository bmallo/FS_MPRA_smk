"""Tests for tools/r01_figure/reorient.py (run: PYTHONPATH=. pytest -q tools/r01_figure)."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

_spec = importlib.util.spec_from_file_location('reorient', Path(__file__).with_name('reorient.py'))
ro = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ro)


def test_span_distance():
    d = ro.span_distance(np.array([50, 60, 70, 86, 90]), np.array([60] * 5), np.array([86] * 5))
    assert d.tolist() == [10, 0, 0, 0, 4]


@pytest.mark.parametrize('lo,hi,cls', [(10, 13, 'sub-TF (10-19)'), (24, 30, 'TF (~20-40)'), (17, 23, 'TF (~20-40)'),
                                       (35, 45, 'TF (~20-40)'), (69, 75, 'PIC-size (41-80)'), (82, 88, 'long (81+)')])
def test_size_class(lo, hi, cls):
    assert ro.size_class(SimpleNamespace(len_lo=lo, len_hi=hi)) == cls


@pytest.mark.skipif(not ro.MPRA_DIR.exists(), reason='published MPRA files not on this machine')
def test_activity_matches_function_table():
    """Activity parsed here equals A11's act/act2 for every designed SNV (same source, same mapping)."""
    import pandas as pd
    fn = pd.read_csv('results/analysis_v2_20261028/a11_motifs_function/function_snvs.tsv', sep='\t').set_index('snv')
    act = ro.activity().reindex(fn.index)
    assert act.act.notna().all() and act.act2.notna().all()
    assert np.allclose(act.act, fn.act) and np.allclose(act.act2, fn.act2)


_spec2 = importlib.util.spec_from_file_location('clinical', Path(__file__).with_name('clinical.py'))
cl = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(cl)
REF = 'ACGTTTACGT' + 'A' * 201          # toy element, oligo 21-231


def test_spdi_trims_and_maps():
    # oligo 25 = chr19:11089289 (0-based 11089288); SNV T>G
    assert cl.from_spdi('NC_000019.10:11089288:T:G', REF) == (25, 25, 'G')
    # repeat-expanded deletion of one T in TTT (oligo 24-26) trims to a 1-bp deletion
    s, e, i = cl.from_spdi('NC_000019.10:11089287:TTT:TT', REF)
    assert cl.kind(s, e, i) == '1-bp deletion'
    assert cl.from_spdi('NC_000019.10:11089000:A:G', REF) is None


def test_hgvs_g_forms():
    assert cl.from_hgvs_g('chr19:g.11089289T>G', REF) == (25, 25, 'G')
    assert cl.from_hgvs_g('chr19:g.11089289del', REF) == (25, 25, '')
    assert cl.from_hgvs_g('chr19:g.11089289dup', REF) == (26, 25, 'T')
    assert cl.from_hgvs_g('chr19:g.11089289C>G', REF) is None          # wrong reference base


def test_indel_designs_match_by_sequence():
    """A deletion anywhere in a homopolymer gives the same sequence as the design's deletion."""
    assert cl.edit(REF, 24, 24, '') == cl.edit(REF, 26, 26, '')

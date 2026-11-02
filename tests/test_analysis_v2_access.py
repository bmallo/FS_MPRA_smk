import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('access', ROOT / 'tools/analysis_v2/access.py')
access = importlib.util.module_from_spec(spec)
spec.loader.exec_module(access)
RESULTS = ROOT / access.ROOT
EFFECTS = RESULTS / 'a10_effects/HepG2_run3/effects.parquet'


def test_oligo_coordinates():
    assert access.oligo_of(np.array([1234, 1444]), (1234, 1445), False).tolist() == [21, 231]
    assert access.oligo_of(np.array([3446, 3236]), (3236, 3447), True).tolist() == [21, 231]


@pytest.mark.skipif(not EFFECTS.exists(), reason='analysis v2 results not present')
def test_reproduces_a10_e1(monkeypatch):
    monkeypatch.chdir(ROOT)
    d = access.load_element('HepG2')
    assert d.X.shape[1] == 101 and np.all(np.diff(d.oligo) > 0)
    masks = access.a7_masks('HepG2')
    assert masks['61A>T'].tolist() == list(range(59, 66))
    ref = access.E1Reference(d)
    eff = pd.read_parquet(EFFECTS, filters=[('endpoint', '==', 'E1'), ('snv', 'in', ['79C>A', '146C>T'])])
    for v in ('79C>A', '146C>T'):
        t = ref.test(d, v, mask=masks[v])
        e = eff[eff.snv == v].sort_values('oligo')
        assert np.allclose(t.D, e.D, equal_nan=True) and np.allclose(t.se, e.se, equal_nan=True)
        assert (t.testable.to_numpy() == e.testable.to_numpy()).all()

"""Tests for tools/r01_figure/repeat3.py pooling (run: PYTHONPATH=. pytest -q tools/r01_figure)."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

_spec = importlib.util.spec_from_file_location('repeat3', Path(__file__).with_name('repeat3.py'))
r3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(r3)


def test_dl_homogeneous_equals_fixed_effect():
    r = r3.dersimonian_laird([1.0, 1.0, 1.0], [0.5, 0.5, 0.5])
    assert r['tau2'] == 0 and r['est'] == pytest.approx(1.0) and r['se'] == pytest.approx(0.5 / np.sqrt(3))


def test_dl_heterogeneous_matches_hand_calculation():
    est, se = np.array([0.0, 2.0]), np.array([0.5, 0.5])
    w = 1 / se ** 2
    Q = (w * (est - 1) ** 2).sum()                       # fixed mean 1
    tau2 = (Q - 1) / (w.sum() - (w ** 2).sum() / w.sum())
    r = r3.dersimonian_laird(est, se)
    assert r['tau2'] == pytest.approx(tau2) and r['est'] == pytest.approx(1.0)
    assert r['se'] == pytest.approx(np.sqrt((se[0] ** 2 + tau2) / 2))


def test_dl_drops_missing():
    assert r3.dersimonian_laird([1.0, np.nan], [0.1, 0.1])['k'] == 1
    assert r3.dersimonian_laird([], [])['k'] == 0

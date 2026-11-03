"""Tests for tools/r01_figure/gnomad_check.py (run: PYTHONPATH=. pytest -q tools/r01_figure)."""
import importlib.util
import json
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location('gnomad_check', Path(__file__).with_name('gnomad_check.py'))
gc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gc)


def test_gnomad_table_sums_datasets_and_skips_indels(tmp_path):
    v = [dict(variant_id='19-11089335-C-T', pos=11089335, ref='C', alt='T', rsids=['rs1'],
              genome=dict(ac=2, an=100, af=0.02, homozygote_count=1, filters=[]),
              exome=dict(ac=1, an=300, af=1 / 300, homozygote_count=0, filters=['AS_VQSR'])),
         dict(variant_id='19-11089336-CT-C', pos=11089336, ref='CT', alt='C', rsids=[], genome=None,
              exome=dict(ac=5, an=10, af=0.5, homozygote_count=0, filters=[]))]
    p = tmp_path / 'g.json'
    p.write_text(json.dumps({'data': {'region': {'variants': v}}}))
    t = gc.gnomad_table(p)
    assert len(t) == 1
    r = t.iloc[0]
    assert r['name'] == '71C>T' and r.ac == 3 and r.an == 400 and r.af == pytest.approx(3 / 400)
    assert r.hom == 1 and bool(r.filtered)


def test_wilson_bounds():
    lo, hi = gc.wilson(0, 10)
    assert lo == pytest.approx(0, abs=1e-12) and hi == pytest.approx(0.2775, abs=1e-3)
    lo, hi = gc.wilson(5, 10)
    assert lo < 0.5 < hi

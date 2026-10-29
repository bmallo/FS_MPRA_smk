import importlib.util
from pathlib import Path

import pandas as pd

spec = importlib.util.spec_from_file_location('a3', Path(__file__).parents[1] / 'tools/analysis_v2/a3_genotype_audit.py')
a3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a3)


def test_relation():
    cls = dict(wt='WT', s1='SNV', s1b='SNV', s2='SNV', d='DEL')
    eq = dict(wt=frozenset(), s1=frozenset({30}), s1b=frozenset({30}), s2=frozenset({50}), d=frozenset({30, 31}))
    assert a3.relation('wt', 's1', cls, eq) == 'one edit apart'
    assert a3.relation('s1', 's1b', cls, eq) == 'one edit apart'
    assert a3.relation('s1', 'd', cls, eq) == 'one edit apart'
    assert a3.relation('s1', 's2', cls, eq) == 'different positions'


def test_barcode_table_majority():
    q = pd.DataFrame(dict(barcode=['A'] * 4 + ['B'] * 2, observed_variant_id=['x', 'x', 'x', 'y', 'z', 'z']))
    t = a3.barcode_table(q)
    assert t.loc['A', 'majority'] == 'x' and t.loc['A', 'reads'] == 4 and t.loc['A', 'share'] == 0.75 and t.loc['A', 'genotypes'] == 2
    assert t.loc['B', 'share'] == 1.0 and t.loc['B', 'genotypes'] == 1

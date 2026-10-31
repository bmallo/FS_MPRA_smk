import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location('a11', Path(__file__).parents[1] / 'tools/analysis_v2/a11_motifs_function.py')
a11 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a11)


def lod_for(consensus):
    """log2-odds matrix scoring +2 for the consensus base and -1 otherwise."""
    lod = np.full((len(consensus), 4), -1.0)
    for i, b in enumerate(consensus):
        lod[i, 'ACGT'.index(b)] = 2.0
    return lod


def test_best_site_finds_both_strands():
    wt = 'T' * 30 + 'GGGCGG' + 'T' * 30
    lod = lod_for('GGGCGG')
    sc, start, strand = a11.best_site(lod, wt, 21 + 28, 21 + 40)
    assert (sc, start, strand) == (12.0, 21 + 30, '+')
    rc = wt.translate(a11.COMP)[::-1]
    sc, start, strand = a11.best_site(lod, rc, 21 + 28, 21 + 40)
    assert sc == 12.0 and strand == '-'


def test_delta_pwm_of_a_consensus_break():
    wt = 'T' * 30 + 'GGGCGG' + 'T' * 30
    lod = lod_for('GGGCGG')
    sc, start, strand = a11.best_site(lod, wt, 21 + 28, 21 + 40)
    k = 32
    var = wt[:k] + 'A' + wt[k + 1:]
    assert a11.site_score(lod, var, start, strand) - sc == -3.0


def test_bh():
    q = a11.bh(np.array([0.01, 0.04, np.nan, 0.03]))
    assert np.allclose(q[[0, 1, 3]], [0.03, 0.04, 0.04]) and np.isnan(q[2])


def test_activity_join(tmp_path):
    p = tmp_path / 'a.tsv'
    pd.DataFrame(dict(Position=[10, 10], Ref=['C', 'C'], Alt=['A', 'G'], Value=[0.5, -0.2])).to_csv(p, sep='\t', index=False)
    snvs = pd.DataFrame(dict(source_variant=['x', 'y', 'z'], genomic_pos1=[10, 10, 11], ref=['C', 'C', 'A'], alt=['G', 'A', 'T']))
    act = a11.activity(p, '\t', snvs)
    assert act['x'] == -0.2 and act['y'] == 0.5 and np.isnan(act['z'])

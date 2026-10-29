import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location('a2', Path(__file__).parents[1] / 'tools/analysis_v2/a2_read_accounting.py')
a2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a2)

WT = 'ACGTTTGCAGCTAGCATCGATCGGATCCAGT'
AL = a2.aligner()


def apply(seq, i0, i1, rep):
    return seq[:i0] + rep + seq[i1:]


def test_two_designed_edits():
    s1, s2 = apply(WT, 2, 3, 'A'), apply(WT, 20, 21, 'G')
    designed = {WT, s1, s2}
    dbl = apply(apply(WT, 20, 21, 'G'), 2, 3, 'A')
    cls, kinds = a2.classify(WT, designed, a2.edits(AL, WT, dbl))
    assert cls == 'two designed edits' and kinds == ['designed', 'designed']


def test_single_nondesigned_sub_and_homopolymer_deletion():
    designed = {WT}
    cls, _ = a2.classify(WT, designed, a2.edits(AL, WT, apply(WT, 10, 11, 'T')))
    assert cls == 'one non-designed substitution'
    cls, _ = a2.classify(WT, designed, a2.edits(AL, WT, apply(WT, 4, 5, '')))  # TTT -> TT
    assert cls == 'one non-designed indel (homopolymer)'


def test_large_deletion_and_designed_plus_other():
    designed = {WT, apply(WT, 2, 3, 'A')}
    cls, _ = a2.classify(WT, designed, a2.edits(AL, WT, apply(WT, 8, 22, '')))
    assert cls.startswith('large')
    cls, _ = a2.classify(WT, designed, a2.edits(AL, WT, apply(apply(WT, 25, 26, 'A'), 2, 3, 'A')))
    assert cls == 'designed edit + one other edit'


def test_designed_homopolymer_deletion_is_placement_invariant():
    d = apply(WT, 3, 4, '')  # delete one T of TTT
    cls, kinds = a2.classify(WT, {WT, d}, a2.edits(AL, WT, apply(apply(WT, 5, 6, ''), 20, 21, 'G')))
    assert kinds.count('designed') == 1 and cls == 'designed edit + one other edit'

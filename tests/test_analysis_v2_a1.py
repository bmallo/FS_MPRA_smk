import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location('a1', Path(__file__).parents[1] / 'tools/analysis_v2/a1_inventory.py')
a1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a1)


def test_software_and_caller_jasmine_flag():
    pg = [dict(ID='ccs', VN='8.7.0 (commit v8.7.0)', CL='ccs ...'),
          dict(ID='jasmine', VN='2.7.99 (commit x)', CL='jasmine --6ma --trt', DS='Call select base modifications'),
          dict(ID='ft.1', VN='0.8.0', CL='/bin/ft add-nucleosome -t 12 in.bam -'),
          dict(ID='ft.2', VN='0.8.0', CL='/bin/ft clear-kinetics -t 12 - out.bam')]
    sw = a1.software(pg)
    assert sw['ccs'] == '8.7.0' and sw['jasmine'] == '2.7.99' and sw['ft_steps'] == ['add-nucleosome', 'clear-kinetics']
    assert a1.m6a_caller(sw, 'R/P2-C2/5.0-25M') == 'Jasmine 2.7.99 (on instrument)'


def test_caller_fibertools_predict():
    pg = [dict(ID='jasmine', VN='2.2.1 (commit v2.2.1)', CL='jasmine --streamed', DS='Predict 5mC in PacBio HiFi reads.'),
          dict(ID='ft.1', VN='0.5.4', CL='/bin/ft predict-m6a -t 12 a.bam b.bam')]
    assert a1.m6a_caller(a1.software(pg), 'R/P1-C1/5.0-25M') == 'fibertools predict-m6a 0.5.4'


def test_caller_sprq_without_flag():
    pg = [dict(ID='jasmine', VN='26.1.3 (commit v26.1.3)', CL='jasmine --streamed --trt')]
    assert a1.m6a_caller(a1.software(pg), 'R/P2-C3/5.0-25M') == 'Jasmine 26.1.3 (on instrument)'

"""Library-wide per-site m6A with same-reaction proxy controls (U/M/C by distal methylation)."""
import json
from pathlib import Path
REPO=Path(workflow.basedir).parent
CFG=str(REPO/config.get('site_config','config/ldlr_site_methylation.json'))
with open(CFG) as h: SPEC=json.load(h)
OUT=str(REPO/SPEC['output_dir'])
PY=SPEC['processing_python']
ENV='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLCONFIGDIR=/tmp/fsmpra-matplotlib'
SENS=[f'u{u}_m{m}' for u,m in SPEC['analysis']['sensitivity_thresholds']]

rule all:
    input: OUT+'/effects', expand(OUT+'/sensitivity/{s}', s=SENS), OUT+'/effects_supported_barcodes'

rule site_extract:
    input:
        cfg=CFG,
        features=str(REPO/SPEC['features']),
        feature_audit=str(REPO/SPEC['feature_audit']),
        catalog=str(REPO/SPEC['catalog']),
        contract=str(REPO/SPEC['footprint_contract']),
        code=[str(REPO/'fsmpra'/p) for p in ['site_extract.py','raw_candidate_extract.py','heterogeneity.py','footprint_join.py','io.py','processing.py']]
    output: directory(OUT+'/extract')
    threads: 1
    resources: mem_mb=8192, runtime=180
    params: py=PY,repo=str(REPO)
    shell: 'cd {params.repo:q} && '+ENV+' {params.py:q} -m fsmpra.site_extract --config {input.cfg:q} --output {output:q}'

rule site_effects:
    input:
        cfg=CFG,extract=OUT+'/extract',
        code=[str(REPO/'fsmpra'/p) for p in ['site_effects.py','library_stats.py','io.py']]
    output: directory(OUT+'/effects')
    threads: 1
    resources: mem_mb=16384, runtime=60
    params: py=PY,repo=str(REPO)
    shell: 'cd {params.repo:q} && '+ENV+' {params.py:q} -m fsmpra.site_effects --config {input.cfg:q} --extract {input.extract:q} --output {output:q}'

rule site_effects_sensitivity:
    input:
        cfg=CFG,extract=OUT+'/extract',
        code=[str(REPO/'fsmpra'/p) for p in ['site_effects.py','library_stats.py','io.py']]
    output: directory(OUT+'/sensitivity/u{u}_m{m}')
    threads: 1
    resources: mem_mb=16384, runtime=60
    params: py=PY,repo=str(REPO)
    shell: 'cd {params.repo:q} && '+ENV+' {params.py:q} -m fsmpra.site_effects --config {input.cfg:q} --extract {input.extract:q} --output {output:q} --u-max {wildcards.u} --m-min {wildcards.m} --aa-draws 100'

rule site_effects_supported_barcodes:
    input:
        cfg=CFG,extract=OUT+'/extract',
        code=[str(REPO/'fsmpra'/p) for p in ['site_effects.py','library_stats.py','io.py']]
    output: directory(OUT+'/effects_supported_barcodes')
    threads: 1
    resources: mem_mb=16384, runtime=60
    params: py=PY,repo=str(REPO)
    shell: 'cd {params.repo:q} && '+ENV+' {params.py:q} -m fsmpra.site_effects --config {input.cfg:q} --extract {input.extract:q} --output {output:q} --cohort eligible_with_barcode --aa-draws 100'

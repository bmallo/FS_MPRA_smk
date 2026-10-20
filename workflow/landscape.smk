"""Descriptive WT/SNV footprint landscape and local atlas."""
import json
from pathlib import Path
REPO=Path(workflow.basedir).parent
CFG=str(REPO/config.get('landscape_config','config/ldlr_landscape.json'))
with open(CFG) as h: SPEC=json.load(h)
OUT=str(REPO/SPEC['output_dir'])
PY=SPEC['processing_python']

rule all:
    input: OUT+'/atlas', OUT+'/patterns'

rule landscape:
    input:
        cfg=CFG,
        joined=str(REPO/SPEC['joined_dir']/'molecule_footprints.tsv.gz'),
        join_audit=str(REPO/SPEC['joined_dir']/'audit.json'),
        catalog=str(REPO/SPEC['catalog']),
        aggregates=[str(REPO/SPEC['aggregates']/p) for p in ['audit.json','barcode_counts.npy','barcode_molecules.npy','barcode_metadata.tsv.gz','variant_metadata.tsv']],
        code=[str(REPO/'fsmpra'/p) for p in ['landscape.py','io.py','processing.py']]
    output: directory(OUT+'/landscape')
    threads: 1
    resources: mem_mb=4096, runtime=30
    params: py=PY,repo=str(REPO)
    shell: 'cd {params.repo:q} && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 {params.py:q} -m fsmpra.landscape --config {input.cfg:q} --output {output:q}'

rule landscape_atlas:
    input:
        cfg=CFG,landscape=OUT+'/landscape',
        heterogeneity=str(REPO/SPEC['heterogeneity']),
        mpra=str(REPO/'results/pilot_20261014/catalog/catalog.json'),
        coordinates=str(REPO/config.get('coordinate_audit','results/rep3_landscape_20261018/genomic_coordinates_20261019/audit.json')),
        code=[str(REPO/'fsmpra'/p) for p in ['landscape_atlas.py','landscape.py','io.py','processing.py','genomic_coordinates.py']],
        assets=[str(REPO/'dashboard'/p) for p in ['landscape.html','landscape.js']]
    output: directory(OUT+'/atlas')
    threads: 1
    resources: mem_mb=4096, runtime=15
    params: py=PY,repo=str(REPO)
    shell: 'cd {params.repo:q} && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLCONFIGDIR=/tmp/fsmpra-matplotlib {params.py:q} -m fsmpra.landscape_atlas --config {input.cfg:q} --mpra-catalog {input.mpra:q} --coordinate-audit {input.coordinates:q} --output {output:q}'

rule landscape_patterns:
    input:
        atlas=OUT+'/atlas',
        code=str(REPO/'tools/summarize_landscape_patterns.py'),
        helpers=[str(REPO/'fsmpra'/p) for p in ['io.py','processing.py']]
    output: directory(OUT+'/patterns')
    threads: 1
    resources: mem_mb=2048, runtime=5
    params: py=PY,repo=str(REPO)
    shell: 'cd {params.repo:q} && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLCONFIGDIR=/tmp/fsmpra-matplotlib {params.py:q} -m tools.summarize_landscape_patterns --atlas {input.atlas:q} --output {output:q}'

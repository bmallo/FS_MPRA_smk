"""Replicate-3 fixed-endpoint barcode statistics; run with --snakefile."""
import json
import shlex
import sys
from pathlib import Path

REPO = Path(workflow.basedir).parent
CFG = str(Path(config.get("statistics_config", REPO / "config/ldlr_statistics.json")).resolve())
with open(CFG) as handle:
    SPEC = json.load(handle)
JOIN = str(Path(config.get("statistics_joined_dir", REPO / SPEC["joined_dir"])).resolve())
CATALOG = str(Path(config.get("statistics_catalog", REPO / SPEC["catalog"])).resolve())
OUT = str(Path(config.get("statistics_output", REPO / "results/rep3_statistics_20261017")).resolve())
AGG = OUT + "/aggregates"
INF = OUT + "/inference"
CAL = OUT + "/calibration"
RUN_CALIBRATION = config.get("statistics_run_calibration", True)
PY = config.get("statistics_python", "/mmfs1/home/bmallo/miniconda3/bin/python")
CODE = str(REPO / "fsmpra/statistics.py")
CAL_CODE = str(REPO / "fsmpra/statistics_calibration.py")

rule all:
    input: INF

rule aggregate_statistics:
    input: cfg=CFG, catalog=CATALOG, joined=JOIN + "/molecule_footprints.tsv.gz",
           audit=JOIN + "/audit.json", code=CODE
    output: directory(AGG)
    threads: 1
    resources: mem_mb=8192, runtime=240
    shell:
        "cd {REPO:q} && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "
        "NUMEXPR_MAX_THREADS=1 {PY:q} -m fsmpra.statistics aggregate "
        "--config {input.cfg:q} --catalog {input.catalog:q} "
        "--joined-dir {JOIN:q} --output {output:q}"

if RUN_CALIBRATION:
    rule calibrate_statistics:
        input: cfg=CFG, aggregates=AGG, code=CAL_CODE, statistics_code=CODE
        output: directory(CAL)
        threads: 1
        resources: mem_mb=8192, runtime=180
        shell:
            "cd {REPO:q} && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "
            "NUMEXPR_MAX_THREADS=1 {PY:q} -m fsmpra.statistics_calibration "
            "--config {input.cfg:q} --aggregated {input.aggregates:q} --output {output:q}"

rule infer_statistics:
    input: cfg=CFG, catalog=CATALOG, aggregates=AGG, code=CODE,
           calibration=CAL if RUN_CALIBRATION else []
    output: directory(INF)
    threads: 1
    resources: mem_mb=8192, runtime=120
    params: calibration_option=("--calibration " + shlex.quote(CAL + "/report.json")) if RUN_CALIBRATION else ""
    shell:
        "cd {REPO:q} && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "
        "NUMEXPR_MAX_THREADS=1 {PY:q} -m fsmpra.statistics infer "
        "--config {input.cfg:q} --catalog {input.catalog:q} "
        "--aggregates {input.aggregates:q} --output {output:q} {params.calibration_option}"

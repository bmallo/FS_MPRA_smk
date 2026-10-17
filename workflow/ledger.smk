"""Existing indexed alignment -> complete barcode/genotype ledger.

This entrypoint explicitly reuses alignment; it does not claim a new raw-BAM run.
Run with --cores 4, processing_config, ledger_bam, and ledger_run_dir.
"""
import json
import sys
from pathlib import Path

REPO = Path(workflow.basedir).parent
CFG = str(Path(config.get("processing_config", "config/ldlr_processing.json")).resolve())
with open(CFG) as handle:
    C = json.load(handle)
BAM = str(Path(config["ledger_bam"]).resolve())
RUN = str(Path(config["ledger_run_dir"]).resolve())
PY = config.get("processing_python", sys.executable)
CODE = [str(p) for p in (REPO / "fsmpra").glob("*.py")]
def source(p):
    return str(Path(p) if Path(p).is_absolute() else Path(C["project_root"]) / p)

rule all:
    input: RUN + "/ledger"

rule catalog:
    input: CODE + [CFG, source(C["library"]["reference_fasta"]), source(C["library"]["design_oligos"]), C["alignment"]["reference_fasta"]]
    output: directory(RUN + "/catalog")
    params: py=PY, cfg=CFG, repo=str(REPO)
    threads: 1
    resources: mem_mb=512
    shell: "cd {params.repo:q} && {params.py:q} -m fsmpra process --config {params.cfg:q} --stage prepare --output {output:q}"

rule extract:
    input: catalog=rules.catalog.output, bam=BAM, bai=BAM + ".bai", code=CODE, cfg=CFG
    output: directory(RUN + "/extracted")
    params: py=PY, cfg=CFG, repo=str(REPO), workers=lambda wildcards, threads: threads - 1
    threads: 4
    resources: mem_mb=4096, runtime=180
    shell: "cd {params.repo:q} && {params.py:q} -m fsmpra.ledger extract --config {params.cfg:q} --catalog {input.catalog:q}/catalog.json --bam {input.bam:q} --workers {params.workers} --allocated-cpus {threads} --output {output:q}"

rule barcodes:
    input: catalog=rules.catalog.output, extracted=rules.extract.output, code=CODE, cfg=CFG
    output: directory(RUN + "/barcodes")
    params: py=PY, cfg=CFG, repo=str(REPO)
    threads: 1
    resources: mem_mb=4096, runtime=180
    shell: "cd {params.repo:q} && {params.py:q} -m fsmpra process --config {params.cfg:q} --stage barcodes --catalog {input.catalog:q}/catalog.json --extracted {input.extracted:q} --output {output:q}"

rule ledger:
    input: catalog=rules.catalog.output, extracted=rules.extract.output, barcodes=rules.barcodes.output, code=CODE, cfg=CFG
    output: directory(RUN + "/ledger")
    params: py=PY, cfg=CFG, repo=str(REPO)
    threads: 1
    resources: mem_mb=1024, runtime=120
    shell: "cd {params.repo:q} && {params.py:q} -m fsmpra.ledger join --config {params.cfg:q} --catalog {input.catalog:q}/catalog.json --extracted {input.extracted:q} --barcodes {input.barcodes:q} --output {output:q}"

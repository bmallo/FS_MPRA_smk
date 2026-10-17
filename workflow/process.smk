"""Raw HiFi BAM -> molecule features and MPRAsnakeflow assignment export.

Run from the repository root with --config processing_config=config/ldlr_processing.json.
Snakemake owns resource scheduling; caller subprocesses never submit jobs.
"""
import sys
from pathlib import Path

REPO = Path(workflow.basedir).parent
CFG = str(Path(config.get("processing_config", "config/ldlr_processing.json")).resolve())
# Configuration is passed by path so exactly the same file is audited by every CLI.
import json
with open(CFG) as handle:
    C = json.load(handle)
RUN = str(Path(config.get("processing_run_dir", C["output_dir"])).resolve())
PY = config.get("processing_python", sys.executable)
CODE = [str(p) for p in (REPO / "fsmpra").glob("*.py")]
BASE = CODE
SOURCE = C["input"]["bam"]
LIB = C["library"]
def source(p):
    return str(Path(p) if Path(p).is_absolute() else Path(C["project_root"]) / p)

def signature(*keys):
    return json.dumps({key: C[key] for key in keys}, sort_keys=True)

STRUCTURE = {k: C["qc"][k] for k in ("end_slop", "max_softclip", "concatemer_length_fraction", "ml_threshold")}

rule all:
    input:
        RUN + "/features",
        RUN + "/barcodes",

rule preflight:
    input: BASE + [SOURCE]
    output: directory(RUN + "/preflight")
    params: py=PY, cfg=CFG, repo=str(REPO), settings=signature("input")
    threads: 1
    resources: mem_mb=512
    shell: "cd {params.repo:q} && {params.py:q} -m fsmpra process --config {params.cfg:q} --stage preflight --output {output:q}"

rule catalog:
    input: BASE + [source(LIB["reference_fasta"]), source(LIB["design_oligos"]), C["alignment"]["reference_fasta"], C["alignment"]["reference_fasta"] + ".fai"]
    output: directory(RUN + "/catalog")
    params: py=PY, cfg=CFG, repo=str(REPO), settings=signature("library", "project_root")
    threads: 1
    resources: mem_mb=512
    shell: "cd {params.repo:q} && {params.py:q} -m fsmpra process --config {params.cfg:q} --stage prepare --output {output:q}"

rule align:
    input:
        checked=rules.preflight.output, catalog=rules.catalog.output,
        raw=SOURCE, reference=C["alignment"].get("index") or C["alignment"]["reference_fasta"], code=BASE,
    output: directory(RUN + "/alignment")
    params: py=PY, cfg=CFG, repo=str(REPO), settings=signature("input", "alignment", "tools", "tool_versions")
    threads: 4
    resources: mem_mb=40000, runtime=240
    shell: "cd {params.repo:q} && {params.py:q} -m fsmpra process --config {params.cfg:q} --stage align --threads {threads} --output {output:q}"

rule extract:
    input: aligned=rules.align.output, catalog=rules.catalog.output, code=BASE
    output: directory(RUN + "/extracted")
    params: py=PY, cfg=CFG, repo=str(REPO), settings=signature("qc", "sample")
    threads: 1
    resources: mem_mb=2048, runtime=180
    shell: "cd {params.repo:q} && {params.py:q} -m fsmpra process --config {params.cfg:q} --stage extract --catalog {input.catalog:q}/catalog.json --bam {input.aligned:q}/aligned.bam --output {output:q}"

rule barcodes:
    input: extracted=rules.extract.output, catalog=rules.catalog.output, code=BASE
    output: directory(RUN + "/barcodes")
    params: py=PY, cfg=CFG, repo=str(REPO), settings=signature("barcodes", "qc")
    threads: 1
    resources: mem_mb=2048, runtime=180
    shell: "cd {params.repo:q} && {params.py:q} -m fsmpra process --config {params.cfg:q} --stage barcodes --catalog {input.catalog:q}/catalog.json --extracted {input.extracted:q} --output {output:q}"

rule partition:
    input: aligned=rules.align.output, catalog=rules.catalog.output, code=BASE
    output: directory(RUN + "/partition")
    params: py=PY, cfg=CFG, repo=str(REPO), settings=json.dumps(STRUCTURE, sort_keys=True)
    threads: 1
    resources: mem_mb=1024, runtime=120
    shell: "cd {params.repo:q} && {params.py:q} -m fsmpra process --config {params.cfg:q} --stage partition --catalog {input.catalog:q}/catalog.json --bam {input.aligned:q}/aligned.bam --output {output:q}"

rule footprints:
    input: extracted=rules.partition.output, model=C["footprints"]["model"], code=BASE
    output: directory(RUN + "/footprints")
    params: py=PY, cfg=CFG, repo=str(REPO), settings=signature("footprints") + json.dumps({"ml": C["qc"]["ml_threshold"], "tool": C["tools"]["fiberhmm"], "version": C["tool_versions"]["fiberhmm"]})
    threads: 4
    resources: mem_mb=16000, runtime=240
    shell: "cd {params.repo:q} && {params.py:q} -m fsmpra process --config {params.cfg:q} --stage footprints --extracted {input.extracted:q} --threads {threads} --output {output:q}"

rule features:
    input: extracted=rules.extract.output, barcodes=rules.barcodes.output, footprints=rules.footprints.output, catalog=rules.catalog.output, code=BASE
    output: directory(RUN + "/features")
    params: py=PY, cfg=CFG, repo=str(REPO)
    threads: 1
    resources: mem_mb=2048, runtime=180
    shell: "cd {params.repo:q} && {params.py:q} -m fsmpra process --config {params.cfg:q} --stage features --catalog {input.catalog:q}/catalog.json --extracted {input.extracted:q} --barcodes {input.barcodes:q} --footprints {input.footprints:q} --output {output:q}"

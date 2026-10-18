"""Replicate-3 molecule features followed by descriptive heterogeneity summaries."""
import json
from pathlib import Path

REPO = Path(workflow.basedir).parent
CFG = str((REPO / config.get("heterogeneity_config", "config/ldlr_heterogeneity.json")).resolve())
with open(CFG) as handle:
    SPEC = json.load(handle)
if SPEC.get("schema_version") != "fsmpra.heterogeneity_config.v1":
    raise ValueError("Expected heterogeneity configuration")
for key in ("ledger_run_dir", "footprint_run_dir", "footprint_contract", "output_dir"):
    if key in config and config[key] != SPEC[key]:
        raise ValueError(f"{key} must come from the heterogeneity configuration, not a Snakemake override")
RUN = (REPO / SPEC["ledger_run_dir"]).resolve()
JOIN = (REPO / SPEC["footprint_run_dir"]).resolve()
CONTRACT = (REPO / SPEC["footprint_contract"]).resolve()
with CONTRACT.open() as handle:
    BAM = json.load(handle)["bam"]
OUT = (REPO / config.get("output_override", SPEC["output_dir"])).resolve()
FEATURES = str(OUT / "features")
SUMMARY = str(OUT / "summary")
PY = SPEC.get("processing_python", "/mmfs1/home/bmallo/miniconda3/bin/python")

rule all:
    input: SUMMARY

rule extract_heterogeneity:
    input:
        cfg=CFG,
        joined=str(JOIN / "molecule_footprints.tsv.gz"),
        join_audit=str(JOIN / "audit.json"),
        ledger_cfg=str(RUN / "config.json"),
        catalog=str(RUN / "catalog/catalog.json"),
        variant_coverage=str(RUN / "barcodes/variant_barcode_coverage.tsv"),
        contract=str(CONTRACT),
        bam=BAM,
        bai=BAM + ".bai",
        code=[str(REPO / "fsmpra" / name) for name in
              ("heterogeneity.py", "footprint_join.py", "processing.py", "io.py")]
    output: directory(FEATURES)
    threads: 1
    resources: mem_mb=8192, runtime=720
    params: py=PY, repo=str(REPO)
    shell:
        "cd {params.repo:q} && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 {params.py:q} -m fsmpra.heterogeneity "
        "--config {input.cfg:q} --output {output:q}"

rule summarize_heterogeneity:
    input:
        cfg=CFG,
        features=FEATURES,
        code=str(REPO / "fsmpra/heterogeneity_summary.py")
    output: directory(SUMMARY)
    threads: 1
    resources: mem_mb=8192, runtime=240
    params: py=PY, repo=str(REPO)
    shell:
        "cd {params.repo:q} && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 {params.py:q} -m fsmpra.heterogeneity_summary "
        "--config {input.cfg:q} --features {input.features:q} --output {output:q}"

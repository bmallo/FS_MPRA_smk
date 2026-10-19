"""Frozen WT barcode sample, descriptive call-quality summary, and benchmark."""
import json
from pathlib import Path

REPO = Path(workflow.basedir).parent
CFG = str((REPO / config.get("footprint_quality_config", "config/ldlr_footprint_quality.json")).resolve())
with open(CFG) as handle:
    SPEC = json.load(handle)
if SPEC.get("schema_version") != "fsmpra.footprint_quality_config.v1":
    raise ValueError("Expected footprint-quality configuration")
for key in ("features", "feature_audit", "barcode_metadata", "barcode_audit", "catalog",
            "ledger_config", "footprint_contract", "seed", "sample_barcodes"):
    if key in config and config[key] != SPEC[key]:
        raise ValueError(f"{key} must come from the footprint-quality configuration")
with open((REPO / SPEC["footprint_contract"]).resolve()) as handle:
    BAM = json.load(handle)["bam"]
OUT = (REPO / config.get("output_override", SPEC["output_dir"])).resolve()
SAMPLE = str(OUT / "sample")
SUMMARY = str(OUT / "summary")
BENCHMARK = str(OUT / "benchmark")
PY = SPEC.get("processing_python", "/mmfs1/home/bmallo/miniconda3/bin/python")

rule all:
    input: BENCHMARK

rule extract_footprint_quality:
    input:
        cfg=CFG,
        features=str((REPO / SPEC["features"]).resolve()),
        feature_audit=str((REPO / SPEC["feature_audit"]).resolve()),
        barcode_metadata=str((REPO / SPEC["barcode_metadata"]).resolve()),
        barcode_audit=str((REPO / SPEC["barcode_audit"]).resolve()),
        catalog=str((REPO / SPEC["catalog"]).resolve()),
        ledger_config=str((REPO / SPEC["ledger_config"]).resolve()),
        contract=str((REPO / SPEC["footprint_contract"]).resolve()),
        bam=BAM,
        bai=BAM + ".bai",
        code=[str(REPO / "fsmpra" / name) for name in
              ("footprint_quality.py", "heterogeneity.py", "footprint_join.py",
               "processing.py", "io.py")]
    output: directory(SAMPLE)
    threads: 1
    resources: mem_mb=4096, runtime=30
    params: py=PY, repo=str(REPO)
    shell:
        "cd {params.repo:q} && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 {params.py:q} -m fsmpra.footprint_quality "
        "--config {input.cfg:q} --output {output:q}"

rule summarize_footprint_quality:
    input:
        cfg=CFG,
        sample=SAMPLE,
        code=[str(REPO / "fsmpra" / name) for name in
              ("footprint_quality_summary.py", "statistics.py", "io.py")]
    output: directory(SUMMARY)
    threads: 1
    resources: mem_mb=4096, runtime=30
    params: py=PY, repo=str(REPO)
    shell:
        "cd {params.repo:q} && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLCONFIGDIR=/tmp/fsmpra-matplotlib "
        "{params.py:q} -m fsmpra.footprint_quality_summary "
        "--config {input.cfg:q} --sample {input.sample:q} --output {output:q}"

rule benchmark_footprint_quality:
    input:
        cfg=CFG,
        summary=SUMMARY,
        code=str(REPO / "fsmpra/footprint_quality_benchmark.py")
    output: directory(BENCHMARK)
    threads: 1
    resources: mem_mb=4096, runtime=30
    params: py=PY, repo=str(REPO)
    shell:
        "cd {params.repo:q} && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLCONFIGDIR=/tmp/fsmpra-matplotlib "
        "{params.py:q} -m fsmpra.footprint_quality_benchmark "
        "--config {input.cfg:q} --summary {input.summary:q} --output {output:q}"

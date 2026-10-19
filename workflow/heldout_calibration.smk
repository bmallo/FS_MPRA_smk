"""Frozen held-out WT endpoint outcomes, calibration, and figures."""
import json
from pathlib import Path

REPO = Path(workflow.basedir).parent
CFG = str((REPO / config.get("heldout_config", "config/ldlr_heldout_calibration.json")).resolve())
with open(CFG) as handle:
    SPEC = json.load(handle)
if SPEC.get("schema_version") != "fsmpra.heldout_calibration_config.v1":
    raise ValueError("Expected held-out calibration configuration")
for key in ("aggregates", "joined_dir", "excluded_sample", "catalog", "seed"):
    if key in config and config[key] != SPEC[key]:
        raise ValueError(f"{key} must come from the held-out configuration")

OUT = (REPO / config.get("output_override", SPEC["output_dir"])).resolve()
OUTCOMES = str(OUT / "outcomes")
CALIBRATION = str(OUT / "calibration")
PLOTS = str(OUT / "plots")
AGGREGATES = (REPO / SPEC["aggregates"]).resolve()
JOINED = (REPO / SPEC["joined_dir"]).resolve()
SAMPLE = (REPO / SPEC["excluded_sample"]).resolve()
CATALOG = (REPO / SPEC["catalog"]).resolve()
PY = SPEC.get("processing_python", "/mmfs1/home/bmallo/miniconda3/bin/python")


rule all:
    input: PLOTS


rule heldout_outcomes:
    input:
        cfg=CFG,
        catalog=str(CATALOG),
        joined=str(JOINED / "molecule_footprints.tsv.gz"),
        joined_audit=str(JOINED / "audit.json"),
        barcode_metadata=str(AGGREGATES / "barcode_metadata.tsv.gz"),
        barcode_counts=str(AGGREGATES / "barcode_counts.npy"),
        barcode_molecules=str(AGGREGATES / "barcode_molecules.npy"),
        aggregate_audit=str(AGGREGATES / "audit.json"),
        sample=str(SAMPLE / "sampled_barcodes.tsv"),
        sample_audit=str(SAMPLE / "audit.json"),
        sample_summary=str(SAMPLE.parent / "summary" / "summary.json"),
        code=[str(REPO / "fsmpra" / name) for name in
              ("heldout_outcomes.py", "statistics.py", "io.py")],
    output: directory(OUTCOMES)
    threads: 1
    resources: mem_mb=4096, runtime=30
    params: py=PY, repo=str(REPO)
    shell:
        "cd {params.repo:q} && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 {params.py:q} -m fsmpra.heldout_outcomes "
        "--config {input.cfg:q} --output {output:q}"


rule heldout_calibration:
    input:
        cfg=CFG,
        outcomes=OUTCOMES,
        code=[str(REPO / "fsmpra" / name) for name in
              ("heldout_calibration.py", "statistics.py", "processing.py", "io.py")],
    output: directory(CALIBRATION)
    threads: 1
    resources: mem_mb=4096, runtime=30
    params: py=PY, repo=str(REPO)
    shell:
        "cd {params.repo:q} && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 {params.py:q} -m fsmpra.heldout_calibration "
        "--config {input.cfg:q} --outcomes {input.outcomes:q} --output {output:q}"


rule plot_heldout_calibration:
    input:
        calibration=CALIBRATION,
        code=str(REPO / "tools/plot_heldout_calibration.py"),
    output: directory(PLOTS)
    threads: 1
    resources: mem_mb=4096, runtime=10
    params: py=PY, repo=str(REPO)
    shell:
        "cd {params.repo:q} && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLCONFIGDIR=/tmp/fsmpra-matplotlib "
        "{params.py:q} tools/plot_heldout_calibration.py "
        "--calibration {input.calibration:q} --output {output:q}"

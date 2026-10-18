"""Join existing annotations without rebuilding the validated sequence ledger."""
import json
import sys
from pathlib import Path

REPO = Path(workflow.basedir).parent
RUN = str(Path(config.get("ledger_run_dir", "results/rep3_ledger_20261016")).resolve())
CONTRACT = str(Path(config.get("footprint_contract", "config/rep3_existing_footprints.json")).resolve())
with open(CONTRACT) as handle:
    SPEC = json.load(handle)
BAM = SPEC["bam"]
OUT = str(Path(config.get("footprint_output", "results/rep3_footprints_20261017")).resolve())
PY = config.get("processing_python", sys.executable)
CODE = [str(REPO / "fsmpra" / name) for name in ("footprint_join.py", "processing.py", "storage.py", "io.py")]

rule all:
    input: OUT

rule join_footprints:
    input:
        ledger=RUN + "/ledger/molecule_ledger.tsv.gz",
        validation=RUN + "/ledger/validation.json",
        catalog=RUN + "/catalog/catalog.json",
        cfg=RUN + "/config.json",
        contract=CONTRACT, bam=BAM, bai=BAM + ".bai", code=CODE
    output: directory(OUT)
    threads: 1
    resources: mem_mb=2048, runtime=180
    params: py=PY, repo=str(REPO)
    shell:
        "cd {params.repo:q} && {params.py:q} -m fsmpra.footprint_join "
        "--config {input.cfg:q} --catalog {input.catalog:q} --ledger {input.ledger:q} "
        "--bam {input.bam:q} --contract {input.contract:q} --output {output:q}"

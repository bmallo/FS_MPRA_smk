"""Compare existing measurements; run workflow/pilot.smk to produce new ones."""
import sys
from pathlib import Path

REPO = Path(workflow.basedir).parent
configfile: str(REPO / "config/ldlr_comparison.json")


def relative_to_repo(value):
    path = Path(value)
    return str(path if path.is_absolute() else REPO / path)


PILOTS = [relative_to_repo(p) for p in config["comparison_inputs"]]
MANIFEST = relative_to_repo(config["comparison_manifest"])
RUN = relative_to_repo(config["comparison_output"])
PYTHON = str(config.get("pilot_python", sys.executable))
CODE = [str(p) for p in sorted((REPO / "fsmpra").glob("*.py"))]
ASSETS = [str(REPO / "dashboard" / name) for name in ("index.html", "app.js", "style.css")]


rule all:
    input:
        RUN + "/dashboard",


rule compare_conditions:
    input:
        pilots=PILOTS,
        manifest=MANIFEST,
        code=CODE,
    output:
        bundle=directory(RUN + "/comparison"),
    params:
        python=PYTHON, repo=str(REPO), minimum=int(config.get("comparison_min_reads", 10)),
    threads: 1
    resources:
        mem_mb=1024,
    shell:
        "cd {params.repo:q} && {params.python:q} -m fsmpra compare --pilots {input.pilots:q} "
        "--manifest {input.manifest:q} --min-reads {params.minimum} --output {output.bundle:q}"


rule comparison_dashboard:
    input:
        bundle=rules.compare_conditions.output.bundle,
        code=CODE,
        assets=ASSETS,
    output:
        bundle=directory(RUN + "/dashboard"),
    params:
        python=PYTHON, repo=str(REPO),
    threads: 1
    resources:
        mem_mb=1024,
    shell:
        "cd {params.repo:q} && {params.python:q} -m fsmpra dashboard "
        "--pilot {input.bundle:q}/collection.json --output {output.bundle:q}"

"""Independent redesign workflow; the legacy Snakefile is unchanged.

Run from repository root:
  snakemake -s workflow/pilot.smk --cores 1
Use a new pilot_run_dir when changing inputs or measurement settings.
"""

import json
import sys
from pathlib import Path

REPO = Path(workflow.basedir).parent
MANIFEST = str(Path(config.get("pilot_manifest", REPO / "config/ldlr_pilot.json")).resolve())
with open(MANIFEST) as handle:
    MAN = json.load(handle)
RUN = str(Path(config.get("pilot_run_dir", REPO / "results/pilot_run")).resolve())
SAMPLE = config.get("pilot_sample", "HepG2_rep3")
FRACTION = float(config.get("pilot_fraction", 0.01))
SEED = int(config.get("pilot_seed", 42))
CAP = int(config.get("pilot_max_molecules", 200000))
PYTHON = str(config.get("pilot_python", sys.executable))
BACKEND = config.get("pilot_backend", "memory")
SCRATCH = str(config.get("pilot_scratch_dir", "/tmp"))
CACHE = int(config.get("pilot_cache_mib", 64))
DISK = float(config.get("pilot_max_disk_gib", 10))
BATCH = int(config.get("pilot_batch_size", 10000))
SOURCE_SAMPLE = next(s for s in MAN["samples"] if s["id"] == SAMPLE)


def source(path):
    p = Path(path)
    return str(p if p.is_absolute() else Path(MAN["project_root"]) / p)


SMALL_INPUTS = [MANIFEST]
for library in MAN["libraries"]:
    SMALL_INPUTS.extend(source(library[k]) for k in ("reference_fasta", "design_oligos"))
SMALL_INPUTS.extend(source(m["path"]) for m in MAN["mpra_sources"])
SMALL_INPUTS.extend(source(s["inputs"]["design_coverage"]) for s in MAN["samples"] if "design_coverage" in s["inputs"])
SOURCE_INPUTS = sorted({source(value) for s in MAN["samples"] for value in s["inputs"].values()})
CODE = [str(p) for p in sorted((REPO / "fsmpra").glob("*.py"))]
VIEWER = [str(REPO / "dashboard" / name) for name in ("index.html", "app.js", "style.css")]


rule all:
    input:
        RUN + "/" + SAMPLE + "/measurements",
        RUN + "/" + SAMPLE + "/dashboard",


rule pilot_catalog:
    input:
        SMALL_INPUTS + SOURCE_INPUTS + CODE,
    output:
        bundle=directory(RUN + "/catalog"),
    params:
        repo=str(REPO), python=PYTHON, manifest=MANIFEST,
    threads: 1
    resources:
        mem_mb=512,
    shell:
        "cd {params.repo:q} && {params.python:q} -m fsmpra catalog "
        "--manifest {params.manifest:q} --output {output.bundle:q}"


rule pilot_molecules:
    input:
        catalog=rules.pilot_catalog.output.bundle,
        ledger=source(SOURCE_SAMPLE["inputs"]["ledger"]),
        metrics=source(SOURCE_SAMPLE["inputs"]["metrics"]),
        code=CODE,
    output:
        bundle=directory(RUN + "/" + SAMPLE + "/molecules"),
    params:
        repo=str(REPO), python=PYTHON, sample=SAMPLE,
        fraction=FRACTION, seed=SEED, cap=CAP,
        backend=BACKEND, scratch=SCRATCH, cache=CACHE, disk=DISK, batch=BATCH,
    threads: 1
    resources:
        mem_mb=2048,
    shell:
        "cd {params.repo:q} && {params.python:q} -m fsmpra molecules "
        "--catalog {input.catalog:q}/catalog.json --sample {params.sample:q} "
        "--fraction {params.fraction} --seed {params.seed} --max-molecules {params.cap} "
        "--backend {params.backend:q} --scratch-dir {params.scratch:q} --cache-mib {params.cache} "
        "--max-disk-gib {params.disk} --batch-size {params.batch} "
        "--output {output.bundle:q}"


rule pilot_measure:
    input:
        catalog=rules.pilot_catalog.output.bundle,
        molecules=rules.pilot_molecules.output.bundle,
        bam=source(SOURCE_SAMPLE["inputs"]["footprints_bam"]),
        bai=source(SOURCE_SAMPLE["inputs"]["footprints_bam"]) + ".bai",
        code=CODE,
    output:
        bundle=directory(RUN + "/" + SAMPLE + "/measurements"),
    params:
        repo=str(REPO), python=PYTHON,
        backend=BACKEND, scratch=SCRATCH, cache=CACHE, disk=DISK,
        min_passes=int(config.get("pilot_min_passes", 10)),
        tf_min=int(config.get("pilot_tf_min", 10)), tf_max=int(config.get("pilot_tf_max", 60)),
        nuc_min=int(config.get("pilot_nuc_min", 85)), mask=int(config.get("pilot_mask_radius", 2)),
    threads: 1
    resources:
        mem_mb=2048,
    shell:
        "cd {params.repo:q} && {params.python:q} -m fsmpra measure "
        "--catalog {input.catalog:q}/catalog.json --molecules {input.molecules:q} "
        "--min-passes {params.min_passes} --tf-min {params.tf_min} --tf-max {params.tf_max} "
        "--backend {params.backend:q} --scratch-dir {params.scratch:q} --cache-mib {params.cache} --max-disk-gib {params.disk} "
        "--nuc-min {params.nuc_min} --mask-radius {params.mask} --output {output.bundle:q}"


rule pilot_dashboard:
    input:
        measurements=rules.pilot_measure.output.bundle,
        code=CODE,
        assets=VIEWER,
    output:
        bundle=directory(RUN + "/" + SAMPLE + "/dashboard"),
    params:
        repo=str(REPO), python=PYTHON,
    threads: 1
    resources:
        mem_mb=512,
    shell:
        "cd {params.repo:q} && {params.python:q} -m fsmpra dashboard "
        "--pilot {input.measurements:q}/pilot.json --output {output.bundle:q}"

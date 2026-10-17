"""Small, strict IO helpers; raw inputs are always opened read-only."""

import csv
import gzip
import hashlib
import json
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path


def open_text(path, mode="rt"):
    return gzip.open(path, mode, encoding="utf-8", newline="") if str(path).endswith(".gz") else open(path, mode, encoding="utf-8", newline="")


def read_json(path):
    with open(path) as handle:
        return json.load(handle)


def write_json(path, value):
    with open(path, "w") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def write_table(path, rows, fields):
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: ("" if row.get(key) is None else row.get(key)) for key in fields})


def source_info(path, checksum=False):
    path = Path(path).resolve()
    stat = path.stat()
    info = {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if checksum:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        info["sha256"] = digest.hexdigest()
    return info


def require_columns(reader, columns, path):
    missing = set(columns) - set(reader.fieldnames or [])
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")


def read_fasta(path, contig):
    records = {}
    current = None
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(f"Duplicate FASTA contig {current}")
                records[current] = []
            elif line.strip():
                if current is None:
                    raise ValueError("FASTA sequence precedes its header")
                records[current].append(line.strip().upper())
    if contig not in records:
        raise ValueError(f"Missing contig {contig} in {path}")
    seq = "".join(records[contig])
    if not seq or set(seq) - set("ACGTN"):
        raise ValueError(f"Invalid reference sequence for {contig}")
    return seq


def resolve_input(manifest, path):
    path = Path(path)
    return path if path.is_absolute() else Path(manifest["project_root"]) / path


def load_manifest(path):
    manifest = read_json(path)
    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported manifest schema_version")
    for group in ("libraries", "samples", "mpra_sources"):
        ids = [item["id"] for item in manifest[group]]
        if len(ids) != len(set(ids)):
            raise ValueError(f"Duplicate IDs in {group}")
    libraries = {item["id"] for item in manifest["libraries"]}
    for sample in manifest["samples"]:
        if sample["library_id"] not in libraries:
            raise ValueError(f"Unknown library for {sample['id']}")
    return manifest


def new_output(path):
    """Fail on existing results; callers must choose a new run directory."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=False)
    return path


@contextmanager
def atomic_output(path):
    """Publish a complete output directory only after every stage succeeds."""
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{path.name}-", dir=path.parent))
    try:
        yield stage
        if path.exists():
            raise FileExistsError(path)
        stage.rename(path)
    finally:
        if stage.exists():
            shutil.rmtree(stage)

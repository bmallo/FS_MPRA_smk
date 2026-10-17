"""Deterministic read sampling and auditable ledger/metrics joins.

Samples across the entire ledger, never the head of a coordinate-sorted BAM.
No methylation-based exclusion is performed. Memory is bounded by an explicit
selected-molecule cap, not by available system RAM.
"""

import csv
import hashlib
import logging
from collections import Counter

from . import SCHEMA_VERSION, __version__
from .catalog import source_variant_sequence
from .io import (new_output, open_text, read_json, require_columns, source_info,
                 write_json, write_table)

LEDGER_COLUMNS = ["read", "len", "np", "rq", "strand", "ref_start", "ref_end", "structure", "class", "variant", "design", "n_edits", "max_del_run", "tier", "barcode"]
METRIC_COLUMNS = ["read", "contig", "pc", "mf", "nc", "rf", "fs", "at", "m6a"]
FIELDS = ["sample_id", "read_id", "variant_id", "genotype_class", "source_variant",
          "source_design", "barcode", "passes", "read_quality", "read_length", "strand",
          "n_edits", "max_del_run", "genotype_tier", "ledger_structure", "plasmid_class",
          "global_methylation_fraction", "nucleosome_count", "unresolved_fraction",
          "source_filter_status", "metrics_present", "source_at_count", "source_m6a_count"]


def selected(read_id, fraction, seed):
    value = int.from_bytes(hashlib.blake2b(f"{seed}\0{read_id}".encode(), digest_size=8).digest(), "big")
    return value < int(fraction * (1 << 64))


def finite_float(value, name):
    result = float(value)
    if not -float("inf") < result < float("inf"):
        raise ValueError(f"Nonfinite {name}")
    return result



def ledger_molecule(row, sample_id, lib, alias, by_seq):
    label = row["variant"]
    vid = None
    if row["class"] in ("WT", "SNV", "INS", "DEL"):
        if label not in alias:
            seq = source_variant_sequence(lib["mutable_sequence"], label, lib["oligo_mutable_start1"])
            alias[label] = by_seq[seq]
        vid = alias[label]
    return {
        "sample_id": sample_id, "read_id": row["read"], "variant_id": vid,
        "genotype_class": row["class"], "source_variant": label,
        "source_design": row["design"], "barcode": row["barcode"],
        "passes": int(row["np"]), "read_quality": finite_float(row["rq"], "rq"),
        "read_length": int(row["len"]), "strand": row["strand"],
        "n_edits": int(row["n_edits"]), "max_del_run": int(row["max_del_run"]),
        "genotype_tier": row["tier"], "ledger_structure": row["structure"],
        "metrics_present": 0}


def join_metric(mol, row, lib):
    mf = finite_float(row["mf"], "mf")
    if not 0 <= mf <= 1:
        raise ValueError(f"Invalid methylation fraction for {row['read']}")
    if row["contig"] != lib["contig"]:
        raise ValueError(f"Ledger/metrics contig disagreement for {row['read']}")
    mol.update(plasmid_class=row["pc"], global_methylation_fraction=mf,
               nucleosome_count=int(row["nc"]), unresolved_fraction=finite_float(row["rf"], "rf"),
               source_filter_status=row["fs"], metrics_present=1,
               source_at_count=int(row["at"]), source_m6a_count=int(row["m6a"]))

def import_molecules(catalog_path, sample_id, output, fraction=0.01, seed=42, max_molecules=200000,
                     backend="memory", scratch_dir=None, cache_mib=64, max_disk_gib=10, batch_size=10000):
    if backend == "sqlite":
        from .disk_import import import_disk
        return import_disk(catalog_path, sample_id, output, fraction, seed, scratch_dir, cache_mib, max_disk_gib, batch_size)
    if backend != "memory":
        raise ValueError("Unknown storage backend")
    if not 0 < fraction <= 1 or max_molecules < 1:
        raise ValueError("fraction must be in (0,1]; max_molecules must be positive")
    catalog = read_json(catalog_path)
    sample = next(s for s in catalog["samples"] if s["id"] == sample_id)
    lib = next(l for l in catalog["libraries"] if l["id"] == sample["library_id"])
    by_seq = {v["sequence"]: v["variant_id"] for v in lib["variants"]}
    alias = {v["source_variant"]: v["variant_id"] for v in lib["variants"]}
    ledger_path, metrics_path = (sample["inputs"][k]["path"] for k in ("ledger", "metrics"))
    sources = [source_info(ledger_path), source_info(metrics_path)]
    for key, current in zip(("ledger", "metrics"), sources):
        if current != sample["inputs"][key]:
            raise ValueError(f"{key} changed since catalog creation; build a new catalog")
    molecules = {}
    all_classes, selected_classes = Counter(), Counter()
    with open_text(ledger_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader, LEDGER_COLUMNS, ledger_path)
        for number, row in enumerate(reader, 1):
            all_classes[row["class"]] += 1
            if number % 1000000 == 0:
                logging.info("Ledger: %s rows, %s selected", f"{number:,}", f"{len(molecules):,}")
            if not selected(row["read"], fraction, seed):
                continue
            if row["read"] in molecules:
                raise ValueError(f"Duplicate selected ledger read: {row['read']}")
            if len(molecules) >= max_molecules:
                raise ValueError("Selected molecule cap exceeded; reduce fraction or explicitly raise --max-molecules")
            selected_classes[row["class"]] += 1
            molecules[row["read"]] = ledger_molecule(row, sample_id, lib, alias, by_seq)
    n_metrics = 0
    with open_text(metrics_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader, METRIC_COLUMNS, metrics_path)
        for number, row in enumerate(reader, 1):
            n_metrics = number
            if number % 2000000 == 0:
                logging.info("Metrics: %s rows scanned", f"{number:,}")
            mol = molecules.get(row["read"])
            if mol is None:
                continue
            if mol["metrics_present"]:
                raise ValueError(f"Duplicate selected metric read: {row['read']}")
            join_metric(mol, row, lib)
    for source in sources:
        if source_info(source["path"]) != source:
            raise ValueError("An input changed during import")
    audit = {"schema_version": SCHEMA_VERSION, "software_version": __version__,
             "sample_id": sample_id, "catalog": source_info(catalog_path, True),
             "sampling": {"method": "blake2b-64 threshold on seed + NUL + read_id", "fraction": fraction, "seed": seed,
                          "max_molecules": max_molecules, "population": "entire primary-plasmid genotype ledger"},
             "sources": sources, "ledger_rows": sum(all_classes.values()), "metrics_rows": n_metrics,
             "ledger_classes": dict(all_classes), "selected_classes": dict(selected_classes),
             "selected_reads": len(molecules), "missing_metrics": sum(not m["metrics_present"] for m in molecules.values()),
             "methylation_filter": None, "limitations": ["Genotypes are imported exploratory assignments, not independently recalled.",
                                                        "Missing metric rows are retained explicitly.",
                                                        "Large-input provenance uses size and mtime; small catalog files use SHA256."]}
    out = new_output(output)
    write_table(out / "molecules.tsv.gz", molecules.values(), FIELDS)
    write_json(out / "audit.json", audit)
    return audit

"""Stream a ledger/metrics join through bounded SQLite scratch storage."""

import csv
import itertools
import json
import logging
import sqlite3
from collections import Counter

from . import SCHEMA_VERSION, __version__
from .io import atomic_output, open_text, read_json, require_columns, source_info, write_json, write_table
from .molecules import FIELDS, LEDGER_COLUMNS, METRIC_COLUMNS, join_metric, ledger_molecule, selected
from .storage import scratch_database


def batches(rows, size):
    while True:
        batch = list(itertools.islice(rows, size))
        if not batch:
            return
        yield batch


def encode(row):
    return json.dumps(row, separators=(",", ":"), allow_nan=False)


def pass_bin(value):
    return "<5" if value < 5 else "5-9" if value < 10 else "10-19" if value < 20 else "20-49" if value < 50 else "50+"


def length_bin(value, reference_length):
    ratio = value / reference_length
    return "<0.5x" if ratio < .5 else "0.5-0.9x" if ratio < .9 else "0.9-1.1x" if ratio < 1.1 else "1.1-2x" if ratio < 2 else ">=2x"


def import_disk(catalog_path, sample_id, output, fraction=1, seed=42, scratch_dir=None,
                cache_mib=64, max_disk_gib=10, batch_size=10000):
    if not 0 < fraction <= 1 or not 1 <= batch_size <= 100000:
        raise ValueError("fraction must be in (0,1]; batch size must be 1–100000")
    catalog_identity = source_info(catalog_path, True)
    catalog = read_json(catalog_path)
    sample = next(s for s in catalog["samples"] if s["id"] == sample_id)
    lib = next(l for l in catalog["libraries"] if l["id"] == sample["library_id"])
    by_seq = {v["sequence"]: v["variant_id"] for v in lib["variants"]}
    alias = {v["source_variant"]: v["variant_id"] for v in lib["variants"]}
    sources = [source_info(sample["inputs"][k]["path"]) for k in ("ledger", "metrics")]
    for key, identity in zip(("ledger", "metrics"), sources):
        if identity != sample["inputs"][key]:
            raise ValueError(f"{key} changed since catalog creation; build a new catalog")
    all_classes, chosen_classes, metric_classes = Counter(), Counter(), Counter()
    n_selected, n_metrics = 0, 0
    with atomic_output(output) as out, scratch_database(scratch_dir, cache_mib, max_disk_gib) as db:
        db.execute("CREATE TABLE molecules (ordinal INTEGER PRIMARY KEY, read_id TEXT NOT NULL UNIQUE, payload TEXT NOT NULL, matched INTEGER NOT NULL DEFAULT 0)")
        with open_text(sources[0]["path"]) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            require_columns(reader, LEDGER_COLUMNS, sources[0]["path"])
            for batch in batches(reader, batch_size):
                inserts = []
                for row in batch:
                    all_classes[row["class"]] += 1
                    if fraction < 1 and not selected(row["read"], fraction, seed):
                        continue
                    mol = ledger_molecule(row, sample_id, lib, alias, by_seq)
                    chosen_classes[row["class"]] += 1
                    n_selected += 1
                    inserts.append((n_selected, row["read"], encode(mol)))
                try:
                    db.executemany("INSERT INTO molecules(ordinal,read_id,payload) VALUES (?,?,?)", inserts)
                    db.commit()
                except sqlite3.IntegrityError as error:
                    raise ValueError("Duplicate selected ledger read") from error
                if sum(all_classes.values()) % 500000 < batch_size:
                    logging.info("Ledger: %s rows; %s selected on disk", f"{sum(all_classes.values()):,}", f"{n_selected:,}")
        with open_text(sources[1]["path"]) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            require_columns(reader, METRIC_COLUMNS, sources[1]["path"])
            for batch in batches(reader, batch_size):
                matched = {}
                # 500 placeholders also works with SQLite's older 999-variable limit.
                for part in batches(iter(batch), 500):
                    ids = [row["read"] for row in part]
                    query = "SELECT read_id,payload,matched FROM molecules WHERE read_id IN (" + ",".join("?" for _ in ids) + ")"
                    for rid, payload, present in db.execute(query, ids):
                        matched[rid] = (json.loads(payload), present)
                updates = []
                for row in batch:
                    n_metrics += 1
                    metric_classes[row["pc"]] += 1
                    if row["read"] not in matched:
                        continue
                    mol, present = matched[row["read"]]
                    if present or mol["metrics_present"]:
                        raise ValueError(f"Duplicate selected metric read: {row['read']}")
                    join_metric(mol, row, lib)
                    updates.append((encode(mol), row["read"]))
                db.executemany("UPDATE molecules SET payload=?,matched=1 WHERE read_id=?", updates)
                db.commit()
                if n_metrics % 1000000 < batch_size:
                    logging.info("Metrics: %s rows scanned", f"{n_metrics:,}")
        for source in sources:
            if source_info(source["path"]) != source:
                raise ValueError("An input changed during import")
        if source_info(catalog_path, True) != catalog_identity:
            raise ValueError("Catalog changed during import")
        missing = db.execute("SELECT count(*) FROM molecules WHERE matched=0").fetchone()[0]
        detailed, recovery = Counter(), {}
        dimensions = ["plasmid_class", "genotype_class", "ledger_structure", "genotype_tier", "pass_bin", "read_length_bin", "methylation_bin"]

        def export_rows():
            for exported, (payload,) in enumerate(db.execute("SELECT payload FROM molecules ORDER BY ordinal"), 1):
                mol = json.loads(payload)
                pc = mol.get("plasmid_class", "missing_metrics")
                cell = recovery.setdefault(pc, {"molecules": 0, "read_length_sum": 0, "passes_below_10": 0,
                                                "genotype_classes": Counter(), "ledger_structures": Counter(),
                                                "genotype_tiers": Counter(), "length_bins": Counter()})
                cell["molecules"] += 1
                cell["read_length_sum"] += mol["read_length"]
                cell["passes_below_10"] += mol["passes"] < 10
                for source_key, key in (("genotype_class", "genotype_classes"), ("ledger_structure", "ledger_structures"), ("genotype_tier", "genotype_tiers")):
                    cell[key][mol[source_key]] += 1
                length_group = length_bin(mol["read_length"], lib["reference_length"])
                cell["length_bins"][length_group] += 1
                methylation = str(min(19, int(mol["global_methylation_fraction"] * 20))) if mol["metrics_present"] else "missing"
                key = (pc, mol["genotype_class"], mol["ledger_structure"], mol["genotype_tier"], pass_bin(mol["passes"]), length_group, methylation)
                detailed[key] += 1
                if len(detailed) > 100000:
                    raise ValueError("Recovery summary exceeds 100000 categorical groups")
                if exported % 500000 == 0:
                    logging.info("Exported %s joined molecules", f"{exported:,}")
                yield mol

        write_table(out / "molecules.tsv.gz", export_rows(), FIELDS)
        recovery_rows = [dict(zip(dimensions + ["molecules"], list(key) + [count])) for key, count in sorted(detailed.items())]
        write_table(out / "recovery.tsv", recovery_rows, dimensions + ["molecules"])
        for cell in recovery.values():
            cell["mean_read_length"] = cell.pop("read_length_sum") / cell["molecules"]
        write_json(out / "recovery.json", {"sample_id": sample_id, "reference_length": lib["reference_length"],
                   "population": "selected genotype-ledger molecules, including missing metrics", "fraction": fraction,
                   "by_plasmid_class": recovery, "all_metrics_plasmid_classes": dict(metric_classes),
                   "methylation_bins": "0..19: floor(mf*20), with mf=1 in bin19; missing is separate",
                   "interpretation": "Imported structural/genotype classifications, not causal explanations or methylation-based exclusions."})
        database_bytes = db.execute("PRAGMA page_count").fetchone()[0] * db.execute("PRAGMA page_size").fetchone()[0]
        audit = {"schema_version": SCHEMA_VERSION, "software_version": __version__, "sample_id": sample_id,
                 "catalog": catalog_identity, "sources": sources,
                 "sampling": {"method": "blake2b-64 threshold on seed + NUL + read_id", "fraction": fraction,
                              "seed": seed, "max_molecules": None, "population": "entire primary-plasmid genotype ledger"},
                 "storage": {"backend": "sqlite", "cache_mib": cache_mib, "batch_size": batch_size,
                             "max_database_gib": max_disk_gib, "database_bytes": database_bytes,
                             "temporary_database_retained": False},
                 "ledger_rows": sum(all_classes.values()), "metrics_rows": n_metrics, "ledger_classes": dict(all_classes),
                 "selected_classes": dict(chosen_classes), "selected_reads": n_selected, "missing_metrics": missing,
                 "methylation_filter": None,
                 "limitations": ["Genotypes and structural classes are imported exploratory assignments.",
                                 "Missing metric rows are retained explicitly; methylation levels are not filtered.",
                                 "Large-input provenance uses size and mtime; small catalog files use SHA256."]}
        write_json(out / "audit.json", audit)
    return audit

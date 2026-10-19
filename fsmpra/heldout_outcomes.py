"""Build barcode outcomes for WT barcodes outside the frozen quality sample.

The joined molecule table is the only molecule input. All observations are
computed before averaging within barcode; no BAM or caller is opened here.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np

from .io import atomic_output, open_text, read_json, source_info, write_json
from .statistics import fixed_endpoints

SCHEMA = "fsmpra.heldout_outcomes.v1"
CONFIG_SCHEMA = "fsmpra.heldout_calibration_config.v1"
WT_BARCODES = 268283
WT_MOLECULES = 985264
SAMPLED_BARCODES = 10000
SAMPLED_MOLECULES = 37075
REQUIRED = {"read_id", "sample_id", "library_id", "corrected_barcode",
            "assigned_variant_id", "genotype_class", "footprint_status",
            "footprint_analysis_eligible_with_barcode", "eligible_ref_intervals0",
            "tf_small_ref_intervals0", "msp_ref_intervals0"}


def _path(value):
    return Path(value).resolve()


def _frozen_source(audit, path, checksum):
    current = source_info(path, checksum)
    if current not in audit.get("sources", []):
        raise ValueError(f"Frozen provenance disagrees: {path}")
    return current


def _mask(value, start, end, label):
    """Validate disjoint reference intervals and return their element bitmap."""
    try:
        intervals = json.loads(value)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid {label} intervals") from exc
    if not isinstance(intervals, list):
        raise ValueError(f"Invalid {label} intervals")
    result = 0
    previous = start
    for pair in intervals:
        if (not isinstance(pair, list) or len(pair) != 2 or
                any(type(x) is not int for x in pair)):
            raise ValueError(f"Invalid {label} interval")
        a, b = pair
        if not start <= a < b <= end or a < previous:
            raise ValueError(f"Overlapping or out-of-range {label} interval")
        result |= ((1 << (b-a)) - 1) << (a-start)
        previous = b
    return result


def _endpoint_observations(row, endpoints, start, end):
    aligned = _mask(row["eligible_ref_intervals0"], start, end, "aligned")
    if aligned != (1 << (end-start)) - 1:
        raise ValueError("Eligible WT molecule does not align across full element")
    tf = _mask(row["tf_small_ref_intervals0"], start, end, "TF")
    msp = _mask(row["msp_ref_intervals0"], start, end, "MSP")
    if tf & ~aligned or msp & ~aligned:
        raise ValueError("Call outside aligned reference bases")
    presence = np.empty(len(endpoints), dtype=np.float64)
    coverage = np.empty(len(endpoints), dtype=np.float64)
    for j, endpoint in enumerate(endpoints):
        a, b = endpoint["relative_start"], endpoint["relative_end"]
        window = ((1 << (b-a)) - 1) << a
        mask = tf if endpoint["track"] == "tf_small" else msp
        covered = bin(mask & window).count("1")
        presence[j] = float(covered > 0)
        coverage[j] = covered / (b-a)
    return presence, coverage


def _sample(sample_path, sample_audit, summary, wt_id):
    if (sample_audit.get("schema_version") != "fsmpra.footprint_quality.v1" or
            sample_audit.get("complete") is not True or
            sample_audit.get("barcode_limit") is not None or
            sample_audit.get("counts", {}).get("wt_barcodes") != WT_BARCODES or
            sample_audit["counts"].get("wt_eligible_molecules") != WT_MOLECULES or
            sample_audit.get("counts", {}).get("sampled_barcodes") != SAMPLED_BARCODES or
            sample_audit["counts"].get("sampled_expected_molecules") != SAMPLED_MOLECULES):
        raise ValueError("Expected complete frozen 10,000-WT-barcode extraction")
    if (summary.get("schema_version") != "fsmpra.footprint_quality_summary.v1" or
            summary.get("complete") is not True or
            summary.get("sampled_barcodes") != SAMPLED_BARCODES or
            summary.get("sampled_molecules") != SAMPLED_MOLECULES):
        raise ValueError("Expected complete frozen sample summary")
    if _path(summary["source_audit"]) != _path(sample_path.parent / "audit.json"):
        raise ValueError("Sample summary refers to different extraction audit")
    _frozen_source(summary, sample_path, True)
    _frozen_source(summary, sample_path.parent / "audit.json", True)
    result = {}
    with open_text(sample_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"barcode", "variant_id", "eligible_depth", "split", "weight"} <= set(reader.fieldnames or []):
            raise ValueError("Incomplete sampled barcode table")
        for row in reader:
            barcode = row["barcode"]
            depth = int(row["eligible_depth"])
            if (barcode in result or row["variant_id"] != wt_id or depth < 1 or
                    row["split"] not in ("A", "B") or
                    abs(float(row["weight"]) - 1 / depth) > 1e-12):
                raise ValueError("Invalid sampled WT barcode")
            result[barcode] = depth
    if len(result) != SAMPLED_BARCODES or sum(result.values()) != SAMPLED_MOLECULES:
        raise ValueError("Sampled WT barcode totals disagree with audit")
    return result


def _metadata(path, wt_id, excluded, molecules):
    wt = {}
    all_rows = 0
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"barcode_index", "barcode", "variant_id", "molecules"} <= set(reader.fieldnames or []):
            raise ValueError("Incomplete barcode metadata")
        for row in reader:
            index = int(row["barcode_index"])
            depth = int(row["molecules"])
            if index != all_rows or depth < 1 or depth != int(molecules[index]):
                raise ValueError("Barcode metadata/order/count disagreement")
            all_rows += 1
            if row["variant_id"] == wt_id:
                barcode = row["barcode"]
                if barcode in wt or len(barcode) != 15 or set(barcode) - set("ACGT"):
                    raise ValueError("Invalid/duplicate WT barcode metadata")
                wt[barcode] = (index, depth)
    if all_rows != len(molecules) or len(wt) != WT_BARCODES:
        raise ValueError("Barcode metadata cohort changed")
    if sum(depth for _, depth in wt.values()) != WT_MOLECULES:
        raise ValueError("WT molecule cohort changed")
    if not set(excluded) <= set(wt):
        raise ValueError("Sample contains barcode outside WT cohort")
    for barcode, depth in excluded.items():
        if wt[barcode][1] != depth:
            raise ValueError("Sampled WT depth disagrees with aggregate")
    return wt


def _reconcile(coverage, depths, indices, counts, endpoints):
    max_error = 0.0
    for first in range(0, len(indices), 16_384):
        last = min(first + 16_384, len(indices))
        # Read the old per-base count array once per bounded barcode batch.
        block = np.asarray(counts[indices[first:last]], dtype=np.uint64)
        cumulative = np.concatenate((np.zeros((last-first, 2, 1), dtype=np.uint64),
                                     np.cumsum(block, axis=2)), axis=2)
        for j, endpoint in enumerate(endpoints):
            track = 0 if endpoint["track"] == "tf_small" else 1
            a, b = endpoint["relative_start"], endpoint["relative_end"]
            expected = ((cumulative[:, track, b] - cumulative[:, track, a]) /
                        ((b-a) * depths[first:last]))
            error = float(np.max(np.abs(coverage[first:last, j] - expected)))
            max_error = max(max_error, error)
    if max_error > 1e-12:
        raise ValueError(f"Held-out coverage disagrees with aggregate counts: {max_error}")
    return max_error


def run(config, output=None, progress_every=1_000_000):
    config = Path(config)
    config_source = source_info(config, True)
    cfg = read_json(config)
    if cfg.get("schema_version") != CONFIG_SCHEMA:
        raise ValueError("Invalid held-out calibration config schema")
    for key in ("aggregates", "joined_dir", "excluded_sample", "catalog", "output_dir", "seed"):
        if key not in cfg:
            raise ValueError(f"Missing held-out config key: {key}")
    if type(cfg["seed"]) is not int or cfg["seed"] != 2026101802:
        raise ValueError("Held-out calibration seed changed")
    aggregates = _path(cfg["aggregates"])
    joined = _path(cfg["joined_dir"])
    sample = _path(cfg["excluded_sample"])
    if sample.is_dir():
        sample = sample / "sampled_barcodes.tsv"
    catalog_path = _path(cfg["catalog"])
    output = _path(output or Path(cfg["output_dir"]) / "outcomes")
    joined_table = joined / "molecule_footprints.tsv.gz"
    paths = [Path(config), catalog_path, aggregates / "audit.json",
             aggregates / "barcode_metadata.tsv.gz", aggregates / "barcode_counts.npy",
             aggregates / "barcode_molecules.npy", joined / "audit.json", joined_table,
             sample, sample.parent / "audit.json", sample.parent.parent / "summary" / "summary.json",
             Path(__file__), Path(__file__).with_name("statistics.py"),
             Path(__file__).with_name("io.py")]
    checksums = {str(Path(config)), str(catalog_path), str(aggregates / "audit.json"),
                 str(aggregates / "barcode_metadata.tsv.gz"), str(joined / "audit.json"),
                 str(sample), str(sample.parent / "audit.json"),
                 str(sample.parent.parent / "summary" / "summary.json"), str(Path(__file__)),
                 str(Path(__file__).with_name("statistics.py")),
                 str(Path(__file__).with_name("io.py"))}
    sources = [source_info(path, str(path) in checksums) for path in paths]
    if sources[0] != config_source:
        raise ValueError("Held-out config changed while loading")
    aggregate_audit = read_json(aggregates / "audit.json")
    join_audit = read_json(joined / "audit.json")
    sample_audit = read_json(sample.parent / "audit.json")
    summary = read_json(sample.parent.parent / "summary" / "summary.json")
    catalog = read_json(catalog_path)["library"]
    wt = [v for v in catalog["variants"] if v["class"] == "WT"]
    if len(wt) != 1 or len(wt[0]["sequence"]) != 211:
        raise ValueError("Expected one 211-base catalog WT")
    wt_id = wt[0]["variant_id"]
    start, end = catalog["mutable_ref_start0"], catalog["mutable_ref_end0"]
    if (end-start != 211 or aggregate_audit.get("schema_version") != "fsmpra.statistics_aggregate.v1" or
            aggregate_audit.get("interval") != [start, end] or
            aggregate_audit.get("counts", {}).get("barcodes") != 568354 or
            aggregate_audit["counts"].get("footprint_analysis_eligible_with_barcode") != 2102906 or
            join_audit.get("interval") != [start, end] or
            join_audit.get("counts", {}).get("output_molecules") != 8715827 or
            join_audit["counts"].get("footprint_analysis_eligible_with_barcode") != 2102906 or
            sample_audit.get("element") != [start, end] or
            summary.get("element") != [start, end]):
        raise ValueError("Expected complete existing join and aggregate cohort")
    if aggregate_audit["counts"]["output_molecules"] != join_audit["counts"]["output_molecules"]:
        raise ValueError("Join and aggregate molecule counts differ")
    _frozen_source(aggregate_audit, catalog_path, True)
    _frozen_source(aggregate_audit, joined / "audit.json", True)
    _frozen_source(aggregate_audit, joined_table, False)
    _frozen_source(sample_audit, catalog_path, True)
    _frozen_source(sample_audit, aggregates / "audit.json", True)
    _frozen_source(sample_audit, aggregates / "barcode_metadata.tsv.gz", True)
    excluded = _sample(sample, sample_audit, summary, wt_id)
    counts = np.load(aggregates / "barcode_counts.npy", mmap_mode="r", allow_pickle=False)
    molecules = np.load(aggregates / "barcode_molecules.npy", mmap_mode="r", allow_pickle=False)
    if counts.shape != (568354, 2, 211) or molecules.shape != (568354,):
        raise ValueError("Unexpected barcode count array dimensions")
    metadata = _metadata(aggregates / "barcode_metadata.tsv.gz", wt_id, excluded, molecules)
    barcodes = sorted(set(metadata) - set(excluded))
    if len(barcodes) != WT_BARCODES-SAMPLED_BARCODES:
        raise ValueError("Held-out WT barcode total changed")
    depths = np.array([metadata[b][1] for b in barcodes], dtype=np.int32)
    indices = np.array([metadata[b][0] for b in barcodes], dtype=np.int32)
    if int(depths.sum(dtype=np.int64)) != WT_MOLECULES-SAMPLED_MOLECULES:
        raise ValueError("Held-out WT molecule total changed")
    endpoints = fixed_endpoints(start, end)
    barcode_index = {barcode: i for i, barcode in enumerate(barcodes)}
    presence_sum = np.zeros((len(barcodes), len(endpoints)), dtype=np.float64)
    coverage_sum = np.zeros_like(presence_sum)
    valid = np.zeros((len(barcodes), len(endpoints)), dtype=np.int32)
    seen_ids = set()
    tally = Counter()
    with open_text(joined_table) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not REQUIRED <= set(reader.fieldnames or []):
            raise ValueError("Joined molecule table lacks required columns")
        for row in reader:
            tally["joined_rows"] += 1
            if progress_every and tally["joined_rows"] % progress_every == 0:
                print(f"Held-out outcomes: {tally['joined_rows']:,} joined rows", flush=True)
            if row["footprint_analysis_eligible_with_barcode"] not in ("0", "1"):
                raise ValueError("Invalid eligible-with-barcode flag")
            if row["footprint_analysis_eligible_with_barcode"] != "1" or row["assigned_variant_id"] != wt_id:
                continue
            barcode = row["corrected_barcode"]
            if barcode not in metadata or row["genotype_class"] != "WT":
                raise ValueError("Eligible WT row has invalid barcode/genotype")
            read_id = row["read_id"]
            if not read_id or read_id in seen_ids:
                raise ValueError("Duplicate or missing original eligible WT read ID")
            seen_ids.add(read_id)
            tally["wt_rows"] += 1
            if barcode in excluded:
                tally["excluded_rows"] += 1
                continue
            if row["footprint_status"] != "evaluated" or row["sample_id"] != join_audit["sample_id"] or row["library_id"] != join_audit["library_id"]:
                raise ValueError("Held-out eligible WT identity/status mismatch")
            i = barcode_index[barcode]
            p, c = _endpoint_observations(row, endpoints, start, end)
            presence_sum[i] += p
            coverage_sum[i] += c
            valid[i] += 1
            tally["heldout_rows"] += 1
    if (tally["joined_rows"] != join_audit["counts"]["output_molecules"] or
            tally["wt_rows"] != WT_MOLECULES or tally["excluded_rows"] != SAMPLED_MOLECULES or
            tally["heldout_rows"] != WT_MOLECULES-SAMPLED_MOLECULES or
            len(seen_ids) != WT_MOLECULES or np.any(valid != depths[:, None])):
        raise ValueError("Held-out WT molecule depth or row conservation failed")
    presence = presence_sum / valid
    coverage = coverage_sum / valid
    if (not np.all(np.isfinite(presence)) or not np.all(np.isfinite(coverage)) or
            np.any(presence < 0) or np.any(presence > 1) or
            np.any(coverage < 0) or np.any(coverage > 1)):
        raise ValueError("Invalid held-out endpoint values")
    max_error = _reconcile(coverage, depths, indices, counts, endpoints)
    with atomic_output(output) as out:
        np.savez_compressed(out / "endpoint_matrix.npz",
            barcodes=np.array(barcodes, dtype="U15"),
            endpoint_ids=np.array([x["endpoint_id"] for x in endpoints], dtype="U"),
            presence=presence, coverage=coverage, valid_molecules=valid, depths=depths,
            relative_starts=np.array([x["relative_start"] for x in endpoints], dtype=np.int16),
            relative_ends=np.array([x["relative_end"] for x in endpoints], dtype=np.int16),
            reference_starts0=np.array([x["reference_start0"] for x in endpoints], dtype=np.int32),
            reference_ends0=np.array([x["reference_end0"] for x in endpoints], dtype=np.int32))
        matrix_source = source_info(out / "endpoint_matrix.npz", True)
        matrix_source["path"] = str(output / "endpoint_matrix.npz")
        write_json(out / "audit.json", dict(schema_version=SCHEMA, complete=True,
            input_cohort={"wt_barcodes": WT_BARCODES, "wt_molecules": WT_MOLECULES,
                          "excluded_barcodes": SAMPLED_BARCODES, "excluded_molecules": SAMPLED_MOLECULES},
            heldout_barcodes=len(barcodes), heldout_molecules=tally["heldout_rows"],
            counts={"joined_rows": tally["joined_rows"], "wt_rows": tally["wt_rows"],
                    "excluded_rows": tally["excluded_rows"], "heldout_rows": tally["heldout_rows"],
                    "heldout_barcodes": len(barcodes)},
            joined_rows=tally["joined_rows"], endpoint_count=len(endpoints),
            coverage_reconciliation_max_abs_error=max_error,
            seed=cfg["seed"], sources=sources + [matrix_source],
            interpretation="WT barcodes outside the prior 10,000-barcode audit sample, not globally untouched. One observation per barcode; original eligible molecules only."))
        for old, path in zip(sources, paths):
            if source_info(path, "sha256" in old) != old:
                raise ValueError(f"Input changed during held-out extraction: {path}")
    return {"heldout_barcodes": len(barcodes), "heldout_molecules": tally["heldout_rows"],
            "coverage_reconciliation_max_abs_error": max_error}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    print(json.dumps(run(args.config, args.output), sort_keys=True))


if __name__ == "__main__":
    main()

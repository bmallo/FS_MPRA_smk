#!/usr/bin/env python3
"""Independently reconcile the frozen WT sample with earlier barcode arrays.

This validator streams metadata and molecules, and reads only selected rows of
the old per-base count array. It does not import the extraction or summary code.
"""

import argparse
import bisect
import csv
import gzip
import hashlib
import heapq
import json
import os
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np


SEED = 20261018
SAMPLE_N = 10_000
TOLERANCE = 1e-12


def _json(path):
    with Path(path).open() as handle:
        return json.load(handle)


def _rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty table: {path}")
        yield from reader


def _source(path, checksum):
    path = Path(path).resolve()
    stat = path.stat()
    record = {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if checksum:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        record["sha256"] = digest.hexdigest()
    return record


def _check_sources(sources):
    for record in sources:
        if _source(record["path"], "sha256" in record) != record:
            raise ValueError("Input changed during validation: " + record["path"])


def _require(record, fields, name):
    missing = set(fields) - set(record)
    if missing:
        raise ValueError(f"{name} missing columns: {sorted(missing)}")


def _hash(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _selected_metadata(path, wt_id, sample_n=SAMPLE_N, seed=SEED):
    """Independently rank every WT row, retaining only the selected 10,000."""
    heap = []
    seen = set()
    total = 0
    total_depth = 0
    metadata_rows = 0
    for row in _rows(path):
        if metadata_rows == 0:
            _require(row, ("barcode_index", "barcode", "variant_id", "molecules"), "old metadata")
        index = int(row["barcode_index"])
        if index != metadata_rows:
            raise ValueError("Old barcode indices are not contiguous")
        metadata_rows += 1
        if row["variant_id"] != wt_id:
            continue
        barcode = row["barcode"]
        depth = int(row["molecules"])
        if not barcode or barcode in seen or depth <= 0:
            raise ValueError("Invalid or duplicate WT metadata barcode")
        seen.add(barcode)
        total += 1
        total_depth += depth
        score = _hash(f"{seed}|{barcode}")
        # Worst selected hash sorts first in the min-heap.
        candidate = (-int(score, 16), barcode, index, depth, score)
        if len(heap) < sample_n:
            heapq.heappush(heap, candidate)
        elif candidate > heap[0]:
            heapq.heapreplace(heap, candidate)
    if total < sample_n:
        raise ValueError("Too few WT barcodes in old metadata")
    selected = {barcode: (index, depth, score) for _, barcode, index, depth, score in heap}
    return selected, dict(metadata_rows=metadata_rows, wt_barcodes=total,
                          wt_eligible_molecules=total_depth,
                          sampled_expected_molecules=sum(item[1] for item in selected.values()))


def _offsets(value, name, required=True):
    if value == "" and not required:
        return None
    values = json.loads(value)
    if (not isinstance(values, list) or
            any(type(pos) is not int or not 0 <= pos < 211 for pos in values) or
            values != sorted(set(values))):
        raise ValueError("Invalid " + name + " bitmap")
    return values


def _windows(start):
    begins = list(range(0, 192, 10))
    if begins[-1] != 191:
        begins.append(191)
    result = [(f"tf_{lo:03d}_{lo+20:03d}", "tf_small", lo, lo+20, start+lo, start+lo+20)
              for lo in begins]
    result.append(("msp_promoter", "msp", 0, 211, start, start+211))
    return result


def _window_count(sorted_positions, lo, hi):
    return bisect.bisect_left(sorted_positions, hi) - bisect.bisect_left(sorted_positions, lo)


def _compare(actual, expected, label):
    if actual.shape != expected.shape or not np.all(np.isfinite(actual)) or not np.all(np.isfinite(expected)):
        raise ValueError(label + " shape or finite-value mismatch")
    largest = float(np.max(np.abs(actual - expected)))
    if largest > TOLERANCE:
        raise ValueError(f"{label} maximum absolute error {largest:.6g} exceeds {TOLERANCE}")
    return largest


def validate(sample, summary, aggregates):
    sample, summary, aggregates = map(Path, (sample, summary, aggregates))
    audit_path = sample / "audit.json"
    summary_path = summary / "summary.json"
    barcodes_path = sample / "sampled_barcodes.tsv"
    molecules_path = sample / "sampled_molecules.tsv.gz"
    matrix_path = summary / "barcode_endpoint_matrix.npz"
    metadata_path = aggregates / "barcode_metadata.tsv.gz"
    counts_path = aggregates / "barcode_counts.npy"
    depths_path = aggregates / "barcode_molecules.npy"
    aggregate_audit_path = aggregates / "audit.json"
    source_paths = [(audit_path, True), (summary_path, True), (barcodes_path, True),
                    (molecules_path, False), (matrix_path, True), (metadata_path, True),
                    (counts_path, False), (depths_path, False), (aggregate_audit_path, True),
                    (Path(__file__), True)]
    sources = [_source(path, checksum) for path, checksum in source_paths]
    audit, report, aggregate_audit = map(_json, (audit_path, summary_path, aggregate_audit_path))
    if (audit.get("schema_version") != "fsmpra.footprint_quality.v1" or audit.get("complete") is not True or
            audit.get("seed") != SEED or audit.get("counts", {}).get("sampled_barcodes") != SAMPLE_N):
        raise ValueError("Expected complete frozen 10,000-barcode sample")
    if (report.get("schema_version") != "fsmpra.footprint_quality_summary.v1" or
            report.get("complete") is not True or report.get("pilot") is not False or
            report.get("sampled_barcodes") != SAMPLE_N or report.get("endpoint_count") != 22 or
            report.get("endpoint_matrix") != matrix_path.name):
        raise ValueError("Expected complete 22-endpoint summary")
    if (aggregate_audit.get("schema_version") != "fsmpra.statistics_aggregate.v1" or
            aggregate_audit.get("sample_id") != report.get("sample_id") or
            aggregate_audit.get("tracks") != ["tf_small", "msp"]):
        raise ValueError("Old aggregate contract differs from summary")
    start, end = audit.get("element", [None, None])
    if (start, end) != (1234, 1445) or aggregate_audit.get("interval") != [start, end]:
        raise ValueError("Element coordinates disagree")
    source_records = {item["path"]: item for item in audit.get("sources", [])}
    if source_records.get(str(metadata_path.resolve())) != _source(metadata_path, True):
        raise ValueError("Sample was selected from different old barcode metadata")
    summary_sources = {item["path"]: item for item in report.get("sources", [])}
    for path in (audit_path, barcodes_path, molecules_path, matrix_path):
        record = summary_sources.get(str(path.resolve()))
        if record is None or record != _source(path, "sha256" in record):
            raise ValueError("Summary source provenance mismatch: " + str(path))
    sampled_rows = {}
    for row in _rows(barcodes_path):
        if not sampled_rows:
            _require(row, ("barcode", "variant_id", "eligible_depth", "split", "weight", "selection_hash"),
                     "sampled barcode")
        barcode = row["barcode"]
        if not barcode or barcode in sampled_rows:
            raise ValueError("Invalid or duplicate sampled barcode")
        sampled_rows[barcode] = row
    if len(sampled_rows) != SAMPLE_N:
        raise ValueError("Sampled barcode table does not contain 10,000 unique rows")
    wt_ids = {row["variant_id"] for row in sampled_rows.values()}
    if len(wt_ids) != 1:
        raise ValueError("Sample contains multiple variants")
    wt_id = next(iter(wt_ids))
    selected, selection = _selected_metadata(metadata_path, wt_id, sample_n=SAMPLE_N)
    if set(sampled_rows) != set(selected):
        raise ValueError("Sampled barcodes differ from lowest SHA256 WT selection")
    if selection["metadata_rows"] != aggregate_audit.get("counts", {}).get("barcodes"):
        raise ValueError("Old aggregate barcode count differs from metadata")
    for key in ("wt_barcodes", "wt_eligible_molecules", "sampled_expected_molecules"):
        if selection[key] != audit["counts"].get(key):
            raise ValueError("WT selection count differs from extraction audit: " + key)
    sorted_barcodes = sorted(sampled_rows)
    expected_depths = np.empty(SAMPLE_N, dtype=np.int64)
    old_indices = np.empty(SAMPLE_N, dtype=np.int64)
    expected_splits = []
    for i, barcode in enumerate(sorted_barcodes):
        row = sampled_rows[barcode]
        old_index, depth, score = selected[barcode]
        split = "A" if int(_hash(f"split|{SEED}|{barcode}"), 16) % 2 == 0 else "B"
        if (int(row["eligible_depth"]) != depth or row["selection_hash"] != score or
                row["split"] != split or abs(float(row["weight"]) - 1/depth) > TOLERANCE):
            raise ValueError("Sample barcode depth, selection hash, split, or weight differs: " + barcode)
        expected_depths[i], old_indices[i] = depth, old_index
        expected_splits.append(split)
    split_counts = dict(Counter(expected_splits))
    if split_counts != report.get("split_barcodes"):
        raise ValueError("Summary A/B split counts differ from independent parity")
    counts = np.load(counts_path, mmap_mode="r", allow_pickle=False)
    old_depths = np.load(depths_path, mmap_mode="r", allow_pickle=False)
    if (counts.shape != (selection["metadata_rows"], 2, 211) or
            old_depths.shape != (selection["metadata_rows"],) or
            not np.issubdtype(counts.dtype, np.unsignedinteger) or
            not np.issubdtype(old_depths.dtype, np.unsignedinteger)):
        raise ValueError("Unexpected old aggregate array shape or dtype")
    if not np.array_equal(old_depths[old_indices], expected_depths):
        raise ValueError("Old aggregate depths differ from sampled barcode metadata")
    windows = _windows(start)
    with np.load(matrix_path, allow_pickle=False) as matrix:
        required = {"barcodes", "endpoint_ids", "relative_starts", "relative_ends",
                    "reference_starts0", "reference_ends0", "depths", "splits",
                    "valid_molecules", "presence", "coverage"}
        if not required <= set(matrix.files):
            raise ValueError("Summary matrix missing endpoint arrays")
        if list(matrix["barcodes"]) != sorted_barcodes or list(matrix["splits"]) != expected_splits:
            raise ValueError("Matrix barcode order or A/B split differs from sampled table")
        for key, column in (("endpoint_ids", 0), ("relative_starts", 2), ("relative_ends", 3),
                            ("reference_starts0", 4), ("reference_ends0", 5)):
            if list(matrix[key]) != [window[column] for window in windows]:
                raise ValueError("Matrix endpoint grid differs: " + key)
        matrix_depths = matrix["depths"]
        valid = matrix["valid_molecules"]
        presence = matrix["presence"]
        coverage = matrix["coverage"]
    if (matrix_depths.shape != (SAMPLE_N,) or valid.shape != (SAMPLE_N, 22) or
            presence.shape != (SAMPLE_N, 22) or coverage.shape != (SAMPLE_N, 22)):
        raise ValueError("Summary matrix array shapes differ")
    if not np.array_equal(matrix_depths, expected_depths):
        raise ValueError("Matrix barcode depths differ from sampled metadata")
    # Read one old barcode slice at a time; no full (568k,2,211) copy.
    old_expected = np.empty((SAMPLE_N, 22), dtype=float)
    for i, (index, depth) in enumerate(zip(old_indices, expected_depths)):
        perbase = counts[index]
        for j, (_, track, lo, hi, _, _) in enumerate(windows):
            old_expected[i, j] = int(perbase[0 if track == "tf_small" else 1, lo:hi].sum(dtype=np.uint64)) / (int(depth) * (hi-lo))
    old_error = _compare(coverage, old_expected, "old aggregate vs summary coverage")
    barcode_to_index = {barcode: i for i, barcode in enumerate(sorted_barcodes)}
    molecule_depths = np.zeros(SAMPLE_N, dtype=np.int64)
    molecule_valid = np.zeros((SAMPLE_N, 22), dtype=np.int32)
    presence_sum = np.zeros((SAMPLE_N, 22), dtype=np.int64)
    coverage_sum = np.zeros((SAMPLE_N, 22), dtype=float)
    molecule_count = 0
    for row in _rows(molecules_path):
        if molecule_count == 0:
            _require(row, ("corrected_barcode", "assigned_variant_id", "sample_split",
                           "element_aligned_offsets0", "element_tf_small_offsets0", "element_msp_offsets0",
                           "element_tf_small_fraction", "element_msp_fraction"), "sampled molecule")
        barcode = row["corrected_barcode"]
        i = barcode_to_index.get(barcode)
        if i is None or row["assigned_variant_id"] != wt_id or row["sample_split"] != expected_splits[i]:
            raise ValueError("Molecule barcode, WT identity, or split differs")
        aligned = _offsets(row["element_aligned_offsets0"], "aligned")
        tf = _offsets(row["element_tf_small_offsets0"], "TF", required=False)
        msp = _offsets(row["element_msp_offsets0"], "MSP", required=False)
        if aligned != list(range(211)) or tf is None or msp is None or not set(tf).issubset(aligned) or not set(msp).issubset(aligned):
            raise ValueError("Old cohort molecule lacks full aligned/callable element")
        if row["element_tf_small_fraction"] == "" or row["element_msp_fraction"] == "":
            raise ValueError("Old cohort molecule has missing TF/MSP measurement")
        if abs(float(row["element_tf_small_fraction"]) - len(tf)/211) > TOLERANCE or abs(float(row["element_msp_fraction"]) - len(msp)/211) > TOLERANCE:
            raise ValueError("Sampled molecule fraction differs from bitmap")
        molecule_depths[i] += 1
        molecule_count += 1
        for j, (_, track, lo, hi, _, _) in enumerate(windows):
            positions = tf if track == "tf_small" else msp
            denominator = _window_count(aligned, lo, hi)
            covered = _window_count(positions, lo, hi)
            molecule_valid[i, j] += 1
            presence_sum[i, j] += covered > 0
            coverage_sum[i, j] += covered / denominator
    if not np.array_equal(molecule_depths, expected_depths) or molecule_count != selection["sampled_expected_molecules"]:
        raise ValueError("Sampled molecule depth conservation failed")
    if (molecule_count != audit["counts"].get("matched_molecules") or
            molecule_count != report.get("sampled_molecules") or
            not np.array_equal(molecule_valid, valid) or not np.all(valid == expected_depths[:, None])):
        raise ValueError("Molecule counts or endpoint denominators differ")
    expected_presence = presence_sum / molecule_valid
    expected_coverage = coverage_sum / molecule_valid
    presence_error = _compare(presence, expected_presence, "bitmap vs summary presence")
    coverage_error = _compare(coverage, expected_coverage, "bitmap vs summary coverage")
    _check_sources(sources)
    return {"schema_version": "fsmpra.footprint_quality_validation.v1", "status": "passed",
            "sample_id": report["sample_id"], "sampled_barcodes": SAMPLE_N,
            "sampled_molecules": molecule_count, "wt_barcodes_ranked": selection["wt_barcodes"],
            "wt_eligible_molecules": selection["wt_eligible_molecules"],
            "split_barcodes": split_counts, "endpoints": len(windows),
            "max_abs_coverage_error_vs_old_aggregate": old_error,
            "max_abs_presence_error_vs_molecule_bitmaps": presence_error,
            "max_abs_coverage_error_vs_molecule_bitmaps": coverage_error,
            "full_element_aligned_molecules": molecule_count,
            "sources": sources,
            "interpretation": "Independent consistency check of frozen sampled WT barcodes, original aggregate coverage, and sampled molecule endpoint fractions."}


def _write_atomic(path, report):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
    fd, temporary = tempfile.mkstemp(prefix="." + path.name + "-", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(report, handle, indent=2, allow_nan=False)
            handle.write("\n")
        _check_sources(report["sources"])
        os.link(temporary, path)
    finally:
        os.unlink(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--aggregates", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = validate(args.sample, args.summary, args.aggregates)
    _write_atomic(args.output, report)
    print(json.dumps({key: report[key] for key in ("status", "sampled_barcodes", "sampled_molecules",
                                                 "max_abs_coverage_error_vs_old_aggregate")}, sort_keys=True))


if __name__ == "__main__":
    main()

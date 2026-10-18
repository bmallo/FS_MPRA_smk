"""Stream replicate-3 Fiber-seq heterogeneity features from a frozen joined BAM.

This module does not revise genotype, quality control, or barcode assignments.
Reference intervals are zero-based, half-open. Missing MM/ML or MA yields blank
outcome features; sparse but present MA tracks are evaluated as empty.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from .footprint_join import merge_intervals, region_blocks, region_calls
from .io import atomic_output, read_json, source_info, write_json, write_table
from .processing import check_sources

META_FIELDS = [
    "read_id", "sample_id", "library_id", "observed_variant_id", "genotype_class",
    "genotype_status", "corrected_barcode", "assigned_variant_id",
    "barcode_assignment_status", "passes", "read_quality", "mapq", "read_length",
    "strand", "plasmid_class", "circular", "m6a_status", "at_count", "m6a_count",
    "methylation_fraction", "sequence_qc", "length_qc", "methylation_qc",
    "analysis_eligible", "footprint_status", "region_complete", "tf_small_fraction",
    "nucleosome_fraction", "msp_fraction", "footprint_analysis_eligible",
    "footprint_analysis_eligible_with_barcode",
]
REGION_SUFFIXES = [
    "aligned_bases", "coverage_fraction", "at_count", "m6a_count", "m6a_fraction",
    "m6a_count_204", "m6a_fraction_204", "tf_small_fraction",
    "nucleosome_fraction", "msp_fraction",
]
EXTRA_FIELDS = [
    "whole_nucleosome_segments", "whole_nucleosome_length_median",
    "whole_nucleosome_spacing_median", "whole_tf_calls_all",
    "whole_tf_calls_small", "distal_tile_status", "distal_tiles",
    "distal_tile_methylation_q10", "distal_tile_methylation_median",
    "distal_tile_methylation_q90",
]
FEATURE_FIELDS = META_FIELDS + [f"{name}_{suffix}" for name in ("whole", "element", "distal")
                              for suffix in REGION_SUFFIXES] + EXTRA_FIELDS
RETENTION_FIELDS = [
    "observed_variant_id", "genotype_class", "genotype_status", "sequence_qc",
    "length_qc", "methylation_qc", "barcode_assignment_status", "total_primary",
    "sequence_qualified", "analysis_eligible", "footprint_analysis_eligible",
    "footprint_analysis_eligible_with_barcode", "failed_sequence_qc", "failed_length_qc",
    "failed_methylation_qc",
]
MM_HEADER = re.compile(r"([ACGTUN])([+-])([a-z]+|[0-9]+)([.?]?)\Z")


def selected(row):
    return (row["genotype_status"] == "exact_design" and row["sequence_qc"] == "1"
            and row["genotype_class"] in ("WT", "SNV"))


def decode_m6a(read):
    """Return status and stored-query positive positions at ML 128 and 204.

    Decode all MM groups against original forward SEQ, consuming every ML byte.
    Repeated target positions use the largest ML probability. A '?' on either
    target channel makes fractions unknown; observed positive calls remain usable.
    """
    seq = read.query_sequence or ""
    n = len(seq)
    at = seq.count("A") + seq.count("T")
    empty = np.empty(0, dtype=np.int32)
    if not read.has_tag("MM") or not read.has_tag("ML"):
        return "missing", at, empty, empty
    if read.has_tag("MN") and read.get_tag("MN") != n:
        raise ValueError("MN/SEQ length mismatch: " + read.query_name)
    forward = read.get_forward_sequence() or ""
    if len(forward) != n:
        raise ValueError("Forward SEQ length mismatch: " + read.query_name)
    forward_bytes = np.frombuffer(forward.encode("ascii"), dtype=np.uint8)
    ml_tag = read.get_tag("ML")
    if getattr(ml_tag, "typecode", None) != "B":
        raise ValueError("Expected unsigned byte ML array: " + read.query_name)
    ml = np.frombuffer(ml_tag, dtype=np.uint8)
    scores = np.zeros(n, dtype=np.uint8)
    cursor = 0
    target_seen = unknown = False
    for section in read.get_tag("MM").split(";"):
        if not section:
            continue
        header, *tokens = section.split(",")
        match = MM_HEADER.fullmatch(header)
        if match is None:
            raise ValueError("Unsupported MM group: " + header)
        base, strand, codes, omitted = match.groups()
        modifications = [codes] if codes.isdigit() else list(codes)
        count = len(tokens)
        end = cursor + count * len(modifications)
        if end > len(ml):
            raise ValueError("ML shorter than MM: " + read.query_name)
        try:
            skips = np.asarray([int(token) for token in tokens], dtype=np.int64)
        except ValueError as exc:
            raise ValueError("Invalid MM delta: " + read.query_name) from exc
        if np.any(skips < 0):
            raise ValueError("Negative MM delta: " + read.query_name)
        opportunities = n if base == "N" else forward.count(base)
        if count and int(skips.sum()) + count > opportunities:
            raise ValueError("MM/SEQ length disagreement: " + read.query_name)
        target_columns = [j for j, code in enumerate(modifications)
                          if (base, strand, code) in (("A", "+", "a"), ("T", "-", "a"))]
        if target_columns:
            target_seen = True
            unknown |= omitted == "?"
            if count:
                opportunity_positions = (np.arange(n, dtype=np.int32) if base == "N"
                                         else np.flatnonzero(forward_bytes == ord(base)).astype(np.int32))
                indices = np.cumsum(skips + 1) - 1
                forward_positions = opportunity_positions[indices]
                positions = n - 1 - forward_positions if read.is_reverse else forward_positions
                probabilities = ml[cursor:end].reshape(count, len(modifications))
                for j in target_columns:
                    np.maximum.at(scores, positions, probabilities[:, j])
        cursor = end
    if cursor != len(ml):
        raise ValueError("Extra ML probabilities: " + read.query_name)
    if not target_seen:
        return "missing_m6a", at, empty, empty
    status = "unknown_opportunities" if unknown or not at else "ok"
    return status, at, np.flatnonzero(scores >= 128), np.flatnonzero(scores >= 204)


def _length(intervals):
    return sum(b - a for a, b in intervals)


def _clip(spans, intervals):
    return merge_intervals([[max(a, lo), min(b, hi)] for a, b in spans
                            for lo, hi in intervals if max(a, lo) < min(b, hi)])


def _fraction(numerator, denominator):
    return numerator / denominator if denominator else ""


def _region_boolean(ref_positions, intervals):
    mask = np.zeros(len(ref_positions), dtype=bool)
    for lo, hi in intervals:
        mask |= (ref_positions >= lo) & (ref_positions < hi)
    return mask


def extract_features(read, row, spec):
    """Extract one selected molecule, independent of cohort/QC classification."""
    regions = spec["regions"]
    ref_length = int(spec["reference_length"])
    seq = read.query_sequence or ""
    qlen = len(seq)
    blocks = region_blocks(read, 0, ref_length)
    qref = np.full(qlen, -1, dtype=np.int32)
    for qa, qb, ra in blocks:
        qref[qa:qb] = np.arange(ra, ra + qb - qa, dtype=np.int32)
    bases = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    at_mask = (bases == ord("A")) | (bases == ord("T"))
    methyl_status, global_at, pos128, pos204 = decode_m6a(read)
    if row["m6a_status"] == "ok":
        if methyl_status != "ok" or global_at != int(row["at_count"]) or len(pos128) != int(row["m6a_count"]):
            raise ValueError("Global m6A count/AT disagreement: " + read.query_name)
    elif row["m6a_status"] != methyl_status:
        raise ValueError("m6A status disagreement: " + read.query_name)
    methyl_mask128 = np.zeros(qlen, dtype=bool)
    methyl_mask204 = np.zeros(qlen, dtype=bool)
    methyl_mask128[pos128] = True
    methyl_mask204[pos204] = True
    if np.any(methyl_mask128 & ~at_mask) or np.any(methyl_mask204 & ~at_mask):
        raise ValueError("m6A call not on stored A/T: " + read.query_name)
    has_ma = read.has_tag("MA")
    if has_ma:
        calls, _ = region_calls(read, blocks)
        if row["footprint_status"] != "evaluated":
            raise ValueError("MA status disagrees with joined table: " + read.query_name)
    else:
        calls = None
        if row["footprint_status"] != "missing_MA":
            raise ValueError("MA status disagrees with joined table: " + read.query_name)
    result = {name: row[name] for name in META_FIELDS}
    for name, intervals in regions.items():
        prefix = name + "_"
        region_mask = _region_boolean(qref, intervals)
        aligned = int(region_mask.sum())
        target_length = _length(intervals)
        if aligned > target_length:
            raise ValueError("Regional aligned bases exceed interval length: " + read.query_name)
        at_here = int(np.count_nonzero(region_mask & at_mask))
        result[prefix + "aligned_bases"] = aligned
        result[prefix + "coverage_fraction"] = aligned / target_length
        result[prefix + "at_count"] = at_here
        if methyl_status in ("ok", "unknown_opportunities"):
            n128 = int(np.count_nonzero(region_mask & methyl_mask128))
            n204 = int(np.count_nonzero(region_mask & methyl_mask204))
            if n204 > n128 or n128 > at_here:
                raise ValueError("Regional m6A count exceeds aligned A/T: " + read.query_name)
            result[prefix + "m6a_count"] = n128
            result[prefix + "m6a_count_204"] = n204
            result[prefix + "m6a_fraction"] = _fraction(n128, at_here) if methyl_status == "ok" else ""
            result[prefix + "m6a_fraction_204"] = _fraction(n204, at_here) if methyl_status == "ok" else ""
        else:
            for suffix in ("m6a_count", "m6a_count_204", "m6a_fraction", "m6a_fraction_204"):
                result[prefix + suffix] = ""
        for label, track, minimum, maximum in (("tf_small", "tf", spec["tf_min"], spec["tf_max"]),
                                                ("nucleosome", "nuc", spec["nuc_min"], None),
                                                ("msp", "msp", 1, None)):
            if calls is None:
                result[prefix + label + "_fraction"] = ""
                continue
            spans = _clip([span for call in calls[track]
                           if call["query_length"] >= minimum
                           and (maximum is None or call["query_length"] <= maximum)
                           for span in call["ref_intervals0"]], intervals)
            numerator = _length(spans)
            if numerator > aligned:
                raise ValueError("Footprint union exceeds aligned bases: " + read.query_name)
            result[prefix + label + "_fraction"] = _fraction(numerator, aligned)
    if calls is None:
        for field in EXTRA_FIELDS[:5]:
            result[field] = ""
    else:
        nuc_calls = [call for call in calls["nuc"]
                     if call["query_length"] >= spec["nuc_min"]]
        # An indel can split one MA call into several projected spans. Count
        # the call once, using its outer projected bounds for the midpoint.
        midpoints = sorted((call["ref_intervals0"][0][0] +
                            call["ref_intervals0"][-1][1]) / 2 for call in nuc_calls)
        result["whole_nucleosome_segments"] = len(nuc_calls)
        result["whole_nucleosome_length_median"] = (
            float(np.median([call["query_length"] for call in nuc_calls])) if nuc_calls else "")
        result["whole_nucleosome_spacing_median"] = (
            float(np.median(np.diff(midpoints))) if len(midpoints) > 1 else "")
        result["whole_tf_calls_all"] = len(calls["tf"])
        result["whole_tf_calls_small"] = sum(spec["tf_min"] <= call["query_length"] <= spec["tf_max"]
                                                  for call in calls["tf"])
    if methyl_status != "ok":
        result["distal_tile_status"] = methyl_status
        result["distal_tiles"] = ""
        for field in EXTRA_FIELDS[7:]:
            result[field] = ""
    else:
        tile_size = spec["tile_size"]
        min_at = spec["tile_min_at"]
        tile_ids = spec.get("distal_tile_ids")
        if tile_ids is None:
            tile_ids = [tile_start // tile_size for lo, hi in regions["distal"]
                        for tile_start in range(((lo + tile_size - 1) // tile_size) * tile_size,
                                                hi - tile_size + 1, tile_size)]
        n_bins = (ref_length + tile_size - 1) // tile_size
        aligned_at = at_mask & (qref >= 0)
        at_by_tile = np.bincount(qref[aligned_at] // tile_size, minlength=n_bins)
        methyl_by_tile = np.bincount(qref[methyl_mask128 & (qref >= 0)] // tile_size,
                                     minlength=n_bins)
        valid_tiles = np.asarray(tile_ids, dtype=np.int32)
        valid_tiles = valid_tiles[at_by_tile[valid_tiles] >= min_at]
        fractions = methyl_by_tile[valid_tiles] / at_by_tile[valid_tiles]
        result["distal_tiles"] = len(fractions)
        result["distal_tile_status"] = "evaluated" if len(fractions) else "no_qualifying_tiles"
        if len(fractions):
            q10, median, q90 = np.quantile(fractions, [0.1, 0.5, 0.9])
            result["distal_tile_methylation_q10"] = float(q10)
            result["distal_tile_methylation_median"] = float(median)
            result["distal_tile_methylation_q90"] = float(q90)
        else:
            for field in EXTRA_FIELDS[7:]:
                result[field] = ""
    return result


def _read_config(path):
    cfg = read_json(path)
    if cfg.get("schema_version") != "fsmpra.heterogeneity_config.v1":
        raise ValueError("Unsupported heterogeneity configuration")
    regions = cfg["regions"]
    if set(regions) != {"whole", "element", "distal"}:
        raise ValueError("Expected whole, element, distal regions")
    catalog = read_json(Path(cfg["ledger_run_dir"]) / "catalog" / "catalog.json")
    reference_length = catalog["library"]["reference_length"]
    normalized = {}
    for name, intervals in regions.items():
        items = sorted((int(lo), int(hi)) for lo, hi in intervals)
        if not items or any(lo < 0 or lo >= hi or hi > reference_length for lo, hi in items):
            raise ValueError("Invalid reference interval for " + name)
        if any(items[i][1] > items[i + 1][0] for i in range(len(items) - 1)):
            raise ValueError("Overlapping intervals for " + name)
        normalized[name] = items
    if normalized["whole"] != [(0, reference_length)]:
        raise ValueError("Whole region must cover the complete plasmid")
    if normalized["element"] != [(catalog["library"]["mutable_ref_start0"],
                                  catalog["library"]["mutable_ref_end0"])]:
        raise ValueError("Element region disagrees with catalog")
    if any(max(a, c) < min(b, d) for a, b in normalized["distal"]
           for c, d in normalized["element"]):
        raise ValueError("Distal and element regions overlap")
    excluded = merge_intervals([
        [max(0, normalized["element"][0][0] - 100),
         min(reference_length, normalized["element"][0][1] + 100)],
        [max(0, catalog["library"]["barcode_start0"] - 25),
         min(reference_length, catalog["library"]["barcode_end0"] + 25)],
    ])
    expected_distal, cursor = [], 0
    for lo, hi in excluded:
        if cursor < lo:
            expected_distal.append((cursor, lo))
        cursor = hi
    if cursor < reference_length:
        expected_distal.append((cursor, reference_length))
    if normalized["distal"] != expected_distal:
        raise ValueError("Distal intervals disagree with element/barcode flanks")
    cfg["regions"] = normalized
    cfg["reference_length"] = reference_length
    cfg["tile_size"] = int(cfg.get("tile_size", 100))
    cfg["tile_min_at"] = int(cfg.get("tile_min_at", 20))
    cfg["tf_min"] = int(cfg.get("tf_min", 10))
    cfg["tf_max"] = int(cfg.get("tf_max", 60))
    cfg["nuc_min"] = int(cfg.get("nuc_min", 85))
    if cfg.get("thresholds", [128, 204]) != [128, 204] or cfg["tile_size"] < 1 or cfg["tile_min_at"] < 1:
        raise ValueError("Invalid methylation tile parameters")
    if (cfg["tf_min"], cfg["tf_max"], cfg["nuc_min"]) != (10, 60, 85):
        raise ValueError("Footprint thresholds disagree with established region features")
    cfg["distal_tile_ids"] = [tile_start // cfg["tile_size"]
                              for lo, hi in normalized["distal"]
                              for tile_start in range(((lo + cfg["tile_size"] - 1) // cfg["tile_size"])
                                                      * cfg["tile_size"], hi - cfg["tile_size"] + 1,
                                                      cfg["tile_size"])]
    return cfg, catalog


def extract(config, output, limit_primary=None, progress_every=100000):
    """Publish selected features and retention audit after a strict ordered scan."""
    import pysam

    cfg, catalog = _read_config(config)
    run_dir = Path(cfg["ledger_run_dir"])
    footprint_dir = Path(cfg["footprint_run_dir"])
    joined = footprint_dir / "molecule_footprints.tsv.gz"
    join_audit_path = footprint_dir / "audit.json"
    join_audit = read_json(join_audit_path)
    contract_path = Path(cfg["footprint_contract"])
    contract = read_json(contract_path)
    bam_path = Path(contract["bam"])
    if join_audit.get("schema_version") != "fsmpra.footprint_join.v1":
        raise ValueError("Expected complete existing footprint join")
    if join_audit["footprint_contract"]["bam"] != str(bam_path):
        raise ValueError("Footprint contract/BAM disagreement")
    frozen_contract = next((s for s in join_audit["sources"]
                            if s["path"] == str(contract_path.resolve())), None)
    if frozen_contract is None or source_info(contract_path, True) != frozen_contract:
        raise ValueError("Footprint contract differs from frozen join source")
    frozen = next((s for s in join_audit["sources"] if s["path"] == str(bam_path.resolve())), None)
    if frozen is None or source_info(bam_path) != frozen:
        raise ValueError("BAM differs from frozen footprint join source")
    bai_path = Path(str(bam_path) + ".bai")
    frozen_bai = next((s for s in join_audit["sources"] if s["path"] == str(bai_path.resolve())), None)
    if frozen_bai is None or source_info(bai_path) != frozen_bai:
        raise ValueError("BAM index differs from frozen footprint join source")
    for path in (run_dir / "config.json", run_dir / "catalog" / "catalog.json"):
        original = next((s for s in join_audit["sources"] if s["path"] == str(path.resolve())), None)
        if original is None or source_info(path, True) != original:
            raise ValueError("Ledger config/catalog differs from footprint join")
    if contract["sample_id"] != join_audit["sample_id"]:
        raise ValueError("Sample disagrees with footprint join")
    coverage_path = run_dir / "barcodes" / "variant_barcode_coverage.tsv"
    variants = {item["variant_id"]: item["class"] for item in catalog["library"]["variants"]}
    expected_selected = 0
    with coverage_path.open(newline="") as handle:
        coverage_reader = csv.DictReader(handle, delimiter="\t")
        if not {"variant_id", "cohort", "molecules"}.issubset(coverage_reader.fieldnames or []):
            raise ValueError("Invalid variant barcode coverage table")
        seen_coverage = set()
        for entry in coverage_reader:
            if entry["cohort"] == "sequence_qualified" and variants.get(entry["variant_id"]) in ("WT", "SNV"):
                if entry["variant_id"] in seen_coverage:
                    raise ValueError("Repeated sequence-qualified coverage variant")
                seen_coverage.add(entry["variant_id"])
                expected_selected += int(entry["molecules"])
    sources = [source_info(config, True), source_info(joined), source_info(join_audit_path, True),
               source_info(contract_path, True), source_info(bam_path), source_info(str(bam_path) + ".bai"),
               source_info(run_dir / "config.json", True),
               source_info(run_dir / "catalog" / "catalog.json", True), source_info(coverage_path, True)]
    sources += [source_info(Path(__file__).with_name(name), True) for name in
                ("heterogeneity.py", "footprint_join.py", "processing.py", "io.py")]
    if limit_primary is not None and limit_primary < 1:
        raise ValueError("--limit-primary must be positive")
    counts = Counter()
    retention = defaultdict(Counter)
    started = time.monotonic()
    with atomic_output(output) as out, gzip.open(joined, "rt", newline="") as source, \
            gzip.open(out / "molecule_heterogeneity.tsv.gz", "wt", compresslevel=1, newline="") as dest, \
            pysam.AlignmentFile(str(bam_path), "rb", threads=1) as bam:
        reader = csv.DictReader(source, delimiter="\t")
        if not set(META_FIELDS + ["genotype_status"]).issubset(reader.fieldnames or []):
            raise ValueError("Joined table lacks required columns")
        writer = csv.DictWriter(dest, fieldnames=FEATURE_FIELDS, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        if bam.get_reference_length(catalog["library"]["contig"]) != cfg["reference_length"]:
            raise ValueError("BAM/plasmid length disagreement")
        for read in bam.fetch(catalog["library"]["contig"]):
            counts["bam_records"] += 1
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                counts["nonprimary_records"] += 1
                continue
            row = next(reader, None)
            if row is None:
                raise ValueError("Joined table exhausted before BAM primary reads")
            counts["primary_rows"] += 1
            if (row["read_id"] != read.query_name or int(row["read_length"]) != read.query_length
                    or row["strand"] != ("-" if read.is_reverse else "+")):
                raise ValueError("Ordered joined/BAM identity mismatch at primary " + str(counts["primary_rows"]))
            if row["sample_id"] != join_audit["sample_id"] or row["library_id"] != join_audit["library_id"]:
                raise ValueError("Joined sample/library disagreement")
            retained = selected(row)
            key = tuple(row.get(name, "") for name in RETENTION_FIELDS[:7])
            group = retention[key]
            group["total_primary"] += 1
            group["sequence_qualified"] += int(retained)
            for flag in ("analysis_eligible", "footprint_analysis_eligible",
                         "footprint_analysis_eligible_with_barcode"):
                group[flag] += int(row[flag] == "1")
            for flag in ("sequence_qc", "length_qc", "methylation_qc"):
                group["failed_" + flag] += int(row[flag] != "1")
            if retained:
                if variants.get(row["observed_variant_id"]) != row["genotype_class"]:
                    raise ValueError("Selected genotype disagrees with catalog: " + read.query_name)
                feature = extract_features(read, row, cfg)
                if row["region_complete"] == "1" and row["footprint_status"] == "evaluated":
                    for label in ("tf_small", "nucleosome", "msp"):
                        original = row[label + "_fraction"]
                        observed = feature["element_" + label + "_fraction"]
                        if original and (observed == "" or abs(float(original) - observed) > 1e-10):
                            raise ValueError("Element footprint fraction disagrees with join: " + read.query_name)
                writer.writerow(feature)
                counts["selected_rows"] += 1
                counts["selected_analysis_eligible"] += int(row["analysis_eligible"] == "1")
                counts["selected_footprint_analysis_eligible"] += int(row["footprint_analysis_eligible"] == "1")
                counts["selected_footprint_analysis_eligible_with_barcode"] += int(
                    row["footprint_analysis_eligible_with_barcode"] == "1")
                counts["selected_missing_MA"] += int(row["footprint_status"] == "missing_MA")
                counts["selected_missing_m6a"] += int(row["m6a_status"] != "ok")
            if progress_every and counts["primary_rows"] % progress_every == 0:
                print(f"Heterogeneity: {counts['primary_rows']:,} primary, {counts['selected_rows']:,} selected", flush=True)
            if limit_primary is not None and counts["primary_rows"] >= limit_primary:
                break
        if limit_primary is None:
            if next(reader, None) is not None:
                raise ValueError("Joined table has rows after BAM traversal")
            if counts["primary_rows"] != join_audit["counts"]["matched_primary_molecules"]:
                raise ValueError("Primary count differs from prior join")
            if counts["primary_rows"] != join_audit["counts"]["output_molecules"]:
                raise ValueError("Joined table unexpectedly contains missing BAM rows")
            if counts["selected_rows"] != expected_selected:
                raise ValueError("Selected sequence-qualified WT/SNV count differs from independent coverage")
            if counts["selected_analysis_eligible"] != join_audit["counts"]["ledger_analysis_eligible"]:
                raise ValueError("Existing analysis cohort not conserved")
            for flag, expected in (("selected_footprint_analysis_eligible", "footprint_analysis_eligible"),
                                   ("selected_footprint_analysis_eligible_with_barcode",
                                    "footprint_analysis_eligible_with_barcode")):
                if counts[flag] != join_audit["counts"][expected]:
                    raise ValueError("Existing footprint cohort not conserved: " + flag)
        rows = [dict(zip(RETENTION_FIELDS[:7], key),
                     **{name: counter[name] for name in RETENTION_FIELDS[7:]})
                for key, counter in sorted(retention.items())]
        write_table(out / "retention.tsv", rows, RETENTION_FIELDS)
        check_sources(sources)
        report = {
            "schema_version": "fsmpra.heterogeneity.v1", "complete": limit_primary is None,
            "limit_primary": limit_primary, "counts": dict(counts), "sources": sources,
            "expected_selected_from_coverage": expected_selected,
            "regions": cfg["regions"], "tile_size": cfg["tile_size"],
            "tile_min_at": cfg["tile_min_at"], "thresholds": [128, 204],
            "footprint_size_thresholds": {"tf_min": cfg["tf_min"], "tf_max": cfg["tf_max"],
                                          "nuc_min": cfg["nuc_min"]},
            "seconds": time.monotonic() - started,
            "denominators": "Regional methylation uses aligned A/T; footprint fractions use aligned bases.",
            "missingness": "Missing MM/ML or MA is blank. Unknown omitted MM opportunities have observed counts but blank fractions. Present sparse MA tracks can be zero.",
            "segment_interpretation": "Nucleosome counts are projected MA segments, not biological nucleosome counts; circular split segments are not reunited and spacing excludes the circular wrap.",
            "tile_interpretation": "100-bp fixed-reference bins wholly inside distal intervals; tile fractions summarize observed calls, not labeling efficiency.",
        }
        write_json(out / "audit.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--limit-primary", type=int)
    parser.add_argument("--progress-every", type=int, default=100000)
    args = parser.parse_args()
    report = extract(args.config, args.output, args.limit_primary, args.progress_every)
    print(json.dumps(report["counts"], sort_keys=True))


if __name__ == "__main__":
    main()

"""Descriptive WT sample audit of original FiberHMM calls and callable sequence.

The input is a deterministic barcode sample, never a full WT population. A/B
comparisons describe recurrence on the same WT sequence, without p-values.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

from .io import atomic_output, open_text, read_json, source_info, write_json, write_table
from .processing import check_sources
from .statistics import fixed_endpoints
from .footprint_quality import band_for

SCHEMA = "fsmpra.footprint_quality_summary.v1"
REGIONS = ("whole", "element")
TRACKS = ("tf", "nuc", "msp")
AQ_TRACKS = ("tf", "nuc")
REP_GROUPS = ("no_nuc_with_localTF", "no_nuc_without_localTF",
              "exactly_one_nuc_longer1000", "other")
BANDS = {"tf": ("lt10", "10_60", "gt60"),
         "nuc": ("lt85", "85_250", "251_500", "501_1000", "gt1000")}
FEATURES = ("whole_tf_small_fraction", "whole_nucleosome_fraction", "whole_msp_fraction",
            "element_tf_small_fraction", "element_nucleosome_fraction", "element_msp_fraction",
            "whole_m6a_fraction", "element_m6a_fraction", "whole_nucleosome_segments",
            "whole_tf_calls_small") + tuple(
                f"{region}_{track}_{band}_{suffix}"
                for region in REGIONS for track in ("tf", "nuc")
                for band in BANDS[track] for suffix in ("count", "fraction"))


def _num(value):
    if value in ("", None):
        return math.nan
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("Nonfinite sampled feature")
    return number


def _offsets(value, length, name, required=True):
    if value in ("", None):
        if required:
            raise ValueError("Missing offset array: " + name)
        return None
    raw = json.loads(value)
    if not isinstance(raw, list) or any(not isinstance(x, int) or isinstance(x, bool) or x < 0 or x >= length for x in raw):
        raise ValueError("Invalid offset array: " + name)
    if len(raw) != len(set(raw)):
        raise ValueError("Duplicate offset: " + name)
    return np.array(raw, dtype=np.int16)


def _groups(group):
    if group not in REP_GROUPS:
        raise ValueError("Unknown representative group: " + group)
    broad = "no_nuc" if group.startswith("no_nuc_") else "other_nuc"
    return ("all", broad, group)


def _read_barcodes(path, audit):
    barcodes = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"barcode", "variant_id", "eligible_depth", "split", "weight", "selection_hash"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Sampled barcode table lacks fields")
        for row in reader:
            barcode = row["barcode"]
            if not barcode or barcode in barcodes or row["split"] not in ("A", "B"):
                raise ValueError("Duplicate/invalid sampled barcode")
            depth = int(row["eligible_depth"])
            if depth < 1 or not math.isclose(float(row["weight"]), 1 / depth, abs_tol=1e-12):
                raise ValueError("Sample barcode weight/depth disagreement")
            row["eligible_depth"] = depth
            barcodes[barcode] = row
    if len(barcodes) != audit["counts"]["sampled_barcodes"]:
        raise ValueError("Sampled barcode count disagrees with audit")
    if sum(item["eligible_depth"] for item in barcodes.values()) != audit["counts"]["sampled_expected_molecules"]:
        raise ValueError("Sampled barcode depths disagree with audit")
    if len({item["variant_id"] for item in barcodes.values()}) != 1:
        raise ValueError("Expected a single sampled WT variant")
    return barcodes


def _rows(path, required):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(required) <= set(reader.fieldnames or []):
            raise ValueError("Missing input columns: " + str(sorted(set(required) - set(reader.fieldnames or []))))
        for row in reader:
            yield row


class Calls:
    """Stream calls in sampled-molecule order; do not retain the full call table."""

    def __init__(self, path):
        self.handle = open_text(path)
        self.reader = csv.DictReader(self.handle, delimiter="\t")
        required = {"read_id", "barcode", "split", "weight", "track", "band",
                    "original_length", "ref_intervals0", "quality_spec", "quality",
                    "query_at_count", "query_m6a128_count", "query_m6a204_count"}
        if not required <= set(self.reader.fieldnames or []):
            raise ValueError("Sampled call table lacks fields")
        self.next = next(self.reader, None)
        self.total = 0
        self.hist = defaultdict(lambda: [0, 0.0])
        self.aq_hist = defaultdict(lambda: [0, 0.0, 0, 0.0, 0, 0.0])
        self.observability = defaultdict(lambda: [0, 0.0])
        self.observability_bands = defaultdict(lambda: [0, 0, 0, 0])
        self.past = set()

    def for_molecule(self, row):
        read_id = row["read_id"]
        self.past.add(read_id)
        found = []
        while self.next is not None and self.next["read_id"] == read_id:
            call = self.next
            if (call["barcode"] != row["corrected_barcode"] or call["split"] != row["sample_split"]
                    or not math.isclose(float(call["weight"]), float(row["barcode_weight"]), abs_tol=1e-12)):
                raise ValueError("Sampled call/molecule identity disagreement")
            if call["track"] not in TRACKS or not call["band"] or int(call["original_length"]) < 1:
                raise ValueError("Invalid sampled call track/length")
            if call["band"] != band_for(call["track"], int(call["original_length"])):
                raise ValueError("Sampled call length/band disagreement")
            quality = json.loads(call["quality"])
            if (not isinstance(quality, (list, type(None))) or
                    (quality is not None and (len(quality) != len(call["quality_spec"]) or
                     any(type(score) is not int or not 0 <= score <= 255 for score in quality)))):
                raise ValueError("Invalid sampled call AQ values")
            call["aq100"] = (quality is not None and len(quality) > 0 and
                             quality[0] >= 100)
            call["aq_known"] = quality is not None and len(quality) > 0
            call["aq_applicable"] = call["track"] in AQ_TRACKS and bool(call["quality_spec"])
            at_count = int(call["query_at_count"])
            m6a128_count = int(call["query_m6a128_count"])
            m6a204_count = int(call["query_m6a204_count"])
            if not 0 <= m6a204_count <= m6a128_count <= at_count <= int(call["original_length"]):
                raise ValueError("Invalid call query A/T or m6A observability")
            spans = json.loads(call["ref_intervals0"])
            if not isinstance(spans, list) or any(not isinstance(span, list) or len(span) != 2 or
                    not all(isinstance(x, int) for x in span) or span[0] < 0 or span[1] <= span[0]
                    for span in spans):
                raise ValueError("Invalid call reference intervals")
            call["spans"] = spans
            key = (call["track"], call["band"], int(call["original_length"]))
            self.hist[key][0] += 1
            self.hist[key][1] += float(call["weight"])
            for metric, value in (("query_at_count", at_count),
                                  ("query_m6a128_count", m6a128_count),
                                  ("query_m6a204_count", m6a204_count)):
                self.observability[(call["track"], call["band"], metric, value)][0] += 1
                self.observability[(call["track"], call["band"], metric, value)][1] += float(call["weight"])
            band_obs = self.observability_bands[(call["track"], call["band"])]
            band_obs[0] += at_count
            band_obs[1] += m6a128_count
            band_obs[2] += m6a204_count
            band_obs[3] += int(at_count == 0)
            if call["aq_applicable"] and call["aq100"]:
                self.aq_hist[(call["track"], call["band"])][0] += 1
                self.aq_hist[(call["track"], call["band"])][1] += float(call["weight"])
            if call["aq_applicable"] and not call["aq_known"]:
                self.aq_hist[(call["track"], call["band"])][2] += 1
                self.aq_hist[(call["track"], call["band"])][3] += float(call["weight"])
            if not call["aq_applicable"]:
                self.aq_hist[(call["track"], call["band"])][4] += 1
                self.aq_hist[(call["track"], call["band"])][5] += float(call["weight"])
            self.total += 1
            found.append(call)
            self.next = next(self.reader, None)
        if self.next is not None and self.next["read_id"] in self.past:
            raise ValueError("Sampled call rows are out of molecule order")
        return found

    def close(self):
        try:
            if self.next is not None:
                raise ValueError("Trailing sampled call rows")
        finally:
            self.handle.close()


def _call_offsets(calls, start, length, track, exclude_band=None, aq100=False):
    covered = set()
    for call in calls:
        if call["track"] != track or call["band"] == exclude_band or (aq100 and not call["aq100"]):
            continue
        for lo, hi in call["spans"]:
            first, last = max(lo, start), min(hi, start + length)
            if first < last:
                covered.update(range(first - start, last - start))
    return np.array(sorted(covered), dtype=np.int16)


def _profile_bucket():
    return {name: np.zeros(211, dtype=np.float64) for name in
            ("aligned", "at", "at_observed", "m6a", "unmodified_opportunity", "tf_callable",
             "tf_covered", "nuc_callable", "nuc_covered", "msp_callable", "msp_covered",
             "tf_aq100_callable", "tf_aq100_covered", "nuc_aq100_callable",
             "nuc_aq100_covered")}


def _profile_add(bucket, aligned, at, m6a, tf, nuc, msp, aq100, aq100_valid,
                 weight, methylation_valid, tf_valid, nuc_valid, msp_valid):
    bucket["aligned"][aligned] += weight
    bucket["at"][at] += weight
    if methylation_valid:
        bucket["at_observed"][at] += weight
        bucket["m6a"][m6a] += weight
        bucket["unmodified_opportunity"][np.setdiff1d(at, m6a, assume_unique=True)] += weight
    for valid, covered, name in ((tf_valid, tf, "tf"), (nuc_valid, nuc, "nuc"), (msp_valid, msp, "msp")):
        if valid:
            bucket[name + "_callable"][aligned] += weight
            bucket[name + "_covered"][covered] += weight
        if name in AQ_TRACKS and valid and aq100_valid[name]:
            bucket[name + "_aq100_callable"][aligned] += weight
            bucket[name + "_aq100_covered"][aq100[name]] += weight


def _ratio(top, bottom):
    return float(top / bottom) if bottom else None


def _correlation(a, b):
    valid = np.isfinite(a) & np.isfinite(b)
    x, y = a[valid], b[valid]
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return {"positions": int(len(x)), "pearson": None, "spearman": None}
    return {"positions": int(len(x)), "pearson": float(np.corrcoef(x, y)[0, 1]),
            "spearman": float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])}


def _nearest_at(sequence):
    at = np.flatnonzero(np.fromiter((base in "AT" for base in sequence), dtype=np.bool_, count=len(sequence)))
    if not len(at):
        raise ValueError("WT element has no A/T bases")
    return np.min(np.abs(np.arange(len(sequence))[:, None] - at[None, :]), axis=1)


def summarize(config, sample, output, allow_pilot=False):
    config, sample, output = Path(config), Path(sample), Path(output)
    cfg = read_json(config)
    if cfg.get("schema_version") != "fsmpra.footprint_quality_config.v1":
        raise ValueError("Unexpected footprint-quality configuration")
    audit_path = sample / "audit.json"
    audit = read_json(audit_path)
    if audit.get("schema_version") != "fsmpra.footprint_quality.v1":
        raise ValueError("Unexpected sample extraction audit")
    if not audit.get("complete") and not allow_pilot:
        raise ValueError("Pilot sample requires --allow-pilot")
    if audit["complete"] and audit["counts"]["sampled_barcodes"] != cfg["sample_barcodes"]:
        raise ValueError("Complete sample barcode count disagrees with config")
    current_config = source_info(config, True)
    frozen_config = next((item for item in audit.get("sources", []) if item["path"] == current_config["path"]), None)
    if frozen_config != current_config:
        raise ValueError("Sample extraction used different config")
    check_sources(audit["sources"])
    catalog_path = Path(cfg["catalog"])
    catalog = read_json(catalog_path)["library"]
    frozen_catalog = next((item for item in audit["sources"] if item["path"] == str(catalog_path.resolve())), None)
    if frozen_catalog != source_info(catalog_path, True):
        raise ValueError("Sample extraction used different catalog")
    wt = [variant for variant in catalog["variants"] if variant["class"] == "WT"]
    if len(wt) != 1 or len(wt[0]["sequence"]) != 211:
        raise ValueError("Expected one 211-base WT element")
    sequence = wt[0]["sequence"]
    start, end = audit["element"]
    if end - start != 211:
        raise ValueError("Sample element length differs from 211")
    endpoints = fixed_endpoints(start, end)
    if len(endpoints) != 22:
        raise ValueError("Expected 21 TF windows plus MSP endpoint")
    barcode_path = sample / "sampled_barcodes.tsv"
    molecule_path = sample / "sampled_molecules.tsv.gz"
    call_path = sample / "sampled_calls.tsv.gz"
    length_path = sample / "call_length_histograms.tsv"
    rep_path = sample / "representatives.jsonl.gz"
    source_records = [source_info(path, checksum) for path, checksum in (
        (config, True), (audit_path, True), (barcode_path, True),
        (molecule_path, False), (call_path, False), (length_path, True),
        (rep_path, False), (catalog_path, True), (Path(__file__), True),
        (Path(__file__).with_name("footprint_quality.py"), True),
        (Path(__file__).with_name("statistics.py"), True),
        (Path(__file__).with_name("processing.py"), True),
        (Path(__file__).with_name("io.py"), True))]
    barcodes = _read_barcodes(barcode_path, audit)
    barcode_ids = sorted(barcodes)
    barcode_index = {barcode: i for i, barcode in enumerate(barcode_ids)}
    n = len(barcode_ids)
    e = len(endpoints)
    presence_sum = np.zeros((n, e), dtype=np.float64)
    coverage_sum = np.zeros((n, e), dtype=np.float64)
    valid = np.zeros((n, e), dtype=np.int32)
    seen_depth = Counter()
    split_molecules = Counter()
    group_molecules = Counter()
    aq100_missing_molecules = Counter()
    group_barcodes = defaultdict(set)
    sample_id = library_id = None
    profiles = defaultdict(_profile_bucket)
    moments = defaultdict(lambda: [0, 0.0, 0.0])
    calls = Calls(call_path)
    required = {"read_id", "sample_id", "library_id", "corrected_barcode", "assigned_variant_id",
                "sample_split", "barcode_weight", "representative_group", "element_aligned_offsets0",
                "element_at_offsets0", "element_m6a_offsets0", "element_tf_small_offsets0",
                "element_msp_offsets0", "element_m6a_fraction", "element_tf_small_fraction",
                "element_nucleosome_fraction", "element_msp_fraction"} | set(FEATURES)
    rows = 0
    try:
        for row in _rows(molecule_path, required):
            rows += 1
            if rows % 25000 == 0:
                print(f"Footprint-quality summary: {rows:,} sampled molecules", flush=True)
            barcode = row["corrected_barcode"]
            if barcode not in barcodes or row["assigned_variant_id"] != wt[0]["variant_id"]:
                raise ValueError("Sampled molecule has non-WT/unknown barcode")
            info = barcodes[barcode]
            if row["sample_split"] != info["split"] or not math.isclose(float(row["barcode_weight"]), 1 / info["eligible_depth"], abs_tol=1e-12):
                raise ValueError("Molecule split/weight disagrees with barcode")
            if sample_id is None:
                sample_id, library_id = row["sample_id"], row["library_id"]
            elif (row["sample_id"], row["library_id"]) != (sample_id, library_id):
                raise ValueError("Mixed sample/library IDs")
            if not sample_id or not library_id:
                raise ValueError("Missing sample/library ID")
            group = row["representative_group"]
            groups = _groups(group)
            seen_depth[barcode] += 1
            split_molecules[info["split"]] += 1
            group_molecules[group] += 1
            for label in groups:
                group_barcodes[label].add(barcode)
            weight = 1 / info["eligible_depth"]
            aligned = _offsets(row["element_aligned_offsets0"], 211, "aligned")
            at = _offsets(row["element_at_offsets0"], 211, "A/T")
            m6a = _offsets(row["element_m6a_offsets0"], 211, "m6A", required=False)
            tf = _offsets(row["element_tf_small_offsets0"], 211, "TF", required=False)
            msp = _offsets(row["element_msp_offsets0"], 211, "MSP", required=False)
            if not set(at).issubset(aligned) or (m6a is not None and not set(m6a).issubset(at)):
                raise ValueError("Methylation offsets exceed callable opportunities")
            observed_calls = calls.for_molecule(row)
            projected_tf = _call_offsets([call for call in observed_calls if call["band"] == "10_60"], start, 211, "tf")
            projected_msp = _call_offsets(observed_calls, start, 211, "msp")
            if tf is not None and not np.array_equal(tf, projected_tf):
                raise ValueError("TF bitmap disagrees with original calls")
            if msp is not None and not np.array_equal(msp, projected_msp):
                raise ValueError("MSP bitmap disagrees with original calls")
            nuc = _call_offsets(observed_calls, start, 211, "nuc", exclude_band="lt85")
            if (tf is not None and not set(tf).issubset(aligned)) or not set(nuc).issubset(aligned) or (msp is not None and not set(msp).issubset(aligned)):
                raise ValueError("Projected call offsets exceed aligned positions")
            methylation_valid = row["element_m6a_fraction"] != "" and m6a is not None
            tf_valid = row["element_tf_small_fraction"] != "" and tf is not None
            nuc_valid = row["element_nucleosome_fraction"] != ""
            msp_valid = row["element_msp_fraction"] != "" and msp is not None
            aq100 = {}
            aq100_valid = {}
            for track, excluded in (("tf", {"lt10", "gt60"}), ("nuc", {"lt85"})):
                relevant = [call for call in observed_calls if call["track"] == track and
                            call["band"] not in excluded]
                aq100_valid[track] = all(call["aq_applicable"] and call["aq_known"] for call in relevant)
                if not aq100_valid[track]:
                    aq100_missing_molecules[track] += 1
                aq100[track] = _call_offsets(relevant, start, 211, track, aq100=True)
            if row["element_m6a_fraction"] != "" and m6a is None:
                raise ValueError("Missing m6A positions despite valid fraction")
            feature_values = {feature: _num(row[feature]) for feature in FEATURES}
            for feature, value in feature_values.items():
                if feature.endswith("_fraction") and math.isfinite(value) and not 0 <= value <= 1:
                    raise ValueError("Sampled fraction outside [0,1]: " + feature)
            for split in ("all", info["split"]):
                for label in groups:
                    for weighting, factor in (("molecule", 1.0), ("equal_barcode", weight)):
                        _profile_add(profiles[(split, label, weighting)], aligned, at,
                                     m6a if m6a is not None else np.empty(0, np.int16),
                                     tf if tf is not None else np.empty(0, np.int16), nuc,
                                     msp if msp is not None else np.empty(0, np.int16),
                                     aq100, aq100_valid, factor,
                                     methylation_valid, tf_valid, nuc_valid, msp_valid)
                        for feature in FEATURES:
                            value = feature_values[feature]
                            bucket = moments[(split, label, weighting, feature)]
                            bucket[0] += factor
                            if math.isfinite(value):
                                bucket[1] += factor
                                bucket[2] += factor * value
            aligned_set = set(aligned)
            tf_set = set(tf) if tf is not None else set()
            msp_set = set(msp) if msp is not None else set()
            i = barcode_index[barcode]
            for j, endpoint in enumerate(endpoints):
                lo, hi = endpoint["relative_start"], endpoint["relative_end"]
                callable_bases = sum(lo <= x < hi for x in aligned_set)
                if endpoint["track"] == "tf_small":
                    if not tf_valid or not callable_bases:
                        continue
                    covered = sum(lo <= x < hi for x in tf_set)
                else:
                    if not msp_valid or not callable_bases:
                        continue
                    covered = sum(lo <= x < hi for x in msp_set)
                presence_sum[i, j] += int(covered > 0)
                coverage_sum[i, j] += covered / callable_bases
                valid[i, j] += 1
    except BaseException:
        calls.handle.close()
        raise
    else:
        calls.close()
    if rows != audit["counts"]["matched_molecules"] or calls.total != audit["counts"]["sampled_calls"]:
        raise ValueError("Sample molecule/call count disagrees with extraction audit")
    if rows != sum(seen_depth.values()) or rows != audit["counts"]["sampled_expected_molecules"]:
        raise ValueError("Sample molecule counts do not conserve")
    for barcode, info in barcodes.items():
        if seen_depth[barcode] != info["eligible_depth"]:
            raise ValueError("Barcode molecule depth disagrees with sample table")
    if sum(split_molecules.values()) != rows or sum(group_molecules.values()) != rows:
        raise ValueError("Split/group molecule counts do not conserve")
    if len(barcodes) != sum(item["split"] == "A" for item in barcodes.values()) + sum(item["split"] == "B" for item in barcodes.values()):
        raise ValueError("Split barcode counts do not conserve")
    if not math.isclose(moments[("all", "all", "equal_barcode", "whole_nucleosome_segments")][0], n,
                        rel_tol=0, abs_tol=1e-7):
        raise ValueError("Equal-barcode weights do not conserve")
    length_hist = {}
    for row in _rows(length_path, ("track", "band", "original_length", "call_count", "weighted_call_count")):
        key = (row["track"], row["band"], int(row["original_length"]))
        if key in length_hist:
            raise ValueError("Duplicate call-length histogram row")
        length_hist[key] = (int(row["call_count"]), float(row["weighted_call_count"]))
    if set(length_hist) != set(calls.hist):
        raise ValueError("Call-length histogram keys differ from sampled calls")
    for key, (count, weighted) in length_hist.items():
        observed = calls.hist[key]
        if count != observed[0] or not math.isclose(weighted, observed[1], rel_tol=1e-9, abs_tol=1e-7):
            raise ValueError("Call-length histogram disagrees with sampled calls")
    depths = np.array([barcodes[b]["eligible_depth"] for b in barcode_ids], dtype=np.int32)
    splits = np.array([barcodes[b]["split"] for b in barcode_ids], dtype="U1")
    presence = np.full((n, e), np.nan, dtype=np.float64)
    coverage = np.full((n, e), np.nan, dtype=np.float64)
    np.divide(presence_sum, valid, out=presence, where=valid > 0)
    np.divide(coverage_sum, valid, out=coverage, where=valid > 0)
    if np.any(valid > depths[:, None]) or np.any((presence < 0) | (presence > 1)) or np.any((coverage < 0) | (coverage > 1)):
        raise ValueError("Endpoint matrix invalid")
    return _write_outputs(config, sample, output, cfg, audit, barcodes, barcode_ids, depths, splits,
                          endpoints, presence, coverage, valid, profiles, moments, sequence,
                          split_molecules, group_molecules, group_barcodes, length_hist,
                          calls.aq_hist, calls.observability, calls.observability_bands,
                          aq100_missing_molecules, rep_path,
                          sample_id, library_id, rows, source_records)


def _write_outputs(config, sample, output, cfg, audit, barcodes, barcode_ids, depths, splits,
                   endpoints, presence, coverage, valid, profiles, moments, sequence,
                   split_molecules, group_molecules, group_barcodes, length_hist,
                   aq_hist, observability, observability_bands,
                   aq100_missing_molecules, rep_path,
                   sample_id, library_id, rows, source_records):
    start, end = audit["element"]
    distances = _nearest_at(sequence)
    with atomic_output(output) as out:
        npz_path = out / "barcode_endpoint_matrix.npz"
        np.savez_compressed(npz_path,
            barcodes=np.array(barcode_ids, dtype="U"),
            endpoint_ids=np.array([x["endpoint_id"] for x in endpoints], dtype="U"),
            presence=presence, coverage=coverage, valid_molecules=valid,
            depths=depths, splits=splits,
            relative_starts=np.array([x["relative_start"] for x in endpoints], dtype=np.int16),
            relative_ends=np.array([x["relative_end"] for x in endpoints], dtype=np.int16),
            reference_starts0=np.array([x["reference_start0"] for x in endpoints], dtype=np.int32),
            reference_ends0=np.array([x["reference_end0"] for x in endpoints], dtype=np.int32))
        endpoint_rows = []
        for i, barcode in enumerate(barcode_ids):
            for j, endpoint in enumerate(endpoints):
                endpoint_rows.append(dict(barcode=barcode, split=splits[i], eligible_depth=int(depths[i]),
                    endpoint_id=endpoint["endpoint_id"], track=endpoint["track"],
                    relative_start=endpoint["relative_start"], relative_end=endpoint["relative_end"],
                    reference_start0=endpoint["reference_start0"], reference_end0=endpoint["reference_end0"],
                    valid_molecules=int(valid[i, j]), presence_mean=presence[i, j] if valid[i, j] else None,
                    coverage_mean=coverage[i, j] if valid[i, j] else None,
                    weighted_presence_sum=presence[i, j] * valid[i, j] / depths[i] if valid[i, j] else 0,
                    weighted_coverage_sum=coverage[i, j] * valid[i, j] / depths[i] if valid[i, j] else 0))
        write_table(out / "barcode_endpoints.tsv.gz", endpoint_rows, list(endpoint_rows[0]))
        for split in ("all", "A", "B"):
            for group in ("all", "no_nuc", "other_nuc") + REP_GROUPS:
                for weighting in ("molecule", "equal_barcode"):
                    profiles[(split, group, weighting)]
        profile_rows = []
        for (split, group, weighting), bucket in sorted(profiles.items()):
            for pos in range(211):
                at_den = bucket["at_observed"][pos]
                profile_rows.append(dict(split=split, group=group, weighting=weighting,
                    relative_pos0=pos, reference_pos0=start + pos, reference_base=sequence[pos],
                    distance_to_nearest_AT=int(distances[pos]),
                    aligned_weight=bucket["aligned"][pos], at_weight=bucket["at"][pos],
                    at_observed_weight=at_den, m6a_weight=bucket["m6a"][pos],
                    unmethylated_opportunity_weight=bucket["unmodified_opportunity"][pos],
                    observed_m6a_fraction=_ratio(bucket["m6a"][pos], at_den),
                    unmethylated_opportunity_fraction=_ratio(bucket["unmodified_opportunity"][pos], at_den),
                    tf_callable_weight=bucket["tf_callable"][pos], tf_covered_weight=bucket["tf_covered"][pos],
                    tf_coverage_fraction=_ratio(bucket["tf_covered"][pos], bucket["tf_callable"][pos]),
                    nuc_callable_weight=bucket["nuc_callable"][pos], nuc_covered_weight=bucket["nuc_covered"][pos],
                    nuc_coverage_fraction=_ratio(bucket["nuc_covered"][pos], bucket["nuc_callable"][pos]),
                    msp_callable_weight=bucket["msp_callable"][pos], msp_covered_weight=bucket["msp_covered"][pos],
                    msp_coverage_fraction=_ratio(bucket["msp_covered"][pos], bucket["msp_callable"][pos]),
                    **{f"{track}_aq100_{suffix}": value for track in AQ_TRACKS for suffix, value in (
                        ("callable_weight", bucket[f"{track}_aq100_callable"][pos]),
                        ("covered_weight", bucket[f"{track}_aq100_covered"][pos]),
                        ("coverage_fraction", _ratio(bucket[f"{track}_aq100_covered"][pos],
                                                     bucket[f"{track}_aq100_callable"][pos]))) }))
        write_table(out / "element_profiles.tsv.gz", profile_rows, list(profile_rows[0]))
        context_rows = []
        for row in profile_rows:
            pos = row["relative_pos0"]
            context_rows.append(dict(split=row["split"], group=row["group"],
                weighting=row["weighting"], relative_pos0=pos, reference_pos0=start + pos,
                trinucleotide=(sequence[pos - 1] if pos else "N") + sequence[pos] +
                              (sequence[pos + 1] if pos + 1 < len(sequence) else "N"),
                reference_base=sequence[pos], distance_to_nearest_AT=row["distance_to_nearest_AT"],
                observed_at_weight=row["at_observed_weight"], m6a_weight=row["m6a_weight"],
                unmethylated_opportunity_weight=row["unmethylated_opportunity_weight"]))
        write_table(out / "sequence_context_opportunities.tsv.gz", context_rows, list(context_rows[0]))
        moment_rows = []
        for (split, group, weighting, feature), (total_weight, valid_weight, sum_value) in sorted(moments.items()):
            moment_rows.append(dict(split=split, group=group, weighting=weighting, feature=feature,
                                    total_weight=total_weight, valid_weight=valid_weight,
                                    missing_weight=total_weight - valid_weight,
                                    mean=_ratio(sum_value, valid_weight)))
        write_table(out / "group_feature_summaries.tsv.gz", moment_rows, list(moment_rows[0]))
        band_rows = []
        band_counts = defaultdict(lambda: [0, 0.0])
        for (track, band, length), (count, weighted) in length_hist.items():
            band_counts[(track, band)][0] += count
            band_counts[(track, band)][1] += weighted
        for (track, band), (count, weighted) in sorted(band_counts.items()):
            (aq100_calls, aq100_weighted, aq_missing, aq_missing_weighted,
             aq_unscored, aq_unscored_weighted) = aq_hist[(track, band)]
            at_sites, m6a128_sites, m6a204_sites, zero_at_calls = observability_bands[(track, band)]
            band_rows.append(dict(track=track, band=band, calls=count,
                                  equal_barcode_weighted_calls=weighted,
                                  aq100_calls=aq100_calls,
                                  equal_barcode_weighted_aq100_calls=aq100_weighted,
                                  missing_aq_calls=aq_missing,
                                  equal_barcode_weighted_missing_aq_calls=aq_missing_weighted,
                                  aq_unscored_calls=aq_unscored,
                                  equal_barcode_weighted_aq_unscored_calls=aq_unscored_weighted,
                                  query_at_opportunities=at_sites,
                                  query_m6a128_calls=m6a128_sites,
                                  query_m6a204_calls=m6a204_sites,
                                  zero_query_at_calls=zero_at_calls))
        write_table(out / "call_length_bands.tsv", band_rows,
                    ["track", "band", "calls", "equal_barcode_weighted_calls",
                     "aq100_calls", "equal_barcode_weighted_aq100_calls", "missing_aq_calls",
                     "equal_barcode_weighted_missing_aq_calls", "aq_unscored_calls",
                     "equal_barcode_weighted_aq_unscored_calls", "query_at_opportunities",
                     "query_m6a128_calls", "query_m6a204_calls", "zero_query_at_calls"])
        observability_rows = [dict(track=track, band=band, metric=metric, value=value,
                                   calls=counts[0], equal_barcode_weighted_calls=counts[1])
                              for (track, band, metric, value), counts in sorted(observability.items())]
        write_table(out / "call_observability_histograms.tsv.gz", observability_rows,
                    ["track", "band", "metric", "value", "calls", "equal_barcode_weighted_calls"])
        _plots(out, profile_rows, band_counts, rep_path, start, end, audit["complete"])
        check_sources(source_records)
        check_sources(audit["sources"])
        matrix_source = source_info(npz_path, True)
        matrix_source["path"] = str((output / npz_path.name).resolve())
        a = np.array([profiles[("A", "all", "equal_barcode")]["tf_covered"][i] /
                      profiles[("A", "all", "equal_barcode")]["tf_callable"][i]
                      if profiles[("A", "all", "equal_barcode")]["tf_callable"][i] else np.nan
                      for i in range(211)])
        b = np.array([profiles[("B", "all", "equal_barcode")]["tf_covered"][i] /
                      profiles[("B", "all", "equal_barcode")]["tf_callable"][i]
                      if profiles[("B", "all", "equal_barcode")]["tf_callable"][i] else np.nan
                      for i in range(211)])
        report = {"schema_version": SCHEMA, "complete": bool(audit["complete"]),
            "pilot": not bool(audit["complete"]), "sample_id": sample_id,
            "library_id": library_id, "sampled_barcodes": len(barcode_ids),
            "sampled_molecules": rows, "split_barcodes": dict(Counter(splits)),
            "split_molecules": dict(split_molecules), "representative_group_molecules": dict(group_molecules),
            "group_unique_barcodes": {key: len(value) for key, value in group_barcodes.items()},
            "aq100_missing_quality_molecules": dict(aq100_missing_molecules),
            "element": [start, end], "max_distance_to_nearest_AT": int(distances.max()),
            "tf_A_B_211base_concordance": _correlation(a, b),
            "endpoint_matrix": "barcode_endpoint_matrix.npz", "endpoint_count": len(endpoints),
            "source_audit": str((sample / "audit.json").resolve()),
            "sources": source_records + [matrix_source],
            "interpretation": "Descriptive deterministic WT barcode sample. Same WT sequence underlies A/B profiles; concordance is not independent motif or TF validation. No p-values or compartment assignment.",
            "missingness": "Blank measurements remain missing; evaluated empty calls are zero. Unmethylated opportunities are observed A/T sites without an ML-thresholded m6A call, not probabilities or true absence of methylation. TF/nucleosome AQ>=100 sensitivity uses the first raw AQ byte of each scored call; a molecule with any relevant call missing AQ is excluded from that track's AQ>=100 per-base denominator. MSP has no defined first-byte score."}
        write_json(out / "summary.json", report)
    return report


def _plots(out, profile_rows, band_counts, rep_path, start, end, complete):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    prefix = "WT SAMPLE" if complete else "WT PILOT SAMPLE"
    profiles = {(row["split"], row["group"], row["weighting"], row["relative_pos0"]): row
                for row in profile_rows}
    x = np.arange(211)
    fig, axes = plt.subplots(4, 1, figsize=(12, 11), sharex=True, constrained_layout=True)
    for split in ("A", "B"):
        values = [profiles[(split, "all", "equal_barcode", i)]["tf_coverage_fraction"] for i in x]
        axes[0].plot(x, [np.nan if v is None else v for v in values], label=f"Split {split}")
    axes[0].set(ylabel="TF-call covered / aligned weight", title="A/B local TF-call recurrence")
    axes[0].legend()
    for group in ("no_nuc", "other_nuc"):
        values = [profiles[("all", group, "equal_barcode", i)]["tf_coverage_fraction"] for i in x]
        axes[1].plot(x, [np.nan if v is None else v for v in values], label=group.replace("_", " "))
    axes[1].set(ylabel="TF-call covered / aligned weight", title="Separate observed nucleosome-call groups")
    axes[1].legend()
    for field, label in (("tf_coverage_fraction", "All original 10–60 bp calls"),
                         ("tf_aq100_coverage_fraction", "First AQ byte ≥100")):
        values = [profiles[("all", "all", "equal_barcode", i)][field] for i in x]
        axes[2].plot(x, [np.nan if v is None else v for v in values], label=label)
    axes[2].set(ylabel="TF covered / callable weight", title="Raw AQ sensitivity; missing AQ excluded")
    axes[2].legend()
    values = [profiles[("all", "all", "equal_barcode", i)]["observed_m6a_fraction"] for i in x]
    axes[3].plot(x, [np.nan if v is None else v for v in values], color="tab:blue", label="Observed m6A / A/T")
    values = [profiles[("all", "all", "equal_barcode", i)]["unmethylated_opportunity_fraction"] for i in x]
    axes[3].plot(x, [np.nan if v is None else v for v in values], color="tab:orange", label="No called m6A / observed A/T")
    axes[3].set(xlabel="Position within 211-bp WT element (0-based)", ylabel="Observed A/T opportunity fraction",
                title="A/T context and ML-thresholded marks")
    axes[3].legend()
    fig.suptitle(prefix + " · original call annotations and observable opportunities")
    fig.savefig(out / "element_quality_profiles.png", dpi=170)
    fig.savefig(out / "element_quality_profiles.svg")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    nuc = [(band, value[0]) for (track, band), value in band_counts.items() if track == "nuc"]
    order = ["lt85", "85_250", "251_500", "501_1000", "gt1000"]
    ax.bar(order, [dict(nuc).get(band, 0) for band in order])
    ax.set(xlabel="Original nucleosome-call query-length band (bp)", ylabel="Calls in sampled WT molecules",
           title=prefix + " · original MA nucleosome-call lengths")
    fig.savefig(out / "nucleosome_length_bands.png", dpi=170)
    fig.savefig(out / "nucleosome_length_bands.svg")
    plt.close(fig)

    rep_dir = out / "representatives"
    rep_dir.mkdir()
    with open_text(rep_path) as handle:
        for number, line in enumerate(handle):
            rep = json.loads(line)
            fig, axes = plt.subplots(2, 1, figsize=(13, 6), constrained_layout=True)
            for ax, left, right, title in ((axes[0], 0, max(end, 4647), "Whole plasmid"),
                                           (axes[1], start, end, "Mutable WT element, 211 bp")):
                for field, y, color, label in (("whole_at_positions0", 0, "0.7", "Aligned A/T"),
                                                ("whole_m6a_positions0", .25, "tab:blue", "Called m6A, ML≥128")):
                    points = [pos for pos in rep[field] if left <= pos < right]
                    ax.scatter(points, [y] * len(points), s=2 if left == 0 else 7,
                               color=color, label=label, rasterized=True)
                for call in rep["calls"]:
                    y = {"tf": 1, "nuc": 1.5, "msp": 2}[call["track"]]
                    color = {"tf": "tab:green", "nuc": "tab:red", "msp": "tab:purple"}[call["track"]]
                    for lo, hi in call["ref_intervals0"]:
                        a, b = max(left, lo), min(right, hi)
                        if a < b:
                            ax.add_patch(Rectangle((a, y - .09), b - a, .18, color=color, alpha=.6))
                ax.set(xlim=(left, right), ylim=(-.15, 2.3), yticks=[0, .25, 1, 1.5, 2],
                       yticklabels=["A/T", "m6A", "TF", "Nuc", "MSP"], xlabel="Reference position, 0-based", title=title)
            fig.suptitle(prefix + f" · {rep['representative_group']} · read {rep['read_id']}\n"
                         "Original MA calls; projected reference spans; colors denote track, not calibrated probability", fontsize=10)
            safe = str(number).zfill(2) + "_" + rep["representative_group"]
            fig.savefig(rep_dir / (safe + ".png"), dpi=140)
            fig.savefig(rep_dir / (safe + ".svg"))
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--allow-pilot", action="store_true")
    args = parser.parse_args()
    result = summarize(args.config, args.sample, args.output, args.allow_pilot)
    print(json.dumps({key: result[key] for key in ("complete", "sampled_barcodes", "sampled_molecules")}, sort_keys=True))


if __name__ == "__main__":
    main()

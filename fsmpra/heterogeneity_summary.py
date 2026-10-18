"""Descriptive, streaming summaries of replicate-3 heterogeneity features.

Numerical distal bins are observational display bins, not inferred cell states.
Primary-cohort molecular distributions give each barcode total weight one while
retaining every molecule's measured value.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from .io import atomic_output, open_text, read_json, source_info, write_json, write_table
from .processing import check_sources

COHORTS = ("sequence_qualified", "footprint_analysis_eligible",
           "footprint_analysis_eligible_with_barcode")
GROUPS = ("WT", "pooled_SNV")
FEATURES = (
    "whole_m6a_fraction", "whole_m6a_fraction_204", "whole_tf_small_fraction",
    "whole_nucleosome_fraction", "whole_msp_fraction",
    "element_m6a_fraction", "element_m6a_fraction_204", "element_tf_small_fraction",
    "element_nucleosome_fraction", "element_msp_fraction",
    "distal_m6a_fraction", "distal_m6a_fraction_204", "distal_tf_small_fraction",
    "distal_nucleosome_fraction", "distal_msp_fraction",
    "distal_tile_methylation_q10", "distal_tile_methylation_median",
    "distal_tile_methylation_q90", "distal_tile_methylation_range",
    "whole_nucleosome_segments", "whole_nucleosome_length_median",
    "whole_nucleosome_spacing_median", "whole_tf_calls_all", "whole_tf_calls_small",
    "distal_tiles",
)
FRACTIONS = FEATURES[:19]
COUNT_FEATURES = ("whole_nucleosome_segments", "whole_tf_calls_small")
COVERAGE_FEATURES = {region: region + "_coverage_fraction" for region in ("whole", "element", "distal")}
PAIRS = (
    ("distal_m6a_fraction", "distal_nucleosome_fraction"),
    ("whole_m6a_fraction", "whole_nucleosome_fraction"),
    ("element_msp_fraction", "distal_nucleosome_fraction"),
    ("element_tf_small_fraction", "distal_m6a_fraction"),
)
PRIMARY = COHORTS[-1]
TF = "element_tf_small_fraction"
MSP = "element_msp_fraction"
DISTAL_M6A = "distal_m6a_fraction"
DISTAL_NUC = "distal_nucleosome_fraction"


def _number(raw):
    if raw in (None, ""):
        return math.nan
    value = float(raw)
    if not math.isfinite(value):
        raise ValueError("Nonfinite feature value")
    return value


def _values(row):
    values = np.array([_number(row.get(f)) for f in FEATURES], dtype=np.float64)
    lo = values[FEATURES.index("distal_tile_methylation_q10")]
    hi = values[FEATURES.index("distal_tile_methylation_q90")]
    values[FEATURES.index("distal_tile_methylation_range")] = hi - lo if np.isfinite(lo + hi) else math.nan
    for value in values[:len(FRACTIONS)]:
        if math.isfinite(value) and not 0 <= value <= 1:
            raise ValueError("Fraction outside [0,1]")
    return values


def _bin(value, edges):
    if not math.isfinite(value):
        return None
    if not 0 <= value <= 1:
        raise ValueError("Fraction outside [0,1]")
    return min(int(np.searchsorted(edges, value, side="right") - 1), len(edges) - 2)


class Moments:
    def __init__(self):
        n = len(FEATURES)
        self.rows = 0
        self.valid = np.zeros(n, dtype=np.int64)
        self.total = np.zeros(n, dtype=np.float64)
        self.square = np.zeros(n, dtype=np.float64)

    def add(self, values):
        self.rows += 1
        valid = np.isfinite(values)
        self.valid += valid
        self.total[valid] += values[valid]
        self.square[valid] += values[valid] ** 2


class Barcodes:
    """Explicit bounded arrays; keys are assigned sample/library/barcode tuples."""

    def __init__(self, maximum):
        self.maximum = maximum
        self.capacity = min(8192, maximum)
        self.keys = []
        self.index = {}
        self.variants = []
        self.depth = np.zeros(self.capacity, dtype=np.int32)
        self.valid = np.zeros((self.capacity, len(FEATURES)), dtype=np.int32)
        self.total = np.zeros((self.capacity, len(FEATURES)), dtype=np.float64)
        self.square = np.zeros((self.capacity, len(FEATURES)), dtype=np.float64)
        self.state_mask = np.zeros(self.capacity, dtype=np.uint32)

    def _grow(self):
        if len(self.keys) < self.capacity:
            return
        if self.capacity >= self.maximum:
            raise ValueError("Primary barcode count exceeds summary.max_barcodes")
        new = min(self.maximum, max(self.capacity + 1, int(self.capacity * 1.5)))
        for name in ("depth", "state_mask", "valid", "total", "square"):
            array = getattr(self, name)
            shape = (new,) + array.shape[1:]
            expanded = np.zeros(shape, dtype=array.dtype)
            expanded[:self.capacity] = array
            setattr(self, name, expanded)
        self.capacity = new

    def add(self, row, values):
        barcode = row["corrected_barcode"]
        variant = row["assigned_variant_id"]
        if not barcode or not variant or variant != row["observed_variant_id"]:
            raise ValueError("Primary molecule lacks matching supported barcode/variant")
        key = (row["sample_id"], row["library_id"], barcode)
        index = self.index.get(key)
        if index is None:
            self._grow()
            index = len(self.keys)
            self.keys.append(key)
            self.variants.append(variant)
            self.index[key] = index
        elif self.variants[index] != variant:
            raise ValueError("Barcode assigned to multiple variants")
        self.depth[index] += 1
        valid = np.isfinite(values)
        self.valid[index, valid] += 1
        self.total[index, valid] += values[valid]
        self.square[index, valid] += values[valid] ** 2
        return index

    def find(self, row):
        key = (row["sample_id"], row["library_id"], row["corrected_barcode"])
        index = self.index.get(key)
        if index is None or self.variants[index] != row["assigned_variant_id"]:
            raise ValueError("Primary barcode changed between passes")
        return index


def _rows(path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = (set(FEATURES) - {"distal_tile_methylation_range"}) | set(COVERAGE_FEATURES.values())
        required |= {"sample_id", "library_id", "observed_variant_id", "genotype_class",
                     "corrected_barcode", "assigned_variant_id", "footprint_analysis_eligible",
                     "footprint_analysis_eligible_with_barcode"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Feature table missing fields: " + str(sorted(required - set(reader.fieldnames or []))))
        for row in reader:
            yield row


def _hist_key(cohort, group, weighting, feature):
    return (cohort, group, weighting, feature)


def _hist_add(hist, missing, boundary, key, value, weight, bins):
    if not math.isfinite(value):
        missing[key] += weight
    else:
        if not 0 <= value <= 1:
            raise ValueError("Fraction outside [0,1]")
        index = min(int(value * bins), bins - 1)
        hist[key][index] += weight
        if value == 0:
            boundary[key][0] += weight
        if value == 1:
            boundary[key][1] += weight


def _joint_add(joint, missing, key, x, y, weight, bins):
    if not math.isfinite(x + y):
        missing[key] += weight
    else:
        if not 0 <= x <= 1 or not 0 <= y <= 1:
            raise ValueError("Fraction outside [0,1]")
        joint[key][min(int(x * bins), bins - 1), min(int(y * bins), bins - 1)] += weight


def _count_add(hist, missing, key, value, weight):
    if not math.isfinite(value):
        missing[key] += weight
    else:
        if value < 0 or value > 100000 or value != int(value):
            raise ValueError("Invalid discrete call count")
        hist[key][int(value)] += weight


def _coverage_add(summary, cohort, group, row):
    for region, feature in COVERAGE_FEATURES.items():
        value = _number(row[feature])
        if math.isfinite(value) and not 0 <= value <= 1:
            raise ValueError("Coverage fraction outside [0,1]")
        for label in ("all", group):
            counts = summary[(cohort, label, region)]
            counts["rows"] += 1
            counts["valid"] += int(math.isfinite(value))
            counts["complete"] += int(value == 1)


def _mean(total, valid):
    return float(total / valid) if valid else None


def _source_label(catalog, variant):
    return catalog.get(variant, {}).get("source_variant", variant)


def summarize(config, features_dir, output, allow_partial=False):
    """Summarize one immutable feature directory into a new atomic directory."""
    config = Path(config)
    features_dir = Path(features_dir)
    output = Path(output)
    cfg = read_json(config)
    summary_cfg = cfg["summary"]
    maximum = int(summary_cfg["max_barcodes"])
    bins = int(summary_cfg["histogram_bins"])
    edges = np.array(summary_cfg["state_edges"], dtype=np.float64)
    if maximum < 1 or bins < 2 or len(edges) < 3 or edges[0] != 0 or edges[-1] != 1 or np.any(np.diff(edges) <= 0):
        raise ValueError("Invalid summary bins or barcode capacity")
    if (len(edges) - 1) ** 2 > 32:
        raise ValueError("State bins exceed uint32 barcode mask")
    audit_path = features_dir / "audit.json"
    feature_path = features_dir / "molecule_heterogeneity.tsv.gz"
    audit = read_json(audit_path)
    if audit.get("schema_version") != "fsmpra.heterogeneity.v1":
        raise ValueError("Expected heterogeneity extraction audit")
    if not audit.get("complete") and not allow_partial:
        raise ValueError("Partial extraction requires --allow-partial")
    if audit.get("complete") and audit.get("limit_primary") is not None:
        raise ValueError("Complete extraction has a primary-row limit")
    source_records = [source_info(config, True), source_info(audit_path, True),
                      source_info(feature_path), source_info(__file__, True)]
    if audit.get("complete"):
        check_sources(audit.get("sources", []))
    catalog_path = Path(cfg["ledger_run_dir"]) / "catalog" / "catalog.json"
    catalog_data = read_json(catalog_path)
    catalog = {v["variant_id"]: v for v in catalog_data["library"]["variants"]}
    source_records.append(source_info(catalog_path, True))
    for current in (source_records[0], source_records[-1]):
        frozen = next((item for item in audit.get("sources", []) if item["path"] == current["path"]), None)
        if frozen != current:
            raise ValueError("Summary config/catalog differs from extraction audit: " + current["path"])
    barcodes = Barcodes(maximum)
    molecule = defaultdict(Moments)
    cohort_counts = Counter()
    variant_class = {}
    hist = defaultdict(lambda: np.zeros(bins, dtype=np.float64))
    hist_missing = defaultdict(float)
    hist_boundary = defaultdict(lambda: np.zeros(2, dtype=np.float64))
    joint = defaultdict(lambda: np.zeros((bins, bins), dtype=np.float64))
    joint_missing = defaultdict(float)
    count_hist = defaultdict(Counter)
    count_missing = defaultdict(float)
    coverage_summary = defaultdict(Counter)
    n_rows = 0
    for row in _rows(feature_path):
        n_rows += 1
        if n_rows % 250000 == 0:
            print(f"Heterogeneity summary pass 1: {n_rows:,} feature rows", flush=True)
        variant = row["observed_variant_id"]
        group = "WT" if row["genotype_class"] == "WT" else "pooled_SNV"
        if row["genotype_class"] not in ("WT", "SNV") or variant not in catalog:
            raise ValueError("Unexpected selected genotype or variant")
        if variant in variant_class and variant_class[variant] != group:
            raise ValueError("Variant class changed within table")
        variant_class[variant] = group
        values = _values(row)
        cohorts = [COHORTS[0]]
        if row["footprint_analysis_eligible"] == "1":
            cohorts.append(COHORTS[1])
        if row["footprint_analysis_eligible_with_barcode"] == "1":
            if COHORTS[1] not in cohorts:
                raise ValueError("Barcode cohort not nested in footprint cohort")
            cohorts.append(PRIMARY)
            barcodes.add(row, values)
        for cohort in cohorts:
            cohort_counts[cohort] += 1
            _coverage_add(coverage_summary, cohort, group, row)
            for item in (("variant", cohort, variant), ("group", cohort, group)):
                molecule[item].add(values)
            for i, feature in enumerate(FRACTIONS):
                _hist_add(hist, hist_missing, hist_boundary,
                          _hist_key(cohort, group, "molecule", feature), values[i], 1, bins)
            for x_name, y_name in PAIRS:
                _joint_add(joint, joint_missing, (cohort, group, "molecule", x_name, y_name),
                           values[FEATURES.index(x_name)], values[FEATURES.index(y_name)], 1, bins)
            for feature in COUNT_FEATURES:
                _count_add(count_hist, count_missing, _hist_key(cohort, group, "molecule", feature),
                           values[FEATURES.index(feature)], 1)
    expected = audit["counts"]
    if n_rows != expected["selected_rows"]:
        raise ValueError("Feature row count differs from extraction audit")
    for cohort, key in ((COHORTS[1], "selected_footprint_analysis_eligible"),
                        (PRIMARY, "selected_footprint_analysis_eligible_with_barcode")):
        if cohort_counts[cohort] != expected[key]:
            raise ValueError("Cohort count differs from extraction audit: " + cohort)
    if sum(int(v) for v in barcodes.depth[:len(barcodes.keys)]) != cohort_counts[PRIMARY]:
        raise ValueError("Primary barcode depths do not conserve molecules")
    if audit["complete"] and (expected["primary_rows"] != 8715827 or
                              n_rows != 2660013 or
                              cohort_counts[COHORTS[1]] != 2594368 or
                              cohort_counts[PRIMARY] != 2102906 or len(barcodes.keys) != 568354):
        raise ValueError("Complete replicate-3 cohort does not reconcile")

    # Equal-barcode variant moments, with each feature's own valid denominator.
    # Second pass applies molecular weights 1/depth; barcode means are never histogrammed.
    state = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0, 0.0])
    state_missing_weight = Counter()
    state_missing_molecules = Counter()
    binned = Counter()
    primary_weight = 0.0
    second_rows = 0
    for row in _rows(feature_path):
        second_rows += 1
        if second_rows % 250000 == 0:
            print(f"Heterogeneity summary pass 2: {second_rows:,} feature rows", flush=True)
        if row["footprint_analysis_eligible_with_barcode"] != "1":
            continue
        i = barcodes.find(row)
        variant = row["observed_variant_id"]
        group = variant_class[variant]
        weight = 1 / int(barcodes.depth[i])
        primary_weight += weight
        values = _values(row)
        for j, feature in enumerate(FRACTIONS):
            _hist_add(hist, hist_missing, hist_boundary,
                      _hist_key(PRIMARY, group, "equal_barcode", feature), values[j], weight, bins)
        for x_name, y_name in PAIRS:
            _joint_add(joint, joint_missing, (PRIMARY, group, "equal_barcode", x_name, y_name),
                       values[FEATURES.index(x_name)], values[FEATURES.index(y_name)], weight, bins)
        for feature in COUNT_FEATURES:
            _count_add(count_hist, count_missing, _hist_key(PRIMARY, group, "equal_barcode", feature),
                       values[FEATURES.index(feature)], weight)
        x = values[FEATURES.index(DISTAL_M6A)]
        y = values[FEATURES.index(DISTAL_NUC)]
        if not math.isfinite(x + y):
            state_missing_weight[variant] += weight
            state_missing_molecules[variant] += 1
            continue
        x_bin, y_bin = _bin(x, edges), _bin(y, edges)
        state_index = x_bin * (len(edges) - 1) + y_bin
        barcodes.state_mask[i] |= np.uint32(1 << state_index)
        bucket = state[(variant, state_index)]
        bucket[0] += 1
        bucket[1] += weight
        tf = values[FEATURES.index(TF)]
        msp = values[FEATURES.index(MSP)]
        if math.isfinite(tf):
            bucket[2] += weight * tf
            bucket[3] += weight
        if math.isfinite(msp):
            bucket[4] += weight * msp
            bucket[5] += weight
        binned[variant] += 1
    state_unique = Counter()
    if second_rows != n_rows:
        raise ValueError("Feature rows changed between passes")
    for i, mask in enumerate(barcodes.state_mask[:len(barcodes.keys)]):
        mask = int(mask)
        while mask:
            bit = (mask & -mask).bit_length() - 1
            state_unique[(barcodes.variants[i], bit)] += 1
            mask &= mask - 1

    with atomic_output(output) as out:
        barcode_variant = defaultdict(Moments)
        barcode_sd_sum = {group: np.zeros(len(FEATURES), dtype=np.float64) for group in GROUPS}
        barcode_sd_count = {group: np.zeros(len(FEATURES), dtype=np.int64) for group in GROUPS}
        def barcode_rows():
            for i, key in enumerate(barcodes.keys):
                valid = barcodes.valid[i]
                means = np.full(len(FEATURES), np.nan)
                sd = np.full(len(FEATURES), np.nan)
                present = valid > 0
                spread_valid = valid >= 2
                means[present] = barcodes.total[i, present] / valid[present]
                sd[spread_valid] = np.sqrt(np.maximum(0, barcodes.square[i, spread_valid] / valid[spread_valid]
                                                 - means[spread_valid] ** 2))
                barcode_variant[barcodes.variants[i]].add(means)
                group = variant_class[barcodes.variants[i]]
                barcode_sd_sum[group][spread_valid] += sd[spread_valid]
                barcode_sd_count[group][spread_valid] += 1
                out_row = {"sample_id": key[0], "library_id": key[1], "corrected_barcode": key[2],
                           "variant_id": barcodes.variants[i], "source_variant": _source_label(catalog, barcodes.variants[i]),
                           "molecules": int(barcodes.depth[i])}
                for j, feature in enumerate(FEATURES):
                    out_row[feature + "_valid_molecules"] = int(valid[j])
                    out_row[feature + "_mean"] = float(means[j]) if present[j] else None
                    out_row[feature + "_within_sd"] = float(sd[j]) if spread_valid[j] else None
                yield out_row
        fields = ["sample_id", "library_id", "corrected_barcode", "variant_id", "source_variant", "molecules"]
        for feature in FEATURES:
            fields += [feature + suffix for suffix in ("_valid_molecules", "_mean", "_within_sd")]
        write_table(out / "barcode_heterogeneity.tsv.gz", barcode_rows(), fields)

        wt_ids = [v for v, cls in variant_class.items() if cls == "WT"]
        if len(wt_ids) != 1:
            raise ValueError("Expected one WT variant")
        wt_id = wt_ids[0]
        variant_rows = []
        for cohort in COHORTS:
            for variant in sorted(variant_class):
                m = molecule.get(("variant", cohort, variant))
                if not m:
                    continue
                b = barcode_variant.get(variant) if cohort == PRIMARY else None
                wt_m = molecule.get(("variant", cohort, wt_id))
                wt_b = barcode_variant.get(wt_id) if cohort == PRIMARY else None
                for j, feature in enumerate(FEATURES):
                    molecule_mean = _mean(m.total[j], m.valid[j])
                    wt_molecule = _mean(wt_m.total[j], wt_m.valid[j]) if wt_m else None
                    barcode_mean = _mean(b.total[j], b.valid[j]) if b else None
                    wt_barcode = _mean(wt_b.total[j], wt_b.valid[j]) if wt_b else None
                    variant_rows.append({
                        "cohort": cohort, "variant_id": variant, "source_variant": _source_label(catalog, variant),
                        "genotype_group": variant_class[variant], "feature": feature,
                        "molecules": m.rows, "valid_molecules": int(m.valid[j]),
                        "missing_molecules": int(m.rows - m.valid[j]), "molecule_mean": molecule_mean,
                        "molecule_delta_vs_WT": molecule_mean - wt_molecule if molecule_mean is not None and wt_molecule is not None else None,
                        "barcodes": b.rows if b else None, "valid_barcodes": int(b.valid[j]) if b else None,
                        "barcode_mean": barcode_mean,
                        "barcode_delta_vs_WT": barcode_mean - wt_barcode if barcode_mean is not None and wt_barcode is not None else None,
                    })
        variant_fields = list(variant_rows[0]) if variant_rows else []
        write_table(out / "variant_heterogeneity.tsv.gz", variant_rows, variant_fields)

        hist_rows = []
        for key in sorted(set(hist) | set(hist_missing)):
            cohort, group, weighting, feature = key
            for i, count in enumerate(hist[key]):
                hist_rows.append(dict(cohort=cohort, genotype_group=group, weighting=weighting, feature=feature,
                                      bin_low=i / bins, bin_high=(i + 1) / bins,
                                      count_weight=float(count), missing_weight=hist_missing[key] if i == 0 else 0))
        write_table(out / "histograms.tsv.gz", hist_rows,
                    ["cohort", "genotype_group", "weighting", "feature", "bin_low", "bin_high", "count_weight", "missing_weight"])
        boundary_rows = []
        for key in sorted(set(hist) | set(hist_missing)):
            cohort, group, weighting, feature = key
            valid_weight = float(hist[key].sum())
            boundary_rows.append(dict(cohort=cohort, genotype_group=group, weighting=weighting,
                feature=feature, valid_weight=valid_weight, missing_weight=hist_missing[key],
                exact_zero_weight=float(hist_boundary[key][0]),
                exact_one_weight=float(hist_boundary[key][1]),
                exact_zero_proportion=_mean(hist_boundary[key][0], valid_weight),
                exact_one_proportion=_mean(hist_boundary[key][1], valid_weight)))
        write_table(out / "histogram_boundaries.tsv.gz", boundary_rows, list(boundary_rows[0]))
        count_rows = []
        for key in sorted(set(count_hist) | set(count_missing)):
            cohort, group, weighting, feature = key
            for count_value, count_weight in sorted(count_hist[key].items()):
                count_rows.append(dict(cohort=cohort, genotype_group=group, weighting=weighting,
                    feature=feature, count_value=count_value, count_weight=count_weight,
                    missing_weight=0))
            if count_missing[key] or not count_hist[key]:
                count_rows.append(dict(cohort=cohort, genotype_group=group, weighting=weighting,
                    feature=feature, count_value=None, count_weight=0,
                    missing_weight=count_missing[key]))
        write_table(out / "count_histograms.tsv.gz", count_rows,
                    ["cohort", "genotype_group", "weighting", "feature", "count_value",
                     "count_weight", "missing_weight"])
        joint_rows = []
        for key in sorted(set(joint) | set(joint_missing)):
            cohort, group, weighting, x_name, y_name = key
            for x_i in range(bins):
                for y_i in range(bins):
                    joint_rows.append(dict(cohort=cohort, genotype_group=group, weighting=weighting,
                                           x_feature=x_name, y_feature=y_name,
                                           x_low=x_i / bins, x_high=(x_i + 1) / bins,
                                           y_low=y_i / bins, y_high=(y_i + 1) / bins,
                                           count_weight=float(joint[key][x_i, y_i]),
                                           missing_weight=joint_missing[key] if x_i == y_i == 0 else 0))
        write_table(out / "joint_histograms.tsv.gz", joint_rows,
                    ["cohort", "genotype_group", "weighting", "x_feature", "y_feature", "x_low", "x_high", "y_low", "y_high", "count_weight", "missing_weight"])

        n_edge = len(edges) - 1
        wt_state = {s: state.get((wt_id, s), [0] * 6) for s in range(n_edge ** 2)}
        state_rows = []
        for variant in sorted(barcode_variant):
            denom = barcode_variant[variant].rows
            variant_bins = [state.get((variant, s), [0] * 6) for s in range(n_edge ** 2)]
            if sum(bucket[0] for bucket in variant_bins) + state_missing_molecules[variant] != molecule[("variant", PRIMARY, variant)].rows:
                raise ValueError("State molecules do not conserve variant: " + variant)
            if abs(sum(bucket[1] for bucket in variant_bins) + state_missing_weight[variant] - denom) > 1e-7:
                raise ValueError("State barcode weight does not conserve variant: " + variant)
            for s in range(n_edge ** 2):
                bucket = state.get((variant, s), [0] * 6)
                wt_bucket = wt_state[s]
                tf_mean = _mean(bucket[2], bucket[3])
                msp_mean = _mean(bucket[4], bucket[5])
                wt_tf = _mean(wt_bucket[2], wt_bucket[3])
                wt_msp = _mean(wt_bucket[4], wt_bucket[5])
                state_rows.append(dict(variant_id=variant, source_variant=_source_label(catalog, variant),
                    genotype_group=variant_class[variant], distal_m6a_low=float(edges[s // n_edge]),
                    distal_m6a_high=float(edges[s // n_edge + 1]),
                    distal_nucleosome_low=float(edges[s % n_edge]),
                    distal_nucleosome_high=float(edges[s % n_edge + 1]),
                    n_molecules=bucket[0], unique_barcodes=state_unique[(variant, s)],
                    sum_weight=bucket[1], bin_proportion=bucket[1] / denom,
                    local_tf_valid_weight=bucket[3], local_tf_mean=tf_mean,
                    local_tf_delta_vs_WT=tf_mean - wt_tf if tf_mean is not None and wt_tf is not None else None,
                    local_msp_valid_weight=bucket[5], local_msp_mean=msp_mean,
                    local_msp_delta_vs_WT=msp_mean - wt_msp if msp_mean is not None and wt_msp is not None else None))
        write_table(out / "variant_state_summaries.tsv.gz", state_rows, list(state_rows[0]))

        _plot(out, hist, joint, state, state_unique, wt_id, variant_rows, barcode_variant, bins, edges,
              int(summary_cfg["plot_min_barcodes"]), partial=not audit["complete"])
        check_sources(source_records)
        if audit.get("complete"):
            check_sources(audit.get("sources", []))
        audit_summary = {
            "schema_version": "fsmpra.heterogeneity_summary.v1", "complete": bool(audit["complete"]),
            "partial": not bool(audit["complete"]), "extraction_limit_primary": audit.get("limit_primary"),
            "extraction_primary_rows": expected.get("primary_rows"), "sources": source_records,
            "cohort_molecules": dict(cohort_counts), "primary_barcodes": len(barcodes.keys),
            "primary_molecules_from_barcode_depths": int(sum(barcodes.depth[:len(barcodes.keys)])),
            "primary_equal_barcode_weight": primary_weight,
            "histogram_bins": bins, "state_edges": edges.tolist(),
            "state_binned_molecules": int(sum(binned.values())),
            "state_missing_joint_molecules": int(cohort_counts[PRIMARY] - sum(binned.values())),
            "state_binned_weight": float(sum(v[1] for v in state.values())),
            "state_missing_joint_weight": float(sum(state_missing_weight.values())),
            "SNV_variants_at_plot_threshold": sum(
                variant_class[variant] == "pooled_SNV" and value.rows >= int(summary_cfg["plot_min_barcodes"])
                for variant, value in barcode_variant.items()),
            "WT_means": {feature: {
                "molecule": _mean(molecule[("variant", PRIMARY, wt_id)].total[j], molecule[("variant", PRIMARY, wt_id)].valid[j]),
                "equal_barcode": _mean(barcode_variant[wt_id].total[j], barcode_variant[wt_id].valid[j]),
                "valid_molecules": int(molecule[("variant", PRIMARY, wt_id)].valid[j]),
                "valid_barcodes": int(barcode_variant[wt_id].valid[j]),
            } for j, feature in enumerate(FEATURES)},
            "WT_approximate_histogram_quantiles": {
                weighting: {feature: _hist_quantiles(hist[(PRIMARY, "WT", weighting, feature)])
                            for feature in FRACTIONS}
                for weighting in ("molecule", "equal_barcode")},
            "WT_exact_boundary_proportions": {
                weighting: {feature: {
                    "zero": _mean(hist_boundary[(PRIMARY, "WT", weighting, feature)][0],
                                  hist[(PRIMARY, "WT", weighting, feature)].sum()),
                    "one": _mean(hist_boundary[(PRIMARY, "WT", weighting, feature)][1],
                                 hist[(PRIMARY, "WT", weighting, feature)].sum())}
                    for feature in FRACTIONS}
                for weighting in ("molecule", "equal_barcode")},
            "WT_discrete_call_count_distribution": {
                weighting: {feature: _count_summary(
                    count_hist[(PRIMARY, "WT", weighting, feature)],
                    count_missing[(PRIMARY, "WT", weighting, feature)])
                    for feature in COUNT_FEATURES}
                for weighting in ("molecule", "equal_barcode")},
            "cohort_region_coverage": {
                cohort: {group: {region: {
                    "rows": coverage_summary[(cohort, group, region)]["rows"],
                    "valid": coverage_summary[(cohort, group, region)]["valid"],
                    "missing": (coverage_summary[(cohort, group, region)]["rows"]
                                - coverage_summary[(cohort, group, region)]["valid"]),
                    "complete": coverage_summary[(cohort, group, region)]["complete"],
                } for region in COVERAGE_FEATURES}
                    for group in ("all",) + GROUPS}
                for cohort in COHORTS},
            "cohort_feature_missingness": {
                cohort: {group: {feature: {
                    "rows": molecule[("group", cohort, group)].rows,
                    "valid": int(molecule[("group", cohort, group)].valid[j]),
                    "missing": int(molecule[("group", cohort, group)].rows - molecule[("group", cohort, group)].valid[j])}
                    for j, feature in enumerate(FEATURES)}
                    for group in GROUPS if ("group", cohort, group) in molecule}
                for cohort in COHORTS},
            "within_barcode_sd": {group: {feature: {
                "mean": _mean(barcode_sd_sum[group][j], barcode_sd_count[group][j]),
                "barcodes_with_at_least_two_valid_molecules": int(barcode_sd_count[group][j])}
                for j, feature in enumerate(FEATURES)} for group in GROUPS},
            "barcode_depth": _barcode_depth_summary(barcodes, variant_class),
            "per_variant_equal_barcode_mean_range": {feature: [
                min((m.total[j] / m.valid[j] for m in barcode_variant.values() if m.valid[j]), default=None),
                max((m.total[j] / m.valid[j] for m in barcode_variant.values() if m.valid[j]), default=None)]
                for j, feature in enumerate(FEATURES)},
            "SNV_equal_barcode_delta_range": {
                scope: {feature: _delta_range(barcode_variant, variant_class, wt_id, j,
                    min_barcodes=0 if scope == "all_SNV" else int(summary_cfg["plot_min_barcodes"]))
                    for j, feature in enumerate(FEATURES)}
                for scope in ("all_SNV", "plotted_SNV")},
            "interpretation": summary_cfg["interpretation"],
            "weighting": "Each primary molecule receives 1 / all eligible molecules in its supported barcode. Endpoint missingness does not change barcode depth.",
            "within_barcode_sd_definition": "Population SD of observed molecules within each barcode (divisor n), reported only when at least two molecules have a valid feature; not a variance component.",
            "state_interpretation": "Fixed observed numerical bins; conditional comparisons are descriptive, not causal or compartment assignments.",
        }
        if abs(audit_summary["primary_equal_barcode_weight"] - len(barcodes.keys)) > 1e-5:
            raise ValueError("Equal-barcode weights do not conserve barcodes")
        if abs(audit_summary["state_binned_weight"] + audit_summary["state_missing_joint_weight"] - len(barcodes.keys)) > 1e-5:
            raise ValueError("State-bin weights do not conserve barcodes")
        write_json(out / "summary.json", audit_summary)
    return audit_summary


def _hist_quantiles(counts):
    """Approximate quantiles at the midpoint of fixed-width histogram bins."""
    total = float(counts.sum())
    if total == 0:
        return {"q10": None, "median": None, "q90": None, "method": "histogram_bin_midpoint"}
    cumulative = np.cumsum(counts)
    n = len(counts)
    return {"q10": (int(np.searchsorted(cumulative, .1 * total)) + .5) / n,
            "median": (int(np.searchsorted(cumulative, .5 * total)) + .5) / n,
            "q90": (int(np.searchsorted(cumulative, .9 * total)) + .5) / n,
            "method": "histogram_bin_midpoint"}


def _count_summary(counts, missing_weight):
    valid_weight = float(sum(counts.values()))
    if not valid_weight:
        return {"valid_weight": 0.0, "missing_weight": float(missing_weight),
                "exact_zero_weight": 0.0, "exact_zero_proportion": None,
                "q10": None, "median": None, "q90": None,
                "quantile_method": "exact_weighted_discrete_cdf"}
    ordered = sorted(counts.items())
    quantiles = {}
    for label, probability in (("q10", .1), ("median", .5), ("q90", .9)):
        cumulative = 0.0
        for count, weight in ordered:
            cumulative += weight
            if cumulative >= probability * valid_weight - 1e-10:
                quantiles[label] = count
                break
    zero_weight = float(counts.get(0, 0))
    return {"valid_weight": valid_weight, "missing_weight": float(missing_weight),
            "exact_zero_weight": zero_weight, "exact_zero_proportion": zero_weight / valid_weight,
            "quantile_method": "exact_weighted_discrete_cdf", **quantiles}


def _delta_range(barcode_variant, variant_class, wt_id, feature_index, min_barcodes):
    wt = barcode_variant.get(wt_id)
    if wt is None or not wt.valid[feature_index]:
        return [None, None]
    reference = wt.total[feature_index] / wt.valid[feature_index]
    values = [m.total[feature_index] / m.valid[feature_index] - reference
              for variant, m in barcode_variant.items()
              if variant_class[variant] == "pooled_SNV" and m.rows >= min_barcodes and m.valid[feature_index]]
    return [float(min(values)), float(max(values))] if values else [None, None]


def _barcode_depth_summary(barcodes, variant_class):
    depths = barcodes.depth[:len(barcodes.keys)]
    result = {}
    for group in GROUPS:
        mask = np.fromiter((variant_class[variant] == group for variant in barcodes.variants),
                           dtype=np.bool_, count=len(barcodes.variants))
        selected = depths[mask]
        if not len(selected):
            result[group] = {"barcodes": 0}
            continue
        q10, median, q90 = np.quantile(selected, [.1, .5, .9])
        result[group] = {"barcodes": int(len(selected)), "n1": int(np.count_nonzero(selected == 1)),
                         "n_at_least_2": int(np.count_nonzero(selected >= 2)),
                         "n_at_least_5": int(np.count_nonzero(selected >= 5)),
                         "min": int(selected.min()), "q10": float(q10), "median": float(median),
                         "q90": float(q90), "max": int(selected.max()), "mean": float(selected.mean())}
    return result


def _plot(out, hist, joint, state, state_unique, wt_id, variant_rows, barcode_variant, bins, state_edges, min_barcodes, partial):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
    centers = (np.arange(bins) + .5) / bins
    for feature, label in ((DISTAL_M6A, "Distal observed m6A fraction"),
                           (DISTAL_NUC, "Distal nucleosome-call coverage fraction")):
        counts = hist[(PRIMARY, "WT", "equal_barcode", feature)]
        if counts.sum():
            axes[0, 0].step(centers, counts / counts.sum(), where="mid", label=label)
    axes[0, 0].set(xlabel="Observed fraction", ylabel="Proportion of valid WT barcode weight", title="WT molecular variation")
    axes[0, 0].legend(fontsize=8)
    counts = joint[(PRIMARY, "WT", "equal_barcode", DISTAL_M6A, DISTAL_NUC)]
    positive = counts[counts > 0]
    if len(positive) and counts.max() / positive.min() >= 100:
        image = axes[0, 1].imshow(np.ma.masked_less_equal(counts.T, 0), origin="lower",
                                   extent=(0, 1, 0, 1), aspect="auto", cmap="viridis",
                                   norm=LogNorm(vmin=float(positive.min()), vmax=float(counts.max())))
    else:
        image = axes[0, 1].imshow(counts.T, origin="lower", extent=(0, 1, 0, 1), aspect="auto", cmap="viridis")
    fig.colorbar(image, ax=axes[0, 1], label="Equal-barcode weight", shrink=.8)
    axes[0, 1].set(xlabel="Distal observed m6A fraction", ylabel="Distal nucleosome-call coverage fraction",
                    title="WT joint molecular distribution")
    n = len(state_edges) - 1
    for ax, sum_i, weight_i, title in ((axes[1, 0], 2, 3, "WT mean local TF-call coverage"),
                                        (axes[1, 1], 4, 5, "WT mean local MSP coverage")):
        grid = np.full((n, n), np.nan)
        for s in range(n * n):
            bucket = state.get((wt_id, s))
            if bucket and bucket[weight_i] and state_unique[(wt_id, s)] >= min_barcodes:
                grid[s % n, s // n] = bucket[sum_i] / bucket[weight_i]
        finite = grid[np.isfinite(grid)]
        if sum_i == 2 and len(finite):
            vmin, vmax = 0, max(.001, float(finite.max()))
        else:
            vmin, vmax = 0, 1
        image = ax.pcolormesh(state_edges, state_edges, np.ma.masked_invalid(grid), vmin=vmin,
                              vmax=vmax, cmap="viridis", shading="flat")
        fig.colorbar(image, ax=ax, label="Mean coverage fraction", shrink=.8)
        if not len(finite):
            ax.text(.5, .5, f"No bins with ≥{min_barcodes} WT barcodes", ha="center", va="center",
                    transform=ax.transAxes, fontsize=9)
        ax.set(xlabel="Distal observed m6A fraction (fixed bins)",
               ylabel="Distal nucleosome-call coverage fraction (fixed bins)",
               title=title + f" (≥{min_barcodes} barcodes/bin)")
    prefix = "PARTIAL PILOT — " if partial else ""
    fig.suptitle(prefix + "Observed WT heterogeneity; each barcode has total weight one", fontsize=12)
    fig.savefig(out / "wt_heterogeneity.png", dpi=170)
    fig.savefig(out / "wt_heterogeneity.svg")
    plt.close(fig)

    lookup = {(r["variant_id"], r["feature"]): r for r in variant_rows if r["cohort"] == PRIMARY}
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    comparisons = ((TF, DISTAL_NUC, "Local TF-call coverage Δ vs WT", "Distal nucleosome-call coverage Δ vs WT"),
                   (MSP, DISTAL_M6A, "Local MSP coverage Δ vs WT", "Distal observed m6A fraction Δ vs WT"))
    for ax, (x_name, y_name, x_label, y_label) in zip(axes, comparisons):
        points = []
        for variant, m in barcode_variant.items():
            if variant == wt_id or m.rows < min_barcodes:
                continue
            x = lookup.get((variant, x_name), {}).get("barcode_delta_vs_WT")
            y = lookup.get((variant, y_name), {}).get("barcode_delta_vs_WT")
            if x is not None and y is not None:
                points.append((x, y))
        if points:
            x, y = np.array(points).T
            ax.scatter(x, y, s=10, alpha=.55)
        ax.axhline(0, color="0.7", linewidth=.7)
        ax.axvline(0, color="0.7", linewidth=.7)
        ax.set(xlabel=x_label, ylabel=y_label, title=f"SNVs with ≥{min_barcodes} barcodes")
    fig.suptitle(prefix + "Observed per-variant differences; equal-barcode means", fontsize=12)
    fig.savefig(out / "variant_heterogeneity.png", dpi=170)
    fig.savefig(out / "variant_heterogeneity.svg")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args()
    result = summarize(args.config, args.features, args.output, args.allow_partial)
    print(json.dumps({k: result[k] for k in ("complete", "cohort_molecules", "primary_barcodes")}, sort_keys=True))


if __name__ == "__main__":
    main()

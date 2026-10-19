"""Deterministic WT barcode sample to audit existing FiberHMM call quality.

This is a descriptive inspection of frozen annotations. It does not recall
footprints, change eligible cohorts, or infer plasmid compartment.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from .footprint_join import merge_intervals, region_blocks, region_calls
from .heterogeneity import decode_m6a
from .io import atomic_output, open_text, read_json, source_info, write_json, write_table
from .processing import check_sources

BANDS = {
    "tf": ("lt10", "10_60", "gt60"),
    "nuc": ("lt85", "85_250", "251_500", "501_1000", "gt1000"),
    "msp": ("all",),
}
GROUPS = ("no_nuc_with_localTF", "no_nuc_without_localTF",
          "exactly_one_nuc_longer1000", "other")
BARCODE_FIELDS = ["barcode", "variant_id", "eligible_depth", "split", "weight", "selection_hash"]
CALL_FIELDS = ["read_id", "barcode", "split", "weight", "track", "band",
               "query_start0", "original_length", "ref_intervals0", "quality_spec",
               "quality", "call_name", "query_at_count", "query_m6a128_count",
               "query_m6a204_count"]
BITMAP_FIELDS = ["element_aligned_offsets0", "element_at_offsets0",
                 "element_m6a_offsets0", "element_m6a_204_offsets0",
                 "element_tf_small_offsets0", "element_msp_offsets0"]
BAND_FIELDS = [f"{region}_{track}_{band}_{suffix}"
               for region in ("whole", "element") for track in ("tf", "nuc")
               for band in BANDS[track] for suffix in ("count", "union_bases", "fraction")]
EXTRA_FIELDS = ["sample_split", "barcode_weight", "representative_group", "is_representative"] + BAND_FIELDS + BITMAP_FIELDS
LENGTH_FIELDS = ["track", "band", "original_length", "call_count", "weighted_call_count"]


def digest(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def select_barcodes(metadata_path, wt_variant_id, sample_n, seed):
    """Keep the lowest SHA256(seed|barcode) WT barcodes using bounded memory."""
    heap = []
    total_barcodes = total_depth = 0
    with open_text(metadata_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"barcode", "variant_id", "molecules"}.issubset(reader.fieldnames or []):
            raise ValueError("Invalid barcode metadata columns")
        seen = set()
        for row in reader:
            if row["variant_id"] != wt_variant_id:
                continue
            barcode = row["barcode"]
            if barcode in seen:
                raise ValueError("Duplicate WT barcode metadata: " + barcode)
            seen.add(barcode)
            depth = int(row["molecules"])
            if depth < 1:
                raise ValueError("WT barcode has no eligible molecules")
            total_barcodes += 1
            total_depth += depth
            score = digest(f"{seed}|{barcode}")
            # Negated integer makes heap[0] the current worst selected score.
            item = (-int(score, 16), barcode, depth, score)
            if len(heap) < sample_n:
                heapq.heappush(heap, item)
            elif item > heap[0]:
                heapq.heapreplace(heap, item)
    if total_barcodes < sample_n:
        raise ValueError("Fewer WT barcodes than sample size")
    sampled = {}
    for _, barcode, depth, score in heap:
        split_score = digest(f"split|{seed}|{barcode}")
        sampled[barcode] = {"barcode": barcode, "variant_id": wt_variant_id,
                            "eligible_depth": depth,
                            "split": "A" if int(split_score, 16) % 2 == 0 else "B",
                            "weight": 1 / depth, "selection_hash": score}
    return sampled, {"wt_barcodes": total_barcodes, "wt_eligible_molecules": total_depth,
                     "sampled_barcodes": len(sampled),
                     "sampled_expected_molecules": sum(row["eligible_depth"] for row in sampled.values())}


def representative_group(row):
    count = row.get("whole_nucleosome_segments", "")
    tf = row.get("element_tf_small_fraction", "")
    if count == "0" and tf:
        return "no_nuc_with_localTF" if float(tf) > 0 else "no_nuc_without_localTF"
    if count == "1" and row.get("whole_nucleosome_length_median", ""):
        if float(row["whole_nucleosome_length_median"]) > 1000:
            return "exactly_one_nuc_longer1000"
    return "other"


def choose_representatives(rows, seed, n_per_group):
    """Select stable read-ID minima from frozen features, before BAM inspection."""
    heaps = {group: [] for group in GROUPS}
    for read_id, row in rows.items():
        group = row["representative_group"]
        score = digest(f"representative|{seed}|{read_id}")
        item = (-int(score, 16), read_id)
        heap = heaps[group]
        if len(heap) < n_per_group:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)
    return {read_id for heap in heaps.values() for _, read_id in heap}


def band_for(track, length):
    if track == "tf":
        return "lt10" if length < 10 else "10_60" if length <= 60 else "gt60"
    if track == "nuc":
        if length < 85:
            return "lt85"
        if length <= 250:
            return "85_250"
        if length <= 500:
            return "251_500"
        if length <= 1000:
            return "501_1000"
        return "gt1000"
    return "all"


def interval_length(intervals):
    return sum(end - start for start, end in intervals)


def clip_intervals(spans, start, end):
    return merge_intervals([[max(a, start), min(b, end)] for a, b in spans
                            if max(a, start) < min(b, end)])


def positions(spans, start=0):
    return [pos - start for lo, hi in merge_intervals(spans) for pos in range(lo, hi)]


def query_call_counts(start, length, at_prefix, pos128, pos204):
    """Count bases and decoded m6A calls over the original stored-query span."""
    end = start + length
    if start < 0 or length <= 0 or end > len(at_prefix) - 1:
        raise ValueError(f"MA query interval outside stored SEQ: {start}, {length}")
    at = int(at_prefix[end] - at_prefix[start])
    m128 = int(np.searchsorted(pos128, end) - np.searchsorted(pos128, start))
    m204 = int(np.searchsorted(pos204, end) - np.searchsorted(pos204, start))
    return at, m128, m204


def _close(actual, expected, read_id, label):
    if expected == "":
        if actual != "":
            raise ValueError(f"{label} missingness disagrees with features: {read_id}")
    elif actual == "" or abs(float(actual) - float(expected)) > 1e-10:
        raise ValueError(f"{label} disagrees with features: {read_id}")


def audit_read(read, feature, barcode_info, spec, representative=False):
    """Return one molecule, all projected call rows, and optional plot record."""
    read_id = read.query_name
    if (read_id != feature["read_id"] or read.query_length != int(feature["read_length"])
            or ("-" if read.is_reverse else "+") != feature["strand"]):
        raise ValueError("BAM/features read identity disagreement: " + read_id)
    if feature["footprint_status"] != "evaluated" or not read.has_tag("MA"):
        raise ValueError("Expected evaluated MA on eligible molecule: " + read_id)
    ref_length = spec["reference_length"]
    element_start, element_end = spec["element"]
    blocks = region_blocks(read, 0, ref_length)
    calls, _ = region_calls(read, blocks)
    status, global_at, pos128, pos204 = decode_m6a(read)
    if (status != "ok" or feature["m6a_status"] != "ok" or global_at != int(feature["at_count"])
            or len(pos128) != int(feature["m6a_count"])):
        raise ValueError("Global MM/ML disagreement: " + read_id)
    seq = read.query_sequence or ""
    qref = np.full(len(seq), -1, dtype=np.int32)
    for qa, qb, ra in blocks:
        qref[qa:qb] = np.arange(ra, ra + qb - qa, dtype=np.int32)
    base = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    at_mask = (base == ord("A")) | (base == ord("T"))
    at_prefix = np.empty(len(seq) + 1, dtype=np.int32)
    at_prefix[0] = 0
    np.cumsum(at_mask, dtype=np.int32, out=at_prefix[1:])
    element_mask = (qref >= element_start) & (qref < element_end)
    aligned_offsets = (qref[element_mask] - element_start).tolist()
    at_offsets = (qref[element_mask & at_mask] - element_start).tolist()
    meth_offsets = (qref[pos128[(qref[pos128] >= element_start) &
                                 (qref[pos128] < element_end)]] - element_start).tolist()
    meth204_offsets = (qref[pos204[(qref[pos204] >= element_start) &
                                   (qref[pos204] < element_end)]] - element_start).tolist()
    for label, observed, field in (("aligned", len(aligned_offsets), "element_aligned_bases"),
                                   ("AT", len(at_offsets), "element_at_count"),
                                   ("m6A128", len(meth_offsets), "element_m6a_count"),
                                   ("m6A204", len(meth204_offsets), "element_m6a_count_204")):
        if observed != int(feature[field]):
            raise ValueError(f"Element {label} count disagrees with features: {read_id}")
    molecule = dict(feature)
    molecule.update(sample_split=barcode_info["split"], barcode_weight=barcode_info["weight"],
                    is_representative=int(representative))
    call_rows = []
    spans_by_band = defaultdict(list)
    counts = Counter()
    length_hist = Counter()
    tf_small_element_spans = []
    msp_element_spans = []
    msp_whole_spans = []
    for track in ("tf", "nuc", "msp"):
        for call in calls[track]:
            length = call["query_length"]
            query_at, query_m128, query_m204 = query_call_counts(
                call["query_start0"], length, at_prefix, pos128, pos204)
            band = band_for(track, length)
            spans = call["ref_intervals0"]
            element_spans = clip_intervals(spans, element_start, element_end)
            call_rows.append({"read_id": read_id, "barcode": barcode_info["barcode"],
                              "split": barcode_info["split"], "weight": barcode_info["weight"],
                              "track": track, "band": band, "query_start0": call["query_start0"],
                              "original_length": length,
                              "ref_intervals0": json.dumps(spans, separators=(",", ":")),
                              "quality_spec": call["quality_spec"],
                              "quality": json.dumps(call["quality"], separators=(",", ":")),
                              "call_name": call["name"] or "",
                              "query_at_count": query_at,
                              "query_m6a128_count": query_m128,
                              "query_m6a204_count": query_m204})
            length_hist[(track, band, length)] += 1
            if track in ("tf", "nuc"):
                counts[("whole", track, band)] += 1
                spans_by_band[("whole", track, band)].extend(spans)
                if element_spans:
                    counts[("element", track, band)] += 1
                    spans_by_band[("element", track, band)].extend(element_spans)
            if track == "tf" and band == "10_60":
                tf_small_element_spans.extend(element_spans)
            if track == "msp":
                msp_element_spans.extend(element_spans)
                msp_whole_spans.extend(spans)
    whole_aligned = int(feature["whole_aligned_bases"])
    element_aligned = int(feature["element_aligned_bases"])
    for region, denominator in (("whole", whole_aligned), ("element", element_aligned)):
        for track in ("tf", "nuc"):
            for band in BANDS[track]:
                prefix = f"{region}_{track}_{band}_"
                molecule[prefix + "count"] = counts[(region, track, band)]
                union = interval_length(merge_intervals(spans_by_band[(region, track, band)]))
                if union > denominator:
                    raise ValueError("Band union exceeds aligned bases: " + read_id)
                molecule[prefix + "union_bases"] = union
                molecule[prefix + "fraction"] = union / denominator if denominator else ""
    tf_all = sum(counts[("whole", "tf", band)] for band in BANDS["tf"])
    tf_small = counts[("whole", "tf", "10_60")]
    nuc85 = sum(counts[("whole", "nuc", band)] for band in BANDS["nuc"] if band != "lt85")
    if (tf_all != int(feature["whole_tf_calls_all"]) or
            tf_small != int(feature["whole_tf_calls_small"]) or
            nuc85 != int(feature["whole_nucleosome_segments"])):
        raise ValueError("Baseline whole-plasmid MA counts disagree with features: " + read_id)
    for region in ("whole", "element"):
        tf_fraction = molecule[f"{region}_tf_10_60_fraction"]
        nuc_spans = [span for band in BANDS["nuc"] if band != "lt85"
                     for span in spans_by_band[(region, "nuc", band)]]
        denominator = whole_aligned if region == "whole" else element_aligned
        nuc_fraction = interval_length(merge_intervals(nuc_spans)) / denominator if denominator else ""
        _close(tf_fraction, feature[f"{region}_tf_small_fraction"], read_id, region + " TF fraction")
        _close(nuc_fraction, feature[f"{region}_nucleosome_fraction"], read_id, region + " nuc fraction")
    molecule.update(element_aligned_offsets0=json.dumps(sorted(aligned_offsets), separators=(",", ":")),
                    element_at_offsets0=json.dumps(sorted(at_offsets), separators=(",", ":")),
                    element_m6a_offsets0=json.dumps(sorted(meth_offsets), separators=(",", ":")),
                    element_m6a_204_offsets0=json.dumps(sorted(meth204_offsets), separators=(",", ":")),
                    element_tf_small_offsets0=json.dumps(positions(tf_small_element_spans, element_start),
                                                         separators=(",", ":")),
                    element_msp_offsets0=json.dumps(positions(msp_element_spans, element_start),
                                                    separators=(",", ":")))
    _close(len(positions(msp_element_spans)) / element_aligned if element_aligned else "",
           feature["element_msp_fraction"], read_id, "element MSP fraction")
    _close(interval_length(merge_intervals(msp_whole_spans)) / whole_aligned if whole_aligned else "",
           feature["whole_msp_fraction"], read_id, "whole MSP fraction")
    representative_record = None
    if representative:
        ml = read.get_tag("ML")
        representative_record = {
            "read_id": read_id, "barcode": barcode_info["barcode"],
            "split": barcode_info["split"], "weight": barcode_info["weight"],
            "representative_group": feature["representative_group"],
            "strand": feature["strand"], "circular": feature["circular"],
            "feature_metadata": {key: feature.get(key, "") for key in
                                 ("passes", "read_quality", "mapq", "read_length", "plasmid_class",
                                  "m6a_status", "methylation_fraction", "whole_m6a_fraction",
                                  "element_m6a_fraction", "whole_nucleosome_segments",
                                  "whole_nucleosome_length_median", "element_tf_small_fraction")},
            "whole_aligned_intervals0": merge_intervals([[ra, ra + qb - qa] for qa, qb, ra in blocks]),
            "whole_at_positions0": sorted(qref[(qref >= 0) & at_mask].tolist()),
            "whole_m6a_positions0": sorted(qref[pos128[qref[pos128] >= 0]].tolist()),
            "whole_m6a_204_positions0": sorted(qref[pos204[qref[pos204] >= 0]].tolist()),
            "raw_ml_histogram": np.bincount(np.frombuffer(ml, dtype=np.uint8), minlength=256).tolist(),
            "raw_mm": read.get_tag("MM"), "raw_ma": read.get_tag("MA"),
            "aq_present": read.has_tag("AQ"),
            "aq_length": len(read.get_tag("AQ")) if read.has_tag("AQ") else None,
            "an_present": read.has_tag("AN"),
            "calls": [{"track": row["track"], "band": row["band"],
                       "query_start0": row["query_start0"], "original_length": row["original_length"],
                       "ref_intervals0": json.loads(row["ref_intervals0"]),
                       "quality_spec": row["quality_spec"], "quality": json.loads(row["quality"]),
                       "call_name": row["call_name"],
                       "query_at_count": row["query_at_count"],
                       "query_m6a128_count": row["query_m6a128_count"],
                       "query_m6a204_count": row["query_m6a204_count"]} for row in call_rows],
        }
    return molecule, call_rows, length_hist, representative_record


def _path(value):
    return Path(value).resolve()


def load_config(path):
    cfg = read_json(path)
    if cfg.get("schema_version") != "fsmpra.footprint_quality_config.v1":
        raise ValueError("Unsupported footprint quality configuration")
    required = ("features", "feature_audit", "barcode_metadata", "barcode_audit", "catalog",
                "ledger_config", "footprint_contract", "seed", "sample_barcodes",
                "max_molecules", "representative_n", "tf_min", "tf_max", "nuc_min", "nuc_band_edges")
    if any(key not in cfg for key in required):
        raise ValueError("Incomplete footprint quality configuration")
    for key in ("features", "feature_audit", "barcode_metadata", "barcode_audit", "catalog",
                "ledger_config", "footprint_contract"):
        cfg[key] = str(_path(cfg[key]))
    for key in ("seed", "sample_barcodes", "max_molecules", "representative_n",
                "tf_min", "tf_max", "nuc_min"):
        cfg[key] = int(cfg[key])
    if cfg["sample_barcodes"] < 1 or cfg["max_molecules"] < 1 or cfg["representative_n"] < 1:
        raise ValueError("Invalid sample or memory cap")
    if cfg["sample_barcodes"] != 10000 or cfg["seed"] != 20261018:
        raise ValueError("This audit requires the frozen 10,000-barcode sample and seed 20261018")
    if ((cfg["tf_min"], cfg["tf_max"], cfg["nuc_min"]) != (10, 60, 85)
            or cfg["nuc_band_edges"] != [85, 251, 501, 1001]):
        raise ValueError("Footprint band boundaries disagree with frozen design")
    return cfg


def run(config, output, barcode_limit=None, progress_every=1000000):
    """Stream frozen evidence, publishing only after all sample counts reconcile."""
    import pysam

    started = time.monotonic()
    config_source = source_info(config, True)
    cfg = load_config(config)
    full_n = cfg["sample_barcodes"]
    if barcode_limit is not None and not 1 <= barcode_limit < full_n:
        raise ValueError("--barcode-limit is for pilots smaller than the full sample")
    sample_n = barcode_limit or full_n
    features_path = Path(cfg["features"])
    feature_audit_path = Path(cfg["feature_audit"])
    # Snapshot all known inputs before reading audits, metadata, or features.
    sources = [config_source, source_info(cfg["catalog"], True),
               source_info(cfg["ledger_config"], True), source_info(cfg["footprint_contract"], True),
               source_info(cfg["barcode_metadata"], True), source_info(cfg["barcode_audit"], True),
               source_info(features_path), source_info(feature_audit_path, True)]
    sources += [source_info(Path(__file__).with_name(name), True) for name in
                ("footprint_quality.py", "heterogeneity.py", "footprint_join.py", "processing.py", "io.py")]
    contract = read_json(cfg["footprint_contract"])
    bam_path = _path(contract["bam"])
    sources += [source_info(bam_path), source_info(str(bam_path) + ".bai")]
    feature_audit = read_json(feature_audit_path)
    if feature_audit.get("schema_version") != "fsmpra.heterogeneity.v1" or not feature_audit.get("complete"):
        raise ValueError("Expected complete heterogeneity features")
    barcode_audit = read_json(cfg["barcode_audit"])
    if barcode_audit.get("schema_version") != "fsmpra.statistics_aggregate.v1":
        raise ValueError("Expected complete existing barcode aggregates")
    if barcode_audit["counts"]["barcodes"] != 568354 or barcode_audit["counts"]["footprint_analysis_eligible_with_barcode"] != 2102906:
        raise ValueError("Primary barcode aggregate counts changed")
    if (barcode_audit["counts"]["footprint_analysis_eligible_with_barcode"] !=
            feature_audit["counts"]["selected_footprint_analysis_eligible_with_barcode"]):
        raise ValueError("Barcode and heterogeneity eligible cohorts disagree")
    catalog = read_json(cfg["catalog"])["library"]
    wt_variants = [item["variant_id"] for item in catalog["variants"] if item["class"] == "WT"]
    if len(wt_variants) != 1:
        raise ValueError("Expected exactly one catalog WT genotype")
    wt_variant = wt_variants[0]
    frozen_bam = next((source for source in feature_audit["sources"]
                       if source["path"] == str(bam_path)), None)
    if frozen_bam is None or sources[-2] != frozen_bam:
        raise ValueError("Annotation BAM differs from frozen feature input")
    for path in (cfg["catalog"], cfg["ledger_config"], cfg["footprint_contract"]):
        frozen = next((source for source in feature_audit["sources"]
                       if source["path"] == str(_path(path))), None)
        captured = next(source for source in sources if source["path"] == str(_path(path)))
        if frozen is None or captured != frozen:
            raise ValueError("Input differs from frozen feature source: " + path)
    if feature_audit["regions"]["whole"] != [[0, catalog["reference_length"]]]:
        raise ValueError("Feature/reference length disagreement")
    element = feature_audit["regions"]["element"]
    if element != [[catalog["mutable_ref_start0"], catalog["mutable_ref_end0"]]]:
        raise ValueError("Element coordinates disagree with catalog")
    if element[0][1] - element[0][0] != 211:
        raise ValueError("Expected 211-base element")
    if barcode_audit["interval"] != element[0]:
        raise ValueError("Barcode aggregate and heterogeneity intervals disagree")
    spec = {"reference_length": catalog["reference_length"], "element": element[0]}
    sampled, selection_counts = select_barcodes(cfg["barcode_metadata"], wt_variant,
                                                 sample_n, cfg["seed"])
    if selection_counts["wt_barcodes"] != 268283 or selection_counts["wt_eligible_molecules"] != 985264:
        raise ValueError("WT primary barcode cohort changed")
    if selection_counts["sampled_expected_molecules"] > cfg["max_molecules"]:
        raise ValueError("Sample exceeds configured molecule memory cap")
    selected = {}
    selected_depth = Counter()
    with gzip.open(features_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        necessary = {"read_id", "corrected_barcode", "assigned_variant_id",
                     "footprint_analysis_eligible_with_barcode", "whole_nucleosome_segments",
                     "element_tf_small_fraction"}
        if not necessary <= set(reader.fieldnames or []):
            raise ValueError("Heterogeneity features lack required columns")
        feature_fields = list(reader.fieldnames)
        for row in reader:
            barcode = row["corrected_barcode"]
            if barcode not in sampled:
                continue
            if row["footprint_analysis_eligible_with_barcode"] != "1":
                continue
            if row["assigned_variant_id"] != wt_variant or row["genotype_class"] != "WT":
                raise ValueError("Sampled barcode assignment disagrees with WT genotype")
            read_id = row["read_id"]
            if read_id in selected:
                raise ValueError("Duplicate sampled feature read ID: " + read_id)
            row["representative_group"] = representative_group(row)
            selected[read_id] = row
            selected_depth[barcode] += 1
            if len(selected) > cfg["max_molecules"]:
                raise ValueError("Selected molecule memory cap exceeded")
    for barcode, info in sampled.items():
        if selected_depth[barcode] != info["eligible_depth"]:
            raise ValueError("Sampled barcode depth disagrees with frozen metadata: " + barcode)
    if len(selected) != selection_counts["sampled_expected_molecules"]:
        raise ValueError("Sampled molecule total disagrees with barcode metadata")
    representatives = choose_representatives(selected, cfg["seed"], cfg["representative_n"])
    counts = Counter()
    histogram = defaultdict(lambda: [0, 0.0])
    matched = set()
    molecule_fields = feature_fields + EXTRA_FIELDS
    with atomic_output(output) as out, \
            gzip.open(out / "sampled_molecules.tsv.gz", "wt", compresslevel=1, newline="") as mol_dest, \
            gzip.open(out / "sampled_calls.tsv.gz", "wt", compresslevel=1, newline="") as call_dest, \
            gzip.open(out / "representatives.jsonl.gz", "wt", compresslevel=1) as rep_dest, \
            pysam.AlignmentFile(str(bam_path), "rb", threads=1) as bam:
        if bam.get_reference_length(catalog["contig"]) != spec["reference_length"]:
            raise ValueError("Annotation BAM/reference length disagreement")
        mol_writer = csv.DictWriter(mol_dest, fieldnames=molecule_fields, delimiter="\t", extrasaction="raise")
        call_writer = csv.DictWriter(call_dest, fieldnames=CALL_FIELDS, delimiter="\t", extrasaction="raise")
        mol_writer.writeheader()
        call_writer.writeheader()
        for read in bam.fetch(catalog["contig"]):
            counts["bam_records"] += 1
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                counts["nonprimary_records"] += 1
                continue
            counts["bam_primary"] += 1
            feature = selected.get(read.query_name)
            if feature is None:
                if progress_every and counts["bam_primary"] % progress_every == 0:
                    print(f"Footprint-quality BAM: {counts['bam_primary']:,} primary, {counts['matched_molecules']:,} sampled", flush=True)
                continue
            if read.query_name in matched:
                raise ValueError("Duplicate sampled primary BAM read: " + read.query_name)
            matched.add(read.query_name)
            barcode = feature["corrected_barcode"]
            info = sampled[barcode]
            molecule, calls, lengths, representative_record = audit_read(
                read, feature, info, spec, read.query_name in representatives)
            mol_writer.writerow(molecule)
            call_writer.writerows(calls)
            for key, count in lengths.items():
                histogram[key][0] += count
                histogram[key][1] += count * info["weight"]
            if representative_record is not None:
                rep_dest.write(json.dumps(representative_record, separators=(",", ":")) + "\n")
                counts["representatives"] += 1
            counts["matched_molecules"] += 1
            counts["sampled_calls"] += len(calls)
            if progress_every and counts["bam_primary"] % progress_every == 0:
                print(f"Footprint-quality BAM: {counts['bam_primary']:,} primary, {counts['matched_molecules']:,} sampled", flush=True)
        if len(matched) != len(selected):
            raise ValueError("Sampled feature reads missing from annotation BAM")
        if counts["bam_primary"] != feature_audit["counts"]["primary_rows"]:
            raise ValueError("Annotation BAM primary count differs from complete feature audit")
        if counts["representatives"] != len(representatives):
            raise ValueError("Representative reads missing from BAM")
        write_table(out / "sampled_barcodes.tsv",
                    (sampled[barcode] for barcode in sorted(sampled)), BARCODE_FIELDS)
        write_table(out / "call_length_histograms.tsv",
                    (dict(track=track, band=band, original_length=length,
                          call_count=histogram[(track, band, length)][0],
                          weighted_call_count=histogram[(track, band, length)][1])
                     for track, band, length in sorted(histogram)), LENGTH_FIELDS)
        check_sources(sources)
        report = {"schema_version": "fsmpra.footprint_quality.v1",
                  "complete": barcode_limit is None, "barcode_limit": barcode_limit,
                  "seed": cfg["seed"], "counts": {**selection_counts, **dict(counts)},
                  "representatives_per_group": cfg["representative_n"],
                  "representative_groups": dict(Counter(selected[read_id]["representative_group"]
                                                        for read_id in representatives)),
                  "element": spec["element"], "sources": sources,
                  "seconds": time.monotonic() - started,
                  "definitions": {
                      "sampling": "Lowest SHA256(seed|barcode) among supported WT barcodes; every eligible molecule retained.",
                      "split": "SHA256(split|seed|barcode) parity; barcode-disjoint A/B.",
                      "weights": "Each sampled eligible molecule has weight 1 / full eligible barcode depth.",
                      "calls": "Original MA query lengths and query-coordinate AT/m6A counts, including unaligned query bases; band union coverage projected to aligned reference bases. Circular split MA calls remain separate.",
                      "representatives": "Lowest SHA256(representative|seed|read_id) by frozen-feature group, selected before BAM inspection.",
                      "compartment": "No compartment or naked-plasmid assignment.",
                  }}
        write_json(out / "audit.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--barcode-limit", type=int)
    parser.add_argument("--progress-every", type=int, default=1000000)
    args = parser.parse_args()
    report = run(args.config, args.output, args.barcode_limit, args.progress_every)
    print(json.dumps(report["counts"], sort_keys=True))


if __name__ == "__main__":
    main()

"""Single-preparation barcode-summary statistics for replicate 3.

Candidate marginal p-values are approximate; calibration controls release.
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

from .io import atomic_output, open_text, read_json, source_info, write_json, write_table


def _validate_config(cfg):
    frozen = dict(schema_version="fsmpra.statistics_config.v1", interval=[1234, 1445],
        tracks=["tf_small", "msp"], primary_cohort="footprint_analysis_eligible_with_barcode",
        primary_weighting="equal_barcode", tf_window_width=20, tf_window_step=10,
        include_terminal_window=True, min_barcodes=50, min_molecules=100,
        min_nonzero_barcodes=10, min_nonone_barcodes=10, confidence_level=.95,
        fdr_method="BY", fdr_alpha=.05, sensitivity_mask_radius=2, seed=17012026)
    for key, value in frozen.items():
        if cfg.get(key) != value:
            raise ValueError(f"Frozen statistical configuration differs at {key}")
    if not 1 <= cfg.get("max_barcodes", 0) <= 1_000_000:
        raise ValueError("Invalid barcode capacity")


def _calibration_report(calibration, cfg, config_path, aggregates):
    if not calibration:
        return {}
    report = read_json(calibration)
    if report.get("schema_version") != "fsmpra.statistics_calibration.v1" or report.get("sample_id") != cfg["sample_id"]:
        raise ValueError("Calibration report schema or sample mismatch")
    recorded = {item.get("path"): item for item in report.get("sources", [])}
    for path, checksum in ((config_path, True), (aggregates / "audit.json", True),
                           (Path(__file__), True), (aggregates / "barcode_counts.npy", False),
                           (aggregates / "barcode_molecules.npy", False),
                           (aggregates / "barcode_metadata.tsv.gz", True)):
        current = source_info(path, checksum)
        if recorded.get(current["path"]) != current:
            raise ValueError(f"Calibration provenance mismatch: {path}")
    return report


def fixed_endpoints(start=1234, end=1445, width=20, step=10):
    length = end - start
    if length != 211 or width != 20 or step != 10:
        raise ValueError("The frozen 211-bp/20-bp/10-bp endpoint grid is required")
    begins = list(range(0, length - width + 1, step))
    if begins[-1] != length - width:
        begins.append(length - width)
    result = [dict(endpoint_id=f"tf_{a:03d}_{a+width:03d}", track="tf_small",
                   relative_start=a, relative_end=a + width,
                   reference_start0=start + a, reference_end0=start + a + width)
              for a in begins]
    result.append(dict(endpoint_id="msp_promoter", track="msp", relative_start=0,
                       relative_end=length, reference_start0=start, reference_end0=end))
    return result


def endpoint_matrix(counts, molecules, endpoints=None):
    """Return barcode fractions, preserving one barcode as one observation."""
    endpoints = endpoints or fixed_endpoints()
    counts = np.asarray(counts)
    molecules = np.asarray(molecules)
    if counts.ndim != 3 or counts.shape[1:] != (2, 211) or molecules.shape != (counts.shape[0],):
        raise ValueError("Expected counts (barcode,track,211) and aligned molecules")
    if np.any(molecules <= 0):
        raise ValueError("Every stored barcode must have a molecule")
    result = np.empty((len(counts), len(endpoints)), dtype=np.float64)
    for first in range(0, len(counts), 16_384):
        last = min(first + 16_384, len(counts))
        cum = np.concatenate((np.zeros((last-first, 2, 1), dtype=np.uint64),
                              np.cumsum(counts[first:last], axis=2, dtype=np.uint64)), axis=2)
        for j, endpoint in enumerate(endpoints):
            track = 0 if endpoint["track"] == "tf_small" else 1
            a, b = endpoint["relative_start"], endpoint["relative_end"]
            result[first:last, j] = (cum[:, track, b] - cum[:, track, a]) / ((b - a) * molecules[first:last])
    return result


def welch_from_summary(n_variant, mean_variant, var_variant, nonzero_variant, nonone_variant,
                       n_wt, mean_wt, var_wt, nonzero_wt, nonone_wt,
                       min_barcodes=50, min_nonzero=10, min_nonone=10):
    """Welch on barcode summary statistics; sample variances use ddof=1.

    Scalars or aligned endpoint arrays are accepted. Status takes one of ok,
    insufficient_barcodes, sparse_endpoint, or degenerate_variance.
    """
    arrays = np.broadcast_arrays(*[np.asarray(x, dtype=float) for x in
        (n_variant, mean_variant, var_variant, nonzero_variant, nonone_variant,
         n_wt, mean_wt, var_wt, nonzero_wt, nonone_wt)])
    nv, mv, vv, zv, ov, nw, mw, vw, zw, ow = arrays
    delta = mv - mw
    shape = delta.shape
    status = np.full(shape, "ok", dtype="U24")
    low_n = (nv < min_barcodes) | (nw < min_barcodes)
    sparse = (zv < min_nonzero) | (ov < min_nonone) | (zw < min_nonzero) | (ow < min_nonone)
    # Floating summation can report ~1e-32 for a mathematically constant arm.
    degenerate = (~np.isfinite(vv)) | (~np.isfinite(vw)) | (vv <= 1e-30) | (vw <= 1e-30)
    status[degenerate] = "degenerate_variance"
    status[sparse] = "sparse_endpoint"
    status[low_n] = "insufficient_barcodes"
    se = np.full(shape, np.nan)
    df = np.full(shape, np.nan)
    lo = np.full(shape, np.nan)
    hi = np.full(shape, np.nan)
    p = np.full(shape, np.nan)
    valid = status == "ok"
    if np.any(valid):
        a, b = vv[valid] / nv[valid], vw[valid] / nw[valid]
        se[valid] = np.sqrt(a + b)
        df[valid] = (a + b)**2 / (a*a/(nv[valid]-1) + b*b/(nw[valid]-1))
        p[valid] = 2 * student_t.sf(np.abs(delta[valid] / se[valid]), df[valid])
        half = student_t.ppf(0.975, df[valid]) * se[valid]
        lo[valid] = np.maximum(-1., delta[valid] - half)
        hi[valid] = np.minimum(1., delta[valid] + half)
    return dict(mean_variant=mv, mean_wt=mw, delta=delta, se=se, df=df,
                ci_low=lo, ci_high=hi, p_value=p, status=status)


def welch_contrast(x, y, min_barcodes=50, min_nonzero=10, min_nonone=10):
    """Variant and WT barcode-by-endpoint fractions; no molecule-level testing."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 2 or y.ndim != 2 or x.shape[1] != y.shape[1] or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("Expected finite, aligned 2D barcode endpoint matrices")
    if np.any((x < 0) | (x > 1)) or np.any((y < 0) | (y > 1)):
        raise ValueError("Endpoint fractions must be within [0,1]")
    def summarize(a):
        n = len(a)
        return n, a.mean(axis=0) if n else np.full(a.shape[1], np.nan), \
            a.var(axis=0, ddof=1) if n > 1 else np.full(a.shape[1], np.nan), \
            (a > 0).sum(axis=0), (a < 1).sum(axis=0)
    return welch_from_summary(*summarize(x), *summarize(y), min_barcodes=min_barcodes,
                              min_nonzero=min_nonzero, min_nonone=min_nonone)


def by_adjust(p_values, family_size):
    """Joint BY including omitted members internally as p=1."""
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or len(p) != family_size or np.any((p[np.isfinite(p)] < 0) | (p[np.isfinite(p)] > 1)):
        raise ValueError("Invalid fixed-family p-values")
    internal = np.where(np.isfinite(p), p, 1.)
    order = np.argsort(internal, kind="stable")
    harmonic = np.sum(1 / np.arange(1, family_size + 1, dtype=float))
    ranked = np.minimum.accumulate((internal[order] * family_size * harmonic /
                                    np.arange(1, family_size + 1))[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = np.minimum(1., ranked)
    q[~np.isfinite(p)] = np.nan
    return q


def _intervals(raw, start, end):
    try:
        value = json.loads(raw)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValueError("Malformed coverage interval JSON") from exc
    if not isinstance(value, list):
        raise ValueError("Intervals must be a list")
    result = []
    last = None
    for item in value:
        if not isinstance(item, list) or len(item) != 2 or any(type(v) is not int for v in item):
            raise ValueError("Malformed coverage interval")
        a, b = item
        if not (start <= a < b <= end) or (last is not None and a <= last):
            raise ValueError("Unsorted, overlapping, or out-of-range coverage interval")
        result.append((a - start, b - start))
        last = b
    return result


def aggregate(config_path, catalog_path, joined_dir, output):
    cfg, catalog = read_json(config_path), read_json(catalog_path)["library"]
    _validate_config(cfg)
    start, end = cfg["interval"]
    variants = catalog["variants"]
    ids = {v["variant_id"]: i for i, v in enumerate(variants)}
    if Counter(v["class"] for v in variants)["SNV"] != 633:
        raise ValueError("Expected 633 designed SNVs")
    source = Path(joined_dir) / "molecule_footprints.tsv.gz"
    join_audit = read_json(Path(joined_dir) / "audit.json")
    if join_audit["interval"] != [start, end] or join_audit["sample_id"] != cfg["sample_id"]:
        raise ValueError("Join audit does not match statistics configuration")
    catalog_source = source_info(catalog_path, True)
    if not any(item == catalog_source for item in join_audit.get("sources", [])):
        raise ValueError("Catalog provenance disagrees with joined table")
    capacity = cfg["max_barcodes"]
    # Explicit upper bound: ~1.69 GiB for both uint32 tracks at 1m barcodes.
    counts = np.zeros((capacity, 2, end-start), dtype=np.uint32)
    molecules = np.zeros(capacity, dtype=np.uint32)
    cohort_sums = np.zeros((len(variants), 2, end-start), dtype=np.uint64)
    cohort_molecules = np.zeros(len(variants), dtype=np.uint64)
    barcode_index, barcode_ids, barcode_variants = {}, [], []
    tally = Counter()
    with open_text(source) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"footprint_analysis_eligible", "footprint_analysis_eligible_with_barcode",
                    "observed_variant_id", "assigned_variant_id", "corrected_barcode",
                    "barcode_assignment_status", "tf_small_ref_intervals0", "msp_ref_intervals0"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Missing joined-table fields")
        for row in reader:
            tally["output_molecules"] += 1
            eligible = row["footprint_analysis_eligible"] == "1"
            barcoded = row["footprint_analysis_eligible_with_barcode"] == "1"
            if row["footprint_analysis_eligible"] not in ("0", "1") or row["footprint_analysis_eligible_with_barcode"] not in ("0", "1") or (barcoded and not eligible):
                raise ValueError("Invalid cohort flags")
            if not eligible:
                continue
            vid = row["observed_variant_id"]
            if vid not in ids or variants[ids[vid]]["class"] not in ("WT", "SNV"):
                raise ValueError("Invalid eligible genotype")
            tally["footprint_analysis_eligible"] += 1
            vi = ids[vid]
            cohort_molecules[vi] += 1
            bi = None
            if barcoded:
                bc = row["corrected_barcode"]
                if row["assigned_variant_id"] != vid or row["barcode_assignment_status"] not in ("supported", "corrected_low_quality_substitution") or len(bc) != catalog["barcode_length"] or set(bc) - set("ACGT"):
                    raise ValueError("Unsupported barcode identity or genotype inconsistency")
                bi = barcode_index.get(bc)
                if bi is None:
                    bi = len(barcode_ids)
                    if bi >= capacity:
                        raise MemoryError("Barcode capacity exceeded")
                    barcode_index[bc] = bi
                    barcode_ids.append(bc)
                    barcode_variants.append(vid)
                elif barcode_variants[bi] != vid:
                    raise ValueError("Barcode maps to multiple genotypes")
                if molecules[bi] == np.iinfo(np.uint32).max:
                    raise OverflowError("Barcode molecule count overflow")
                molecules[bi] += 1
                tally["footprint_analysis_eligible_with_barcode"] += 1
            for track, field in enumerate(("tf_small_ref_intervals0", "msp_ref_intervals0")):
                for a, b in _intervals(row[field], start, end):
                    cohort_sums[vi, track, a:b] += 1
                    if bi is not None:
                        counts[bi, track, a:b] += 1
    for key in ("output_molecules", "footprint_analysis_eligible", "footprint_analysis_eligible_with_barcode"):
        if tally[key] != join_audit["counts"][key]:
            raise ValueError(f"Count conservation failed for {key}")
    n = len(barcode_ids)
    if int(molecules[:n].sum(dtype=np.uint64)) != tally["footprint_analysis_eligible_with_barcode"]:
        raise ValueError("Barcode molecule conservation failed")
    with atomic_output(output) as out:
        np.save(out / "barcode_counts.npy", counts[:n])
        np.save(out / "barcode_molecules.npy", molecules[:n])
        np.save(out / "cohort_counts.npy", cohort_sums)
        np.save(out / "cohort_molecules.npy", cohort_molecules)
        write_table(out / "barcode_metadata.tsv.gz",
                    (dict(barcode_index=i, barcode=bc, variant_id=barcode_variants[i], molecules=int(molecules[i]))
                     for i, bc in enumerate(barcode_ids)),
                    ["barcode_index", "barcode", "variant_id", "molecules"])
        write_table(out / "variant_metadata.tsv", (dict(variant_index=i, variant_id=v["variant_id"],
                     source_variant=v["source_variant"], genotype_class=v["class"],
                     ref_pos0=v["ref_pos0"], ref=v["ref"], alt=v["alt"])
                     for i, v in enumerate(variants)),
                    ["variant_index", "variant_id", "source_variant", "genotype_class", "ref_pos0", "ref", "alt"])
        write_json(out / "audit.json", dict(schema_version="fsmpra.statistics_aggregate.v1",
            sample_id=cfg["sample_id"], interval=[start, end], counts={**dict(tally), "barcodes": n},
            sources=[source_info(p, p != source) for p in
                     (config_path, catalog_path, Path(joined_dir) / "audit.json", source, __file__)],
            max_barcodes=capacity, tracks=["tf_small", "msp"],
            interpretation="One count per molecule per covered reference base; no molecule deduplication."))
    return dict(tally)


def _summary(a):
    n = len(a)
    return (n, a.mean(axis=0) if n else np.full(a.shape[1], np.nan),
            a.var(axis=0, ddof=1) if n > 1 else np.full(a.shape[1], np.nan),
            (a > 0).sum(axis=0), (a < 1).sum(axis=0))


def _float(value):
    value = float(value)
    return value if np.isfinite(value) else None


def _at_change(ref, alt):
    was, now = ref in "AT", alt in "AT"
    return "gain" if now and not was else "loss" if was and not now else "unchanged"


def infer(config_path, catalog_path, aggregates, output, calibration=None):
    cfg, catalog = read_json(config_path), read_json(catalog_path)["library"]
    _validate_config(cfg)
    variants = catalog["variants"]
    endpoints = fixed_endpoints(*cfg["interval"], cfg["tf_window_width"], cfg["tf_window_step"])
    snvs = [v for v in variants if v["class"] == "SNV"]
    if len(snvs) != 633 or len(endpoints) != 22:
        raise ValueError("The fixed 633 x 22 family changed")
    aggregates = Path(aggregates)
    audit = read_json(aggregates / "audit.json")
    if audit["sample_id"] != cfg["sample_id"] or audit["interval"] != cfg["interval"]:
        raise ValueError("Aggregate provenance disagrees with configuration")
    for path in (config_path, catalog_path):
        current = source_info(path, True)
        recorded = next((x for x in audit["sources"] if x["path"] == current["path"]), None)
        if recorded != current:
            raise ValueError("Aggregate config/catalog provenance changed")
    calibration_report = _calibration_report(calibration, cfg, config_path, aggregates)
    counts = np.load(aggregates / "barcode_counts.npy", mmap_mode="r")
    molecules = np.load(aggregates / "barcode_molecules.npy", mmap_mode="r")
    cohort_counts = np.load(aggregates / "cohort_counts.npy", mmap_mode="r")
    cohort_molecules = np.load(aggregates / "cohort_molecules.npy", mmap_mode="r")
    if counts.shape != (len(molecules), 2, 211) or cohort_counts.shape != (len(variants), 2, 211):
        raise ValueError("Invalid aggregate array shapes")
    indices = {v["variant_id"]: [] for v in variants}
    with open_text(aggregates / "barcode_metadata.tsv.gz") as handle:
        for i, row in enumerate(csv.DictReader(handle, delimiter="\t")):
            if int(row["barcode_index"]) != i or row["variant_id"] not in indices or int(row["molecules"]) != molecules[i]:
                raise ValueError("Barcode metadata and arrays disagree")
            indices[row["variant_id"]].append(i)
    if sum(map(len, indices.values())) != len(molecules):
        raise ValueError("Barcode metadata incomplete")
    fractions = endpoint_matrix(counts, molecules, endpoints)
    wt = next(v for v in variants if v["class"] == "WT")
    wt_idx = np.asarray(indices[wt["variant_id"]], dtype=np.int64)
    wt_summary = _summary(fractions[wt_idx])
    variant_index = {v["variant_id"]: i for i, v in enumerate(variants)}
    start = cfg["interval"][0]
    profiles = []
    # A per-base mean of barcode means makes the mutation mask a fixed linear
    # diagnostic without repeated WT matrix scans for every alternate allele.
    barcode_base_means = {}
    barcode_track_sums = np.zeros_like(cohort_counts)
    barcode_track_molecules = np.zeros(len(variants), dtype=np.uint64)
    for v in [wt] + snvs:
        vi = variant_index[v["variant_id"]]
        idx = np.asarray(indices[v["variant_id"]], dtype=np.int64)
        if len(idx):
            base_mean = np.zeros((2, 211), dtype=np.float64)
            for first in range(0, len(idx), 16_384):
                block_idx = idx[first:first+16_384]
                block = np.asarray(counts[block_idx])
                base_mean += (block / np.asarray(molecules[block_idx], dtype=np.float64)[:, None, None]).sum(axis=0)
                barcode_track_sums[vi] += block.sum(axis=0, dtype=np.uint64)
            base_mean /= len(idx)
            barcode_track_molecules[vi] = molecules[idx].sum(dtype=np.uint64)
        else:
            base_mean = np.full((2, 211), np.nan)
        barcode_base_means[v["variant_id"]] = base_mean
        for track, name in enumerate(("tf_small", "msp")):
            for pos in range(211):
                profiles.append(dict(variant_id=v["variant_id"], source_variant=v["source_variant"],
                    track=name, reference_pos0=start+pos, relative_pos=pos,
                    n_barcodes=len(idx), n_molecules=int(barcode_track_molecules[vi]),
                    mean_equal_barcode=_float(base_mean[track, pos]),
                    mean_molecule_weighted=_float(barcode_track_sums[vi, track, pos] /
                        barcode_track_molecules[vi]) if barcode_track_molecules[vi] else None,
                    mean_all_eligible=_float(cohort_counts[vi, track, pos] /
                        cohort_molecules[vi]) if cohort_molecules[vi] else None))
    effect_rows, sensitivity_rows, raw_p = [], [], []
    wt_vi = variant_index[wt["variant_id"]]
    wt_molecule_count = int(barcode_track_molecules[wt_vi])
    for v in snvs:
        vi = variant_index[v["variant_id"]]
        idx = np.asarray(indices[v["variant_id"]], dtype=np.int64)
        variant_molecule_count = int(barcode_track_molecules[vi])
        got = welch_from_summary(*_summary(fractions[idx]), *wt_summary,
            min_barcodes=cfg["min_barcodes"], min_nonzero=cfg["min_nonzero_barcodes"],
            min_nonone=cfg["min_nonone_barcodes"])
        pos = v["ref_pos0"] - start
        if not 0 <= pos < 211:
            raise ValueError("SNV outside frozen interval")
        mask_a = max(0, pos - cfg["sensitivity_mask_radius"])
        mask_b = min(211, pos + cfg["sensitivity_mask_radius"] + 1)
        for j, ep in enumerate(endpoints):
            a, b = ep["relative_start"], ep["relative_end"]
            track = 0 if ep["track"] == "tf_small" else 1
            overlap = a <= pos < b
            distance = 0 if overlap else min(abs(pos - a), abs(pos - (b-1)))
            status = str(got["status"][j])
            if min(variant_molecule_count, wt_molecule_count) < cfg["min_molecules"]:
                status = "insufficient_molecules"
            p = got["p_value"][j] if status == "ok" else np.nan
            raw_p.append(p)
            row = dict(variant_id=v["variant_id"], source_variant=v["source_variant"],
                endpoint_id=ep["endpoint_id"], track=ep["track"],
                reference_start0=ep["reference_start0"], reference_end0=ep["reference_end0"],
                relative_start=a, relative_end=b, ref_pos0=v["ref_pos0"],
                ref=v["ref"], alt=v["alt"], at_change=_at_change(v["ref"], v["alt"]),
                mutation_overlaps_window=int(overlap), distance_to_window_bp=distance,
                n_variant_barcodes=len(idx), n_wt_barcodes=len(wt_idx),
                n_variant_molecules=variant_molecule_count,
                n_wt_molecules=wt_molecule_count,
                mean_variant=_float(got["mean_variant"][j]), mean_wt=_float(got["mean_wt"][j]),
                delta=_float(got["delta"][j]), delta_percentage_points=_float(100*got["delta"][j]),
                se=_float(got["se"][j]) if status == "ok" else None,
                df=_float(got["df"][j]) if status == "ok" else None,
                ci_low=_float(got["ci_low"][j]) if status == "ok" else None,
                ci_high=_float(got["ci_high"][j]) if status == "ok" else None,
                candidate_p_value=_float(p), candidate_q_value=None, status=status,
                caller_effect_supported=0, statistical_supported=0,
                inference_label="candidate_approximate")
            effect_rows.append(row)
            # Distinct descriptive estimands on the same barcode cohort and on
            # all eligible molecules, including those lacking supported barcodes.
            for label, sums, totals in (("supported_molecule_weighted", barcode_track_sums, barcode_track_molecules),
                                        ("all_eligible_molecule_weighted", cohort_counts, cohort_molecules)):
                vmean = sums[vi, track, a:b].sum(dtype=np.uint64) / ((b-a)*totals[vi]) if totals[vi] else np.nan
                wmean = sums[wt_vi, track, a:b].sum(dtype=np.uint64) / ((b-a)*totals[wt_vi]) if totals[wt_vi] else np.nan
                sensitivity_rows.append(dict(variant_id=v["variant_id"], source_variant=v["source_variant"],
                    endpoint_id=ep["endpoint_id"], sensitivity=label,
                    mean_variant=_float(vmean), mean_wt=_float(wmean), delta=_float(vmean-wmean),
                    n_variant_molecules=int(totals[vi]), n_wt_molecules=int(totals[wt_vi]),
                    omitted_start0=None, omitted_end0=None))
            omit_a, omit_b = max(a, mask_a), min(b, mask_b)
            kept = (b-a) - max(0, omit_b-omit_a)
            if kept <= 0:
                raise ValueError("Mask removed entire endpoint")
            def masked_mean(base):
                if not np.all(np.isfinite(base)):
                    return np.nan
                return (base[track, a:b].sum() -
                        (base[track, omit_a:omit_b].sum() if omit_a < omit_b else 0)) / kept
            vmean, wmean = masked_mean(barcode_base_means[v["variant_id"]]), masked_mean(barcode_base_means[wt["variant_id"]])
            sensitivity_rows.append(dict(variant_id=v["variant_id"], source_variant=v["source_variant"],
                endpoint_id=ep["endpoint_id"], sensitivity="mutation_pm2_masked_equal_barcode",
                mean_variant=_float(vmean), mean_wt=_float(wmean), delta=_float(vmean-wmean),
                n_variant_molecules=variant_molecule_count,
                n_wt_molecules=wt_molecule_count,
                omitted_start0=start+omit_a if omit_a < omit_b else None,
                omitted_end0=start+omit_b if omit_a < omit_b else None))
    family_size = 633 * 22
    q = by_adjust(np.asarray(raw_p), family_size)
    validated_ids = calibration_report.get("validated_endpoint_ids", [])
    if not isinstance(validated_ids, list) or any(x not in {e["endpoint_id"] for e in endpoints} for x in validated_ids):
        raise ValueError("Invalid calibrated endpoint IDs")
    supported = calibration_report.get("status") == "pass" and bool(validated_ids)
    for row, value in zip(effect_rows, q):
        row["candidate_q_value"] = _float(value)
        released = (supported and row["endpoint_id"] in validated_ids and
                    row["status"] == "ok" and row["candidate_q_value"] is not None and
                    row["candidate_q_value"] <= cfg["fdr_alpha"])
        row["caller_effect_supported"] = int(released)
        row["statistical_supported"] = int(released)
    with atomic_output(output) as out:
        write_table(out / "effects.tsv.gz", effect_rows, list(effect_rows[0]))
        write_table(out / "profiles.tsv.gz", profiles, list(profiles[0]))
        write_table(out / "sensitivities.tsv.gz", sensitivity_rows, list(sensitivity_rows[0]))
        write_json(out / "release_status.json", dict(schema_version="fsmpra.statistics_release.v1",
            statistical_calibration=calibration_report.get("status", "absent"),
            validated_endpoint_ids=validated_ids, statistically_supported=supported,
            supported_effects=sum(x["caller_effect_supported"] for x in effect_rows),
            candidate_inference_withheld=not supported,
            biological_detectability_resolved=False,
            candidate_hypotheses=family_size, testable_hypotheses=sum(x["status"] == "ok" for x in effect_rows),
            interpretation=cfg["interpretation"]))
        write_json(out / "audit.json", dict(schema_version="fsmpra.statistics_infer.v1",
            family_size=family_size, endpoints=endpoints, counts=dict(variants=len(snvs),
            barcode_groups=len(molecules), profiles=len(profiles), effects=len(effect_rows)),
            sources=[source_info(p, Path(p).name != "barcode_counts.npy") for p in
                     (config_path, catalog_path, aggregates / "audit.json", aggregates / "barcode_counts.npy", __file__)] +
                     ([source_info(calibration, True)] if calibration else []),
            inference="Approximate Welch on observed supported barcode means; joint BY includes untestable hypotheses as p=1."))
    return dict(testable=sum(x["status"] == "ok" for x in effect_rows), supported=supported)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    a = sub.add_parser("aggregate")
    for name in ("config", "catalog", "joined-dir", "output"):
        a.add_argument("--" + name, required=True)
    b = sub.add_parser("infer")
    for name in ("config", "catalog", "aggregates", "output"):
        b.add_argument("--" + name, required=True)
    b.add_argument("--calibration")
    args = vars(parser.parse_args())
    command = args.pop("command")
    if command == "aggregate":
        print(json.dumps(aggregate(args["config"], args["catalog"], args["joined_dir"], args["output"])))
    else:
        print(json.dumps(infer(args["config"], args["catalog"], args["aggregates"], args["output"], args["calibration"])))

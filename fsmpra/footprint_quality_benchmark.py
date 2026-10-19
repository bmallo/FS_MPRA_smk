"""Fixed WT barcode A/A and simulation benchmark for sparse caller endpoints.

This diagnostic never updates variant inference, q-values, or release status.
"""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, t as student_t

from .io import atomic_output, read_json, source_info, write_json, write_table
from .processing import check_sources
from .statistics import fixed_endpoints


def welch_diagnostic(x, y):
    """Two-sided unequal-variance t test without sparse-outcome exclusion."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or min(len(x), len(y)) < 2 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("Expected finite 1D barcode vectors with at least two observations per arm")
    mx, my = float(x.mean()), float(y.mean())
    vx = 0. if np.ptp(x) == 0 else float(x.var(ddof=1))
    vy = 0. if np.ptp(y) == 0 else float(y.var(ddof=1))
    ax, ay = vx / len(x), vy / len(y)
    se = float(np.sqrt(ax + ay))
    delta = mx - my
    if se == 0:
        statistic = 0. if delta == 0 else float(np.copysign(np.inf, delta))
        p = 1. if delta == 0 else 0.
        status = "uninformative" if delta == 0 else "zero_se_difference"
        df = None
    else:
        statistic = delta / se
        df = (ax + ay)**2 / (ax*ax/(len(x)-1) + ay*ay/(len(y)-1))
        p = float(2 * student_t.sf(abs(statistic), df))
        status = "ok"
    return dict(delta=delta, se=se, statistic=statistic, df=df, p_value=p, status=status,
                nonzero_variant=int(np.count_nonzero(x)), nonzero_control=int(np.count_nonzero(y)))


def studentized_permutation(x, y, permutations, rng, batch_size=32):
    """Shuffle whole barcode labels, recomputing the Welch denominator each time.

    Two-sided Monte Carlo p=(1+exceedances)/(B+1), including ties within 1e-12.
    This is not asserted to be finite-sample exact for a mean-null with unequal
    barcode variances.
    """
    if permutations < 1 or batch_size < 1 or batch_size > 256:
        raise ValueError("Invalid bounded permutation request")
    observed = welch_diagnostic(x, y)
    if observed["status"] == "uninformative":
        return dict(p_value=1., exceedances=permutations, status="uninformative",
                    statistic=0., permutations=permutations)
    pooled = np.concatenate((np.asarray(x, dtype=float), np.asarray(y, dtype=float)))
    n, m = len(x), len(y)
    center = pooled.mean()
    values = pooled - center
    total = values.sum()
    total_sq = np.square(values).sum()
    threshold = abs(observed["statistic"]) - 1e-12
    exceedances = 0
    for first in range(0, permutations, batch_size):
        size = min(batch_size, permutations-first)
        # Each row is an independent uniform subset of n barcode labels.
        keys = rng.random((size, n+m))
        chosen = np.argpartition(keys, n-1, axis=1)[:, :n]
        selected = values[chosen]
        sx = selected.sum(axis=1)
        qx = np.square(selected).sum(axis=1)
        sy, qy = total-sx, total_sq-qx
        vx = np.maximum(0., (qx-sx*sx/n)/(n-1))
        vy = np.maximum(0., (qy-sy*sy/m)/(m-1))
        se = np.sqrt(vx/n + vy/m)
        delta = sx/n - sy/m
        stat = np.divide(delta, se, out=np.zeros_like(delta), where=se > 0)
        stat[(se == 0) & (delta != 0)] = np.sign(delta[(se == 0) & (delta != 0)]) * np.inf
        exceedances += int(np.count_nonzero(np.abs(stat) >= threshold))
    return dict(p_value=(1+exceedances)/(permutations+1), exceedances=exceedances,
                status=observed["status"], statistic=observed["statistic"],
                permutations=permutations)


def draw_fraction_barcodes(rng, n, mean, rho, depth_p):
    """Independent barcode propensities; depth is independent of propensity."""
    if not (0 <= mean <= 1 and 0 < rho < 1 and 0 < depth_p <= 1):
        raise ValueError("Invalid simulation parameters")
    depth = rng.geometric(depth_p, n)
    if mean == 0 or mean == 1:
        latent = np.full(n, mean)
    else:
        precision = 1/rho - 1
        latent = rng.beta(mean*precision, (1-mean)*precision, n)
    return rng.binomial(depth, latent) / depth


def _validate_config(cfg, pilot_iterations, pilot_permutations):
    b = cfg.get("benchmark", {})
    expected = dict(seed=20261019, iterations=200, permutations=499,
                    variant_barcode_counts=[100, 400], control_barcodes=1000,
                    windows=[[120, 140], [140, 160]], alpha=.05)
    if cfg.get("schema_version") != "fsmpra.footprint_quality_config.v1":
        raise ValueError("Unexpected footprint-quality configuration")
    for key, value in expected.items():
        if b.get(key) != value:
            raise ValueError(f"Frozen benchmark setting changed: {key}")
    iterations = pilot_iterations if pilot_iterations is not None else b["iterations"]
    permutations = pilot_permutations if pilot_permutations is not None else b["permutations"]
    if iterations < 1 or iterations > b["iterations"] or permutations < 1 or permutations > b["permutations"]:
        raise ValueError("Pilot limits must be positive and at most frozen full settings")
    return b, iterations, permutations, iterations != b["iterations"] or permutations != b["permutations"]


def _load_summary(config_path, summary_dir):
    summary_dir = Path(summary_dir)
    summary_path = summary_dir / "summary.json"
    summary_source = source_info(summary_path, True)
    report = read_json(summary_path)
    if report.get("schema_version") != "fsmpra.footprint_quality_summary.v1" or not report.get("sample_id"):
        raise ValueError("Expected a complete footprint-quality summary")
    if (report.get("complete") is not True or report.get("pilot") is not False or
            report.get("sampled_barcodes") != 10_000 or report.get("endpoint_count") != 22 or
            report.get("endpoint_matrix") != "barcode_endpoint_matrix.npz"):
        raise ValueError("Benchmark requires the complete 10,000-barcode summary")
    records = {item["path"]: item for item in report.get("sources", [])}
    if len(records) != len(report.get("sources", [])) or not records:
        raise ValueError("Duplicate or missing summary source records")
    for path in (config_path, summary_dir / "barcode_endpoint_matrix.npz"):
        current = source_info(path, True)
        if records.get(current["path"]) != current:
            raise ValueError("Summary source provenance mismatch")
    check_sources(report["sources"])
    with np.load(summary_dir / "barcode_endpoint_matrix.npz", allow_pickle=False) as data:
        expected_keys = {"barcodes", "endpoint_ids", "presence", "coverage", "valid_molecules", "depths",
                         "splits", "relative_starts", "relative_ends", "reference_starts0", "reference_ends0"}
        if not expected_keys <= set(data.files):
            raise ValueError("Incomplete endpoint matrix")
        arrays = {key: data[key] for key in expected_keys}
    endpoints = fixed_endpoints()
    ids = [x["endpoint_id"] for x in endpoints]
    if list(arrays["endpoint_ids"]) != ids:
        raise ValueError("Endpoint IDs differ from frozen grid")
    for key, field in (("relative_starts", "relative_start"), ("relative_ends", "relative_end"),
                       ("reference_starts0", "reference_start0"), ("reference_ends0", "reference_end0")):
        if list(arrays[key]) != [x[field] for x in endpoints]:
            raise ValueError("Endpoint coordinates differ from frozen grid")
    barcode = arrays["barcodes"]
    if barcode.ndim != 1 or len(barcode) != 10_000 or not np.all(barcode[:-1] < barcode[1:]):
        raise ValueError("Expected exactly 10,000 unique sorted WT barcodes")
    n = len(barcode)
    for key in ("presence", "coverage", "valid_molecules"):
        if arrays[key].shape != (n, 22):
            raise ValueError("Endpoint matrix shape mismatch")
    if arrays["depths"].shape != (n,) or arrays["splits"].shape != (n,) or not set(arrays["splits"]) == {"A", "B"}:
        raise ValueError("Invalid barcode depths or independent splits")
    if dict(Counter(arrays["splits"])) != report.get("split_barcodes"):
        raise ValueError("Matrix and summary split barcode counts differ")
    valid = arrays["valid_molecules"]
    if np.any(valid < 0) or np.any(valid > arrays["depths"][:, None]) or np.any(arrays["depths"] <= 0):
        raise ValueError("Invalid endpoint molecule denominators")
    for key in ("presence", "coverage"):
        x = arrays[key]
        if np.any(np.isfinite(x) != (valid > 0)) or np.any((x[np.isfinite(x)] < 0) | (x[np.isfinite(x)] > 1)):
            raise ValueError("Missingness or fraction range differs from molecule counts")
    if report.get("sampled_barcodes") != n:
        raise ValueError("Summary barcode count mismatch")
    if source_info(summary_path, True) != summary_source:
        raise ValueError("Summary JSON changed while being read")
    return report, arrays


def _draw_wt_groups(rng, pool, n_variant, n_control):
    """Draw unique barcode rows, with disjoint pseudo-variant and control arms."""
    selected = rng.choice(pool, n_variant + n_control, replace=False)
    return selected[:n_variant], selected[n_variant:]


def _targets(endpoints):
    targets = []
    for start, end in ((120, 140), (140, 160)):
        endpoint = next(x for x in endpoints if (x["relative_start"], x["relative_end"], x["track"]) ==
                        (start, end, "tf_small"))
        for metric in ("coverage", "presence"):
            targets.append((endpoint["endpoint_id"], metric))
    targets.append(("msp_promoter", "coverage"))
    return targets


def _row(partition, scenario, iteration, endpoint, metric, method, truth, x, y, result, is_simulation):
    return dict(partition=partition, scenario=scenario, iteration=iteration,
        data_source="simulation" if is_simulation else "wt_aa", endpoint_id=endpoint, metric=metric,
        method=method, true_delta=truth, n_variant=len(x), n_control=len(y),
        nonzero_variant=result["nonzero_variant"], nonzero_control=result["nonzero_control"],
        mean_variant=float(np.mean(x)), mean_control=float(np.mean(y)),
        delta=result["delta"], se=result["se"], statistic=result["statistic"],
        p_value=result["p_value"], status=result["status"],
        exceedances=result.get("exceedances"), permutations=result.get("permutations"))


def _contrast_rows(partition, scenario, iteration, endpoint, metric, truth, x, y, rng, permutations, simulation):
    w = welch_diagnostic(x, y)
    p = studentized_permutation(x, y, permutations, rng)
    return [_row(partition, scenario, iteration, endpoint, metric, "welch", truth, x, y, w, simulation),
            _row(partition, scenario, iteration, endpoint, metric, "studentized_permutation", truth,
                 x, y, {**w, **p}, simulation)]


def _seed(seed, partition, scenario):
    return np.random.default_rng(np.random.SeedSequence([seed, 0 if partition == "A" else 1, scenario]))


def benchmark(config_path, summary_dir, output, pilot_iterations=None, pilot_permutations=None):
    cfg = read_json(config_path)
    b, iterations, permutations, partial = _validate_config(cfg, pilot_iterations, pilot_permutations)
    summary, arrays = _load_summary(config_path, summary_dir)
    source_paths = [config_path, Path(summary_dir) / "summary.json",
                    Path(summary_dir) / "barcode_endpoint_matrix.npz", __file__,
                    Path(__file__).with_name("statistics.py")]
    sources = [source_info(path, True) for path in source_paths]
    endpoints = fixed_endpoints()
    targets = _targets(endpoints)
    endpoint_index = {name: j for j, name in enumerate(arrays["endpoint_ids"])}
    columns = [(metric, endpoint_index[eid]) for eid, metric in targets]
    # One complete-barcode pool per partition keeps all five endpoints on the
    # same disjoint pseudo-variant/control sets in each iteration.
    complete = np.ones(len(arrays["barcodes"]), dtype=bool)
    for metric, j in columns:
        complete &= np.isfinite(arrays[metric][:, j])
    records = []
    available = {}
    missing_by_target = {f"{eid}:{metric}": int(np.count_nonzero(~np.isfinite(arrays[metric][:, j])))
                         for (eid, metric), (_, j) in zip(targets, columns)}
    if any(count == len(arrays["barcodes"]) for count in missing_by_target.values()):
        raise ValueError("A fixed target has no observed WT barcodes")
    complete_total = int(np.count_nonzero(complete))
    for partition in ("A", "B"):
        pool = np.flatnonzero((arrays["splits"] == partition) & complete)
        available[partition] = len(pool)
        if len(pool) < max(b["variant_barcode_counts"]) + b["control_barcodes"]:
            raise ValueError(f"Insufficient complete WT barcodes in split {partition}")
        for scenario, n_variant in enumerate(b["variant_barcode_counts"]):
            rng = _seed(b["seed"], partition, scenario)
            for iteration in range(iterations):
                variant, control = _draw_wt_groups(rng, pool, n_variant, b["control_barcodes"])
                for (eid, metric), (_, j) in zip(targets, columns):
                    x, y = arrays[metric][variant, j], arrays[metric][control, j]
                    records.extend(_contrast_rows(partition, f"wt_aa_n{n_variant}", iteration,
                        eid, metric, 0., x, y, rng, permutations, False))
    scenario_number = 100
    for partition in ("A", "B"):
        for baseline in (.02, .1, .5):
            for n_variant in b["variant_barcode_counts"]:
                rng = _seed(b["seed"], partition, scenario_number)
                scenario = f"null_mean{baseline:g}_n{n_variant}_rhoV.02_rhoC.2"
                for iteration in range(iterations):
                    x = draw_fraction_barcodes(rng, n_variant, baseline, .02, .3)
                    y = draw_fraction_barcodes(rng, b["control_barcodes"], baseline, .2, .5)
                    records.extend(_contrast_rows(partition, scenario, iteration,
                        "synthetic_fraction", "coverage", 0., x, y, rng, permutations, True))
                scenario_number += 1
        for alt, label in ((0., "known_loss"), (.2, "known_gain")):
            rng = _seed(b["seed"], partition, scenario_number)
            for iteration in range(iterations):
                x = draw_fraction_barcodes(rng, 400, alt, .02, .3)
                y = draw_fraction_barcodes(rng, b["control_barcodes"], .1, .2, .5)
                records.extend(_contrast_rows(partition, label, iteration,
                    "synthetic_fraction", "coverage", alt-.1, x, y, rng, permutations, True))
            scenario_number += 1
    groups = defaultdict(list)
    for row in records:
        key = tuple(row[k] for k in ("partition", "scenario", "data_source", "endpoint_id", "metric", "method"))
        groups[key].append(row)
    summaries = []
    for key, rows in groups.items():
        p = np.array([r["p_value"] for r in rows])
        reject = int(np.count_nonzero(p <= b["alpha"]))
        ci = binomtest(reject, len(rows)).proportion_ci(.95, method="exact")
        summaries.append(dict(zip(("partition", "scenario", "data_source", "endpoint_id", "metric", "method"), key)) |
            dict(iterations=len(rows), true_delta=rows[0]["true_delta"],
                 rejection_rate=reject/len(rows), rejected=reject,
                 rejection_ci_low=float(ci.low), rejection_ci_high=float(ci.high),
                 interpretation="null_false_positive_rate" if rows[0]["true_delta"] == 0 else "nominal_power",
                 positive_delta=sum(r["delta"] > 0 for r in rows),
                 negative_delta=sum(r["delta"] < 0 for r in rows),
                 zero_delta=sum(r["delta"] == 0 for r in rows),
                 uninformative=sum(r["status"] == "uninformative" for r in rows),
                 zero_se_difference=sum(r["status"] == "zero_se_difference" for r in rows),
                 median_nonzero_variant=float(np.median([r["nonzero_variant"] for r in rows])),
                 median_nonzero_control=float(np.median([r["nonzero_control"] for r in rows])),
                 **{f"p_q{int(q*100):02d}": float(np.quantile(p, q)) for q in (.01, .05, .5, .95)}))
    audit = dict(schema_version="fsmpra.footprint_quality_benchmark.v1",
        sample_id=summary["sample_id"], partial=partial,
        iterations=iterations, permutations=permutations, alpha=b["alpha"],
        fixed_targets=[dict(endpoint_id=eid, metric=metric) for eid, metric in targets],
        complete_barcodes_by_partition=available,
        complete_barcodes_total=complete_total,
        sampled_barcodes=10_000,
        excluded_barcodes_missing_fixed_target=10_000-complete_total,
        missing_barcodes_by_target=missing_by_target,
        tested_methods=["welch", "studentized_permutation"],
        total_contrasts=len(records)//2, release_withheld=True,
        limitations=["Studentized permutation is not finite-sample exact for mean-null unequal variances.",
                     "Monte Carlo p-value resolution does not assess 13,926-hypothesis tails.",
                     "No variant p/q values or variant-effect confidence intervals are produced."],
        sources=sources)
    check_sources(summary["sources"])
    check_sources(sources)
    with atomic_output(output) as out:
        write_table(out / "iterations.tsv.gz", records, list(records[0]))
        write_table(out / "summary.tsv", summaries, list(summaries[0]))
        write_json(out / "audit.json", audit)
        check_sources(summary["sources"])
        check_sources(sources)
    return dict(contrasts=len(records)//2, partial=partial,
                complete_barcodes_by_partition=available,
                excluded_barcodes_missing_fixed_target=10_000-complete_total)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("config", "summary", "output"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--pilot-iterations", type=int)
    parser.add_argument("--pilot-permutations", type=int)
    args = parser.parse_args()
    print(json.dumps(benchmark(args.config, args.summary, args.output,
                               args.pilot_iterations, args.pilot_permutations)))

"""Frozen sampling checks for barcode-mean statistics; not sequence-bias validation."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

from .io import atomic_output, open_text, read_json, source_info, write_json, write_table
from .statistics import _validate_config, endpoint_matrix, fixed_endpoints, welch_contrast, welch_from_summary


def exact_interval(successes, trials):
    if not trials:
        return None, None
    ci = binomtest(int(successes), int(trials)).proportion_ci(.95, method="exact")
    return float(ci.low), float(ci.high)


def moments(x):
    """Additive summaries allow disjoint WT complements without repeated scans."""
    return (len(x), x.sum(axis=0), np.square(x).sum(axis=0),
            (x > 0).sum(axis=0), (x < 1).sum(axis=0))


def subtract_moments(total, part):
    return tuple(a - b for a, b in zip(total, part))


def welch_summary(raw):
    n, total, squares, nonzero, nonone = raw
    if n < 2:
        raise ValueError("Complement must contain at least two barcodes")
    numerator = squares - total * total / n
    # Constant fractional outcomes can leave a positive roundoff residual.
    tolerance = 32 * np.finfo(float).eps * np.maximum(squares, total * total / n)
    numerator = np.where(np.abs(numerator) <= tolerance, 0., numerator)
    variance = np.maximum(0., numerator / (n - 1))
    return n, total / n, variance, nonzero, nonone


def summarize_stratum(rows, cfg, check_coverage=True):
    attempts = len(rows)
    valid = [r for r in rows if r["status"] == "ok"]
    reject = sum(r["p_value"] < .05 for r in valid)
    covered = sum(r["ci_low"] <= r["truth"] <= r["ci_high"] for r in valid)
    lo, hi = exact_interval(reject, attempts)
    clo, chi = exact_interval(covered, len(valid))
    sufficient = len(valid) >= cfg["minimum_testable_iterations"]
    is_null = all(r["truth"] == 0 for r in rows)
    failures = []
    if sufficient and is_null and hi > cfg["null_rate_upper_limit"]:
        failures.append("null_rejection_upper_bound")
    if sufficient and check_coverage and clo < cfg["ci_coverage_lower_limit"]:
        failures.append("ci_coverage_lower_bound")
    p = np.array([r["p_value"] for r in valid])
    return dict(attempts=attempts, testable=len(valid), rejected=reject,
                unconditional_rejection_rate=reject / attempts,
                rejection_interval95=[lo, hi], ci_covered=covered,
                conditional_ci_coverage=covered / len(valid) if valid else None,
                conditional_ci_interval95=[clo, chi],
                p_quantiles={str(q): float(np.quantile(p, q)) for q in (.01, .05, .5, .95)} if valid else {},
                statuses=dict(Counter(r["status"] for r in rows)),
                status="unresolved" if not sufficient else "fail" if failures else "pass",
                failures=failures, null=is_null)


def draw_barcodes(rng, n, mean, rho, informative, arm):
    """Beta-distributed barcode propensities and binomial molecule observations.

    Depth depends on propensity in the informative scenario, but every barcode
    is observed. This tests weighting/heteroskedasticity, not selection bias.
    """
    if rho:
        precision = 1 / rho - 1
        latent = rng.beta(mean * precision, (1 - mean) * precision, n)
    else:
        latent = np.full(n, mean)
    if informative:
        depth = 1 + rng.poisson(1 + (8 if arm == "variant" else 4) * latent)
    else:
        depth = rng.geometric(.3 if arm == "variant" else .5, n)
    values = rng.binomial(depth, latent) / depth
    return values[:, None], int(depth.sum())


def _record(result, j, truth, label, iteration, enough_molecules=True):
    status = str(result["status"][j]) if enough_molecules else "insufficient_molecules"
    return dict(stratum=label, iteration=iteration, truth=truth, status=status,
                delta=float(result["delta"][j]),
                **{key: float(result[key][j]) if status == "ok" else None
                   for key in ("p_value", "ci_low", "ci_high")})


def simulate(cfg, rng, thresholds):
    records, strata = [], []
    scenarios = []
    for mean in (.02, .1, .5):
        for rho in (0., .2):
            for informative in (False, True):
                for n in (50, 400):
                    scenarios.append((mean, rho, informative, n, 0.))
    for delta in (-.05, .05, .1):
        scenarios.append((.1, .2, True, 400, delta))
    for index, (mean, rho, informative, n, delta) in enumerate(scenarios):
        label = f"simulation_{index:02d}"
        rows = []
        for iteration in range(cfg["iterations"]):
            x, nx = draw_barcodes(rng, n, mean + delta, rho, informative, "variant")
            y, ny = draw_barcodes(rng, 2000, mean, rho, informative, "wt")
            result = welch_contrast(x, y, **thresholds)
            rows.append(_record(result, 0, delta, label, iteration,
                                min(nx, ny) >= cfg.get("min_molecules", 100)))
        strata.append(dict(stratum=label, kind="simulation", baseline=mean,
                           barcode_correlation=rho, informative_depth=informative,
                           n_variant=n, n_wt=2000, true_delta=delta,
                           **summarize_stratum(rows, cfg)))
        records.extend(rows)
    return records, strata


def release_gate(strata, endpoints):
    """Do not silently discard failed cells or call sparse cells validated."""
    assessed = [s for s in strata if s["status"] != "unresolved"]
    failures = [s["stratum"] for s in assessed if s["status"] == "fail"]
    validated = []
    for endpoint in endpoints:
        cells = [s for s in assessed if s.get("endpoint_id") == endpoint["endpoint_id"]]
        if cells and all(s["status"] == "pass" for s in cells):
            validated.append(endpoint["endpoint_id"])
    unresolved_null = [s["stratum"] for s in strata if s["kind"] == "simulation"
                       and s["null"] and s["status"] == "unresolved"]
    enough = (bool(validated) and not unresolved_null and
              any(s["kind"] == "simulation" and s["null"] for s in assessed))
    return dict(status="fail" if failures else "pass" if enough else "unresolved",
                validated_endpoint_ids=validated, failed_strata=failures,
                unresolved_null_simulations=unresolved_null,
                assessed_strata=len(assessed), unresolved_strata=sum(s["status"] == "unresolved" for s in strata))


def calibrate(config_path, aggregated, output):
    spec = read_json(config_path)
    _validate_config(spec)
    cfg = {**spec["calibration"], "min_molecules": spec["min_molecules"]}
    if cfg["iterations"] < cfg["minimum_testable_iterations"]:
        raise ValueError("Iterations cannot meet the frozen minimum")
    thresholds = dict(min_barcodes=spec["min_barcodes"], min_nonzero=spec["min_nonzero_barcodes"],
                      min_nonone=spec["min_nonone_barcodes"])
    aggregated = Path(aggregated)
    audit = read_json(aggregated / "audit.json")
    if audit["sample_id"] != spec["sample_id"] or audit["interval"] != spec["interval"]:
        raise ValueError("Aggregate configuration mismatch")
    config_source = source_info(config_path, True)
    if not any(s.get("path") == config_source["path"] and
               s.get("sha256") == config_source["sha256"] for s in audit["sources"]):
        raise ValueError("Aggregate config SHA does not match calibration")
    with open_text(aggregated / "variant_metadata.tsv") as handle:
        wt = [r["variant_id"] for r in csv.DictReader(handle, delimiter="\t") if r["genotype_class"] == "WT"]
    if len(wt) != 1:
        raise ValueError("Expected exactly one WT genotype")
    counts = np.load(aggregated / "barcode_counts.npy", mmap_mode="r")
    molecules = np.load(aggregated / "barcode_molecules.npy", mmap_mode="r")
    idx = []
    n_metadata = 0
    with open_text(aggregated / "barcode_metadata.tsv.gz") as handle:
        for i, row in enumerate(csv.DictReader(handle, delimiter="\t")):
            if int(row["barcode_index"]) != i or int(row["molecules"]) != molecules[i]:
                raise ValueError("Metadata and arrays disagree")
            if row["variant_id"] == wt[0]:
                idx.append(i)
            n_metadata += 1
    if n_metadata != len(molecules):
        raise ValueError("Incomplete barcode metadata")
    idx = np.asarray(idx, dtype=np.int64)
    if len(idx) <= max(cfg["wt_pseudo_variant_barcodes"]) + spec["min_barcodes"]:
        raise ValueError("Insufficient WT barcodes for disjoint calibration")
    endpoints = fixed_endpoints(*spec["interval"])
    # Chunk conversion bounds temporary memory independently of WT barcode count.
    matrix = np.empty((len(idx), len(endpoints)))
    for a in range(0, len(idx), 10000):
        selected = idx[a:a+10000]
        matrix[a:a+10000] = endpoint_matrix(counts[selected], molecules[selected], endpoints)
    wt_molecules = np.asarray(molecules[idx], dtype=np.uint64)
    total_molecules = int(wt_molecules.sum())
    total = moments(matrix)
    # Independent RNG streams: adding simulation cells cannot change WT splits.
    aa_seed, simulation_seed = np.random.SeedSequence(spec["seed"]).spawn(2)
    rng = np.random.default_rng(aa_seed)
    records, strata = [], []
    for n in cfg["wt_pseudo_variant_barcodes"]:
        rows = [[] for _ in endpoints]
        for iteration in range(cfg["iterations"]):
            selected = rng.choice(len(idx), n, replace=False)
            part = moments(matrix[selected])
            result = welch_from_summary(*welch_summary(part),
                *welch_summary(subtract_moments(total, part)), **thresholds)
            n_molecules = int(wt_molecules[selected].sum())
            enough = min(n_molecules, total_molecules-n_molecules) >= spec["min_molecules"]
            for j, ep in enumerate(endpoints):
                label = f"wt_aa_{n}_{ep['endpoint_id']}"
                rows[j].append(_record(result, j, 0., label, iteration, enough))
        for ep, cell in zip(endpoints, rows):
            strata.append(dict(stratum=cell[0]["stratum"], kind="wt_aa", endpoint_id=ep["endpoint_id"],
                               n_variant=n, n_wt=len(idx)-n,
                               **summarize_stratum(cell, cfg, check_coverage=False)))
            records.extend(cell)
    sim_records, sim_strata = simulate(cfg, np.random.default_rng(simulation_seed), thresholds)
    records.extend(sim_records)
    strata.extend(sim_strata)
    gate = release_gate(strata, endpoints)
    provenance = [source_info(p, Path(p).suffix != ".npy") for p in
        (config_path, aggregated / "audit.json", aggregated / "barcode_counts.npy",
         aggregated / "barcode_molecules.npy", aggregated / "barcode_metadata.tsv.gz",
         __file__, Path(__file__).with_name("statistics.py"))]
    report = dict(schema_version="fsmpra.statistics_calibration.v1", **gate,
                  sample_id=spec["sample_id"], seed=spec["seed"], iterations=cfg["iterations"],
                  wt_barcodes=len(idx), wt_molecules=total_molecules, strata=strata, sources=provenance,
                  biological_detectability_resolved=False,
                  interpretation="Sampling development checks only. Unconditional rejection uses all iterations; CI coverage uses testable iterations. Unresolved sparse cells do not validate inference. Nominal power is not FDR-adjusted power. WT splits reuse a finite empirical population; binomial bounds quantify conditional Monte Carlo precision, not biological replication. Simulation intervals reflect the specified model, not an empirical guarantee.")
    with atomic_output(output) as out:
        write_table(out / "iterations.tsv.gz", records, list(records[0]))
        write_json(out / "report.json", report)
    return gate


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("config", "aggregated", "output"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    print(json.dumps(calibrate(args.config, args.aggregated, args.output)))

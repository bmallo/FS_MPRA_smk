"""Held-out WT A/A and synthetic calibration of fixed barcode-fraction endpoints.

This diagnostic does not infer variant effects or authorize a release.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, t as student_t

from .io import atomic_output, read_json, source_info, write_json, write_table
from .processing import check_sources
from .statistics import fixed_endpoints

SCHEMA = "fsmpra.heldout_calibration_config.v1"
EXPECTED = dict(seed=2026101802, n_variant=[50, 100, 400, 1000],
                wt_iterations=1000, wt_permutations=9999,
                simulation_iterations=200, simulation_permutations=999,
                simulation_control_barcodes=10000, alpha=[.05, .01], confidence_level=.95)


def validate_config(config, pilot_iterations=None, pilot_permutations=None):
    if config.get("schema_version") != SCHEMA:
        raise ValueError("Unexpected held-out calibration configuration")
    settings = config.get("calibration", config)
    for key, expected in EXPECTED.items():
        if settings.get(key) != expected:
            raise ValueError(f"Frozen held-out calibration setting changed: {key}")
    iterations = settings["wt_iterations"] if pilot_iterations is None else pilot_iterations
    simulation_iterations = settings["simulation_iterations"] if pilot_iterations is None else pilot_iterations
    wt_permutations = settings["wt_permutations"] if pilot_permutations is None else pilot_permutations
    simulation_permutations = settings["simulation_permutations"] if pilot_permutations is None else pilot_permutations
    if not 1 <= iterations <= settings["wt_iterations"] or not 1 <= simulation_iterations <= settings["simulation_iterations"]:
        raise ValueError("Pilot iterations exceed the frozen setting")
    if not 1 <= wt_permutations <= settings["wt_permutations"] or not 1 <= simulation_permutations <= settings["simulation_permutations"]:
        raise ValueError("Pilot permutations exceed the frozen setting")
    partial = iterations != settings["wt_iterations"] or wt_permutations != settings["wt_permutations"]
    return settings, iterations, simulation_iterations, wt_permutations, simulation_permutations, partial


def exact_interval(successes, trials, level=.95):
    interval = binomtest(int(successes), int(trials)).proportion_ci(level, method="exact")
    return float(interval.low), float(interval.high)


def target_columns(endpoint_ids):
    expected = [ep["endpoint_id"] for ep in fixed_endpoints()]
    if list(endpoint_ids) != expected:
        raise ValueError("Endpoint IDs differ from frozen grid")
    return ([(ep, "presence", j) for j, ep in enumerate(expected[:-1])] +
            [(ep, "coverage", j) for j, ep in enumerate(expected)])


def load_outcomes(outcomes_dir):
    outcomes_dir = Path(outcomes_dir)
    audit_path = outcomes_dir / "audit.json"
    matrix_path = outcomes_dir / "endpoint_matrix.npz"
    audit_source = source_info(audit_path, True)
    audit = read_json(audit_path)
    if audit.get("schema_version") != "fsmpra.heldout_outcomes.v1" or audit.get("complete") is not True or audit.get("pilot") is True:
        raise ValueError("Held-out outcomes must be complete")
    sources = audit.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ValueError("Held-out outcomes lack source provenance")
    paths = [item.get("path") for item in sources if isinstance(item, dict)]
    if len(paths) != len(sources) or len(paths) != len(set(paths)):
        raise ValueError("Duplicate or malformed outcome source provenance")
    check_sources(sources)
    matrix_source = source_info(matrix_path, True)
    if matrix_source not in sources:
        raise ValueError("Held-out matrix is absent from outcome source provenance")
    with np.load(matrix_path, allow_pickle=False) as data:
        needed = {"barcodes", "endpoint_ids", "presence", "coverage", "depths", "valid_molecules"}
        if not needed <= set(data.files):
            raise ValueError("Incomplete held-out endpoint matrix")
        arrays = {key: data[key] for key in needed}
    targets = target_columns(arrays["endpoint_ids"])
    with np.load(matrix_path, allow_pickle=False) as data:
        for key, field in (("relative_starts", "relative_start"), ("relative_ends", "relative_end"),
                           ("reference_starts0", "reference_start0"), ("reference_ends0", "reference_end0")):
            if key not in data.files or list(data[key]) != [ep[field] for ep in fixed_endpoints()]:
                raise ValueError("Held-out endpoint coordinates differ from frozen grid")
    barcode = arrays["barcodes"]
    n = len(barcode)
    if barcode.ndim != 1 or n != 258283 or not np.all(barcode[:-1] < barcode[1:]):
        raise ValueError("Expected 258283 unique sorted WT barcodes outside the prior audit sample")
    if audit.get("heldout_barcodes") != n or audit.get("endpoint_count") != 22:
        raise ValueError("Held-out audit barcode count differs from matrix")
    if audit.get("input_cohort") != {"wt_barcodes": 268283, "wt_molecules": 985264,
                                     "excluded_barcodes": 10000, "excluded_molecules": 37075}:
        raise ValueError("Held-out WT input cohort or prior audit exclusion changed")
    depths, valid = arrays["depths"], arrays["valid_molecules"]
    if depths.shape != (n,) or valid.shape != (n, 22) or np.any(depths <= 0) or np.any(valid < 0) or np.any(valid > depths[:, None]):
        raise ValueError("Invalid held-out molecule depths")
    if audit.get("heldout_molecules") != 948189 or int(depths.sum()) != 948189:
        raise ValueError("Held-out molecule depth total differs from audit")
    for metric in ("presence", "coverage"):
        values = arrays[metric]
        if values.shape != (n, 22) or values.dtype != np.float64:
            raise ValueError("Held-out fraction shape or dtype mismatch")
        if np.any(np.isfinite(values) != (valid > 0)) or np.any((values[np.isfinite(values)] < 0) | (values[np.isfinite(values)] > 1)):
            raise ValueError("Invalid held-out fractions or missingness")
    matrix = np.column_stack([arrays[metric][:, j] for _, metric, j in targets])
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Fixed held-out targets require all WT barcodes to be observed")
    if source_info(audit_path, True) != audit_source:
        raise ValueError("Held-out audit changed while loading")
    return audit, targets, matrix, [audit_source, matrix_source]


def _stats_from_centered(total, total_sq, selected, n_pool, center):
    """All columns; selected rows versus their entire complement."""
    selected = np.asarray(selected, dtype=np.float64)
    n = len(selected)
    m = n_pool - n
    if n < 2 or m < 2:
        raise ValueError("Each arm needs at least two barcodes")
    sx = selected.sum(axis=0)
    qx = np.square(selected).sum(axis=0)
    sy, qy = total - sx, total_sq - qx
    mx = sx / n + center
    my = sy / m + center
    delta = sx / n - sy / m
    raw_x = qx - sx * sx / n
    raw_y = qy - sy * sy / m
    scale_x = np.maximum(qx, sx * sx / n)
    scale_y = np.maximum(qy, sy * sy / m)
    tol_x = 64 * np.finfo(float).eps * scale_x
    tol_y = 64 * np.finfo(float).eps * scale_y
    vx = np.maximum(0, np.where(np.abs(raw_x) <= tol_x, 0, raw_x) / (n - 1))
    vy = np.maximum(0, np.where(np.abs(raw_y) <= tol_y, 0, raw_y) / (m - 1))
    ax, ay = vx / n, vy / m
    se = np.sqrt(ax + ay)
    statistic = np.divide(delta, se, out=np.zeros_like(delta), where=se > 0)
    statistic[(se == 0) & (delta != 0)] = np.copysign(np.inf, delta[(se == 0) & (delta != 0)])
    df = np.divide(np.square(ax + ay), np.square(ax) / (n - 1) + np.square(ay) / (m - 1),
                   out=np.full_like(se, np.nan), where=se > 0)
    p = np.where(se > 0, 2 * student_t.sf(np.abs(statistic), df), np.where(delta == 0, 1., 0.))
    critical = student_t.ppf(.975, df)
    low = np.clip(np.where(se > 0, delta - critical * se, delta), -1., 1.)
    high = np.clip(np.where(se > 0, delta + critical * se, delta), -1., 1.)
    return dict(mean_variant=mx, mean_control=my, delta=delta, se=se,
                statistic=statistic, df=df, welch_p=p, welch_ci_low=low, welch_ci_high=high)


def centered_pool(matrix):
    matrix = np.asarray(matrix, dtype=np.float64)
    if matrix.ndim != 2 or not np.all(np.isfinite(matrix)):
        raise ValueError("Expected a finite 2D barcode matrix")
    constant = np.all(matrix == matrix[0], axis=0)
    center = np.where(constant, matrix[0], matrix.mean(axis=0))
    values = matrix - center
    values[:, constant] = 0.
    return values, values.sum(axis=0), np.square(values).sum(axis=0), center, np.count_nonzero(matrix > 0, axis=0)


def contrast_from_pool(pool, selected_indices, confidence_level=.95):
    values, total, total_sq, center, total_nonzero = pool
    result = _stats_from_centered(total, total_sq, values[selected_indices], len(values), center)
    n, m = len(selected_indices), len(values) - len(selected_indices)
    radius = math.sqrt(math.log(4 / (1-confidence_level)) / (2*n)) + math.sqrt(math.log(4 / (1-confidence_level)) / (2*m))
    result["hoeffding_ci_low"] = np.maximum(-1., result["delta"] - radius)
    result["hoeffding_ci_high"] = np.minimum(1., result["delta"] + radius)
    result["nonzero_variant"] = np.count_nonzero(values[selected_indices] + center > 0, axis=0)
    result["nonzero_control"] = total_nonzero - result["nonzero_variant"]
    return result


def _reference_statistic(total, total_sq, selected, n_pool):
    n = len(selected)
    m = n_pool - n
    sx = selected.sum(axis=0)
    qx = np.square(selected).sum(axis=0)
    sy, qy = total-sx, total_sq-qx
    rx = qx-sx*sx/n
    ry = qy-sy*sy/m
    rx = np.where(np.abs(rx) <= 64*np.finfo(float).eps*np.maximum(qx, sx*sx/n), 0., rx)
    ry = np.where(np.abs(ry) <= 64*np.finfo(float).eps*np.maximum(qy, sy*sy/m), 0., ry)
    se = np.sqrt(np.maximum(rx, 0.) / (n*(n-1)) + np.maximum(ry, 0.) / (m*(m-1)))
    delta = sx/n-sy/m
    statistic = np.divide(delta, se, out=np.zeros_like(delta), where=se > 0)
    statistic[(se == 0) & (delta != 0)] = np.inf
    return np.abs(statistic)


def reference_bank(pool, n, permutations, rng):
    """One uniform-label Welch reference bank, shared across fixed endpoints and WT draws."""
    values, total, total_sq, _, _ = pool
    bank = np.empty((permutations, values.shape[1]), dtype=float)
    for i in range(permutations):
        selected = rng.choice(len(values), size=n, replace=False, shuffle=False)
        bank[i] = _reference_statistic(total, total_sq, values[selected], len(values))
    bank.sort(axis=0)
    return bank


def bank_p_values(statistic, sorted_bank):
    statistic = np.atleast_1d(np.abs(statistic))
    if sorted_bank.shape[1] != len(statistic):
        raise ValueError("Bank and statistic columns differ")
    b = len(sorted_bank)
    exceeds = np.array([b - np.searchsorted(sorted_bank[:, j], max(0., statistic[j] - 1e-12), side="left")
                        for j in range(len(statistic))], dtype=int)
    return (1 + exceeds) / (b + 1), exceeds


def draw_barcodes(rng, n, mean, rho, depth_mode, arm):
    if not 0 <= mean <= 1 or not 0 < rho < 1:
        raise ValueError("Invalid simulation mean or correlation")
    if mean in (0, 1):
        latent = np.full(n, mean)
    else:
        precision = 1/rho - 1
        latent = rng.beta(mean*precision, (1-mean)*precision, n)
    if depth_mode == "independent":
        depth = rng.geometric(.3 if arm == "variant" else .5, n)
    elif depth_mode == "associated":
        depth = 1 + rng.poisson(1 + 20*latent)
    else:
        raise ValueError("Invalid depth mode")
    return rng.binomial(depth, latent) / depth


def _scalar_result(result):
    return {key: value[0] for key, value in result.items()}


def _value(value):
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _rows(source, scenario, iteration, endpoint_id, metric, truth, result, permutation_p, exceedances, permutations, n, m, confidence_level):
    radius = math.sqrt(math.log(4/(1-confidence_level))/(2*n)) + math.sqrt(math.log(4/(1-confidence_level))/(2*m))
    delta = float(result["delta"])
    hlow, hhigh = max(-1., delta-radius), min(1., delta+radius)
    status = "uninformative" if result["se"] == 0 and delta == 0 else "zero_se_difference" if result["se"] == 0 else "ok"
    base = dict(data_source=source, scenario=scenario, iteration=iteration,
                endpoint_id=endpoint_id, metric=metric, true_delta=truth,
                n_variant=n, n_control=m,
                nonzero_variant=int(result["nonzero_variant"]), nonzero_control=int(result["nonzero_control"]),
                mean_variant=float(result["mean_variant"]), mean_control=float(result["mean_control"]),
                delta=delta, se=float(result["se"]), statistic=float(result["statistic"]),
                df=_value(result["df"]), status=status,
                welch_ci_low=float(result["welch_ci_low"]), welch_ci_high=float(result["welch_ci_high"]),
                hoeffding_ci_low=hlow, hoeffding_ci_high=hhigh)
    return [dict(base, method="welch", p_value=float(result["welch_p"]), exceedances=None, permutations=None),
            dict(base, method="studentized_permutation", p_value=float(permutation_p),
                 exceedances=int(exceedances), permutations=permutations)]


def simulation_scenarios():
    scenarios = []
    for mean in (.002, .02, .1, .5):
        for n in (100, 400):
            for rho_v, rho_c in ((.02, .2), (.2, .02)):
                for mode in ("independent", "associated"):
                    scenarios.append(dict(kind="null", n_variant=n, control_mean=mean, variant_mean=mean,
                                          rho_variant=rho_v, rho_control=rho_c, depth_mode=mode))
    for n in (100, 400, 1000):
        for control_mean, variant_mean in ((.002, 0), (.02, 0), (.02, .015), (.02, .025), (.1, .08), (.1, .12)):
            scenarios.append(dict(kind="alternative", n_variant=n, control_mean=control_mean,
                                  variant_mean=variant_mean, rho_variant=.02, rho_control=.2,
                                  depth_mode="independent"))
    return scenarios


def _summary(rows, alpha, level):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["data_source"], row["scenario"], row["endpoint_id"], row["metric"], row["method"])].append(row)
    summaries = []
    for key, cell in groups.items():
        truth = cell[0]["true_delta"]
        base = dict(zip(("data_source", "scenario", "endpoint_id", "metric", "method"), key))
        base.update(iterations=len(cell), true_delta=truth,
                    nonzero_variant_median=float(np.median([r["nonzero_variant"] for r in cell])),
                    nonzero_control_median=float(np.median([r["nonzero_control"] for r in cell])),
                    statuses=dict(Counter(r["status"] for r in cell)),
                    positive_delta=sum(r["delta"] > 0 for r in cell),
                    negative_delta=sum(r["delta"] < 0 for r in cell),
                    zero_delta=sum(r["delta"] == 0 for r in cell))
        for threshold in alpha:
            code = f"a{int(threshold*100):02d}"
            rejected = sum(r["p_value"] <= threshold for r in cell)
            positive_rejected = sum(r["p_value"] <= threshold and r["delta"] > 0 for r in cell)
            negative_rejected = sum(r["p_value"] <= threshold and r["delta"] < 0 for r in cell)
            zero_rejected = rejected-positive_rejected-negative_rejected
            low, high = exact_interval(rejected, len(cell), level)
            directional = sum(r["p_value"] <= threshold and
                              np.sign(r["delta"]) == np.sign(truth) for r in cell) if truth else None
            base.update({f"{code}_rejected": rejected, f"{code}_rate": rejected/len(cell),
                         f"{code}_rate_ci_low": low, f"{code}_rate_ci_high": high,
                         f"{code}_positive_rejected": positive_rejected,
                         f"{code}_negative_rejected": negative_rejected,
                         f"{code}_zero_rejected": zero_rejected,
                         f"{code}_directional_rejected": directional})
        for name in ("welch", "hoeffding"):
            widths = [r[f"{name}_ci_high"] - r[f"{name}_ci_low"] for r in cell]
            covered = sum(r[f"{name}_ci_low"] <= truth <= r[f"{name}_ci_high"] for r in cell)
            below = sum(r[f"{name}_ci_high"] < truth for r in cell)
            above = sum(r[f"{name}_ci_low"] > truth for r in cell)
            low, high = exact_interval(covered, len(cell), level)
            base.update({f"{name}_covered": covered, f"{name}_coverage": covered/len(cell),
                         f"{name}_coverage_ci_low": low, f"{name}_coverage_ci_high": high,
                         f"{name}_noncoverage_below_truth": below,
                         f"{name}_noncoverage_above_truth": above,
                         f"{name}_mean_width": float(np.mean(widths))})
        summaries.append(base)
    return summaries


def calibrate(config_path, outcomes_dir, output, pilot_iterations=None, pilot_permutations=None):
    config_path, outcomes_dir = Path(config_path), Path(outcomes_dir)
    config_source = source_info(config_path, True)
    cfg = read_json(config_path)
    settings, wt_iterations, sim_iterations, wt_perms, sim_perms, partial = validate_config(cfg, pilot_iterations, pilot_permutations)
    audit, targets, matrix, outcome_sources = load_outcomes(outcomes_dir)
    if config_source not in audit["sources"]:
        raise ValueError("Calibration config source is not recorded in held-out outcomes audit")
    sources = [config_source, *outcome_sources, source_info(__file__, True),
               source_info(Path(__file__).with_name("statistics.py"), True),
               source_info(Path(__file__).with_name("io.py"), True),
               source_info(Path(__file__).with_name("processing.py"), True)]
    pool = centered_pool(matrix)
    seed = settings["seed"]
    records = []
    for n in settings["n_variant"]:
        bank_rng = np.random.default_rng(np.random.SeedSequence([seed, 0, n]))
        draw_rng = np.random.default_rng(np.random.SeedSequence([seed, 1, n]))
        bank = reference_bank(pool, n, wt_perms, bank_rng)
        for iteration in range(wt_iterations):
            chosen = draw_rng.choice(len(matrix), n, replace=False, shuffle=False)
            result = contrast_from_pool(pool, chosen, settings["confidence_level"])
            permutation_p, exceedances = bank_p_values(result["statistic"], bank)
            for j, (endpoint, metric, _) in enumerate(targets):
                scalar = {key: value[j] for key, value in result.items()}
                records.extend(_rows("wt_aa", f"wt_{n}", iteration, endpoint, metric, 0.,
                                     scalar, permutation_p[j], exceedances[j], wt_perms,
                                     n, len(matrix)-n, settings["confidence_level"]))
        print(f"WT A/A n={n}: {wt_iterations} draws, {wt_perms} reference permutations", flush=True)
    for scenario_id, scenario in enumerate(simulation_scenarios()):
        label = f"simulation_{scenario_id:02d}"
        rng = np.random.default_rng(np.random.SeedSequence([seed, 2, scenario_id]))
        n = scenario["n_variant"]
        m = settings["simulation_control_barcodes"]
        truth = scenario["variant_mean"] - scenario["control_mean"]
        for iteration in range(sim_iterations):
            x = draw_barcodes(rng, n, scenario["variant_mean"], scenario["rho_variant"], scenario["depth_mode"], "variant")
            y = draw_barcodes(rng, m, scenario["control_mean"], scenario["rho_control"], scenario["depth_mode"], "control")
            sample = np.concatenate((x, y))[:, None]
            sim_pool = centered_pool(sample)
            result = _scalar_result(contrast_from_pool(sim_pool, np.arange(n), settings["confidence_level"]))
            bank = reference_bank(sim_pool, n, sim_perms, rng)
            permutation_p, exceedances = bank_p_values([result["statistic"]], bank)
            records.extend(_rows("simulation", label, iteration, "synthetic_fraction", "coverage", truth,
                                 result, permutation_p[0], exceedances[0], sim_perms, n, m,
                                 settings["confidence_level"]))
        print(f"{label}: {sim_iterations} draws", flush=True)
    summaries = _summary(records, settings["alpha"], settings["confidence_level"])
    check_sources(audit["sources"])
    check_sources(sources)
    report = dict(schema_version="fsmpra.heldout_calibration.v1", partial=partial,
                  release_withheld=True, seed=seed, sample_id=cfg["sample_id"],
                  wt_barcodes=len(matrix), input_cohort=audit["input_cohort"],
                  cohort_relation="Disjoint from the earlier 10,000-barcode audit sample only; not globally unused or biologically independent.",
                  wt_contrasts=len(settings["n_variant"])*wt_iterations*len(targets),
                  simulation_contrasts=len(simulation_scenarios())*sim_iterations,
                  total_contrasts=len(records)//2,
                  fixed_targets=[dict(endpoint_id=e, metric=m) for e,m,_ in targets],
                  wt_iterations=wt_iterations, wt_permutations=wt_perms,
                  simulation_iterations=sim_iterations, simulation_permutations=sim_perms,
                  simulation_control_barcodes=settings["simulation_control_barcodes"],
                  simulation_scenarios=simulation_scenarios(),
                  reference_bank="Independent uniform subsets for each n, shared across WT iterations and endpoints; conditional interval omits reference-bank Monte Carlo uncertainty.",
                  inference_limitations=["Studentized permutation is not finite-sample exact for an unequal-variance mean null.",
                                         "Binomial intervals condition on the fixed WT barcode pool and fixed permutation reference bank.",
                                         "Simulation intervals assess only the stated models; neither candidate interval is adopted automatically.",
                                         "No variant p/q values or variant-effect confidence intervals are produced."],
                  sources=sources)
    with atomic_output(output) as out:
        write_table(out / "iterations.tsv.gz", records, list(records[0]))
        write_table(out / "summary.tsv", summaries, list(summaries[0]))
        write_json(out / "audit.json", report)
        check_sources(audit["sources"])
        check_sources(sources)
    return dict(partial=partial, wt_barcodes=len(matrix), contrasts=len(records)//2,
                summary_cells=len(summaries), release_withheld=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--outcomes", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--pilot-iterations", type=int)
    parser.add_argument("--pilot-permutations", type=int)
    args = parser.parse_args()
    print(json.dumps(calibrate(args.config, args.outcomes, args.output,
                               args.pilot_iterations, args.pilot_permutations)))

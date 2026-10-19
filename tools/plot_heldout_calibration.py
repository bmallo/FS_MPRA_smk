#!/usr/bin/env python3
"""Figures and compact tables for the completed frozen held-out calibration."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import tempfile
from pathlib import Path

METHODS = ("welch", "studentized_permutation")
METHOD_LABELS = {"welch": "Welch", "studentized_permutation": "Studentized permutation"}
REQUIRED = {"data_source", "scenario", "endpoint_id", "metric", "method", "iterations",
            "true_delta", "a05_rejected", "a05_rate", "a05_rate_ci_low", "a05_rate_ci_high",
            "a01_rejected", "a01_rate", "a01_rate_ci_low", "a01_rate_ci_high",
            "welch_coverage", "welch_coverage_ci_low", "welch_coverage_ci_high",
            "hoeffding_coverage", "hoeffding_coverage_ci_low", "hoeffding_coverage_ci_high"}


def _source(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    stat = path.stat()
    return {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "sha256": digest.hexdigest()}


def _check_sources(records):
    if not isinstance(records, list) or not records:
        raise ValueError("Calibration has no source provenance")
    for item in records:
        current = _source(item["path"]) if "sha256" in item else None
        if current != item:
            raise ValueError("Calibration source changed: " + item["path"])


def _read_rows(path, audit):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not REQUIRED <= set(reader.fieldnames or []):
            raise ValueError("Held-out summary lacks required fields")
        raw_rows = list(reader)
    targets = {(item["endpoint_id"], item["metric"]) for item in audit["fixed_targets"]}
    if len(targets) != 43:
        raise ValueError("Expected 43 fixed TF-presence, TF/MSP-coverage targets")
    scenarios = audit["simulation_scenarios"]
    if (len(scenarios) != 50 or sum(s["kind"] == "null" for s in scenarios) != 32 or
            sum(s["kind"] == "alternative" for s in scenarios) != 18):
        raise ValueError("Expected frozen 32 null and 18 alternative simulations")
    rows = {}
    for raw in raw_rows:
        key = tuple(raw[k] for k in ("data_source", "scenario", "endpoint_id", "metric", "method"))
        if key in rows:
            raise ValueError("Duplicate held-out summary row: " + str(key))
        source, scenario, endpoint, metric, method = key
        if method not in METHODS:
            raise ValueError("Unexpected inference method")
        if source == "wt_aa":
            if scenario not in {f"wt_{n}" for n in (50, 100, 400, 1000)} or (endpoint, metric) not in targets:
                raise ValueError("Unexpected WT A/A target")
            expected_iterations = int(audit["wt_iterations"])
        elif source == "simulation":
            if (endpoint, metric) != ("synthetic_fraction", "coverage"):
                raise ValueError("Unexpected synthetic endpoint")
            index = int(scenario.removeprefix("simulation_")) if scenario.startswith("simulation_") else -1
            if not 0 <= index < 50 or scenario != f"simulation_{index:02d}":
                raise ValueError("Unexpected simulation scenario")
            expected_iterations = int(audit["simulation_iterations"])
            truth = scenarios[index]["variant_mean"] - scenarios[index]["control_mean"]
            if not math.isclose(float(raw["true_delta"]), truth, abs_tol=1e-12):
                raise ValueError("Simulation truth differs from frozen scenario")
        else:
            raise ValueError("Unexpected data source")
        if int(raw["iterations"]) != expected_iterations:
            raise ValueError("Incomplete calibration summary cell")
        for prefix in ("a05", "a01"):
            rate = float(raw[f"{prefix}_rate"])
            rejected = int(raw[f"{prefix}_rejected"])
            low, high = float(raw[f"{prefix}_rate_ci_low"]), float(raw[f"{prefix}_rate_ci_high"])
            if not (0 <= rejected <= expected_iterations and
                    math.isclose(rate, rejected / expected_iterations, abs_tol=1e-12) and
                    0 <= low <= rate <= high <= 1):
                raise ValueError("Invalid rejection rate or exact interval")
        rows[key] = raw
    expected_keys = {
        ("wt_aa", f"wt_{n}", endpoint, metric, method)
        for n in (50, 100, 400, 1000) for endpoint, metric in targets for method in METHODS
    } | {
        ("simulation", f"simulation_{i:02d}", "synthetic_fraction", "coverage", method)
        for i in range(50) for method in METHODS
    }
    if set(rows) != expected_keys:
        raise ValueError("Held-out summary has missing or unexpected cells")
    return rows


def _plot_wt(rows, targets, metric, output, iterations, partial):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    endpoint_ids = [item["endpoint_id"] for item in targets if item["metric"] == metric]
    if metric == "presence":
        endpoint_ids = [item for item in endpoint_ids if item != "msp_promoter"]
    if len(endpoint_ids) != (21 if metric == "presence" else 22):
        raise ValueError("Unexpected WT target order")
    ns = (50, 100, 400, 1000)
    matrices = []
    for method in METHODS:
        matrices.append(np.array([[float(rows[("wt_aa", f"wt_{n}", endpoint, metric, method)]["a05_rate"])
                                   for n in ns] for endpoint in endpoint_ids]))
    vmax = max(.05, *(float(matrix.max()) for matrix in matrices))
    fig, axes = plt.subplots(1, 2, figsize=(11, 11.5), sharey=True)
    fig.subplots_adjust(left=.16, right=.9, top=.9, bottom=.1, wspace=.18)
    for ax, matrix, method in zip(axes, matrices, METHODS):
        image = ax.imshow(matrix, cmap="YlOrRd", vmin=0, vmax=vmax, aspect="auto")
        ax.set_xticks(range(4), labels=[str(n) for n in ns])
        ax.set_yticks(range(len(endpoint_ids)), labels=endpoint_ids, fontsize=8)
        ax.set_title(METHOD_LABELS[method])
        ax.set_xlabel("Pseudo-variant barcode count")
        for i in range(len(endpoint_ids)):
            for j in range(4):
                value = matrix[i, j]
                ax.text(j, i, f"{100*value:.1f}", ha="center", va="center", fontsize=6.5,
                        color="white" if value > vmax*.58 else "black")
    axes[0].set_ylabel("Fixed endpoint")
    cbar = fig.colorbar(image, ax=axes, shrink=.48, pad=.025)
    cbar.set_label("WT A/A rejection rate at α = 0.05")
    cbar.ax.axhline(.05, color="black", linestyle="--", linewidth=1.3)
    fig.suptitle(f"{'PILOT — ' if partial else ''}WT A/A {metric}: {iterations:,} contrasts per endpoint and size", fontsize=13)
    fig.text(.5, .04, "Cells show percentage rejected. Dashed colorbar marker is the nominal 5% target.\n"
             "These repeated draws condition on one WT pool and one permutation reference bank.",
             ha="center", fontsize=9)
    stem = "wt_coverage_rates" if metric == "coverage" else "wt_presence_rates"
    fig.savefig(output / f"{stem}.png", dpi=180)
    fig.savefig(output / f"{stem}.svg")
    plt.close(fig)


def _plot_simulation_null(rows, scenarios, output, iterations, partial):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    nulls = [(i, s) for i, s in enumerate(scenarios) if s["kind"] == "null"]
    fig, ax = plt.subplots(figsize=(11, 13.5))
    fig.subplots_adjust(left=.26, right=.95, top=.93, bottom=.075)
    colors = {"welch": "#2b6cb0", "studentized_permutation": "#c05621"}
    labels = []
    maximum = .05
    for position, (index, scenario) in enumerate(nulls):
        label = (f"μ={scenario['control_mean']:g}, n={scenario['n_variant']}, "
                 f"ρV/C={scenario['rho_variant']:g}/{scenario['rho_control']:g}, "
                 f"{scenario['depth_mode'][:5]}")
        labels.append(label)
        for method, shift in (("welch", -.13), ("studentized_permutation", .13)):
            row = rows[("simulation", f"simulation_{index:02d}", "synthetic_fraction", "coverage", method)]
            rate = float(row["a05_rate"])
            low, high = float(row["a05_rate_ci_low"]), float(row["a05_rate_ci_high"])
            maximum = max(maximum, high)
            ax.errorbar(rate, position + shift, xerr=[[rate-low], [high-rate]], fmt="o",
                        color=colors[method], markersize=3.5, capsize=2, linewidth=.9)
    ax.axvline(.05, color="0.35", linestyle="--", linewidth=1)
    ax.set_xlim(0, min(1, maximum+.025))
    ax.set_yticks(range(len(nulls)), labels=labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Rejection rate at α = 0.05; exact binomial 95% interval")
    ax.grid(axis="x", color="0.9")
    ax.set_axisbelow(True)
    ax.set_title(f"{'PILOT — ' if partial else ''}Synthetic equal-mean stress tests: {iterations:,} contrasts per scenario")
    ax.legend(handles=[Line2D([0], [0], color=colors[m], marker="o", label=METHOD_LABELS[m])
                       for m in METHODS] + [Line2D([0], [0], color="0.35", linestyle="--",
                                                   label="Nominal 5%")], loc="upper right", frameon=False)
    fig.text(.5, .025, "Intervals describe repeated simulations under each stated model; they do not establish family-wide calibration.",
             ha="center", fontsize=9)
    fig.savefig(output / "synthetic_null_rates.png", dpi=180)
    fig.savefig(output / "synthetic_null_rates.svg")
    plt.close(fig)


def _write_effect_table(rows, scenarios, output):
    fields = ["scenario", "n_variant", "control_mean", "variant_mean", "true_delta",
              "method", "iterations", "a05_rejected", "a05_rate", "a05_rate_ci_low",
              "a05_rate_ci_high", "a05_directional_rejected", "a01_rejected", "a01_rate",
              "welch_coverage", "welch_coverage_ci_low", "welch_coverage_ci_high",
              "welch_mean_width", "hoeffding_coverage", "hoeffding_coverage_ci_low",
              "hoeffding_coverage_ci_high", "hoeffding_mean_width"]
    with (output / "effect_and_interval_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for index, scenario in enumerate(scenarios):
            if scenario["kind"] != "alternative":
                continue
            for method in METHODS:
                row = rows[("simulation", f"simulation_{index:02d}", "synthetic_fraction", "coverage", method)]
                writer.writerow({field: (f"simulation_{index:02d}" if field == "scenario" else
                                         scenario[field] if field in scenario else row.get(field, ""))
                                 for field in fields})


def plot(calibration_dir, output, allow_pilot=False):
    calibration_dir, output = Path(calibration_dir), Path(output)
    summary_path, audit_path = calibration_dir / "summary.tsv", calibration_dir / "audit.json"
    sources = [_source(path) for path in (summary_path, audit_path, Path(__file__))]
    with audit_path.open() as handle:
        audit = json.load(handle)
    if (audit.get("schema_version") != "fsmpra.heldout_calibration.v1" or
            audit.get("release_withheld") is not True or
            (audit.get("partial") is not False and not allow_pilot)):
        raise ValueError("Figures require complete held-out calibration; use --allow-pilot for a labeled pilot")
    partial = bool(audit["partial"])
    _check_sources(audit["sources"])
    rows = _read_rows(summary_path, audit)
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{output.name}-", dir=output.parent) as stage:
        stage = Path(stage)
        _plot_wt(rows, audit["fixed_targets"], "coverage", stage,
                 int(audit["wt_iterations"]), partial)
        _plot_wt(rows, audit["fixed_targets"], "presence", stage,
                 int(audit["wt_iterations"]), partial)
        _plot_simulation_null(rows, audit["simulation_scenarios"], stage,
                              int(audit["simulation_iterations"]), partial)
        _write_effect_table(rows, audit["simulation_scenarios"], stage)
        if [_source(item["path"]) for item in sources] != sources:
            raise ValueError("Calibration input changed during plotting")
        _check_sources(audit["sources"])
        with (stage / "audit.json").open("w") as handle:
            json.dump({"schema_version": "fsmpra.heldout_calibration_plots.v1", "partial": partial,
                       "release_withheld": True, "wt_targets": 43, "simulation_null_scenarios": 32,
                       "simulation_alternative_scenarios": 18,
                       "nominal_alpha": .05,
                       "interpretation": "Conditional WT A/A and model-specific synthetic diagnostics; no variant significance released.",
                       "sources": sources}, handle, indent=2, allow_nan=False)
            handle.write("\n")
        if output.exists():
            raise FileExistsError(output)
        stage.rename(output)
    return {"wt_cells": 344, "simulation_cells": 100, "output": str(output)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibration", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--allow-pilot", action="store_true")
    args = parser.parse_args(argv)
    print(json.dumps(plot(args.calibration, args.output, args.allow_pilot), sort_keys=True))


if __name__ == "__main__":
    main()

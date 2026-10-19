#!/usr/bin/env python3
"""Plot fixed false-positive diagnostics from a completed WT benchmark.

The plots display recorded rejection counts and exact binomial intervals. They
do not make release decisions or compute variant-associated statistics.
"""
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
COLORS = {"welch": "#2b6cb0", "studentized_permutation": "#c05621"}
TARGETS = (
    ("tf_120_140", "coverage", "TF [120,140) coverage"),
    ("tf_120_140", "presence", "TF [120,140) presence"),
    ("tf_140_160", "coverage", "TF [140,160) coverage"),
    ("tf_140_160", "presence", "TF [140,160) presence"),
    ("msp_promoter", "coverage", "Element MSP coverage"),
)
NULL_MEANS = (.02, .1, .5)
VARIANT_SIZES = (100, 400)
PARTITIONS = ("A", "B")
SUMMARY_FIELDS = {"partition", "scenario", "data_source", "endpoint_id", "metric",
                  "method", "iterations", "true_delta", "rejection_rate", "rejected",
                  "rejection_ci_low", "rejection_ci_high", "interpretation"}


def _source(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    stat = path.stat()
    return {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "sha256": digest.hexdigest()}


def _read_rows(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not SUMMARY_FIELDS <= set(reader.fieldnames or []):
            raise ValueError("Benchmark summary lacks required columns")
        rows = list(reader)
    if not rows:
        raise ValueError("Empty benchmark summary")
    result = {}
    for raw in rows:
        key = (raw["partition"], raw["scenario"], raw["data_source"],
               raw["endpoint_id"], raw["metric"], raw["method"])
        if key in result:
            raise ValueError("Duplicate benchmark summary row: " + str(key))
        row = dict(raw)
        for name in ("iterations", "rejected"):
            row[name] = int(row[name])
        for name in ("true_delta", "rejection_rate", "rejection_ci_low", "rejection_ci_high"):
            row[name] = float(row[name])
        if (row["iterations"] != 200 or row["rejected"] < 0 or
                row["rejected"] > row["iterations"] or
                not math.isclose(row["rejection_rate"], row["rejected"] / 200,
                                 rel_tol=0, abs_tol=1e-12) or
                not 0 <= row["rejection_ci_low"] <= row["rejection_rate"] <=
                    row["rejection_ci_high"] <= 1):
            raise ValueError("Invalid 200-iteration rejection interval: " + str(key))
        result[key] = row
    return result


def _expect_null(rows, partition, scenario, source, endpoint, metric, method):
    key = (partition, scenario, source, endpoint, metric, method)
    if key not in rows:
        raise ValueError("Missing fixed null comparison: " + str(key))
    row = rows[key]
    if (row["interpretation"] != "null_false_positive_rate" or row["true_delta"] != 0):
        raise ValueError("Expected a true-null false-positive row: " + str(key))
    return row


def _validate(rows, audit):
    if (audit.get("schema_version") != "fsmpra.footprint_quality_benchmark.v1" or
            audit.get("partial") is not False or audit.get("iterations") != 200 or
            audit.get("alpha") != .05 or audit.get("release_withheld") is not True):
        raise ValueError("Plots require the complete frozen benchmark")
    wt = {}
    synthetic = {}
    for partition in PARTITIONS:
        for n in VARIANT_SIZES:
            scenario = f"wt_aa_n{n}"
            for endpoint, metric, _ in TARGETS:
                for method in METHODS:
                    key = (partition, n, endpoint, metric, method)
                    wt[key] = _expect_null(rows, partition, scenario, "wt_aa",
                                           endpoint, metric, method)
        for mean in NULL_MEANS:
            for n in VARIANT_SIZES:
                scenario = f"null_mean{mean:g}_n{n}_rhoV.02_rhoC.2"
                for method in METHODS:
                    key = (partition, mean, n, method)
                    synthetic[key] = _expect_null(rows, partition, scenario, "simulation",
                                                   "synthetic_fraction", "coverage", method)
    if (len(rows) != 72 or
            sum(row["interpretation"] == "null_false_positive_rate" for row in rows.values()) != 64):
        raise ValueError("Unexpected number of fixed null method rows")
    return wt, synthetic


def _check_benchmark_sources(audit):
    recorded = audit.get("sources", [])
    if not isinstance(recorded, list) or not recorded:
        raise ValueError("Benchmark audit lacks source records")
    for item in recorded:
        if _source(item["path"]) != item:
            raise ValueError("Benchmark source changed: " + item["path"])


def _axis(ax, items, title):
    """items is a list of (label, method->summary-row) pairs."""
    for index, (label, method_rows) in enumerate(items):
        y = len(items) - 1 - index
        for method, shift in (("welch", -.14), ("studentized_permutation", .14)):
            row = method_rows[method]
            rate = row["rejection_rate"]
            low = row["rejection_ci_low"]
            high = row["rejection_ci_high"]
            ax.errorbar(rate, y + shift, xerr=[[rate - low], [high - rate]],
                        fmt="o", color=COLORS[method], capsize=3, markersize=5,
                        elinewidth=1.35)
            ax.annotate(f'{row["rejected"]}/200', (high, y + shift), xytext=(5, 0),
                        textcoords="offset points", color=COLORS[method], fontsize=7,
                        va="center", clip_on=True)
    ax.axvline(.05, color="0.35", linestyle="--", linewidth=1.1)
    ax.set_yticks(range(len(items)), labels=[item[0] for item in reversed(items)])
    ax.set_ylim(-.6, len(items) - .4)
    ax.set_xlim(0, max(.22, max(row["rejection_ci_high"] for _, group in items
                             for row in group.values()) + .11))
    ax.set_title(title, fontsize=10)
    ax.grid(axis="x", color="0.9", linewidth=.7)
    ax.set_axisbelow(True)


def _legend(fig):
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", linestyle="-", color=COLORS[method],
                      label=METHOD_LABELS[method]) for method in METHODS]
    handles.append(Line2D([0], [0], linestyle="--", color="0.35", label="Nominal 0.05"))
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(.5, .955), fontsize=9)


def _plot_wt(wt, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.subplots_adjust(left=.19, right=.97, top=.86, bottom=.12, hspace=.38, wspace=.39)
    for i, partition in enumerate(PARTITIONS):
        for j, n in enumerate(VARIANT_SIZES):
            items = []
            for endpoint, metric, label in TARGETS:
                items.append((label, {method: wt[partition, n, endpoint, metric, method]
                                      for method in METHODS}))
            _axis(axes[i, j], items, f"WT split {partition}: {n} versus 1,000 barcodes")
            if j:
                axes[i, j].tick_params(axis="y", labelleft=False)
            axes[i, j].set_xlabel("False-positive rejection fraction at α = 0.05")
    fig.suptitle("WT A/A: repeated random contrasts conditional on one sampled WT pool",
                 fontsize=14, y=.99)
    _legend(fig)
    fig.text(.5, .035, "Each point is rejected / 200 contrasts; bars are exact binomial 95% intervals. "
             "Barcode sets may recur across contrasts. No variant discovery is shown.",
             ha="center", fontsize=9)
    fig.savefig(out / "wt_aa_false_positive.png", dpi=180)
    fig.savefig(out / "wt_aa_false_positive.svg")
    plt.close(fig)


def _plot_synthetic(synthetic, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.subplots_adjust(left=.085, right=.98, top=.84, bottom=.15, hspace=.53, wspace=.28)
    for i, partition in enumerate(PARTITIONS):
        for j, mean in enumerate(NULL_MEANS):
            items = [(f"n = {n}", {method: synthetic[partition, mean, n, method]
                                   for method in METHODS}) for n in VARIANT_SIZES]
            _axis(axes[i, j], items, f"Simulation {partition}: equal mean {mean:g}")
            if j:
                axes[i, j].tick_params(axis="y", labelleft=False)
            axes[i, j].set_xlabel("False-positive rejection fraction at α = 0.05")
    fig.suptitle("Synthetic equal-mean nulls: unequal barcode dispersion and depth",
                 fontsize=14, y=.99)
    _legend(fig)
    fig.text(.5, .035, "Each point is rejected / 200 independent simulations; bars are exact binomial 95% intervals. "
             "Variant dispersion 0.02, control 0.2; no variant discoveries.",
             ha="center", fontsize=9)
    fig.savefig(out / "synthetic_null_false_positive.png", dpi=180)
    fig.savefig(out / "synthetic_null_false_positive.svg")
    plt.close(fig)


def plot(benchmark_dir, output):
    benchmark_dir, output = Path(benchmark_dir), Path(output)
    summary_path = benchmark_dir / "summary.tsv"
    audit_path = benchmark_dir / "audit.json"
    sources = [_source(path) for path in (summary_path, audit_path, Path(__file__))]
    with audit_path.open() as handle:
        audit = json.load(handle)
    _check_benchmark_sources(audit)
    rows = _read_rows(summary_path)
    wt, synthetic = _validate(rows, audit)
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{output.name}-", dir=output.parent) as stage:
        out = Path(stage)
        _plot_wt(wt, out)
        _plot_synthetic(synthetic, out)
        power = [row for row in rows.values() if row["interpretation"] == "nominal_power"]
        expected_power = {(partition, scenario, method)
                          for partition in PARTITIONS for scenario in ("known_loss", "known_gain")
                          for method in METHODS}
        if ({(row["partition"], row["scenario"], row["method"]) for row in power} !=
                expected_power or len(power) != 8):
            raise ValueError("Expected eight known-effect power method rows")
        with (out / "known_effect_power.tsv").open("w", newline="") as handle:
            fields = ["partition", "scenario", "method", "iterations", "true_delta",
                      "rejected", "rejection_rate", "rejection_ci_low", "rejection_ci_high"]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(sorted(power, key=lambda row: (row["partition"], row["scenario"], row["method"])))
        if [_source(record["path"]) for record in sources] != sources:
            raise ValueError("Benchmark inputs changed while plotting")
        _check_benchmark_sources(audit)
        report = {"schema_version": "fsmpra.footprint_quality_benchmark_plots.v1",
                  "sample_id": audit.get("sample_id"), "partial": False,
                  "null_rows": {"wt_aa": len(wt), "synthetic": len(synthetic)},
                  "methods": list(METHODS), "iterations_per_row": 200,
                  "interval": "Exact two-sided binomial 95% interval from benchmark summary",
                  "interpretation": "WT A/A is conditional resampling of one WT pool; synthetic equal-mean nulls are simulations. Figures do not assess variant discoveries or release status.",
                  "sources": sources}
        with (out / "audit.json").open("w") as handle:
            json.dump(report, handle, indent=2, allow_nan=False)
            handle.write("\n")
        if output.exists():
            raise FileExistsError(output)
        out.rename(output)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    report = plot(args.benchmark, args.output)
    print(json.dumps({"null_rows": report["null_rows"], "output": str(Path(args.output).resolve())}))


if __name__ == "__main__":
    main()

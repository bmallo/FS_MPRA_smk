"""Independently replay real WT pilot draws using explicit control complements."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.stats import ttest_ind, t

from fsmpra.io import atomic_output, read_json, source_info, write_json
from fsmpra.processing import check_sources


def validate(outcomes, pilot, output):
    outcomes, pilot = Path(outcomes), Path(pilot)
    paths = [outcomes / "endpoint_matrix.npz", pilot / "audit.json",
             pilot / "iterations.tsv.gz", Path(__file__)]
    sources = [source_info(path, True) for path in paths]
    audit = read_json(pilot / "audit.json")
    if not audit["partial"] or audit["wt_iterations"] != 2 or audit["wt_permutations"] != 19:
        raise ValueError("Expected the frozen two-draw/19-permutation pilot")
    with np.load(outcomes / "endpoint_matrix.npz", allow_pickle=False) as data:
        ids = list(data["endpoint_ids"])
        columns = [(name, "presence", j) for j, name in enumerate(ids[:-1])]
        columns += [(name, "coverage", j) for j, name in enumerate(ids)]
        matrix = np.column_stack([data[metric][:, j] for _, metric, j in columns])
    with gzip.open(pilot / "iterations.tsv.gz", "rt") as handle:
        rows = {(r["scenario"], int(r["iteration"]), r["endpoint_id"], r["metric"], r["method"]): r
                for r in csv.DictReader(handle, delimiter="\t") if r["data_source"] == "wt_aa"}
    errors = {key: 0. for key in ("delta", "se", "welch_p", "ci_low", "ci_high", "permutation_p")}
    for n in (50, 100, 400, 1000):
        bank_rng = np.random.default_rng(np.random.SeedSequence([audit["seed"], 0, n]))
        draw_rng = np.random.default_rng(np.random.SeedSequence([audit["seed"], 1, n]))
        reference = []
        for _ in range(19):
            selected = bank_rng.choice(len(matrix), n, replace=False, shuffle=False)
            controls = np.ones(len(matrix), dtype=bool)
            controls[selected] = False
            reference.append(np.abs(ttest_ind(matrix[selected], matrix[controls], equal_var=False, axis=0).statistic))
        reference = np.asarray(reference)
        for iteration in range(2):
            selected = draw_rng.choice(len(matrix), n, replace=False, shuffle=False)
            controls = np.ones(len(matrix), dtype=bool)
            controls[selected] = False
            x, y = matrix[selected], matrix[controls]
            mx, my = x.mean(axis=0), y.mean(axis=0)
            ax, ay = x.var(axis=0, ddof=1) / len(x), y.var(axis=0, ddof=1) / len(y)
            se = np.sqrt(ax + ay)
            df = (ax + ay)**2 / (ax**2 / (len(x)-1) + ay**2 / (len(y)-1))
            direct = ttest_ind(x, y, axis=0, equal_var=False)
            delta = mx - my
            low, high = np.maximum(-1, delta-t.ppf(.975, df)*se), np.minimum(1, delta+t.ppf(.975, df)*se)
            p = (1 + np.sum(reference >= np.abs(direct.statistic)[None, :] - 1e-12, axis=0)) / 20
            for j, (endpoint, metric, _) in enumerate(columns):
                row = rows[(f"wt_{n}", iteration, endpoint, metric, "welch")]
                perm = rows[(f"wt_{n}", iteration, endpoint, metric, "studentized_permutation")]
                expected = {"delta": (delta[j], row["delta"]), "se": (se[j], row["se"]),
                            "welch_p": (direct.pvalue[j], row["p_value"]),
                            "ci_low": (low[j], row["welch_ci_low"]),
                            "ci_high": (high[j], row["welch_ci_high"]),
                            "permutation_p": (p[j], perm["p_value"])}
                for key, (value, stored) in expected.items():
                    error = abs(float(value) - float(stored))
                    if not np.isfinite(error) or error > 1e-10:
                        raise ValueError(f"Direct complement disagreement: {n}, {endpoint}, {metric}, {key}: {error}")
                    errors[key] = max(errors[key], error)
                if int(row["nonzero_variant"]) != np.count_nonzero(x[:, j]) or int(row["nonzero_control"]) != np.count_nonzero(y[:, j]):
                    raise ValueError("Nonzero support disagreement")
        print(f"Validated explicit complements for n={n}", flush=True)
    check_sources(sources)
    report = {"schema_version": "fsmpra.heldout_pilot_validation.v1", "status": "passed",
              "wt_pool_barcodes": len(matrix), "observed_draws": 8, "endpoints_per_draw": 43,
              "direct_reference_draws": 76, "max_absolute_errors": errors, "sources": sources}
    with atomic_output(output) as out:
        write_json(out / "audit.json", report)
        check_sources(sources)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("outcomes", "pilot", "output"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    print(json.dumps(validate(args.outcomes, args.pilot, args.output)["max_absolute_errors"]))

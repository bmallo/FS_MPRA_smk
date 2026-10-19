"""Numerical and provenance checks for held-out diagnostic kernels."""
import json
import csv
import gzip
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.stats import ttest_ind

from fsmpra.heldout_calibration import (
    EXPECTED, bank_p_values, calibrate, centered_pool, contrast_from_pool, draw_barcodes,
    exact_interval, load_outcomes, reference_bank, simulation_scenarios,
    validate_config,
)
from fsmpra.io import source_info
from fsmpra.statistics import fixed_endpoints


class TestHeldoutCalibration(unittest.TestCase):
    def test_config_is_frozen_and_pilot_is_partial(self):
        cfg = {"schema_version": "fsmpra.heldout_calibration_config.v1", **EXPECTED}
        self.assertFalse(validate_config(cfg)[-1])
        self.assertEqual(validate_config(cfg, 2, 19)[1:5], (2, 2, 19, 19))
        self.assertTrue(validate_config(cfg, 2, 19)[-1])
        cfg["wt_permutations"] = 99
        with self.assertRaisesRegex(ValueError, "Frozen"):
            validate_config(cfg)

    def test_complement_matches_direct_welch_and_hoeffding(self):
        rng = np.random.default_rng(123)
        matrix = rng.binomial(3, [.002, .1, .5], (97, 3)) / 3
        chosen = rng.choice(97, 17, replace=False)
        x = matrix[chosen]
        y = np.delete(matrix, chosen, axis=0)
        result = contrast_from_pool(centered_pool(matrix), chosen)
        for j in range(3):
            direct = ttest_ind(x[:, j], y[:, j], equal_var=False)
            self.assertAlmostEqual(result["delta"][j], x[:, j].mean() - y[:, j].mean(), places=13)
            self.assertAlmostEqual(result["se"][j], np.sqrt(x[:, j].var(ddof=1)/17 + y[:, j].var(ddof=1)/80), places=13)
            if np.isfinite(direct.pvalue):
                self.assertAlmostEqual(result["welch_p"][j], direct.pvalue, places=12)
            self.assertEqual(result["nonzero_variant"][j], np.count_nonzero(x[:, j]))
            self.assertEqual(result["nonzero_control"][j], np.count_nonzero(y[:, j]))
            radius = np.sqrt(np.log(4/.05)/(2*17)) + np.sqrt(np.log(4/.05)/(2*80))
            self.assertAlmostEqual(result["hoeffding_ci_low"][j], max(-1, result["delta"][j]-radius))
            self.assertAlmostEqual(result["hoeffding_ci_high"][j], min(1, result["delta"][j]+radius))

    def test_constant_and_tied_reference(self):
        matrix = np.full((32, 2), .25)
        matrix[:, 1] = 0
        pool = centered_pool(matrix)
        result = contrast_from_pool(pool, np.arange(5))
        np.testing.assert_array_equal(result["welch_p"], [1, 1])
        np.testing.assert_array_equal(result["statistic"], [0, 0])
        bank = reference_bank(pool, 5, 19, np.random.default_rng(44))
        np.testing.assert_array_equal(bank, 0)
        p, exceed = bank_p_values(result["statistic"], bank)
        np.testing.assert_array_equal(p, 1)
        np.testing.assert_array_equal(exceed, 19)
        tie_bank = np.array([[0.], [1.], [1.], [2.]])
        p, exceed = bank_p_values([1.], tie_bank)
        self.assertEqual(exceed[0], 3)
        self.assertEqual(p[0], .8)
        decimal_pool = centered_pool(np.full((10100, 1), .1))
        decimal_result = contrast_from_pool(decimal_pool, np.arange(100))
        self.assertEqual(decimal_result["delta"][0], 0.)
        self.assertEqual(decimal_result["se"][0], 0.)
        self.assertEqual(decimal_result["statistic"][0], 0.)
        self.assertEqual(decimal_result["welch_p"][0], 1.)
        np.testing.assert_array_equal(reference_bank(decimal_pool, 100, 19, np.random.default_rng(1)), 0.)

    def test_welch_interval_is_clipped(self):
        matrix = np.zeros((20, 1))
        matrix[:2] = 1
        result = contrast_from_pool(centered_pool(matrix), [0, 2])
        self.assertGreater(result["delta"][0], 0)
        self.assertEqual(result["welch_ci_high"][0], 1.)
        self.assertGreaterEqual(result["welch_ci_low"][0], -1.)

    def test_bank_is_reused_and_matches_direct_label_t(self):
        matrix = np.array([[0.], [1.], [.2], [.7], [.1], [.8], [.5], [.3]])
        pool = centered_pool(matrix)
        seed = 199
        bank = reference_bank(pool, 3, 9, np.random.default_rng(seed))
        rng = np.random.default_rng(seed)
        direct = []
        for _ in range(9):
            selected = rng.choice(8, size=3, replace=False, shuffle=False)
            x = matrix[selected, 0]
            y = np.delete(matrix[:, 0], selected)
            direct.append(abs(ttest_ind(x, y, equal_var=False).statistic))
        np.testing.assert_allclose(bank[:, 0], sorted(direct), rtol=1e-12, atol=1e-12)
        for selected in ([0, 1, 2], [4, 5, 6]):
            result = contrast_from_pool(pool, selected)
            p, exceed = bank_p_values(result["statistic"], bank)
            self.assertEqual(p[0], (1 + exceed[0])/10)

    def test_scenarios_and_depths(self):
        scenarios = simulation_scenarios()
        self.assertEqual(len(scenarios), 50)
        self.assertEqual(sum(s["kind"] == "null" for s in scenarios), 32)
        self.assertEqual(sum(s["kind"] == "alternative" for s in scenarios), 18)
        rng = np.random.default_rng(2)
        for mode in ("independent", "associated"):
            fractions = draw_barcodes(rng, 100, 0, .02, mode, "variant")
            np.testing.assert_array_equal(fractions, 0)

    def test_exact_interval_and_source_guard(self):
        low, high = exact_interval(0, 19)
        self.assertEqual(low, 0)
        self.assertGreater(high, 0)
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            endpoints = fixed_endpoints()
            np.savez(directory / "endpoint_matrix.npz",
                     barcodes=np.array(["A", "B"]),
                     endpoint_ids=[e["endpoint_id"] for e in endpoints],
                     presence=np.ones((2, 22)), coverage=np.ones((2, 22)),
                     depths=np.ones(2, dtype=np.int32), valid_molecules=np.ones((2, 22), dtype=np.int32),
                     relative_starts=[e["relative_start"] for e in endpoints],
                     relative_ends=[e["relative_end"] for e in endpoints],
                     reference_starts0=[e["reference_start0"] for e in endpoints],
                     reference_ends0=[e["reference_end0"] for e in endpoints])
            matrix_source = source_info(directory / "endpoint_matrix.npz", True)
            (directory / "audit.json").write_text(json.dumps(dict(
                schema_version="fsmpra.heldout_outcomes.v1", complete=True,
                heldout_barcodes=258283, heldout_molecules=2,
                endpoint_count=22, sources=[matrix_source])))
            with self.assertRaisesRegex(ValueError, "258283"):
                load_outcomes(directory)
            (directory / "endpoint_matrix.npz").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "Input changed"):
                load_outcomes(directory)

    def test_pilot_outputs_summary_and_provenance(self):
        import fsmpra.heldout_calibration as calibration
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            config = directory / "config.json"
            config.write_text(json.dumps({"schema_version": "fsmpra.heldout_calibration_config.v1",
                                          "sample_id": "HepG2_rep3", **EXPECTED}))
            input_file = directory / "input.txt"
            input_file.write_text("test input")
            input_source = source_info(input_file, True)
            audit = {"input_cohort": {"wt_barcodes": 268283,
                                      "excluded_barcodes": 10000},
                     "sources": [input_source, source_info(config, True)]}
            targets = [(ep["endpoint_id"], "presence", j) for j, ep in enumerate(fixed_endpoints()[:-1])]
            targets += [(ep["endpoint_id"], "coverage", j) for j, ep in enumerate(fixed_endpoints())]
            matrix = np.random.default_rng(4).binomial(1, .03, (1100, 43)).astype(float)
            scenario = dict(kind="null", n_variant=100, control_mean=.02, variant_mean=.02,
                            rho_variant=.02, rho_control=.2, depth_mode="independent")
            output = directory / "pilot"
            wrong_contract_audit = {**audit, "sources": [input_source]}
            with patch.object(calibration, "load_outcomes", return_value=(wrong_contract_audit, targets, matrix, [input_source])):
                with self.assertRaisesRegex(ValueError, "config source is not recorded"):
                    calibrate(config, directory, directory / "wrong_contract", pilot_iterations=2, pilot_permutations=19)
            self.assertFalse((directory / "wrong_contract").exists())
            with patch.object(calibration, "load_outcomes", return_value=(audit, targets, matrix, [input_source])), \
                 patch.object(calibration, "simulation_scenarios", return_value=[scenario]):
                result = calibrate(config, directory, output, pilot_iterations=2, pilot_permutations=19)
            self.assertTrue(result["partial"])
            self.assertEqual(result["contrasts"], 4*2*43 + 2)
            with gzip.open(output / "iterations.tsv.gz", "rt") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 2*result["contrasts"])
            with (output / "summary.tsv").open() as handle:
                summaries = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(summaries), 4*43*2 + 2)
            for summary in summaries:
                self.assertEqual(int(summary["a05_rejected"]), sum(int(summary[k]) for k in
                    ("a05_positive_rejected", "a05_negative_rejected", "a05_zero_rejected")))
                self.assertEqual(int(summary["iterations"]), sum(int(summary[k]) for k in
                    ("welch_covered", "welch_noncoverage_below_truth", "welch_noncoverage_above_truth")))
            report = json.loads((output / "audit.json").read_text())
            self.assertEqual(report["sample_id"], "HepG2_rep3")
            self.assertEqual(report["input_cohort"], audit["input_cohort"])
            self.assertEqual(report["wt_contrasts"], 4*2*43)
            self.assertEqual(report["simulation_contrasts"], 2)
            self.assertTrue(report["release_withheld"])
            self.assertIn("not globally unused", report["cohort_relation"])
            names = {Path(s["path"]).name for s in report["sources"]}
            self.assertTrue({"config.json", "input.txt", "heldout_calibration.py", "io.py", "processing.py"} <= names)
            with self.assertRaisesRegex(ValueError, "Pilot permutations"):
                validate_config({"schema_version": "fsmpra.heldout_calibration_config.v1", **EXPECTED},
                                pilot_iterations=2, pilot_permutations=10000)


if __name__ == "__main__":
    unittest.main()

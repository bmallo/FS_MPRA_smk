import itertools
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.stats import ttest_ind

from fsmpra import footprint_quality_benchmark as benchmark_module
from fsmpra.footprint_quality_benchmark import (
    _draw_wt_groups, _load_summary, benchmark, draw_fraction_barcodes,
    studentized_permutation, welch_diagnostic,
)
from fsmpra.io import source_info
from fsmpra.statistics import fixed_endpoints


CONFIG = Path(__file__).resolve().parents[1] / "config/ldlr_footprint_quality.json"


def _small_summary(root, missing=0, complete=True):
    """Make a complete-contract summary with simple independently labeled WT rows."""
    root = Path(root)
    config = root / "config.json"
    config.write_text(CONFIG.read_text())
    summary = root / "summary"
    summary.mkdir()
    endpoints = fixed_endpoints()
    n = 10_000
    valid = np.ones((n, len(endpoints)), dtype=np.int32)
    presence = np.zeros_like(valid, dtype=float)
    coverage = np.zeros_like(valid, dtype=float)
    if missing:
        valid[:missing, 12] = 0
        presence[:missing, 12] = np.nan
        coverage[:missing, 12] = np.nan
    matrix = summary / "barcode_endpoint_matrix.npz"
    np.savez_compressed(matrix,
        barcodes=np.array([f"WT{i:05d}" for i in range(n)]),
        endpoint_ids=np.array([e["endpoint_id"] for e in endpoints]),
        presence=presence, coverage=coverage, valid_molecules=valid,
        depths=np.ones(n, dtype=np.int32),
        splits=np.array(["A"] * 5000 + ["B"] * 5000),
        relative_starts=np.array([e["relative_start"] for e in endpoints]),
        relative_ends=np.array([e["relative_end"] for e in endpoints]),
        reference_starts0=np.array([e["reference_start0"] for e in endpoints]),
        reference_ends0=np.array([e["reference_end0"] for e in endpoints]))
    report = dict(schema_version="fsmpra.footprint_quality_summary.v1",
                  complete=complete, pilot=not complete, sample_id="sample",
                  sampled_barcodes=n, endpoint_count=len(endpoints),
                  endpoint_matrix=matrix.name, split_barcodes={"A": 5000, "B": 5000},
                  sources=[source_info(config, True), source_info(matrix, True)])
    (summary / "summary.json").write_text(json.dumps(report))
    return config, summary


class TestStatistics(unittest.TestCase):
    def test_welch_matches_independent_scipy_calculation(self):
        cases = [([0, 0, 1, 0, 0], [0, 1, 1, 0, 1, 1]),
                 ([0, 0, 0, 0], [0, 1, 0, 0, 0]),
                 ([1, 1, 1, 1], [0, 0, 0, 0])]
        for x, y in cases:
            with self.subTest(x=x, y=y):
                got = welch_diagnostic(x, y)
                expected = ttest_ind(x, y, equal_var=False)
                self.assertAlmostEqual(got["statistic"], expected.statistic)
                self.assertAlmostEqual(got["p_value"], expected.pvalue)
        flat = welch_diagnostic([0, 0], [0, 0])
        self.assertEqual((flat["status"], flat["p_value"]), ("uninformative", 1.0))

    def test_permutation_matches_independent_label_draws(self):
        x = np.array([0., 0., 1., .5])
        y = np.array([0., .2, .8, 1., 1.])
        draws = 37
        seed = 313
        got = studentized_permutation(x, y, draws, np.random.default_rng(seed), batch_size=7)
        rng = np.random.default_rng(seed)
        pooled = np.r_[x, y]
        observed = abs(ttest_ind(x, y, equal_var=False).statistic)
        exceeded = 0
        for size in [7] * 5 + [2]:
            keys = rng.random((size, len(pooled)))
            for row in keys:
                selected = np.argpartition(row, len(x)-1)[:len(x)]
                other = np.setdiff1d(np.arange(len(pooled)), selected)
                stat = ttest_ind(pooled[selected], pooled[other], equal_var=False).statistic
                exceeded += abs(stat) >= observed - 1e-12
        self.assertEqual(got["exceedances"], exceeded)
        self.assertEqual(got["p_value"], (exceeded + 1) / (draws + 1))
        self.assertGreaterEqual(got["p_value"], 1 / (draws + 1))

    def test_constant_cases_and_exact_label_enumeration(self):
        x, y = [0, 0], [1, 1]
        got = studentized_permutation(x, y, 30, np.random.default_rng(11))
        self.assertEqual(got["status"], "zero_se_difference")
        self.assertGreaterEqual(got["p_value"], 1 / 31)
        all_flat = studentized_permutation([0, 0], [0, 0], 30, np.random.default_rng(11))
        self.assertEqual((all_flat["p_value"], all_flat["exceedances"]), (1., 30))
        # The four choose two labelings have two infinite and four zero t values.
        pooled = [0, 0, 1, 1]
        stats = [welch_diagnostic([pooled[i] for i in ids],
                                  [pooled[i] for i in range(4) if i not in ids])["statistic"]
                 for ids in itertools.combinations(range(4), 2)]
        self.assertEqual(sum(np.isinf(s) for s in stats), 2)

    def test_simulation_mean_and_bounded_fraction(self):
        rng = np.random.default_rng(42)
        for mean in (.02, .1, .5):
            for rho, depth_p in ((.02, .3), (.2, .5)):
                sample = draw_fraction_barcodes(rng, 20_000, mean, rho, depth_p)
                self.assertTrue(np.all((0 <= sample) & (sample <= 1)))
                self.assertLess(abs(sample.mean() - mean), .025)
        self.assertTrue(np.all(draw_fraction_barcodes(rng, 100, 0, .02, .3) == 0))

    def test_wt_draws_are_disjoint_and_reproducible(self):
        pool = np.arange(5000)
        a, b = _draw_wt_groups(np.random.default_rng(3), pool, 400, 1000)
        aa, bb = _draw_wt_groups(np.random.default_rng(3), pool, 400, 1000)
        np.testing.assert_array_equal(a, aa)
        np.testing.assert_array_equal(b, bb)
        self.assertEqual(len(np.unique(np.r_[a, b])), 1400)


class TestInputContract(unittest.TestCase):
    def test_complete_summary_and_explicit_missing_pool(self):
        with tempfile.TemporaryDirectory() as tmp:
            config, summary = _small_summary(tmp, missing=10)
            report, arrays = _load_summary(config, summary)
            self.assertEqual(report["sampled_barcodes"], 10_000)
            self.assertEqual(len(arrays["barcodes"]), 10_000)
            output = Path(tmp) / "pilot"
            result = benchmark(config, summary, output, pilot_iterations=1, pilot_permutations=1)
            self.assertTrue(result["partial"])
            self.assertEqual(result["excluded_barcodes_missing_fixed_target"], 10)
            audit = json.loads((output / "audit.json").read_text())
            self.assertEqual(audit["complete_barcodes_by_partition"], {"A": 4990, "B": 5000})
            self.assertEqual(audit["sampled_barcodes"], 10_000)
            self.assertEqual(audit["missing_barcodes_by_target"]["tf_120_140:coverage"], 10)

    def test_rejects_pilot_and_changed_upstream_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            config, summary = _small_summary(tmp, complete=False)
            with self.assertRaisesRegex(ValueError, "complete 10,000-barcode"):
                _load_summary(config, summary)
            report_path = summary / "summary.json"
            report = json.loads(report_path.read_text())
            report["complete"], report["pilot"] = True, False
            report_path.write_text(json.dumps(report))
            config.write_text(config.read_text() + " ")
            with self.assertRaisesRegex(ValueError, "provenance mismatch|Input changed"):
                _load_summary(config, summary)

    def test_rejects_source_change_during_benchmark_before_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            config, summary = _small_summary(tmp)
            output = Path(tmp) / "pilot"
            original = studentized_permutation
            changed = False

            def change_once(*args, **kwargs):
                nonlocal changed
                if not changed:
                    config.write_text(config.read_text() + " ")
                    changed = True
                return original(*args, **kwargs)

            with patch.object(benchmark_module, "studentized_permutation", side_effect=change_once):
                with self.assertRaisesRegex(ValueError, "Input changed during processing"):
                    benchmark(config, summary, output, pilot_iterations=1, pilot_permutations=1)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()

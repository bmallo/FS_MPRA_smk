import unittest

import numpy as np

from fsmpra.statistics import welch_contrast, welch_from_summary
from fsmpra.statistics_calibration import (draw_barcodes, exact_interval, moments,
    release_gate, subtract_moments, summarize_stratum, welch_summary)


class CalibrationTests(unittest.TestCase):
    def test_disjoint_complement_matches_direct_contrast(self):
        rng = np.random.default_rng(912)
        x = rng.beta(.3, 2, (500, 3))
        chosen = rng.choice(len(x), 70, replace=False)
        keep = np.ones(len(x), dtype=bool)
        keep[chosen] = False
        part = moments(x[chosen])
        complement = welch_summary(subtract_moments(moments(x), part))
        direct = welch_contrast(x[chosen], x[keep])
        fast = welch_from_summary(*welch_summary(part), *complement)
        for key in ("delta", "se", "df", "ci_low", "ci_high", "p_value"):
            np.testing.assert_allclose(direct[key], fast[key], rtol=1e-10, atol=1e-14)
        self.assertEqual(complement[0], 430)
        constant = np.full((500, 1), .3)
        raw = subtract_moments(moments(constant), moments(constant[:70]))
        self.assertEqual(welch_summary(raw)[2][0], 0.)

    def test_binomial_limits_and_unconditional_denominator(self):
        self.assertEqual(exact_interval(0, 0), (None, None))
        low, high = exact_interval(0, 200)
        self.assertEqual(low, 0)
        self.assertAlmostEqual(high, 1 - .025**(1/200))
        cfg = dict(minimum_testable_iterations=100, null_rate_upper_limit=.1,
                   ci_coverage_lower_limit=.9)
        rows = [dict(status="ok", truth=0., p_value=.5, ci_low=-.1, ci_high=.1)] * 100
        rows += [dict(status="sparse_endpoint", truth=0.)] * 100
        rows[0] = dict(status="ok", truth=0., p_value=.01, ci_low=.01, ci_high=.1)
        summary = summarize_stratum(rows, cfg)
        self.assertEqual(summary["unconditional_rejection_rate"], 1/200)
        self.assertEqual(summary["conditional_ci_coverage"], .99)
        self.assertEqual(summary["status"], "pass")
        self.assertEqual(summarize_stratum(rows[1:], cfg)["status"], "unresolved")
        bad = [dict(status="ok", truth=0., p_value=.01, ci_low=.01, ci_high=.1)] * 200
        self.assertEqual(summarize_stratum(bad, cfg)["status"], "fail")

    def test_release_gate_retains_failures_and_unresolved_cells(self):
        endpoints = [dict(endpoint_id="a"), dict(endpoint_id="b")]
        strata = [dict(stratum="aa", kind="wt_aa", endpoint_id="a", status="pass"),
                  dict(stratum="rare", kind="wt_aa", endpoint_id="b", status="unresolved"),
                  dict(stratum="sim", kind="simulation", null=True, status="pass")]
        got = release_gate(strata, endpoints)
        self.assertEqual(got["status"], "pass")
        self.assertEqual(got["validated_endpoint_ids"], ["a"])
        unresolved = dict(stratum="unresolved_null", kind="simulation", null=True, status="unresolved")
        self.assertEqual(release_gate(strata + [unresolved], endpoints)["status"], "unresolved")
        strata[-1]["status"] = "fail"
        self.assertEqual(release_gate(strata, endpoints)["status"], "fail")
        self.assertEqual(release_gate([], endpoints)["status"], "unresolved")

    def test_generator_reproducibility_and_known_sign(self):
        def draw(seed, mean):
            return draw_barcodes(np.random.default_rng(seed), 10000, mean, .2, True, "variant")
        x, nx = draw(12, .1)
        y, ny = draw(12, .1)
        np.testing.assert_array_equal(x, y)
        self.assertEqual(nx, ny)
        higher, _ = draw(12, .3)
        self.assertGreater(higher.mean() - x.mean(), .15)
        self.assertLess(abs(x.mean() - .1), .02)
        self.assertTrue(np.all((x >= 0) & (x <= 1)))


if __name__ == "__main__":
    unittest.main()

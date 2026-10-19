"""Focused checks for held-out WT molecule-to-barcode outcomes."""
import unittest

import numpy as np

from fsmpra.heldout_outcomes import _endpoint_observations, _mask, _reconcile
from fsmpra.statistics import fixed_endpoints


class HeldoutOutcomesTests(unittest.TestCase):
    def setUp(self):
        self.endpoints = fixed_endpoints()
        self.row = {
            "eligible_ref_intervals0": "[[1234,1445]]",
            "tf_small_ref_intervals0": "[[1354,1364],[1369,1374]]",
            "msp_ref_intervals0": "[[1234,1245]]",
        }

    def test_presence_is_per_molecule_before_barcode_average(self):
        first_presence, first_coverage = _endpoint_observations(self.row, self.endpoints, 1234, 1445)
        empty = dict(self.row, tf_small_ref_intervals0="[]", msp_ref_intervals0="[]")
        second_presence, second_coverage = _endpoint_observations(empty, self.endpoints, 1234, 1445)
        endpoint_ids = [x["endpoint_id"] for x in self.endpoints]
        j = endpoint_ids.index("tf_120_140")
        self.assertEqual((first_presence[j] + second_presence[j]) / 2, 0.5)
        self.assertAlmostEqual((first_coverage[j] + second_coverage[j]) / 2, 15 / 40)
        self.assertAlmostEqual((first_coverage[-1] + second_coverage[-1]) / 2, 11 / 422)
        self.assertEqual((first_presence[-1] + second_presence[-1]) / 2, 0.5)

    def test_bad_intervals_and_incomplete_alignment_fail(self):
        for malformed in ('[[1234,1240],[1239,1241]]', '[[1234,1446]]',
                          '[[1234,1234]]', '[[1234.0,1240]]', 'null'):
            with self.subTest(malformed=malformed), self.assertRaises(ValueError):
                _mask(malformed, 1234, 1445, "TF")
        row = dict(self.row, eligible_ref_intervals0="[[1234,1444]]")
        with self.assertRaisesRegex(ValueError, "full element"):
            _endpoint_observations(row, self.endpoints, 1234, 1445)

    def test_coverage_reconciles_every_endpoint(self):
        counts = np.zeros((2, 2, 211), dtype=np.uint32)
        counts[0, 0, 120:130] = 1
        counts[0, 1, :11] = 1
        counts[1, 0, 120:125] = 2
        counts[1, 1, 90:110] = 2
        depths = np.array([2, 4], dtype=np.int32)
        coverage = np.empty((2, 22), dtype=np.float64)
        for j, endpoint in enumerate(self.endpoints):
            track = 0 if endpoint["track"] == "tf_small" else 1
            a, b = endpoint["relative_start"], endpoint["relative_end"]
            coverage[:, j] = counts[:, track, a:b].sum(axis=1) / ((b-a) * depths)
        self.assertLessEqual(_reconcile(coverage, depths, np.array([0, 1]), counts, self.endpoints), 1e-12)
        coverage[1, -1] += 1e-4
        with self.assertRaisesRegex(ValueError, "disagrees"):
            _reconcile(coverage, depths, np.array([0, 1]), counts, self.endpoints)


if __name__ == "__main__":
    unittest.main()

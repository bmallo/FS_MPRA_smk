import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.stats import false_discovery_control, ttest_ind

from fsmpra import statistics as statistics_module
from fsmpra.statistics import (aggregate, by_adjust, endpoint_matrix,
    fixed_endpoints, infer, welch_contrast, welch_from_summary)
from fsmpra.io import source_info


class StatisticsTests(unittest.TestCase):
    def test_endpoint_grid(self):
        endpoints = fixed_endpoints()
        self.assertEqual(len(endpoints), 22)
        self.assertEqual([(x["relative_start"], x["relative_end"]) for x in endpoints[:3]],
                         [(0, 20), (10, 30), (20, 40)])
        self.assertEqual((endpoints[-2]["relative_start"], endpoints[-2]["relative_end"]), (191, 211))
        self.assertEqual(endpoints[-1]["endpoint_id"], "msp_promoter")

    def test_welch_and_summary_agree_with_scipy(self):
        rng = np.random.default_rng(723)
        x = rng.beta(2, 5, (60, 3))
        y = rng.beta(2, 4, (80, 3))
        result = welch_contrast(x, y)
        reference = ttest_ind(x, y, axis=0, equal_var=False)
        np.testing.assert_allclose(result["p_value"], reference.pvalue, rtol=1e-12)
        np.testing.assert_array_equal(result["status"], "ok")
        def summary(a):
            return len(a), a.mean(axis=0), a.var(axis=0, ddof=1), (a > 0).sum(axis=0), (a < 1).sum(axis=0)
        alternate = welch_from_summary(*summary(x), *summary(y))
        for field in ("delta", "se", "df", "ci_low", "ci_high", "p_value"):
            np.testing.assert_allclose(result[field], alternate[field])
        np.testing.assert_allclose(result["p_value"], welch_contrast(x[::-1], y[::-1])["p_value"])

    def test_sparse_degenerate_and_sign(self):
        sparse = welch_contrast(np.zeros((50, 1)), np.ones((50, 1)))
        self.assertEqual(sparse["status"][0], "sparse_endpoint")
        self.assertTrue(np.isnan(sparse["p_value"][0]))
        self.assertEqual(welch_contrast(np.full((50, 1), .3), np.full((50, 1), .4))["status"][0],
                         "degenerate_variance")
        rng = np.random.default_rng(3)
        x, y = rng.beta(5, 2, (60, 1)), rng.beta(2, 5, (60, 1))
        self.assertGreater(welch_contrast(x, y)["delta"][0], 0)
        self.assertLess(welch_contrast(y, x)["delta"][0], 0)
        self.assertEqual(welch_contrast(x[:49], y)["status"][0], "insufficient_barcodes")

    def test_by_fixed_family(self):
        p = np.array([.01, .03, np.nan, .2])
        q = by_adjust(p, 4)
        h = 1 + 1/2 + 1/3 + 1/4
        expected = np.array([.01*4*h, .03*4*h/2, np.nan, .2*4*h/3])
        np.testing.assert_allclose(q[np.isfinite(q)], expected[np.isfinite(expected)])
        reference = false_discovery_control(np.nan_to_num(p, nan=1.), method="by")
        np.testing.assert_allclose(q[np.isfinite(q)], reference[np.isfinite(q)])
        self.assertTrue(np.isnan(q[2]))

    def test_aggregation_conservation_and_inference(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            variants = [dict(variant_id="wt", source_variant="WT", **{"class": "WT"},
                             ref_pos0=None, ref="", alt="")]
            variants += [dict(variant_id=f"v{i}", source_variant=f"{i}A>G", **{"class": "SNV"},
                              ref_pos0=1234+i%211, ref="A", alt="G") for i in range(633)]
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({"library": {"barcode_length": 4, "variants": variants}}))
            cfg = root / "config.json"
            cfg.write_text(json.dumps(dict(schema_version="fsmpra.statistics_config.v1",
                sample_id="sample", interval=[1234, 1445], max_barcodes=8,
                tracks=["tf_small", "msp"], primary_cohort="footprint_analysis_eligible_with_barcode",
                primary_weighting="equal_barcode", tf_window_width=20, tf_window_step=10,
                include_terminal_window=True, min_barcodes=50, min_molecules=100,
                min_nonzero_barcodes=10, min_nonone_barcodes=10,
                confidence_level=.95, fdr_method="BY", fdr_alpha=.05,
                sensitivity_mask_radius=2, seed=17012026, interpretation="test")))
            joined = root / "joined"
            joined.mkdir()
            rows = []
            for bc, vid, intervals in (("AAAA", "wt", [[1234, 1244]]),
                                       ("AAAA", "wt", [[1239, 1250]]),
                                       ("CCCC", "v0", [[1234, 1238]])):
                rows.append(dict(footprint_analysis_eligible="1",
                    footprint_analysis_eligible_with_barcode="1", observed_variant_id=vid,
                    assigned_variant_id=vid, corrected_barcode=bc, barcode_assignment_status="supported",
                    tf_small_ref_intervals0=json.dumps(intervals), msp_ref_intervals0="[]"))
            rows.append(dict(footprint_analysis_eligible="1",
                footprint_analysis_eligible_with_barcode="0", observed_variant_id="wt",
                assigned_variant_id="", corrected_barcode="", barcode_assignment_status="unsupported",
                tf_small_ref_intervals0="[]", msp_ref_intervals0="[]"))
            rows.append(dict(footprint_analysis_eligible="0",
                footprint_analysis_eligible_with_barcode="0", observed_variant_id="",
                assigned_variant_id="", corrected_barcode="", barcode_assignment_status="unsupported",
                tf_small_ref_intervals0="", msp_ref_intervals0=""))
            with gzip.open(joined / "molecule_footprints.tsv.gz", "wt") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            (joined / "audit.json").write_text(json.dumps(dict(sample_id="sample", interval=[1234,1445],
                sources=[source_info(catalog, True)],
                counts=dict(output_molecules=5, footprint_analysis_eligible=4,
                            footprint_analysis_eligible_with_barcode=3))))
            aggregate(cfg, catalog, joined, root / "aggregates")
            counts = np.load(root / "aggregates/barcode_counts.npy")
            np.testing.assert_array_equal(counts.shape, (2, 2, 211))
            np.testing.assert_array_equal(counts[0, 0, :16],
                [1]*5 + [2]*5 + [1]*6)
            matrix = endpoint_matrix(counts, np.load(root / "aggregates/barcode_molecules.npy"))
            self.assertAlmostEqual(matrix[0, 0], 21 / 40)
            result = infer(cfg, catalog, root / "aggregates", root / "inference")
            self.assertEqual(result, dict(testable=0, supported=False))
            with gzip.open(root / "inference/effects.tsv.gz", "rt") as handle:
                effects = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(effects), 633*22)
            self.assertEqual(effects[0]["candidate_p_value"], "")
            self.assertEqual(effects[0]["status"], "insufficient_molecules")
            with gzip.open(root / "inference/sensitivities.tsv.gz", "rt") as handle:
                sensitivity = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(sensitivity), 633*22*3)
            first = {r["sensitivity"]: r for r in sensitivity
                     if r["variant_id"] == "v0" and r["endpoint_id"] == "tf_000_020"}
            masked = first["mutation_pm2_masked_equal_barcode"]
            self.assertAlmostEqual(float(masked["mean_variant"]), 1/17)
            self.assertAlmostEqual(float(masked["mean_wt"]), 9/17)
            self.assertEqual((masked["omitted_start0"], masked["omitted_end0"]), ("1234", "1237"))
            self.assertAlmostEqual(float(first["all_eligible_molecule_weighted"]["mean_wt"]), 21/60)
            self.assertAlmostEqual(float(first["supported_molecule_weighted"]["mean_wt"]), 21/40)
            report_path = root / "calibration.json"
            sources = [source_info(path, checksum) for path, checksum in (
                (cfg, True), (root / "aggregates/audit.json", True),
                (Path(statistics_module.__file__), True),
                (root / "aggregates/barcode_counts.npy", False),
                (root / "aggregates/barcode_molecules.npy", False),
                (root / "aggregates/barcode_metadata.tsv.gz", True))]
            report = dict(schema_version="fsmpra.statistics_calibration.v1",
                sample_id="sample", status="pass", validated_endpoint_ids=["tf_000_020"],
                sources=sources)
            report_path.write_text(json.dumps(report))
            infer(cfg, catalog, root / "aggregates", root / "calibrated_inference", report_path)
            release = json.loads((root / "calibrated_inference/release_status.json").read_text())
            self.assertEqual(release["statistical_calibration"], "pass")
            self.assertEqual(release["supported_effects"], 0)
            report["sources"][0]["sha256"] = "bad"
            report_path.write_text(json.dumps(report))
            with self.assertRaisesRegex(ValueError, "Calibration provenance mismatch"):
                infer(cfg, catalog, root / "aggregates", root / "stale_inference", report_path)
            with self.assertRaises(FileExistsError):
                aggregate(cfg, catalog, joined, root / "aggregates")


if __name__ == "__main__":
    unittest.main()

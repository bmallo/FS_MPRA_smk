"""Small conservation tests for descriptive barcode weighting and state bins."""
import csv
import json
import math
import tempfile
import unittest
from pathlib import Path

from fsmpra.heterogeneity_summary import FEATURES, summarize
from fsmpra.io import open_text, source_info, write_json, write_table


class HeterogeneitySummaryTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        catalog_dir = self.root / "ledger" / "catalog"
        catalog_dir.mkdir(parents=True)
        write_json(catalog_dir / "catalog.json", {"library": {"variants": [
            {"variant_id": "WT", "source_variant": "WT"},
            {"variant_id": "pooled_SNV", "source_variant": "1A>G"},
        ]}})
        self.cfg = self.root / "config.json"
        write_json(self.cfg, {"ledger_run_dir": str(self.root / "ledger"), "summary": {
            "max_barcodes": 10, "histogram_bins": 4,
            "state_edges": [0, .1, .25, .5, .75, 1],
            "plot_min_barcodes": 2, "interpretation": "Synthetic descriptive test"}})
        self.features = self.root / "features"
        self.features.mkdir()

    def _row(self, barcode, variant, distal_m6a, distal_nuc, tf, msp):
        row = {feature: "" for feature in FEATURES if feature != "distal_tile_methylation_range"}
        row.update(sample_id="sample", library_id="library", observed_variant_id=variant,
                   genotype_class="WT" if variant == "WT" else "SNV", corrected_barcode=barcode,
                   assigned_variant_id=variant, footprint_analysis_eligible="1",
                   footprint_analysis_eligible_with_barcode="1",
                   distal_m6a_fraction=str(distal_m6a), distal_nucleosome_fraction=str(distal_nuc),
                   element_tf_small_fraction="" if tf is None else str(tf),
                   element_msp_fraction=str(msp))
        return row

    def _prepare(self):
        rows = [self._row("A", "WT", 0, 0, 0, 0),
                self._row("A", "WT", 1, 1, 1, 1),
                self._row("B", "WT", 1, 1, None, 0),
                self._row("C", "pooled_SNV", .5, 1, .5, .5)]
        for row in rows:
            row.update(whole_coverage_fraction="1", element_coverage_fraction="1",
                       distal_coverage_fraction="1")
        rows[0].update(distal_coverage_fraction="0.9", whole_nucleosome_segments="0",
                       whole_tf_calls_small="0")
        rows[1].update(whole_nucleosome_segments="2", whole_tf_calls_small="1")
        rows[2].update(whole_coverage_fraction="0.9", whole_nucleosome_segments="2",
                       whole_tf_calls_small="1")
        rows[3].update(whole_coverage_fraction="", element_coverage_fraction="0.9",
                       whole_tf_calls_small="3")
        write_table(self.features / "molecule_heterogeneity.tsv.gz", rows, list(rows[0]))
        write_json(self.features / "audit.json", {"schema_version": "fsmpra.heterogeneity.v1",
            "complete": False, "limit_primary": 4, "counts": {
                "selected_rows": 4, "selected_footprint_analysis_eligible": 4,
                "selected_footprint_analysis_eligible_with_barcode": 4},
            "sources": [source_info(self.cfg, True),
                        source_info(self.root / "ledger" / "catalog" / "catalog.json", True)]})

    @staticmethod
    def _table(path):
        with open_text(path) as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_partial_guard(self):
        self._prepare()
        with self.assertRaisesRegex(ValueError, "--allow-partial"):
            summarize(self.cfg, self.features, self.root / "reject")
        self.assertFalse((self.root / "reject").exists())

    def test_config_identity_must_match_extraction(self):
        self._prepare()
        altered = json.loads(self.cfg.read_text())
        altered["summary"]["histogram_bins"] = 5
        write_json(self.cfg, altered)
        with self.assertRaisesRegex(ValueError, "config/catalog differs"):
            summarize(self.cfg, self.features, self.root / "wrong_config", allow_partial=True)

    def test_molecular_weighting_missingness_and_state_boundaries(self):
        self._prepare()
        output = self.root / "summary"
        result = summarize(self.cfg, self.features, output, allow_partial=True)
        self.assertEqual(result["primary_barcodes"], 3)
        self.assertAlmostEqual(result["primary_equal_barcode_weight"], 3)
        self.assertEqual(result["state_binned_molecules"], 4)
        self.assertAlmostEqual(result["state_binned_weight"], 3)
        self.assertTrue((output / "wt_heterogeneity.svg").exists())
        self.assertTrue((output / "variant_heterogeneity.png").exists())

        variant = self._table(output / "variant_heterogeneity.tsv.gz")
        wt = next(r for r in variant if r["cohort"] == "footprint_analysis_eligible_with_barcode"
                  and r["variant_id"] == "WT" and r["feature"] == "distal_m6a_fraction")
        self.assertEqual(int(wt["molecules"]), 3)  # Literal variant ID cannot collide with group WT.
        self.assertAlmostEqual(float(wt["molecule_mean"]), 2 / 3)
        self.assertAlmostEqual(float(wt["barcode_mean"]), .75)
        tf = next(r for r in variant if r["cohort"] == "footprint_analysis_eligible_with_barcode"
                  and r["variant_id"] == "WT" and r["feature"] == "element_tf_small_fraction")
        self.assertEqual(int(tf["valid_molecules"]), 2)
        self.assertEqual(int(tf["valid_barcodes"]), 1)
        self.assertAlmostEqual(float(tf["barcode_mean"]), .5)

        barcode = self._table(output / "barcode_heterogeneity.tsv.gz")
        a = next(r for r in barcode if r["corrected_barcode"] == "A")
        b = next(r for r in barcode if r["corrected_barcode"] == "B")
        self.assertAlmostEqual(float(a["distal_m6a_fraction_within_sd"]), .5)
        self.assertEqual(b["distal_m6a_fraction_within_sd"], "")  # Singleton is not a spread estimate.
        self.assertEqual(b["element_tf_small_fraction_mean"], "")
        self.assertEqual(int(b["element_tf_small_fraction_valid_molecules"]), 0)

        histogram = self._table(output / "histograms.tsv.gz")
        bins = [r for r in histogram if r["cohort"] == "footprint_analysis_eligible_with_barcode"
                and r["genotype_group"] == "WT" and r["weighting"] == "equal_barcode"
                and r["feature"] == "distal_m6a_fraction"]
        self.assertEqual(len(bins), 4)
        self.assertAlmostEqual(float(bins[0]["count_weight"]), .5)
        self.assertAlmostEqual(float(bins[3]["count_weight"]), 1.5)  # x=1 belongs to final bin.
        boundaries = self._table(output / "histogram_boundaries.tsv.gz")
        boundary = next(r for r in boundaries if r["cohort"] == "footprint_analysis_eligible_with_barcode"
                        and r["genotype_group"] == "WT" and r["weighting"] == "equal_barcode"
                        and r["feature"] == "distal_m6a_fraction")
        self.assertAlmostEqual(float(boundary["exact_zero_weight"]), .5)
        self.assertAlmostEqual(float(boundary["exact_one_weight"]), 1.5)
        self.assertAlmostEqual(result["WT_exact_boundary_proportions"]["equal_barcode"]["distal_m6a_fraction"]["zero"], .25)
        self.assertEqual(result["within_barcode_sd"]["WT"]["distal_m6a_fraction"]["barcodes_with_at_least_two_valid_molecules"], 1)
        self.assertEqual(result["SNV_variants_at_plot_threshold"], 0)
        counts = self._table(output / "count_histograms.tsv.gz")
        wt_nuc = [r for r in counts if r["cohort"] == "footprint_analysis_eligible_with_barcode"
                  and r["genotype_group"] == "WT" and r["weighting"] == "equal_barcode"
                  and r["feature"] == "whole_nucleosome_segments"]
        self.assertEqual({r["count_value"]: float(r["count_weight"]) for r in wt_nuc},
                         {"0": .5, "2": 1.5})
        missing_nuc = next(r for r in counts if r["cohort"] == "footprint_analysis_eligible_with_barcode"
                           and r["genotype_group"] == "pooled_SNV" and r["weighting"] == "equal_barcode"
                           and r["feature"] == "whole_nucleosome_segments" and r["count_value"] == "")
        self.assertEqual(float(missing_nuc["missing_weight"]), 1)
        count_summary = result["WT_discrete_call_count_distribution"]["equal_barcode"]["whole_nucleosome_segments"]
        self.assertEqual(count_summary["median"], 2)
        self.assertAlmostEqual(count_summary["exact_zero_proportion"], .25)
        coverage = result["cohort_region_coverage"]["footprint_analysis_eligible_with_barcode"]["all"]
        self.assertEqual(coverage["whole"], {"rows": 4, "valid": 3, "missing": 1, "complete": 2})
        self.assertEqual(coverage["distal"]["complete"], 3)
        states = self._table(output / "variant_state_summaries.tsv.gz")
        low = next(r for r in states if r["variant_id"] == "WT" and r["distal_m6a_low"] == "0.0"
                   and r["distal_nucleosome_low"] == "0.0")
        high = next(r for r in states if r["variant_id"] == "WT" and r["distal_m6a_low"] == "0.75"
                    and r["distal_nucleosome_low"] == "0.75")
        self.assertAlmostEqual(float(low["sum_weight"]), .5)
        self.assertAlmostEqual(float(high["sum_weight"]), 1.5)
        self.assertEqual(int(high["unique_barcodes"]), 2)
        self.assertAlmostEqual(float(high["bin_proportion"]), .75)


if __name__ == "__main__":
    unittest.main()

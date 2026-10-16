import copy
import json
import shutil
import unittest

from fsmpra.catalog import build_catalog
from fsmpra.compare import check_compatible, compare_pilots, pairwise, pearson
from fsmpra.dashboard import export_dashboard
from fsmpra.features import measure
from fsmpra.io import read_json, source_info
from fsmpra.molecules import import_molecules
import test_pilot


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.fixture = test_pilot.PipelineTests()
        self.fixture.setUp()
        self.root = self.fixture.root
        catalog = self.root / "catalog"
        build_catalog(self.fixture.manifest, catalog)
        import_molecules(catalog / "catalog.json", "rep1", self.root / "molecules", fraction=1)
        measure(catalog / "catalog.json", self.root / "molecules", self.root / "measurement", tf_min=1, tf_max=5)
        self.a = read_json(self.root / "measurement/pilot.json")
        self.b = copy.deepcopy(self.a)
        self.b["sample"]["id"] = "rep2"
        shutil.copyfile(self.root / "input.bam", self.root / "other.bam")
        self.b["provenance"]["bam"] = source_info(self.root / "other.bam")
        self.paths = [self.root / "a.json", self.root / "b.json"]
        for path, data in zip(self.paths, [self.a, self.b]):
            path.write_text(json.dumps(data))
        manifest = read_json(self.fixture.manifest)
        s = manifest["samples"][0]
        s.update(reaction_condition_id="10min", reaction={"duration_minutes": 10, "hia5_volume_uL": 1.5, "cells": 2000000})
        second = copy.deepcopy(s)
        second.update(id="rep2", reaction_condition_id="5min", reaction={"duration_minutes": 5, "hia5_volume_uL": 3, "cells": 2000000})
        second["inputs"]["footprints_bam"] = "other.bam"
        manifest["samples"].append(second)
        self.fixture.manifest.write_text(json.dumps(manifest))

    def tearDown(self):
        self.fixture.tearDown()

    def test_condition_metadata_overlay_does_not_rewrite_inputs(self):
        before = [p.read_bytes() for p in self.paths]
        result = compare_pilots(self.paths, self.fixture.manifest, self.root / "comparison", min_reads=1, include_processing=False)
        self.assertEqual(result["pairwise"][0]["comparison_type"], "different_reaction_conditions")
        self.assertEqual(result["pairwise"][0]["SNVs_above_minimum_in_both"], 2)
        self.assertIsNone(result["pairwise"][0]["msp_common_site_delta_pearson"])
        self.assertEqual(result["sample_qc"][0]["full_plasmid_fraction_of_metrics_matched_sample"], 1)
        self.assertEqual(before, [p.read_bytes() for p in self.paths])
        bundle = self.root / "comparison/collection.json"
        collection = read_json(bundle)
        self.assertEqual(collection["datasets"][0]["sample"]["reaction"]["duration_minutes"], 10)
        self.assertIn("metadata_manifest", collection["datasets"][0]["provenance"])
        out = export_dashboard(bundle, self.root / "dashboard")
        self.assertGreater(out["bytes"], 1000)

    def test_incompatible_inputs_rejected(self):
        changes = [
            (lambda b: b["sample"].update(id="rep1"), "Duplicate sample"),
            (lambda b: b["provenance"].update(bam=self.a["provenance"]["bam"]), "Duplicate BAM"),
            (lambda b: b.update(reference_sequence_sha256="wrong"), "reference mismatch"),
            (lambda b: b["region"].update(start0=5), "interval mismatch"),
            (lambda b: b["method"].update(min_passes=20), "method mismatch"),
            (lambda b: b["method"].update(footprint_contract={"format":"fiberhmm_ma_sparse","evaluated_tracks":["tf"]}), "contracts differ"),
            (lambda b: b["variants"][1].update(alt="X"), "allele mismatch"),
        ]
        for change, error in changes:
            with self.subTest(error=error):
                b = copy.deepcopy(self.b)
                change(b)
                with self.assertRaisesRegex(ValueError, error):
                    check_compatible([self.a, b])

    def test_coverage_intersection_and_reaction_group(self):
        self.a["sample"]["reaction_condition_id"] = self.b["sample"]["reaction_condition_id"] = "same"
        result = pairwise([self.a, self.b], 10)[0]
        self.assertEqual(result["SNVs_above_minimum_in_both"], 0)
        self.assertEqual(result["comparison_type"], "same_stated_reaction")
        self.assertIsNone(result["msp_common_site_delta_pearson"])
        self.assertEqual(pearson([(1, 3), (2, 2), (3, 1)]), -1)
        self.assertIsNone(pearson([(1, 2), (1, 3), (1, 4)]))

    def test_conflicting_reaction_annotation_rejected(self):
        manifest = read_json(self.fixture.manifest)
        manifest["samples"][1]["reaction_condition_id"] = "10min"
        self.fixture.manifest.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "conflicting experimental"):
            compare_pilots(self.paths, self.fixture.manifest, self.root / "bad", include_processing=False)


if __name__ == "__main__":
    unittest.main()

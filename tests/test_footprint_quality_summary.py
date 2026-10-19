import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from fsmpra.footprint_quality_summary import Calls, FEATURES, summarize
from fsmpra.io import source_info


def table(path, rows, fields):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class SummaryTests(unittest.TestCase):
    def test_first_aq_byte_controls_sensitivity(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calls.tsv"
            rows = []
            for read_id, quality in (("r1", [150, 20, 30]), ("r2", [50, 200, 200])):
                rows.append(dict(read_id=read_id, barcode="A", split="A", weight=1,
                    track="tf", band="10_60", original_length=10,
                    ref_intervals0="[[1234,1244]]", quality_spec="QQQ",
                    quality=json.dumps(quality), query_at_count=10,
                    query_m6a128_count=0, query_m6a204_count=0))
            table(path, rows, list(rows[0]))
            calls = Calls(path)
            try:
                selected = [calls.for_molecule(dict(read_id=read_id, corrected_barcode="A",
                    sample_split="A", barcode_weight=1))[0]["aq100"] for read_id in ("r1", "r2")]
                self.assertEqual(selected, [True, False])
            finally:
                calls.close()

    def fixture(self, root):
        root = Path(root)
        sample = root / "sample"
        sample.mkdir()
        start = 1234
        catalog = root / "catalog.json"
        catalog.write_text(json.dumps({"library": {"variants": [
            {"variant_id": "wt", "class": "WT", "sequence": "A" * 211}]}}))
        config = root / "config.json"
        config.write_text(json.dumps({"schema_version": "fsmpra.footprint_quality_config.v1",
                                      "catalog": str(catalog), "sample_barcodes": 2}))
        barcodes = [dict(barcode="A", variant_id="wt", eligible_depth=2, split="A",
                         weight=.5, selection_hash="a"),
                    dict(barcode="B", variant_id="wt", eligible_depth=1, split="B",
                         weight=1, selection_hash="b")]
        table(sample / "sampled_barcodes.tsv", barcodes, list(barcodes[0]))
        molecules = []
        for read, barcode, split, weight, tf, meth in (
                ("r1", "A", "A", .5, list(range(10, 20)), [0]),
                ("r2", "A", "A", .5, [10], None),
                ("r3", "B", "B", 1, [], [])):
            row = {field: "" for field in FEATURES}
            row.update(read_id=read, sample_id="sample", library_id="library",
                       corrected_barcode=barcode, assigned_variant_id="wt",
                       sample_split=split, barcode_weight=weight,
                       representative_group="no_nuc_with_localTF" if tf else "no_nuc_without_localTF",
                       element_aligned_offsets0=json.dumps(list(range(211))),
                       element_at_offsets0=json.dumps(list(range(211))),
                       element_m6a_offsets0=json.dumps(meth or []),
                       element_tf_small_offsets0=json.dumps(tf),
                       element_msp_offsets0="[]",
                       element_m6a_fraction="" if meth is None else len(meth) / 211,
                       element_tf_small_fraction=len(tf) / 211,
                       element_nucleosome_fraction=0, element_msp_fraction=0,
                       whole_nucleosome_segments=0, whole_tf_calls_small=int(bool(tf)))
            molecules.append(row)
        table(sample / "sampled_molecules.tsv.gz", molecules, list(molecules[0]))
        calls = []
        for read, length, projected_length, quality, m6a in (("r1", 10, 10, [120], 1),
                                                              ("r2", 11, 1, None, 0)):
            calls.append(dict(read_id=read, barcode="A", split="A", weight=.5,
                track="tf", band="10_60",
                original_length=length,
                ref_intervals0=json.dumps([[start + 10, start + 10 + projected_length]]),
                quality_spec="Q", quality=json.dumps(quality),
                query_at_count=length, query_m6a128_count=m6a, query_m6a204_count=0))
        table(sample / "sampled_calls.tsv.gz", calls, list(calls[0]))
        hist = [dict(track="tf", band=call["band"], original_length=call["original_length"],
                     call_count=1, weighted_call_count=.5) for call in calls]
        table(sample / "call_length_histograms.tsv", hist, list(hist[0]))
        with gzip.open(sample / "representatives.jsonl.gz", "wt"):
            pass
        audit = {"schema_version": "fsmpra.footprint_quality.v1", "complete": True,
                 "element": [start, start + 211],
                 "counts": {"sampled_barcodes": 2, "sampled_expected_molecules": 3,
                            "matched_molecules": 3, "sampled_calls": 2},
                 "sources": [source_info(config, True), source_info(catalog, True)]}
        (sample / "audit.json").write_text(json.dumps(audit))
        return config, sample

    def test_weights_missingness_endpoints_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            config, sample = self.fixture(tmp)
            output = Path(tmp) / "summary"
            report = summarize(config, sample, output)
            self.assertTrue(report["complete"])
            self.assertEqual(report["split_barcodes"], {"A": 1, "B": 1})
            self.assertEqual(report["sampled_molecules"], 3)
            with np.load(output / "barcode_endpoint_matrix.npz") as matrix:
                self.assertEqual(matrix["presence"].dtype, np.dtype("float64"))
                self.assertEqual(list(matrix["barcodes"]), ["A", "B"])
                self.assertEqual(matrix["presence"].shape, (2, 22))
                self.assertEqual(matrix["presence"][0, 1], 1)
                self.assertAlmostEqual(matrix["coverage"][0, 1], .275)
                self.assertEqual(matrix["presence"][1, 1], 0)
                self.assertEqual(matrix["valid_molecules"][0, 1], 2)
            with gzip.open(output / "element_profiles.tsv.gz", "rt") as handle:
                profiles = list(csv.DictReader(handle, delimiter="\t"))
            p = {(row["split"], row["group"], row["weighting"], int(row["relative_pos0"])): row
                 for row in profiles}
            self.assertAlmostEqual(float(p["all", "all", "molecule", 10]["tf_coverage_fraction"]), 2 / 3)
            self.assertAlmostEqual(float(p["all", "all", "equal_barcode", 10]["tf_coverage_fraction"]), .5)
            self.assertEqual(float(p["all", "all", "molecule", 0]["at_observed_weight"]), 2)
            self.assertEqual(float(p["all", "all", "equal_barcode", 0]["at_observed_weight"]), 1.5)
            self.assertEqual(float(p["all", "all", "molecule", 10]["tf_aq100_covered_weight"]), 1)
            self.assertEqual(report["aq100_missing_quality_molecules"]["tf"], 1)
            with gzip.open(output / "call_observability_histograms.tsv.gz", "rt") as handle:
                observability = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(sum(int(row["calls"]) for row in observability
                                 if row["metric"] == "query_at_count"), 2)
            self.assertEqual(report["sources"][-1]["path"], str((output / "barcode_endpoint_matrix.npz").resolve()))
            with gzip.open(output / "sequence_context_opportunities.tsv.gz", "rt") as handle:
                context = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(context["trinucleotide"], "NAA")

    def test_pilot_guard_and_call_bitmap_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            config, sample = self.fixture(tmp)
            audit_path = sample / "audit.json"
            audit = json.loads(audit_path.read_text())
            audit["complete"] = False
            audit_path.write_text(json.dumps(audit))
            with self.assertRaisesRegex(ValueError, "Pilot sample"):
                summarize(config, sample, Path(tmp) / "blocked")
            report = summarize(config, sample, Path(tmp) / "pilot", allow_pilot=True)
            self.assertTrue(report["pilot"])
            molecules_path = sample / "sampled_molecules.tsv.gz"
            with gzip.open(molecules_path, "rt") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            rows[0]["element_tf_small_offsets0"] = "[]"
            table(molecules_path, rows, list(rows[0]))
            with self.assertRaisesRegex(ValueError, "TF bitmap"):
                summarize(config, sample, Path(tmp) / "mismatch", allow_pilot=True)


if __name__ == "__main__":
    unittest.main()

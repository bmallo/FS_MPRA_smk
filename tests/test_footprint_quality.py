import csv
import gzip
import json
import tempfile
import unittest
from array import array
from pathlib import Path

import pysam

from fsmpra.footprint_quality import (audit_read, band_for, choose_representatives,
                                      digest, query_call_counts, representative_group,
                                      select_barcodes)
from fsmpra.heterogeneity import META_FIELDS, decode_m6a, extract_features


class SamplingTests(unittest.TestCase):
    def test_lowest_hash_sampling_and_barcode_disjoint_split(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "metadata.tsv.gz"
            rows = [{"barcode": f"B{i:02d}", "variant_id": "wt", "molecules": str(i + 1)}
                    for i in range(12)]
            rows.append({"barcode": "OTHER", "variant_id": "snv", "molecules": "99"})
            with gzip.open(path, "wt", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["barcode", "variant_id", "molecules"],
                                        delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            sample, counts = select_barcodes(path, "wt", 4, 20261018)
            expected = sorted(rows[:12], key=lambda row: (digest("20261018|" + row["barcode"]),
                                                           row["barcode"]))[:4]
            self.assertEqual(set(sample), {row["barcode"] for row in expected})
            self.assertEqual(counts["wt_barcodes"], 12)
            self.assertEqual(counts["wt_eligible_molecules"], sum(range(1, 13)))
            self.assertEqual(counts["sampled_expected_molecules"],
                             sum(int(row["molecules"]) for row in expected))
            for row in expected:
                info = sample[row["barcode"]]
                self.assertEqual(info["weight"], 1 / int(row["molecules"]))
                split = "A" if int(digest("split|20261018|" + row["barcode"]), 16) % 2 == 0 else "B"
                self.assertEqual(info["split"], split)

    def test_representative_selection_depends_only_on_frozen_feature_groups(self):
        rows = {f"r{i}": {"representative_group": group}
                for i, group in enumerate(["no_nuc_with_localTF"] * 5 +
                                          ["no_nuc_without_localTF"] * 4 + ["other"] * 3)}
        first = choose_representatives(rows, 20261018, 2)
        reverse = choose_representatives(dict(reversed(list(rows.items()))), 20261018, 2)
        self.assertEqual(first, reverse)
        self.assertEqual(len(first), 6)
        self.assertEqual(sum(rows[read]["representative_group"] == "other" for read in first), 2)
        self.assertEqual(representative_group({"whole_nucleosome_segments": "0",
                                               "element_tf_small_fraction": "0.02"}),
                         "no_nuc_with_localTF")
        self.assertEqual(representative_group({"whole_nucleosome_segments": "1",
                                               "whole_nucleosome_length_median": "1001"}),
                         "exactly_one_nuc_longer1000")

    def test_exact_band_edges(self):
        self.assertEqual([band_for("nuc", n) for n in (84, 85, 250, 251, 500, 501, 1000, 1001)],
                         ["lt85", "85_250", "85_250", "251_500", "251_500",
                          "501_1000", "501_1000", "gt1000"])
        self.assertEqual([band_for("tf", n) for n in (9, 10, 60, 61)],
                         ["lt10", "10_60", "10_60", "gt60"])


class MoleculeAuditTests(unittest.TestCase):
    def make_read_and_feature(self, reverse=False, seq=None, cigar=None, mm=None, ml=None):
        read = pysam.AlignedSegment()
        read.query_name = "r1"
        read.flag = 16 if reverse else 0
        read.reference_start = 0
        read.cigarstring = cigar or "4647M"
        read.query_sequence = seq or "A" * 4647
        read.set_tag("MM", mm if mm is not None else ("T-a.,3412;" if reverse else "A+a.;"))
        read.set_tag("ML", array("B", ml if ml is not None else ([255] if reverse else [])))
        read.set_tag("MA", f"{read.query_length};tf+QQQ:1235-20;nuc+QQQ:1001-1100;msp+:1235-50")
        read.set_tag("AQ", array("B", [10, 20, 30, 40, 50, 60]))
        spec = {"reference_length": 4647, "regions": {"whole": [(0, 4647)],
            "element": [(1234, 1445)], "distal": [(0, 1134), (1545, 3269), (3334, 4647)]},
            "tf_min": 10, "tf_max": 60, "nuc_min": 85, "tile_size": 100, "tile_min_at": 20}
        old = {name: "" for name in META_FIELDS}
        _, at_count, pos128, _ = decode_m6a(read)
        old.update(read_id="r1", read_length="4647", strand="-" if reverse else "+", m6a_status="ok",
                   at_count=str(at_count), m6a_count=str(len(pos128)), footprint_status="evaluated",
                   corrected_barcode="BARCODE", assigned_variant_id="wt",
                   genotype_class="WT", footprint_analysis_eligible_with_barcode="1")
        feature = extract_features(read, old, spec)
        feature["representative_group"] = representative_group(feature)
        return read, feature

    def test_call_bands_bitmaps_and_baseline_reconciliation(self):
        read, feature = self.make_read_and_feature()
        info = {"barcode": "BARCODE", "split": "A", "weight": 0.25}
        molecule, calls, histogram, representative = audit_read(
            read, feature, info, {"reference_length": 4647, "element": [1234, 1445]}, True)
        self.assertEqual(molecule["whole_tf_10_60_count"], 1)
        self.assertEqual(molecule["whole_nuc_gt1000_count"], 1)
        self.assertEqual(molecule["element_nuc_gt1000_count"], 1)
        self.assertEqual(json.loads(molecule["element_aligned_offsets0"]), list(range(211)))
        self.assertEqual(json.loads(molecule["element_at_offsets0"]), list(range(211)))
        self.assertEqual(json.loads(molecule["element_m6a_offsets0"]), [])
        self.assertEqual(json.loads(molecule["element_tf_small_offsets0"]), list(range(20)))
        self.assertEqual(json.loads(molecule["element_msp_offsets0"]), list(range(50)))
        self.assertEqual(len(calls), 3)
        self.assertEqual([call["query_at_count"] for call in calls], [20, 1100, 50])
        self.assertEqual([call["query_m6a128_count"] for call in calls], [0, 0, 0])
        self.assertEqual(histogram[("nuc", "gt1000", 1100)], 1)
        self.assertEqual(representative["raw_ml_histogram"][0], 0)
        self.assertEqual(len(representative["calls"]), 3)
        feature["element_tf_small_fraction"] = "0.5"
        with self.assertRaisesRegex(ValueError, "TF fraction"):
            audit_read(read, feature, info, {"reference_length": 4647, "element": [1234, 1445]})

    def test_reverse_m6a_projection_uses_stored_query_orientation(self):
        read, feature = self.make_read_and_feature(reverse=True)
        info = {"barcode": "BARCODE", "split": "B", "weight": 0.5}
        molecule, _, _, representative = audit_read(
            read, feature, info, {"reference_length": 4647, "element": [1234, 1445]}, True)
        self.assertEqual(json.loads(molecule["element_m6a_offsets0"]), [0])
        self.assertEqual(representative["whole_m6a_positions0"], [1234])
        self.assertEqual([call["query_m6a128_count"] for call in representative["calls"]], [1, 1, 1])
        self.assertEqual([call["query_m6a204_count"] for call in representative["calls"]], [1, 1, 1])

    def test_call_counts_include_unaligned_insertion_and_separate_thresholds(self):
        seq = "A" * 1234 + "CC" + "A" * 3411
        read, feature = self.make_read_and_feature(
            seq=seq, cigar="1234M2I3411M", mm="A+a.,1234,0;", ml=[180, 220])
        info = {"barcode": "BARCODE", "split": "A", "weight": 0.25}
        _, calls, _, representative = audit_read(
            read, feature, info, {"reference_length": 4647, "element": [1234, 1445]}, True)
        by_track = {call["track"]: call for call in calls}
        self.assertEqual(by_track["tf"]["original_length"], 20)
        self.assertEqual(by_track["tf"]["query_at_count"], 18)
        self.assertEqual(by_track["tf"]["query_m6a128_count"], 2)
        self.assertEqual(by_track["tf"]["query_m6a204_count"], 1)
        self.assertEqual(by_track["tf"]["ref_intervals0"], "[[1234,1252]]")
        self.assertEqual(by_track["msp"]["query_at_count"], 48)
        self.assertEqual(representative["calls"][0]["query_m6a204_count"], 1)

    def test_invalid_query_call_ranges_are_explicit(self):
        import numpy as np
        prefix = np.array([0, 1, 2], dtype=np.int32)
        empty = np.array([], dtype=np.int32)
        for start, length in ((-1, 1), (0, 0), (1, 2)):
            with self.subTest(start=start, length=length):
                with self.assertRaisesRegex(ValueError, "MA query interval outside stored SEQ"):
                    query_call_counts(start, length, prefix, empty, empty)


if __name__ == "__main__":
    unittest.main()

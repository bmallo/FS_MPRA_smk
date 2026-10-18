import csv
import json
import unittest
from array import array

import pysam

from fsmpra.footprint_join import annotate, join_footprints
from fsmpra.io import open_text, write_json
from tests import test_processing
from fsmpra.processing import extract
from fsmpra.barcodes import assign_barcodes
from fsmpra.ledger import join_ledger


class AnnotationTests(unittest.TestCase):
    def read(self, cigar="20M", reverse=False):
        r = pysam.AlignedSegment()
        r.query_name = "test"
        r.flag = 16 if reverse else 0
        r.reference_start = 100
        r.cigarstring = cigar
        r.query_sequence = "A" * r.infer_query_length()
        return r

    def test_missing_is_not_empty(self):
        r = self.read()
        row, _, _ = annotate(r, 100, 120)
        self.assertEqual(row["footprint_status"], "missing_MA")
        self.assertEqual(row["msp_fraction"], "")
        r.set_tag("MA", "20")
        row, _, _ = annotate(r, 100, 120)
        self.assertEqual(row["footprint_status"], "evaluated")
        self.assertEqual(row["msp_fraction"], 0)
        self.assertEqual(row["tf_ref_calls0"], "[]")

    def test_indel_projection_and_quality_orientation(self):
        for reverse in (False, True):
            r = self.read("5M2I5M3D5M", reverse)
            r.set_tag("MA", "17;tf+QQQ:4-12")
            r.set_tag("AQ", array("B", [70, 80, 90]))
            r.set_tag("AN", "call1")
            row, _, _ = annotate(r, 100, 118)
            self.assertEqual(json.loads(row["tf_small_ref_intervals0"]), [[103, 110], [113, 116]])
            self.assertEqual(row["aligned_region_bases"], 15)
            self.assertEqual(row["tf_small_fraction"], "")
            call = json.loads(row["tf_ref_calls0"])[0]
            self.assertEqual(call["query_length"], 12)
            self.assertEqual(call["quality"], [70, 80, 90])
            self.assertEqual(call["name"], "call1")

    def test_lengths_apply_before_clipping_and_overlaps_are_unioned(self):
        r = self.read("120M")
        r.set_tag("MA", "120;tf+QQQ:1-80,91-10,96-10")
        row, _, _ = annotate(r, 150, 210)
        self.assertEqual(json.loads(row["tf_small_ref_intervals0"]), [[190, 205]])
        self.assertAlmostEqual(row["tf_small_fraction"], 15 / 60)
        self.assertEqual(len(json.loads(row["tf_ref_calls0"])), 3)

    def test_malformed_annotations_fail(self):
        r = self.read()
        for ma in ("21", "20;tf+QQQ:15-10", "20;tf+QQQ:1-10;tf+QQQ:11-5"):
            r.set_tag("MA", ma)
            with self.assertRaises(ValueError):
                annotate(r, 100, 120)
        r.set_tag("MA", "20;tf+QQQ:1-10")
        r.set_tag("AQ", array("B", [1, 2]))
        with self.assertRaises(ValueError):
            annotate(r, 100, 120)


class JoinTests(unittest.TestCase):
    def test_join_by_id_retains_missing_and_rejects_duplicate_or_changed_sequence(self):
        fixture = test_processing.ProcessingTests()
        fixture.setUp()
        self.addCleanup(fixture.tearDown)
        reads = [fixture.read("a"), fixture.read("b"), fixture.read("c"), fixture.read("d")]
        bam = fixture.write_bam(reads)
        extracted, barcodes, ledger = [fixture.root / n for n in ("extracted", "barcodes", "ledger")]
        catalog = fixture.catalog / "catalog.json"
        extract(fixture.path, catalog, bam, extracted)
        assign_barcodes(fixture.path, catalog, extracted / "observations.tsv.gz", barcodes)
        join_ledger(fixture.path, catalog, extracted, barcodes, ledger)
        reads[0].set_tag("MA", str(reads[0].query_length))
        reads[1].set_tag("MA", str(reads[1].query_length) + ";msp+:9-12")
        fp = fixture.root / "fp.bam"
        contract = fixture.root / "contract.json"
        write_json(contract, {"bam": str(fp), "sample_id": "sample", "format": "fiberhmm_ma_sparse",
                              "evaluated_tracks": ["tf", "nuc", "msp"]})

        def run(records, name):
            with pysam.AlignmentFile(fp, "wb", header=fixture.header) as dest:
                for r in records:
                    dest.write(r)
            pysam.index(str(fp))
            return join_footprints(fixture.path, catalog, ledger / "molecule_ledger.tsv.gz", fp,
                                   contract, fixture.root / name, max_disk_gib=0.1)

        audit = run([reads[2], reads[1], reads[0]], "joined")
        self.assertEqual(audit["counts"]["output_molecules"], 4)
        self.assertEqual(audit["counts"]["missing_bam"], 1)
        self.assertEqual(audit["counts"]["missing_MA"], 1)
        with open_text(fixture.root / "joined/molecule_footprints.tsv.gz") as f:
            rows = {r["read_id"]: r for r in csv.DictReader(f, delimiter="\t")}
        self.assertEqual(rows["a"]["msp_fraction"], "0.0")
        self.assertEqual(rows["b"]["msp_fraction"], "1.0")
        self.assertEqual(rows["d"]["footprint_status"], "missing_bam")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            run([reads[0], reads[0]], "duplicate")
        reads[0].query_sequence = reads[0].query_sequence[:10] + "T" + reads[0].query_sequence[11:]
        with self.assertRaisesRegex(ValueError, "sequence or barcode mismatch"):
            run([reads[0]], "changed")


if __name__ == "__main__":
    unittest.main()

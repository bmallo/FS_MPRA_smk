import csv
import gzip
import tempfile
import unittest
from array import array
from pathlib import Path

import pysam

from fsmpra.heterogeneity import META_FIELDS, decode_m6a, extract, extract_features, selected
from fsmpra.io import read_json, source_info, write_json


def make_read(sequence, reverse=False, cigar=None):
    read = pysam.AlignedSegment()
    read.query_name = "molecule"
    read.flag = 16 if reverse else 0
    read.reference_start = 0
    read.query_sequence = sequence
    read.cigarstring = cigar or f"{len(sequence)}M"
    return read


def row(read, status="ok", count=0, footprint="evaluated"):
    result = {name: "" for name in META_FIELDS}
    result.update(read_id=read.query_name, sample_id="sample", library_id="library",
                  observed_variant_id="variant", genotype_class="SNV", genotype_status="exact_design",
                  read_length=str(read.query_length), strand="-" if read.is_reverse else "+",
                  m6a_status=status, at_count=str(read.query_sequence.count("A") + read.query_sequence.count("T")),
                  m6a_count=str(count), sequence_qc="1", length_qc="0", methylation_qc="0",
                  analysis_eligible="0", footprint_status=footprint, region_complete="0")
    return result


SPEC = {"reference_length": 18, "regions": {"whole": [(0, 18)], "element": [(4, 14)],
                                           "distal": [(0, 4), (14, 18)]},
        "tf_min": 10, "tf_max": 60, "nuc_min": 10, "tile_size": 2, "tile_min_at": 1}


class MethylationDecoderTests(unittest.TestCase):
    def test_reverse_mixed_groups_and_repeated_target_deduplicate(self):
        read = make_read("ATATACGT", reverse=True)
        # Forward SEQ is ACGTATAT. A positions 0,4,6 map to stored 7,3,1;
        # T positions 3,5,7 map to stored 4,2,0.
        read.set_tag("MM", "A+am.,0,0;C+m.,0;T-a.,0,0,0;A+a.,0,0;")
        read.set_tag("ML", array("B", [100, 140, 240, 35, 35, 128, 204, 255, 220, 30]))
        status, at, p128, p204 = decode_m6a(read)
        self.assertEqual((status, at), ("ok", 6))
        self.assertEqual(set(p128), {0, 2, 3, 4, 7})
        self.assertEqual(set(p204), {0, 2, 3, 7})

    def test_unknown_omitted_and_malformed_ml(self):
        read = make_read("AAAA")
        read.set_tag("MM", "A+a?,0;")
        read.set_tag("ML", array("B", [255]))
        status, at, p128, _ = decode_m6a(read)
        self.assertEqual((status, at, list(p128)), ("unknown_opportunities", 4, [0]))
        read.set_tag("ML", array("B", [255, 200]))
        with self.assertRaisesRegex(ValueError, "Extra ML"):
            decode_m6a(read)
        read.set_tag("MM", "A+a.,4;")
        read.set_tag("ML", array("B", [255]))
        with self.assertRaisesRegex(ValueError, "MM/SEQ"):
            decode_m6a(read)

    def test_multiple_codes_consume_interleaved_ml(self):
        read = make_read("AAAA")
        read.set_tag("MM", "A+am.,0,1;")
        read.set_tag("ML", array("B", [130, 255, 220, 10]))
        _, _, p128, p204 = decode_m6a(read)
        self.assertEqual(list(p128), [0, 2])
        self.assertEqual(list(p204), [2])


class RegionalFeatureTests(unittest.TestCase):
    def test_reverse_indel_projection_and_missingness(self):
        read = make_read("A" * 17, reverse=True, cigar="5M2I5M3D5M")
        read.set_tag("MM", "T-a.," + ",".join("0" for _ in range(17)) + ";")
        read.set_tag("ML", array("B", [255] * 17))
        read.set_tag("MA", "17;tf+QQQ:4-12;nuc+QQQ:1-17;msp+:1-17")
        read.set_tag("AQ", array("B", [10, 20, 30, 40, 50, 60]))
        result = extract_features(read, row(read, count=17), SPEC)
        self.assertEqual(result["whole_aligned_bases"], 15)
        self.assertEqual(result["whole_at_count"], 15)
        self.assertEqual(result["whole_m6a_count"], 15)
        self.assertEqual(result["whole_m6a_fraction"], 1)
        self.assertEqual(result["whole_nucleosome_segments"], 1)
        self.assertEqual(result["whole_nucleosome_length_median"], 17)
        self.assertEqual(result["whole_tf_calls_all"], 1)
        self.assertEqual(result["whole_tf_calls_small"], 1)
        self.assertEqual(result["distal_tile_methylation_median"], 1)
        self.assertEqual(result["element_coverage_fraction"], 0.7)
        self.assertEqual(result["element_tf_small_fraction"], 1)
        self.assertEqual(result["element_nucleosome_fraction"], 1)

        read.set_tag("MA", None)
        missing = extract_features(read, row(read, count=17, footprint="missing_MA"), SPEC)
        self.assertEqual(missing["whole_tf_small_fraction"], "")
        self.assertEqual(missing["whole_m6a_fraction"], 1)
        read.set_tag("AQ", None)
        read.set_tag("MA", "17")
        empty = extract_features(read, row(read, count=17), SPEC)
        self.assertEqual(empty["whole_tf_small_fraction"], 0)

    def test_unknown_methylation_does_not_assign_zero_fraction(self):
        read = make_read("A" * 18)
        read.set_tag("MM", "A+a?,0;")
        read.set_tag("ML", array("B", [255]))
        result = extract_features(read, row(read, status="unknown_opportunities", count=1,
                                            footprint="missing_MA"), SPEC)
        self.assertEqual(result["whole_m6a_count"], 1)
        self.assertEqual(result["whole_m6a_fraction"], "")
        self.assertEqual(result["distal_tile_status"], "unknown_opportunities")
        self.assertEqual(result["whole_tf_small_fraction"], "")

    def test_selection_ignores_old_length_methylation_and_barcode_qc(self):
        read = make_read("A" * 18)
        values = row(read)
        self.assertTrue(selected(values))
        values["barcode_assignment_status"] = "unassigned"
        self.assertTrue(selected(values))
        values["sequence_qc"] = "0"
        self.assertFalse(selected(values))


class StreamTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.ledger = self.root / "ledger_run"
        self.join = self.root / "footprint_run"
        for folder in (self.ledger / "catalog", self.ledger / "barcodes", self.join):
            folder.mkdir(parents=True)
        self.bam = self.root / "annotations.bam"
        header = pysam.AlignmentHeader.from_dict({"HD": {"SO": "coordinate"},
                                                   "SQ": [{"SN": "plasmid", "LN": 500}]})
        reads = []
        for name in ("first", "second"):
            read = pysam.AlignedSegment(header)
            read.query_name = name
            read.query_sequence = "A" * 500
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.mapping_quality = 60
            read.cigarstring = "500M"
            read.set_tag("MM", "A+a.;")
            read.set_tag("ML", array("B"))
            read.set_tag("MA", "500")
            reads.append(read)
        with pysam.AlignmentFile(self.bam, "wb", header=header) as bam:
            for read in reads:
                bam.write(read)
        pysam.index(str(self.bam))
        self.contract = self.root / "contract.json"
        write_json(self.contract, {"bam": str(self.bam), "sample_id": "sample"})
        self.ledger_cfg = self.ledger / "config.json"
        write_json(self.ledger_cfg, {"fixture": True})
        self.catalog = self.ledger / "catalog/catalog.json"
        write_json(self.catalog, {"library": {"id": "library", "contig": "plasmid",
            "reference_length": 500, "mutable_ref_start0": 150, "mutable_ref_end0": 160,
            "barcode_start0": 350, "barcode_end0": 365,
            "variants": [{"variant_id": "vWT", "class": "WT"}]}})
        self.coverage = self.ledger / "barcodes/variant_barcode_coverage.tsv"
        self.coverage.write_text("variant_id\tcohort\tmolecules\n"
                                 "vWT\tsequence_qualified\t1\n")
        first = row(reads[0], count=0)
        first.update(observed_variant_id="vWT", genotype_class="WT", region_complete="1",
                     analysis_eligible="1", footprint_analysis_eligible="1",
                     footprint_analysis_eligible_with_barcode="1", tf_small_fraction="0",
                     nucleosome_fraction="0", msp_fraction="0")
        second = row(reads[1], count=0)
        second.update(observed_variant_id="", genotype_class="unresolved",
                      genotype_status="unresolved", sequence_qc="0")
        self.rows = [first, second]
        self.joined = self.join / "molecule_footprints.tsv.gz"
        self.write_joined(self.rows)
        self.audit = self.join / "audit.json"
        write_json(self.audit, {"schema_version": "fsmpra.footprint_join.v1",
            "sample_id": "sample", "library_id": "library",
            "footprint_contract": {"bam": str(self.bam)},
            "sources": [source_info(self.contract, True), source_info(self.bam),
                        source_info(str(self.bam) + ".bai"), source_info(self.ledger_cfg, True),
                        source_info(self.catalog, True)],
            "counts": {"matched_primary_molecules": 2, "output_molecules": 2,
                       "ledger_analysis_eligible": 1, "footprint_analysis_eligible": 1,
                       "footprint_analysis_eligible_with_barcode": 1}})
        self.config = self.root / "heterogeneity.json"
        write_json(self.config, {"schema_version": "fsmpra.heterogeneity_config.v1",
            "ledger_run_dir": str(self.ledger), "footprint_run_dir": str(self.join),
            "footprint_contract": str(self.contract),
            "regions": {"whole": [[0, 500]], "element": [[150, 160]],
                        "distal": [[0, 50], [260, 325], [390, 500]]},
            "thresholds": [128, 204], "tf_min": 10, "tf_max": 60, "nuc_min": 85,
            "tile_size": 100, "tile_min_at": 20})

    def write_joined(self, rows):
        with gzip.open(self.joined, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=META_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def test_complete_and_partial_retention_conservation(self):
        full = extract(self.config, self.root / "full")
        self.assertTrue(full["complete"])
        self.assertEqual(full["counts"]["primary_rows"], 2)
        self.assertEqual(full["counts"]["selected_rows"], 1)
        self.assertEqual(full["expected_selected_from_coverage"], 1)
        with gzip.open(self.root / "full/molecule_heterogeneity.tsv.gz", "rt") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([item["read_id"] for item in rows], ["first"])
        with (self.root / "full/retention.tsv").open() as handle:
            retention = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(sum(int(item["total_primary"]) for item in retention), 2)
        self.assertEqual(sum(int(item["sequence_qualified"]) for item in retention), 1)
        partial = extract(self.config, self.root / "partial", limit_primary=1)
        self.assertFalse(partial["complete"])
        self.assertEqual(partial["counts"]["primary_rows"], 1)

    def test_reordered_or_trailing_rows_fail(self):
        self.write_joined(list(reversed(self.rows)))
        with self.assertRaisesRegex(ValueError, "identity mismatch"):
            extract(self.config, self.root / "reordered")
        self.assertFalse((self.root / "reordered").exists())
        self.write_joined(self.rows + [self.rows[1]])
        with self.assertRaisesRegex(ValueError, "rows after BAM traversal"):
            extract(self.config, self.root / "trailing")
        self.assertFalse((self.root / "trailing").exists())

    def test_old_cohort_count_must_conserve(self):
        audit = read_json(self.audit)
        audit["counts"]["footprint_analysis_eligible_with_barcode"] = 2
        write_json(self.audit, audit)
        with self.assertRaisesRegex(ValueError, "cohort not conserved"):
            extract(self.config, self.root / "bad_cohort")


if __name__ == "__main__":
    unittest.main()

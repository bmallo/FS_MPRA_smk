import csv
from array import array
import hashlib
import json
import random
import sqlite3
import tempfile
import unittest
from pathlib import Path

import pysam

from fsmpra.barcodes import assign_barcodes
from fsmpra.io import open_text, read_json
from fsmpra.processing import (OBS_FIELDS, extract, extract_interval, feature_table,
                               methylation, _methylation_reference, partition, preflight, prepare)
from fsmpra.ledger import extract_complete, join_ledger


def table(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ProcessingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.wt = "ACGTCAGTACGT"
        self.snv = "TCGTCAGTACGT"
        self.ref = "G" * 8 + self.wt + "C" * 12 + "N" * 15 + "G" * 8
        self.fa = self.root / "reference.fa"
        self.fa.write_text(">plasmid\n" + self.ref + "\n")
        pysam.faidx(str(self.fa))
        (self.root / "design.txt").write_text("wt\t" + self.wt + "\nsnv\t" + self.snv + "\n")
        self.cfg = read_json("config/ldlr_processing.json")
        self.cfg.update(project_root=str(self.root), scratch_dir=str(self.root))
        self.cfg["sample"] = {"id": "sample", "library_id": "library"}
        self.cfg["library"] = dict(self.cfg["library"], id="library", reference_fasta=str(self.fa), contig="plasmid",
            reference_sequence_sha256=hashlib.sha256(self.ref.encode()).hexdigest(),
            design_oligos="design.txt", design_trim_left=0, design_trim_right=0,
            mutable_ref_start0=8, mutable_ref_end0=20, oligo_mutable_start1=1,
            barcode_start0=32, barcode_end0=47, expected_unique_classes={"WT": 1, "SNV": 1})
        self.cfg["alignment"] = {"reference_fasta": str(self.fa)}
        self.cfg["qc"]["end_slop"] = 2
        self.path = self.root / "config.json"
        self.save()
        self.catalog = self.root / "catalog"
        prepare(self.path, self.catalog)
        self.header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "plasmid", "LN": len(self.ref)}]}

    def save(self):
        self.path.write_text(json.dumps(self.cfg))

    def tearDown(self):
        self.tmp.cleanup()

    def read(self, name, bc="ACGTACGTACGTACG", snv=False, reverse=False):
        read = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(self.header))
        read.query_name = name
        read.query_sequence = self.ref[:8] + (self.snv if snv else self.wt) + self.ref[20:32] + bc + self.ref[47:]
        read.query_qualities = [40] * len(read.query_sequence)
        read.flag = 16 if reverse else 0
        read.reference_id = 0
        read.reference_start = 0
        read.mapping_quality = 60
        read.cigarstring = str(len(read.query_sequence)) + "M"
        read.set_tag("np", 20)
        read.set_tag("rq", 0.999)
        read.set_tag("MM", "A+a.;T-a.;")
        read.set_tag("ML", array("B"))
        return read

    def write_bam(self, reads):
        path = self.root / "aligned.bam"
        with pysam.AlignmentFile(path, "wb", header=self.header) as bam:
            for read in reads:
                bam.write(read)
        return path

    def test_interval_keeps_indels_and_reverse_orientation(self):
        read = self.read("reverse", reverse=True)
        self.assertEqual(extract_interval(read, 32, 47, self.ref)[0], "ACGTACGTACGTACG")
        seq = read.query_sequence
        read.query_sequence = seq[:35] + "T" + seq[35:]
        read.query_qualities = [40] * len(read.query_sequence)
        read.cigarstring = "35M1I20M"
        barcode, qualities, status = extract_interval(read, 32, 47, self.ref)
        self.assertEqual(status, "ok")
        self.assertEqual(len(barcode), 16)
        self.assertEqual(barcode, "ACGTTACGTACGTACG")

    def test_cigar_anchor_walk_matches_full_alignment_projection(self):
        rng = random.Random(17)
        for trial in range(200):
            read = self.read(str(trial), reverse=bool(trial % 2))
            cigar, sequence, pos = [(5, 2), (4, 3)], "AAA", 0
            while pos < len(self.ref):
                length = min(rng.randint(1, 12), len(self.ref) - pos)
                op = rng.choice([0, 0, 0, 7, 8, 2, 3])
                cigar.append((op, length))
                if op in (0, 7, 8):
                    sequence += self.ref[pos:pos+length]
                pos += length
                if rng.random() < 0.25:
                    cigar.append((1, 1))
                    sequence += "T"
            cigar.append((4, 2))
            sequence += "CC"
            if trial % 3 == 0:
                i = rng.randrange(len(sequence))
                sequence = sequence[:i] + "T" + sequence[i+1:]
            read.query_sequence = sequence
            read.query_qualities = [40] * len(sequence)
            read.cigartuples = cigar
            aligned = {r: q for q, r in read.get_aligned_pairs(matches_only=True)}
            for start, end in ((8, 20), (32, 47)):
                left = [aligned.get(r) for r in range(start-4, start)]
                right = [aligned.get(r) for r in range(end, end+4)]
                if None in left + right:
                    expected = ("", [], "incomplete_flanks")
                elif (left != list(range(left[0], left[0]+4)) or
                      right != list(range(right[0], right[0]+4)) or
                      sequence[left[0]:left[-1]+1] != self.ref[start-4:start] or
                      sequence[right[0]:right[-1]+1] != self.ref[end:end+4]):
                    expected = ("", [], "flank_disagreement")
                else:
                    seq = sequence[left[-1]+1:right[0]]
                    expected = (seq, [40]*len(seq), "ok")
                self.assertEqual(extract_interval(read, start, end, self.ref), expected)

    def test_missing_modifications_are_not_zero(self):
        read = self.read("zero")
        self.assertEqual(methylation(read, 128)[-1], 0)
        read.set_tag("MM", "A+a?;T-a?;")
        self.assertIsNone(methylation(read, 128)[-1])
        read.set_tag("MM", None)
        self.assertEqual(methylation(read, 128)[0], "missing")

    def test_mixed_modification_groups_and_reverse_reads(self):
        for reverse in (False, True):
            read = self.read("mixed", reverse=reverse)
            read.set_tag("MM", "A+a.,0;C+mh?,0;T-a.,0;")
            read.set_tag("ML", array("B", [255, 5, 6, 255]))
            self.assertEqual(methylation(read, 128)[2], 2)
            read.set_tag("ML", array("B", [255]))
            with self.assertRaisesRegex(ValueError, "length disagreement"):
                methylation(read, 128)

    def test_fast_methylation_matches_sparse_multimod_reference(self):
        rng = random.Random(91)
        for i in range(80):
            read = self.read(str(i), reverse=bool(i % 2))
            forward = read.get_forward_sequence()
            headers = ["A+a.", "C+mh?", "T-a."] + (["A+a?"] if i % 3 == 0 else [])
            sections, probabilities = [], []
            for header in headers:
                size = forward.count(header[0])
                selected = sorted(rng.sample(range(size), rng.randrange(size+1)))
                last, deltas = -1, []
                for pos in selected:
                    deltas.append(pos-last-1)
                    last = pos
                    probabilities.extend(rng.choice([0, 127, 128, 255]) for _ in range(2 if header == "C+mh?" else 1))
                sections.append(header + ("," + ",".join(map(str, deltas)) if deltas else ""))
            read.set_tag("MM", ";".join(sections) + ";")
            read.set_tag("ML", array("B", probabilities))
            for threshold in (0, 128, 255):
                self.assertEqual(methylation(read, threshold), _methylation_reference(read, threshold))
        read.set_tag("MM", "A+a.,100000;")
        read.set_tag("ML", array("B", [255]))
        for decoder in (methylation, _methylation_reference):
            with self.assertRaisesRegex(ValueError, "length disagreement"):
                decoder(read, 128)

    def test_chunked_complete_ledger_preserves_rows_and_checks_global_ids(self):
        reads = [self.read("a"), self.read("b"), self.read("c", "T"*15, snv=True)]
        bam = self.write_bam(reads)
        pysam.index(str(bam))
        extracted = self.root / "chunked"
        extract_complete(self.path, self.catalog / "catalog.json", bam, extracted,
                         workers=1, allocated_cpus=2, chunk_reads=1)
        self.assertEqual([r["read_id"] for r in table(extracted / "observations.tsv.gz")], ["a", "b", "c"])
        bars = self.root / "bars"
        assign_barcodes(self.path, self.catalog / "catalog.json", extracted / "observations.tsv.gz", bars)
        result = join_ledger(self.path, self.catalog / "catalog.json", extracted, bars, self.root / "ledger")
        self.assertEqual(result["molecules"], 3)
        self.assertEqual(result["assigned_molecules"], 2)
        self.assertEqual(len(table(self.root / "ledger/molecule_ledger.tsv.gz")), 3)
        # The same ID in different chunks must fail the global barcode load.
        duplicate = self.write_bam([self.read("same"), self.read("same")])
        pysam.index(str(duplicate))
        extract_complete(self.path, self.catalog / "catalog.json", duplicate, self.root / "duplicates",
                         workers=1, allocated_cpus=2, chunk_reads=1)
        with self.assertRaises(sqlite3.IntegrityError):
            assign_barcodes(self.path, self.catalog / "catalog.json", self.root / "duplicates/observations.tsv.gz", self.root / "rejected")
        self.assertFalse((self.root / "rejected").exists())

    def test_barcode_collision_diversity_and_cohort_coverage(self):
        reads = [self.read("wt1"), self.read("wt2", reverse=True),
                 self.read("collision_wt", "C" * 15), self.read("collision_snv", "C" * 15, snv=True),
                 self.read("singleton", "G" * 15)]
        # Only one of two WT observations is in the downstream analysis cohort.
        reads[1].set_tag("np", 1)
        extracted = self.root / "extracted"
        extract(self.path, self.catalog / "catalog.json", self.write_bam(reads), extracted)
        # Low-pass read cannot support barcode identity; add another good WT.
        reads.append(self.read("wt3"))
        extracted = self.root / "extracted2"
        extract(self.path, self.catalog / "catalog.json", self.write_bam(reads), extracted)
        output = self.root / "barcodes"
        assign_barcodes(self.path, self.catalog / "catalog.json", extracted / "observations.tsv.gz", output)
        mapping = {r["observed_barcode"]: r for r in table(output / "barcode_map.tsv.gz")}
        self.assertEqual(mapping["C" * 15]["status"], "collision")
        self.assertEqual(mapping["G" * 15]["status"], "insufficient_support")
        coverage = [r for r in table(output / "variant_barcode_coverage.tsv") if r["source_variant"] == "WT"]
        self.assertEqual(coverage[0]["assigned_molecules"], "3")
        self.assertEqual(coverage[1]["assigned_molecules"], "2")
        with open_text(output / "mprasnakeflow_assignments.tsv.gz") as handle:
            exported = handle.read().strip().splitlines()
        self.assertEqual(len(exported), 1)
        self.assertTrue(exported[0].startswith("ACGTACGTACGTACG\tlibrary:"))
        self.assertEqual(len(table(output / "molecule_assignments.tsv.gz")), len(reads))

    def test_correction_does_not_merge_high_quality_neighbours(self):
        self.cfg["barcodes"]["method"] = "quality_directional"
        self.save()
        parent = "ACGTACGTACGTACG"
        child = "TCGTACGTACGTACG"
        reads = [self.read("seed" + str(i), parent) for i in range(6)]
        error = self.read("error", child)
        qualities = list(error.query_qualities)
        qualities[32] = 10
        error.query_qualities = qualities
        reads.append(error)
        reads.extend([self.read("genuine1", "GCGTACGTACGTACG"), self.read("genuine2", "GCGTACGTACGTACG")])
        extracted = self.root / "extracted"
        extract(self.path, self.catalog / "catalog.json", self.write_bam(reads), extracted)
        output = self.root / "barcodes"
        assign_barcodes(self.path, self.catalog / "catalog.json", extracted / "observations.tsv.gz", output)
        mapping = {r["observed_barcode"]: r for r in table(output / "barcode_map.tsv.gz")}
        self.assertEqual(mapping[child]["corrected_barcode"], parent)
        self.assertEqual(mapping[child]["status"], "corrected_low_quality_substitution")
        self.assertEqual(mapping["GCGTACGTACGTACG"]["corrected_barcode"], "GCGTACGTACGTACG")

    def test_duplicate_ids_fail_atomically(self):
        output = self.root / "extracted"
        with self.assertRaisesRegex(ValueError, "Duplicate primary"):
            extract(self.path, self.catalog / "catalog.json", self.write_bam([self.read("same"), self.read("same")]), output)
        self.assertFalse(output.exists())

    def test_partition_is_independent_of_analysis_quality(self):
        low_pass = self.read("low_pass")
        low_pass.set_tag("np", 1)
        missing = self.read("missing")
        missing.set_tag("MM", None)
        path = self.write_bam([low_pass, missing])
        pysam.index(str(path))
        out = self.root / "partition"
        report = partition(self.path, self.catalog / "catalog.json", path, out)
        self.assertEqual(report["counts"]["circular"], 1)
        self.assertEqual(report["counts"]["uncallable"], 1)
        with pysam.AlignmentFile(out / "circular.bam", "rb") as bam:
            self.assertEqual(next(bam).query_name, "low_pass")

    def test_feature_join_preserves_empty_vs_missing_calls(self):
        reads = [self.read("empty"), self.read("missing")]
        extracted = self.root / "extracted"
        extract(self.path, self.catalog / "catalog.json", self.write_bam(reads), extracted)
        barcodes = self.root / "barcodes"
        assign_barcodes(self.path, self.catalog / "catalog.json", extracted / "observations.tsv.gz", barcodes)
        fp = self.root / "fp"
        fp.mkdir()
        reads[0].set_tag("MA", str(len(self.ref)) + ";msp:")
        with pysam.AlignmentFile(fp / "footprints.bam", "wb", header=self.header) as bam:
            for r in reads:
                bam.write(r)
        out = self.root / "features"
        feature_table(self.path, self.catalog / "catalog.json", extracted, barcodes, fp, out)
        rows = table(out / "molecule_features.tsv.gz")
        self.assertEqual(rows[0]["msp_fraction"], "0.0")
        self.assertEqual(rows[0]["tf_small_ref_intervals0"], "[]")
        self.assertEqual(rows[1]["msp_fraction"], "")
        self.assertEqual(rows[1]["footprint_status"], "missing")

    def test_preflight_missing_calls_fails(self):
        read = self.read("unaligned")
        read.flag = 4
        read.reference_id = -1
        read.reference_start = -1
        read.cigartuples = None
        read.set_tag("MM", None)
        self.cfg["input"]["bam"] = str(self.write_bam([read]))
        self.save()
        with self.assertRaisesRegex(ValueError, "Missing m6A"):
            preflight(self.path, self.root / "preflight")


if __name__ == "__main__":
    unittest.main()

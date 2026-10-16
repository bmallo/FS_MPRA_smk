import csv
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pysam

from fsmpra.catalog import build_catalog, identify_single_edit, source_variant_sequence
from fsmpra.dashboard import export_dashboard
from fsmpra.features import interval_runs, measure, parse_ma, projected_coverage
from fsmpra.io import read_json, write_table
from fsmpra.molecules import LEDGER_COLUMNS, METRIC_COLUMNS, import_molecules, selected


class CoordinateTests(unittest.TestCase):
    def test_equivalent_homopolymer_deletions(self):
        edit = identify_single_edit("CAAAG", "CAAG", 21)
        self.assertEqual(edit["source_variant"], "22delA")
        self.assertEqual(edit["equivalent_positions"], [22, 23, 24])
        self.assertEqual(source_variant_sequence("CAAAG", "24delA", 21), "CAAG")

    def test_equivalent_insertions_and_boundary(self):
        edit = identify_single_edit("CAAAG", "CAAAAG", 21)
        self.assertEqual(edit["source_variant"], "21_22insA")
        self.assertEqual(edit["equivalent_positions"], [22, 23, 24, 25])
        self.assertEqual(source_variant_sequence("CAAAG", "24_25insA", 21), "CAAAAG")
        self.assertEqual(identify_single_edit("ACG", "TACG", 1)["source_variant"], "0_1insT")

    def test_reject_wrong_reference_and_complex_design(self):
        with self.assertRaises(ValueError):
            source_variant_sequence("ACG", "1T>G", 1)
        with self.assertRaises(ValueError):
            identify_single_edit("ACG", "TTG", 1)

    def test_projection_does_not_fill_deletion_or_double_count(self):
        qmap = [10, 11, None, 14, 15]
        coverage, lengths = projected_coverage([(0, 5), (0, 2)], qmap, 10, 16)
        self.assertEqual(coverage, {0, 1, 4, 5})
        self.assertEqual(interval_runs(coverage), [[0, 2], [4, 6]])
        self.assertEqual(lengths, [5, 2])

    def test_empty_track_differs_from_missing(self):
        tracks = parse_ma("10;nuc+Q:1-3;tf+QQQ:", 10)
        self.assertEqual(tracks["tf"], [])
        self.assertNotIn("msp", tracks)
        self.assertEqual(tracks["nuc"], [(0, 3)])
        for bad in ["9;tf:", "10;tf:0-2", "10;tf:9-3", "10;tf:;tf:"]:
            with self.assertRaises(ValueError):
                parse_ma(bad, 10)

    def test_sampling_is_order_independent_and_nested(self):
        ids = [f"read-{i}" for i in range(1000)]
        small = {r for r in ids if selected(r, .1, 42)}
        reverse = {r for r in reversed(ids) if selected(r, .1, 42)}
        larger = {r for r in ids if selected(r, .2, 42)}
        self.assertEqual(small, reverse)
        self.assertTrue(small < larger)
        self.assertTrue(all(selected(r, 1, 42) for r in ids))


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.wt = "ACGTACGTACGT"
        self.ref = "CCCC" + self.wt + "GGGG"
        (self.root / "ref.fa").write_text(">construct\n" + self.ref + "\n")
        variants = [self.wt, "TCGTACGTACGT", "AGGTACGTACGT", self.wt[1:], "T" + self.wt]
        (self.root / "design.txt").write_text("".join(f"d{i}\t{s}\n" for i, s in enumerate(variants)))
        mpra_fields = ["Chromosome", "Position", "Ref", "Alt", "Value", "P-Value", "Element", "Tags", "DNA", "RNA"]
        rows = [dict(zip(mpra_fields, ["1", 100, "A", "T", 1, .2, "external", 5, 20, 40])),
                dict(zip(mpra_fields, ["1", 101, "C", "G", -1, .1, "external", 5, 20, 10]))]
        write_table(self.root / "mpra.tsv", rows, mpra_fields)
        # Both stored sequences align in reference orientation, including reverse read.
        bam_path = self.root / "input.bam"
        with pysam.AlignmentFile(bam_path, "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "construct", "LN": 20}]}) as bam:
            for name, seq, reverse, ma in [
                ("wt_reverse", self.ref, True, "20;tf:;nuc:;msp:5-12"),
                ("snv", "CCCCTCGTACGTACGTGGGG", False, "20;tf:5-4;nuc:;msp:9-8"),
                ("missing_tf", "CCCCAGGTACGTACGTGGGG", False, "20;nuc:;msp:"),
                ("missing_metrics", self.ref, False, "20;tf:;msp:")]:
                read = pysam.AlignedSegment(bam.header)
                read.query_name = name
                read.query_sequence = seq
                read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                read.flag = 16 if reverse else 0
                read.reference_id = 0
                read.reference_start = 0
                read.mapping_quality = 60
                read.cigarstring = "20M"
                read.set_tag("MA", ma)
                bam.write(read)
        pysam.index(str(bam_path))
        ledger = []
        for name, cls, label, mf in [("wt_reverse", "WT", "WT", .95), ("snv", "SNV", "1A>T", .01),
                                     ("missing_tf", "SNV", "2C>G", .4), ("missing_metrics", "WT", "WT", .3)]:
            ledger.append(dict(zip(LEDGER_COLUMNS, [name, 20, 20, .9999, "-" if name == "wt_reverse" else "+", 0, 20,
                                                   "full_length", cls, label, "design", 0 if cls == "WT" else 1, 0, "strict", "ACGT"])))
        write_table(self.root / "ledger.tsv.gz", ledger, LEDGER_COLUMNS)
        metrics = [dict(zip(METRIC_COLUMNS, [name, "construct", "plasmid_full", mf, 0, 0, "old_filter_fail", 8, 3]))
                   for name, mf in [("wt_reverse", .95), ("snv", .01), ("missing_tf", .4)]]
        write_table(self.root / "metrics.tsv.gz", metrics, METRIC_COLUMNS)
        manifest = {"schema_version": 1, "project_root": str(self.root),
                    "libraries": [{"id": "lib", "reference_fasta": "ref.fa", "contig": "construct",
                                   "reference_sequence_sha256": hashlib.sha256(self.ref.encode()).hexdigest(),
                                   "mutable_ref_start0": 4, "mutable_ref_end0": 16, "oligo_mutable_start1": 1,
                                   "genomic_chrom": "chr1", "genomic_mutable_start1": 100, "genomic_strand": "+",
                                   "design_oligos": "design.txt", "design_trim_left": 0, "design_trim_right": 0,
                                   "expected_unique_classes": {"WT": 1, "SNV": 2, "DEL": 1, "INS": 1},
                                   "analysis_interval_status": "synthetic"}],
                    "samples": [{"id": "rep1", "library_id": "lib", "biological_replicate": "1",
                                 "inputs": {"ledger": "ledger.tsv.gz", "metrics": "metrics.tsv.gz", "footprints_bam": "input.bam"}}],
                    "mpra_sources": [{"id": "external", "path": "mpra.tsv", "delimiter": "\t"}]}
        self.manifest = self.root / "manifest.json"
        self.manifest.write_text(json.dumps(manifest))

    def tearDown(self):
        self.temp.cleanup()

    def test_end_to_end_preserves_methylation_extremes_and_missingness(self):
        catalog_dir = self.root / "catalog"
        result = build_catalog(self.manifest, catalog_dir)
        self.assertEqual(result["mpra_audit"][0]["matched_snvs"], 2)
        cat = catalog_dir / "catalog.json"
        mol = self.root / "molecules"
        audit = import_molecules(cat, "rep1", mol, fraction=1)
        self.assertEqual(audit["selected_reads"], 4)
        self.assertEqual(audit["missing_metrics"], 1)
        output = self.root / "features"
        measure(cat, mol, output, tf_min=1, tf_max=5, mask_radius=0)
        pilot = read_json(output / "pilot.json")
        self.assertEqual(pilot["counts"]["measured"], 3)
        v = {r["source_variant"]: r for r in pilot["variants"]}
        self.assertEqual(v["WT"]["tracks"]["msp"]["occupancy"], [1.] * 12)
        self.assertEqual(v["WT"]["n_reverse"], 1)
        self.assertEqual(v["1A>T"]["tracks"]["tf_small"]["occupancy"], [1.] * 4 + [0.] * 8)
        self.assertEqual(v["2C>G"]["tracks"]["tf_small"]["occupancy"], [None] * 12)
        self.assertEqual(v["2C>G"]["tracks"]["msp"]["occupancy"], [0.] * 12)
        self.assertAlmostEqual(v["1A>T"]["summary"]["msp_common_site_delta"], 8 / 11 - 1)
        self.assertEqual(v["WT"]["methylation_histogram"][19], 1)
        self.assertEqual(v["1A>T"]["methylation_histogram"][0], 1)
        self.assertEqual(len(pilot["mpra"]), 2)
        # Embedded provenance must round-trip without being able to terminate
        # the inert JSON script, and the viewer must not need external assets.
        pilot["sample"]["note"] = '</script><script>alert("test")</script>&'
        bundle = output / "pilot.json"
        bundle.write_text(json.dumps(pilot))
        result = export_dashboard(bundle, self.root / "dashboard")
        html = Path(result["html"]).read_text()
        marker = '<script id="pilot-data" type="application/json">'
        embedded = html.split(marker)[1].split('</script>')[0]
        self.assertEqual(json.loads(embedded), pilot)
        self.assertNotIn(pilot["sample"]["note"], html)
        self.assertNotIn('src="http', html)
        self.assertNotIn('PLACEHOLDER', html)
        with self.assertRaises(FileExistsError):
            export_dashboard(bundle, self.root / "dashboard")

    def test_reference_checksum_and_mpra_allele_mismatch_fail(self):
        (self.root / "ref.fa").write_text(">construct\n" + "A" * 20 + "\n")
        with self.assertRaisesRegex(ValueError, "checksum"):
            build_catalog(self.manifest, self.root / "bad")
        (self.root / "ref.fa").write_text(">construct\n" + self.ref + "\n")
        p = self.root / "mpra.tsv"
        p.write_text(p.read_text().replace("100\tA\tT", "100\tC\tT"))
        with self.assertRaisesRegex(ValueError, "reference mismatch"):
            build_catalog(self.manifest, self.root / "bad")

    def test_sparse_writer_omissions_are_zero_only_for_evaluated_tracks(self):
        manifest = read_json(self.manifest)
        manifest["samples"][0]["footprint_contract"] = {
            "format": "fiberhmm_ma_sparse", "evaluated_tracks": ["tf", "msp", "nuc"]}
        self.manifest.write_text(json.dumps(manifest))
        cat_dir = self.root / "catalog"
        build_catalog(self.manifest, cat_dir)
        cat = cat_dir / "catalog.json"
        mol = self.root / "molecules"
        import_molecules(cat, "rep1", mol, fraction=1)
        output = self.root / "features"
        measure(cat, mol, output, tf_min=1, tf_max=5)
        pilot = read_json(output / "pilot.json")
        v = next(v for v in pilot["variants"] if v["source_variant"] == "2C>G")
        self.assertEqual(v["tracks"]["tf_small"]["occupancy"], [0.] * 12)
        self.assertEqual(pilot["counts"]["omitted_empty_track_tf"], 1)

    def test_duplicate_metrics_and_sampling_cap_fail(self):
        cat_dir = self.root / "catalog"
        build_catalog(self.manifest, cat_dir)
        cat = cat_dir / "catalog.json"
        with self.assertRaisesRegex(ValueError, "cap exceeded"):
            import_molecules(cat, "rep1", self.root / "small", fraction=1, max_molecules=1)
        with gzip.open(self.root / "metrics.tsv.gz", "rt") as f:
            lines = f.readlines()
        with gzip.open(self.root / "metrics.tsv.gz", "wt") as f:
            f.writelines(lines + [lines[1]])
        with self.assertRaisesRegex(ValueError, "changed since"):
            import_molecules(cat, "rep1", self.root / "stale", fraction=1)
        fresh = self.root / "fresh"
        build_catalog(self.manifest, fresh)
        with self.assertRaisesRegex(ValueError, "Duplicate selected metric"):
            import_molecules(fresh / "catalog.json", "rep1", self.root / "dup", fraction=1)


if __name__ == "__main__":
    unittest.main()

import csv
import gzip
import sqlite3
import unittest

import pysam

from fsmpra.catalog import build_catalog
from fsmpra.compare import recovery_evidence
from fsmpra.features import measure
from fsmpra.io import read_json, write_table
from fsmpra.molecules import LEDGER_COLUMNS, METRIC_COLUMNS, import_molecules
import test_pilot


class DiskTests(unittest.TestCase):
    def setUp(self):
        self.fixture = test_pilot.PipelineTests()
        self.fixture.setUp()
        self.root = self.fixture.root
        self.scratch = self.root / "scratch"
        self.scratch.mkdir()

    def tearDown(self):
        self.fixture.tearDown()

    def catalog(self, name="catalog"):
        build_catalog(self.fixture.manifest, self.root / name)
        return self.root / name / "catalog.json"

    def rows(self, path):
        with gzip.open(path, "rt") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_full_and_sampled_import_and_measurement_match_memory_backend(self):
        catalog = self.catalog()
        for i, fraction in enumerate([1, .5]):
            with self.subTest(fraction=fraction):
                memory, disk = self.root / f"memory{i}", self.root / f"disk{i}"
                a = import_molecules(catalog, "rep1", memory, fraction=fraction)
                b = import_molecules(catalog, "rep1", disk, fraction=fraction, backend="sqlite", batch_size=1, scratch_dir=self.scratch)
                self.assertEqual(self.rows(memory / "molecules.tsv.gz"), self.rows(disk / "molecules.tsv.gz"))
                for key in ("selected_reads", "missing_metrics", "ledger_classes", "selected_classes", "ledger_rows", "metrics_rows"):
                    self.assertEqual(a[key], b[key])
                qa = read_json(disk / "recovery.json")
                self.assertEqual(sum(g["molecules"] for g in qa["by_plasmid_class"].values()), b["selected_reads"])
                outputs = [self.root / f"features{i}{k}" for k in ("memory", "sqlite")]
                for src, backend, out in zip((memory, disk), ("memory", "sqlite"), outputs):
                    measure(catalog, src, out, tf_min=1, tf_max=5, backend=backend, scratch_dir=self.scratch)
                x, y = [read_json(p / "pilot.json") for p in outputs]
                self.assertIsNone(recovery_evidence(x))
                self.assertEqual(recovery_evidence(y)["summary"]["sample_id"], "rep1")
                for key in ("counts", "variants", "method", "region", "mpra"):
                    self.assertEqual(x[key], y[key])
                self.assertEqual(self.rows(outputs[0] / "molecule_features.tsv.gz"), self.rows(outputs[1] / "molecule_features.tsv.gz"))
        self.assertEqual(list(self.scratch.iterdir()), [])

    def test_duplicate_metrics_across_and_within_batches_fail_atomically(self):
        path = self.root / "metrics.tsv.gz"
        rows = self.rows(path)
        write_table(path, rows + [rows[0]], METRIC_COLUMNS)
        catalog = self.catalog()
        for batch_size in (1, 100):
            out = self.root / f"bad{batch_size}"
            with self.assertRaisesRegex(ValueError, "Duplicate selected metric"):
                import_molecules(catalog, "rep1", out, fraction=1, backend="sqlite", batch_size=batch_size, scratch_dir=self.scratch)
            self.assertFalse(out.exists())
            self.assertEqual(list(self.scratch.iterdir()), [])

    def test_duplicate_ledger_fails_without_partial_output(self):
        path = self.root / "ledger.tsv.gz"
        rows = self.rows(path)
        write_table(path, rows + [rows[0]], LEDGER_COLUMNS)
        with self.assertRaisesRegex(ValueError, "Duplicate selected ledger"):
            import_molecules(self.catalog(), "rep1", self.root / "bad", fraction=1, backend="sqlite", batch_size=1)
        self.assertFalse((self.root / "bad").exists())

    def test_database_limit_and_existing_outputs_are_respected(self):
        catalog = self.catalog()
        with self.assertRaises(sqlite3.OperationalError):
            import_molecules(catalog, "rep1", self.root / "limited", fraction=1, backend="sqlite", max_disk_gib=1e-9, scratch_dir=self.scratch)
        self.assertFalse((self.root / "limited").exists())
        self.assertEqual(list(self.scratch.iterdir()), [])
        existing = self.root / "existing"
        existing.mkdir()
        (existing / "keep").write_text("keep")
        with self.assertRaises(FileExistsError):
            import_molecules(catalog, "rep1", existing, backend="sqlite")
        self.assertEqual((existing / "keep").read_text(), "keep")

    def test_duplicate_primary_alignment_fails_without_feature_file(self):
        path = self.root / "input.bam"
        with pysam.AlignmentFile(str(path), "rb") as bam:
            header, reads = bam.header, list(bam)
        with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
            for read in reads + [reads[0]]:
                bam.write(read)
        pysam.index(str(path))
        catalog = self.catalog()
        import_molecules(catalog, "rep1", self.root / "molecules", fraction=1, backend="sqlite")
        with self.assertRaisesRegex(ValueError, "Duplicate primary alignment"):
            measure(catalog, self.root / "molecules", self.root / "bad", backend="sqlite", scratch_dir=self.scratch)
        self.assertFalse((self.root / "bad").exists())
        self.assertEqual(list(self.scratch.iterdir()), [])

    def test_empty_selection_retains_all_designed_variants(self):
        catalog = self.catalog()
        audit = import_molecules(catalog, "rep1", self.root / "empty", fraction=1e-15, backend="sqlite")
        self.assertEqual(audit["selected_reads"], 0)
        measure(catalog, self.root / "empty", self.root / "features", backend="sqlite")
        data = read_json(self.root / "features/pilot.json")
        self.assertEqual(len(data["variants"]), 3)
        self.assertEqual(data["counts"]["measured"], 0)
        self.assertTrue(all(v["n_reads"] == 0 for v in data["variants"]))
        self.assertTrue(all(x is None for v in data["variants"] for x in v["tracks"]["msp"]["occupancy"]))


if __name__ == "__main__":
    unittest.main()

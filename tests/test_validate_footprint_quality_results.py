import csv
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from tools import validate_footprint_quality_results as validator


def _table(path, rows):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _hash(value):
    return hashlib.sha256(value.encode()).hexdigest()


class ValidationTests(unittest.TestCase):
    def test_independent_selection_arrays_and_molecule_reconciliation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample, summary, aggregates = (root / name for name in ("sample", "summary", "aggregates"))
            for path in (sample, summary, aggregates):
                path.mkdir()
            names = ("AAA", "CCC", "GGG")
            old_rows = [dict(barcode_index=i, barcode=name, variant_id="wt", molecules=1)
                        for i, name in enumerate(names)]
            _table(aggregates / "barcode_metadata.tsv.gz", old_rows)
            selected = sorted(names, key=lambda name: _hash(f"20261018|{name}"))[:2]
            bitmap = {"AAA": ([0, 1, 2], [0, 1]),
                      "CCC": ([10, 11], [1, 2]),
                      "GGG": ([], [])}
            counts = np.zeros((3, 2, 211), dtype=np.uint32)
            for i, name in enumerate(names):
                tf, msp = bitmap[name]
                counts[i, 0, tf] = 1
                counts[i, 1, msp] = 1
            np.save(aggregates / "barcode_counts.npy", counts)
            np.save(aggregates / "barcode_molecules.npy", np.ones(3, dtype=np.uint32))
            aggregate_audit = dict(schema_version="fsmpra.statistics_aggregate.v1", sample_id="s",
                                   tracks=["tf_small", "msp"], interval=[1234, 1445],
                                   counts={"barcodes": 3})
            (aggregates / "audit.json").write_text(json.dumps(aggregate_audit))
            barcodes = []
            molecules = []
            for name in sorted(selected):
                split = "A" if int(_hash(f"split|20261018|{name}"), 16) % 2 == 0 else "B"
                tf, msp = bitmap[name]
                barcodes.append(dict(barcode=name, variant_id="wt", eligible_depth=1,
                                     split=split, weight=1, selection_hash=_hash(f"20261018|{name}")))
                molecules.append(dict(corrected_barcode=name, assigned_variant_id="wt",
                                      sample_split=split,
                                      element_aligned_offsets0=json.dumps(list(range(211))),
                                      element_tf_small_offsets0=json.dumps(tf),
                                      element_msp_offsets0=json.dumps(msp),
                                      element_tf_small_fraction=len(tf)/211,
                                      element_msp_fraction=len(msp)/211))
            _table(sample / "sampled_barcodes.tsv", barcodes)
            _table(sample / "sampled_molecules.tsv.gz", molecules)
            windows = validator._windows(1234)
            coverage = np.zeros((2, 22))
            presence = np.zeros((2, 22))
            for i, name in enumerate(sorted(selected)):
                tf, msp = bitmap[name]
                for j, (_, track, lo, hi, _, _) in enumerate(windows):
                    positions = tf if track == "tf_small" else msp
                    n = sum(lo <= pos < hi for pos in positions)
                    coverage[i, j] = n / (hi-lo)
                    presence[i, j] = n > 0
            matrix_path = summary / "barcode_endpoint_matrix.npz"
            np.savez_compressed(matrix_path,
                barcodes=np.array(sorted(selected)), splits=np.array([row["split"] for row in barcodes]),
                endpoint_ids=np.array([window[0] for window in windows]),
                relative_starts=np.array([window[2] for window in windows]),
                relative_ends=np.array([window[3] for window in windows]),
                reference_starts0=np.array([window[4] for window in windows]),
                reference_ends0=np.array([window[5] for window in windows]),
                depths=np.ones(2, dtype=np.int32), valid_molecules=np.ones((2, 22), dtype=np.int32),
                presence=presence, coverage=coverage)
            audit_path = sample / "audit.json"
            audit = dict(schema_version="fsmpra.footprint_quality.v1", complete=True,
                         seed=20261018, element=[1234, 1445],
                         counts=dict(sampled_barcodes=2, wt_barcodes=3,
                                     wt_eligible_molecules=3, sampled_expected_molecules=2,
                                     matched_molecules=2),
                         sources=[validator._source(aggregates / "barcode_metadata.tsv.gz", True)])
            audit_path.write_text(json.dumps(audit))
            report = dict(schema_version="fsmpra.footprint_quality_summary.v1", complete=True,
                          pilot=False, sample_id="s", sampled_barcodes=2,
                          sampled_molecules=2, endpoint_count=22, endpoint_matrix=matrix_path.name,
                          split_barcodes=dict(__import__("collections").Counter(row["split"] for row in barcodes)),
                          sources=[validator._source(path, checksum) for path, checksum in
                                   ((audit_path, True), (sample / "sampled_barcodes.tsv", True),
                                    (sample / "sampled_molecules.tsv.gz", False), (matrix_path, True))])
            (summary / "summary.json").write_text(json.dumps(report))
            with patch.object(validator, "SAMPLE_N", 2):
                result = validator.validate(sample, summary, aggregates)
            self.assertEqual(result["status"], "passed")
            self.assertEqual(result["sampled_molecules"], 2)
            self.assertEqual(result["max_abs_coverage_error_vs_old_aggregate"], 0)
            self.assertEqual(result["max_abs_coverage_error_vs_molecule_bitmaps"], 0)
            self.assertEqual(result["max_abs_presence_error_vs_molecule_bitmaps"], 0)
            with patch.object(validator, "SAMPLE_N", 2):
                validator._write_atomic(root / "validation.json", result)
            self.assertTrue((root / "validation.json").exists())


if __name__ == "__main__":
    unittest.main()

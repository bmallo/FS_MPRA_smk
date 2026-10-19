import array
from pathlib import Path
import tempfile
import unittest

from tools.audit_rep3_caller_provenance import (
    bam_identity, check_bam_identity, check_small_file, ma_nuc_aq,
    ml_groups, small_file, write_atomic_new,
)


class FakeRead:
    def __init__(self, **tags):
        self.tags = tags

    def has_tag(self, key):
        return key in self.tags

    def get_tag(self, key):
        return self.tags[key]


class ProvenanceParserTests(unittest.TestCase):
    def test_mixed_mm_channels_consume_ml_in_group_order(self):
        read = FakeRead(MM="C+mh?,0,1;A+a.,0;T-a.,0;",
                        ML=array.array("B", [7, 8, 9, 10, 240, 255]))
        target, all_scores, groups = ml_groups(read)
        self.assertEqual(target, [240, 255])
        self.assertEqual(all_scores, [7, 8, 9, 10, 240, 255])
        self.assertEqual([x["calls"] for x in groups], [2, 1, 1])

    def test_ma_aq_cursor_skips_unscored_msp(self):
        read = FakeRead(MA="100;nuc+Q:1-90,91-10;msp+:5-15;tf+QQQ:30-12",
                        AQ=array.array("B", [0, 0, 101, 202, 203]))
        nuc, labels = ma_nuc_aq(read)
        self.assertEqual(nuc, [0, 0])
        self.assertEqual(dict(labels), {"nuc": 2, "msp": 1, "tf": 1})

    def test_source_mutation_is_detected(self):
        with tempfile.TemporaryDirectory() as dirname:
            path = Path(dirname) / "source.txt"
            path.write_text("one")
            bam_before = bam_identity(path)
            small_before = small_file(path)
            path.write_text("two more")
            with self.assertRaisesRegex(RuntimeError, "changed"):
                check_bam_identity(path, bam_before)
            with self.assertRaisesRegex(RuntimeError, "changed"):
                check_small_file(path, small_before)

    def test_atomic_output_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as dirname:
            path = Path(dirname) / "audit.json"
            write_atomic_new(path, '{"complete": true}\n')
            with self.assertRaises(FileExistsError):
                write_atomic_new(path, '{"complete": false}\n')
            self.assertEqual(path.read_text(), '{"complete": true}\n')


if __name__ == "__main__":
    unittest.main()

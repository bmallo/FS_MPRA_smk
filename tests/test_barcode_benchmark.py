import itertools
import os
import tempfile
import unittest
from pathlib import Path

from fsmpra.barcode_benchmark import (edit_distance, genotype_tokens, hash_fraction,
    nearby_pairs, one_edit_neighborhood, run_pacybara, split_holdout, summarize, synthetic_cases)


LIBRARY = {"id": "test", "mutable_sequence": "ACGT", "barcode_length": 15,
    "variants": [{"variant_id": "wt", "class": "WT", "sequence": "ACGT"},
                 {"variant_id": "snv", "class": "SNV", "sequence": "TCGT"}]}


class BarcodeBenchmarkTests(unittest.TestCase):
    def test_edit_candidates_match_exhaustive_including_indels(self):
        sequences = ["".join(x) for n in range(1, 5) for x in itertools.product("AC", repeat=n)]
        for radius in (1, 2):
            expected = {(a, b, edit_distance(a, b)) for a, b in itertools.combinations(sorted(sequences), 2)
                        if edit_distance(a, b) <= radius}
            self.assertEqual(set(nearby_pairs(sequences, radius)), expected)

    def test_sampling_is_reproducible_and_seeded(self):
        self.assertEqual(hash_fraction("read-1", "a"), hash_fraction("read-1", "a"))
        self.assertNotEqual(hash_fraction("read-1", "a"), hash_fraction("read-1", "b"))
        self.assertTrue(0 <= hash_fraction("read-1", "a") < 1)

    def test_holdout_disagreement_not_gated_by_test_genotype(self):
        rows = [{"read_id": "held", "barcode": "AAAA", "variant_id": "snv"}]
        result = summarize([], {"method": ({}, {"AAAA": "wt"})}, heldout=rows)[0]
        self.assertEqual(result["heldout_variant_disagreements"], 1)

    def test_genotype_tokens(self):
        self.assertEqual(genotype_tokens(LIBRARY), {"wt": "=", "snv": "A1T"})

    def test_neighborhood_and_discovery_exclusion(self):
        neighbors = one_edit_neighborhood("ACGT")
        self.assertIn("ACG", neighbors)
        self.assertIn("AACGT", neighbors)
        self.assertIn("ATGT", neighbors)
        self.assertTrue(all(edit_distance("ACGT", bc) <= 1 for bc in neighbors))
        rows = [{"read_id": str(i)} for i in range(100)]
        _, original_test = split_holdout(rows)
        train, test = split_holdout(rows, {r["read_id"] for r in original_test})
        self.assertTrue(original_test)
        self.assertFalse(test)
        self.assertEqual(len(train), 100)

    @unittest.skipUnless(os.environ.get("PACYBARA_CHECKOUT"), "optional pinned upstream integration")
    def test_actual_core_preserves_ids_and_exposes_neighbor_merges(self):
        rows, truth = synthetic_cases(LIBRARY, repeats=2)
        with tempfile.TemporaryDirectory() as tmp:
            methods = {"pacybara": run_pacybara(rows, LIBRARY, os.environ["PACYBARA_CHECKOUT"], Path(tmp) / "run")}
            result = summarize(rows, methods, truth)[0]
            self.assertEqual(set(methods["pacybara"][0]), set(truth))
            self.assertGreater(result["synthetic_cases"]["genuine_neighbor_same_variant"]["incorrect_assignments"], 0)
            self.assertEqual(result["synthetic_cases"]["barcode_collision"].get("assigned_molecules", 0), 0)


if __name__ == "__main__":
    unittest.main()

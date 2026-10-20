import unittest

from fsmpra.genomic_coordinates import interval_fields, reverse_complement, verify_sequence


class CoordinatesTests(unittest.TestCase):
    def test_all_single_bases_and_whole_element_have_consistent_conventions(self):
        mapping = dict(genomic_start1=11089285, length=211)
        for i in range(211):
            row = interval_fields(mapping, i, i+1)
            self.assertEqual(row['genomic_start1'], 11089285+i)
            self.assertEqual(row['genomic_end1'], 11089285+i)
            self.assertEqual(row['genomic_end0']-row['genomic_start0'], 1)
        self.assertEqual(interval_fields(mapping, 0, 211), dict(
            genomic_start0=11089284, genomic_end0=11089495,
            genomic_start1=11089285, genomic_end1=11089495))
        self.assertEqual(interval_fields(mapping, 104, 105)['genomic_start1'], 11089389)

    def test_invalid_intervals_fail_instead_of_mapping_backbone(self):
        for a,b in [(-1,10), (0,212), (10,10), (11,10)]:
            with self.assertRaises(ValueError):
                interval_fields(dict(genomic_start1=11089285, length=211), a, b)

    def test_reversal_offset_mismatch_and_unsupported_strand_fail(self):
        sequence = 'CTTTCGAAGGACTGGA'
        verify_sequence(sequence, sequence, '+')
        for genomic, strand in [(reverse_complement(sequence), '+'),
                                (sequence[1:]+'G', '+'), ('A'+sequence[1:], '+'),
                                (sequence, '-')]:
            with self.assertRaises(ValueError):
                verify_sequence(sequence, genomic, strand)

    def test_ambiguous_orientation_is_not_claimed_as_verified(self):
        with self.assertRaises(ValueError):
            verify_sequence('ACGT', 'ACGT', '+')


if __name__ == '__main__':
    unittest.main()

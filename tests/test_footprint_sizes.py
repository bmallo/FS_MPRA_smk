import unittest

import numpy as np

from fsmpra.footprint_sizes import WIDTH, density, feature_values, find_features, midpoint_index, standardize

MAP = {'smoothing_sigma_bp': {'midpoint': 1.5, 'length': 2.0}, 'peak_neighbourhood_bp': {'midpoint': 5, 'length': 6},
       'box_midpoint_halfwidth_bp': [2, 8], 'box_length_halfwidth_bp': [3, 15], 'min_box_occupancy_C_std': 0.01,
       'max_box_overlap': 0.5, 'near_naked_ratio_flag': 2.0,
       'background_neighbourhood_bp': {'midpoint': 20, 'length': 20}}
LR = (10, 90)


class FootprintSizeTests(unittest.TestCase):
    def test_midpoint_reversal(self):
        # Element [100, 111): a call [100, 103) centres on index 1 forward and index 9 reversed.
        self.assertEqual(int(midpoint_index(np.array([100]), np.array([103]), (100, 111))[0]), 1)
        self.assertEqual(int(midpoint_index(np.array([100]), np.array([103]), (100, 111), reverse=True)[0]), 9)

    def test_two_sizes_at_one_centre_and_a_third_site(self):
        rng = np.random.default_rng(1)
        n = 20000
        groups = np.r_[np.zeros(1000, int), np.ones(1000, int), np.full(n - 2000, 2)]
        rows, mid, length = [], [], []
        for centre, size, frac in ((100, 18, .06), (100, 45, .04), (160, 25, .05)):
            who = np.flatnonzero((groups == 2) & (rng.random(n) < frac))
            rows.append(who)
            mid.append(centre + rng.integers(-1, 2, len(who)))
            length.append(size + rng.integers(-2, 3, len(who)))
        noise = rng.integers(0, n, 20000)
        rows.append(noise)
        mid.append(rng.integers(0, WIDTH, len(noise)))
        length.append(rng.integers(10, 91, len(noise)))
        rows, mid, length = (np.concatenate(x) for x in (rows, mid, length))
        dens, _ = density(rows, mid, length, groups, np.ones(n, bool), 3, LR)
        feats, _, _ = find_features(standardize(dens, np.array([1.])), dens[1], MAP, LR)
        found = sorted((f['mid_peak'], f['len_peak']) for f in feats)
        self.assertEqual(len(found), 3, found)
        for (m, l), (em, el) in zip(found, [(100, 18), (100, 45), (160, 25)]):
            self.assertLessEqual(abs(m - em), 1)
            self.assertLessEqual(abs(l - el), 2)
        vals = feature_values(rows, mid, length, n, feats)
        occ = vals['occupancy'][0][groups == 2].mean(0)
        self.assertTrue(all(o > .03 for o in occ), occ)


if __name__ == '__main__':
    unittest.main()

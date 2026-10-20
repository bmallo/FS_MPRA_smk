import unittest
import numpy as np
from tools.independent_review_20261020 import circular_shift, split_half


class IndependentReviewTests(unittest.TestCase):
    def test_split_half_recovers_signal_and_noise_scales(self):
        rng=np.random.default_rng(0);signal=rng.normal(0,2,20000)
        a,b=signal+rng.normal(0,2,20000),signal+rng.normal(0,2,20000)
        result=split_half(a,b)
        self.assertAlmostEqual(result['pearson_halves'],.5,delta=.02)
        self.assertAlmostEqual(result['reliability_full'],2/3,delta=.02)
        self.assertAlmostEqual(result['signal_sd_pp'],2,delta=.05)
        self.assertAlmostEqual(result['noise_sd_full_pp'],2/np.sqrt(2),delta=.05)

    def test_circular_shift_keeps_sites_together_and_reports_attainable_floor(self):
        rng=np.random.default_rng(1);pos=np.repeat(np.arange(21,232),3)
        profile=np.convolve(rng.normal(size=211),np.ones(9)/9,'same')
        aligned=circular_shift(pos,np.repeat(profile,3)+rng.normal(0,.01,633),np.repeat(profile,3))
        self.assertEqual(aligned['positions'],211)
        self.assertEqual(aligned['rotations'],192)
        self.assertAlmostEqual(aligned['p_circular_shift'],1/193)
        unrelated=circular_shift(pos,np.repeat(profile,3),np.repeat(np.roll(profile,100),3)+rng.normal(0,.5,633))
        self.assertGreater(unrelated['p_circular_shift'],.05)

    def test_missing_measurements_are_dropped_not_zeroed(self):
        pos=np.repeat(np.arange(21,232),3);x=np.tile([1.,2.,np.nan],211)
        result=circular_shift(pos,x+np.repeat(np.sin(np.arange(211)/5),3),np.repeat(np.sin(np.arange(211)/5),3))
        self.assertEqual(result['positions'],211)


if __name__=='__main__':unittest.main()

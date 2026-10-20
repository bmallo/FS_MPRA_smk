import unittest
from array import array
import numpy as np
import pysam
from fsmpra.raw_candidate_extract import decode_scores,project_read,sampled_read
from fsmpra.raw_candidate_summary import barcode_means,scalar,shared,methylation_outcome,resample_one_per_barcode,joint_distribution


def read(seq,mm,ml,reverse=False):
    r=pysam.AlignedSegment();r.query_name='test';r.query_sequence=seq;r.query_qualities=array('B',[40]*len(seq));r.flag=16 if reverse else 0
    r.reference_start=0;r.cigarstring=f'{len(seq)}M';r.set_tag('MM',mm);r.set_tag('ML',array('B',ml));return r


class RawCandidateTests(unittest.TestCase):
    def test_mixed_channels_and_reverse_coordinates(self):
        r=read('ACTA','C+mh.,0;A+a.,1;T-a.,0;',[100,120,250,240])
        scores,known=decode_scores(r)
        np.testing.assert_array_equal(scores,[0,0,240,250]);np.testing.assert_array_equal(known,[True,False,True,True])
        reverse=read('TAGT','C+mh.,0;A+a.,1;T-a.,0;',[100,120,250,240],True)
        rs,rk=decode_scores(reverse)
        np.testing.assert_array_equal(rs,scores[::-1]);np.testing.assert_array_equal(rk,known[::-1])

    def test_unknown_and_absent_channels_stay_missing(self):
        scores,known=decode_scores(read('AATA','A+a?,1;',[255]))
        np.testing.assert_array_equal(known,[False,True,False,False])
        self.assertEqual(scores[1],255)
        _,known=decode_scores(read('AATA','A+a.;',[]))
        np.testing.assert_array_equal(known,[True,True,False,True])

    def test_projection_and_genotype_guard(self):
        r=read('ACTA','A+a.,1;T-a.,0;',[250,240]);r.set_tag('MA','4')
        result=project_read(r,'ACTA',0,4)
        np.testing.assert_array_equal(result[0],[0,0,240,250]);self.assertEqual(result[-1],1)
        with self.assertRaises(ValueError):project_read(r,'ACAA',0,4)

    def test_bad_ml_is_rejected(self):
        with self.assertRaises(ValueError):decode_scores(read('AA','A+a.,0;',[255,255]))
        with self.assertRaises(ValueError):decode_scores(read('AA','A+a.,4;',[255]))
        self.assertEqual(sampled_read('molecule'),sampled_read('molecule'))

    def test_barcode_and_molecule_weights_answer_different_questions(self):
        y=np.array([1.,1.,0.]);bc=np.array(['a','a','b'])
        self.assertEqual(scalar(y,bc,'barcode')['value'],.5)
        self.assertAlmostEqual(scalar(y,bc,'molecule')['value'],2/3)
        y[0]=np.nan
        self.assertEqual(scalar(y,bc,'barcode')['value'],.5)
        self.assertIsNone(scalar(np.array([np.nan]),np.array(['a']),'barcode')['value'])

    def test_raw_window_uses_shared_strand_and_requires_callable_positions(self):
        np.testing.assert_array_equal(shared('ATCGAT',0,6,1,base='A',radius=0),[0,4])
        scores=np.array([[255,0,255],[255,0,0]],dtype=np.uint8)
        known=np.array([[True,False,True],[True,True,True]])
        quality=np.array([[40,40,40],[40,20,40]])
        result=methylation_outcome(scores,known,quality,np.array([0,1]),baseq=30)
        self.assertTrue(np.isnan(result).all())
        result=methylation_outcome(scores,known,quality,np.array([0,2]))
        np.testing.assert_array_equal(result,[1.,.5])

    def test_empty_and_missing_barcode_profiles_do_not_become_zero(self):
        y=np.array([[1.,np.nan],[0.,np.nan],[1.,.5]])
        means,keys=barcode_means(y,np.array(['a','a','b']))
        self.assertEqual(means[0,0],.5);self.assertTrue(np.isnan(means[0,1]))
        self.assertEqual(barcode_means(np.empty((0,211)),np.array([],dtype=str))[0].shape,(0,211))

    def test_one_per_barcode_cannot_be_dominated_by_high_depth(self):
        y=np.r_[np.ones(100),0.];bc=np.array(['a']*100+['b'])
        result=resample_one_per_barcode(y,bc,np.random.default_rng(1),50)
        np.testing.assert_array_equal(result,np.full(50,.5))

    def test_joint_states_are_barcode_weighted_and_conserve_mass(self):
        fractions,nm,nb=joint_distribution(np.array([0.,0.,1.]),np.array([[4],[4],[8]]),np.ones(3,dtype=bool),np.array(['a','a','b']))
        self.assertEqual((nm,nb),(3,2));self.assertEqual(fractions[1],.5);self.assertEqual(fractions[14],.5)
        self.assertEqual(fractions.sum(),1)
        empty,_,_=joint_distribution(np.array([np.nan]),np.array([[4]]),np.ones(1,dtype=bool),np.array(['a']))
        self.assertTrue(np.isnan(empty).all())


if __name__=='__main__':unittest.main()

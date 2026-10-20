import json
import unittest
import numpy as np
from fsmpra.landscape import category, decode_calls, state_code, union, reconcile, fractions, number


class LandscapeTests(unittest.TestCase):
    def test_zero_coverage_is_distinct_from_missing_variant(self):
        result=fractions(np.array([0.,0.,1.]),np.array([2,0,2]))
        self.assertEqual(result[0],0)
        self.assertTrue(np.isnan(result[1]))
        self.assertIsNone(number(result[1]))
        self.assertEqual(result[2],.5)

    def row(self, tf=None, nuc=None, msp=None):
        return {k+'_ref_calls0':json.dumps(v or []) for k,v in [('tf',tf),('nuc',nuc),('msp',msp)]}

    def call(self, length, a=1234, b=1254):
        return dict(query_length=length,ref_intervals0=[[a,b]])

    def test_original_length_survives_clipping(self):
        groups,_=decode_calls(self.row(nuc=[self.call(1200)]))
        self.assertEqual(groups[6],[[1234,1254]])
        self.assertFalse(groups[3])
        self.assertEqual([category('nuc',n) for n in [84,85,250,251,500,501,1000,1001]],[2,3,3,4,4,5,5,6])

    def test_overlap_is_unioned_without_dropping_call_lengths(self):
        groups,calls=decode_calls(self.row(tf=[self.call(20),self.call(20,1244,1264)]))
        self.assertEqual(groups[0],[[1234,1264]])
        self.assertEqual(len(calls),2)

    def test_empty_is_zero_missing_is_error(self):
        groups,_=decode_calls(self.row())
        self.assertEqual(state_code(groups),0)
        row=self.row();row['tf_ref_calls0']=''
        with self.assertRaises(ValueError):decode_calls(row)
        with self.assertRaises(ValueError):union([[1233,1240]])

    def test_states_are_presence_not_coverage_or_short_nucleosomes(self):
        groups,_=decode_calls(self.row(tf=[self.call(20)],nuc=[self.call(150),self.call(900)]))
        self.assertEqual(state_code(groups),7)
        groups,_=decode_calls(self.row(nuc=[self.call(60)],msp=[self.call(1000)]))
        self.assertEqual(state_code(groups),0)

    def test_barcode_weighting_reconciliation(self):
        counts=np.zeros((3,2,211),dtype=int)
        counts[0]=1;counts[1]=2;counts[2]=3
        depth=np.array([1,4,3]);vi=np.array([0,0,1])
        expected=np.array([np.ones((2,211))*1.5,np.ones((2,211))])
        self.assertEqual(reconcile(counts,depth,vi,expected),0)
        with self.assertRaises(ValueError):reconcile(counts,depth,vi,expected+0.1)


if __name__=='__main__':unittest.main()

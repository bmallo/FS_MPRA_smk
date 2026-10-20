import unittest
import numpy as np
from fsmpra.landscape_atlas import regions, lists, at_change


class AtlasTests(unittest.TestCase):
    def test_wt_peak_windows_stay_inside_element(self):
        y=np.zeros(211);y[3:8]=.05;y[190:195]=.1
        result=regions(y)
        peaks=[r for r in result if r['kind']=='WT_peak_neighborhood']
        self.assertEqual(len(peaks),2)
        for r in result:
            self.assertTrue(0<=r['start0']<r['end0']<=211)
            self.assertEqual(r['end0']-r['start0'],20)
        self.assertEqual(len(result)-len(peaks),21)
        self.assertEqual(result[-1]['end0'],211)

    def test_missing_data_survive_json_conversion(self):
        self.assertEqual(lists(np.array([[0,np.nan],[.25,.5]])),[[0,None],[.25,.5]])

    def test_at_count_flags_do_not_claim_context_preserved(self):
        self.assertEqual(at_change('A','G'),'loss')
        self.assertEqual(at_change('C','T'),'gain')
        self.assertEqual(at_change('A','T'),'retained_count')
        self.assertEqual(at_change('C','G'),'retained_count')


if __name__=='__main__':unittest.main()

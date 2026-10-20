import itertools
import unittest
import numpy as np
from fsmpra.biology_motifs import pwm,score_distribution,scan,variant_delta
from fsmpra.biology_audit import shared_positions,standardized


class BiologyTests(unittest.TestCase):
    def test_score_tail_matches_exhaustive_sequence_enumeration(self):
        bg=np.array([.2,.3,.3,.2]);weights=pwm(np.array([[10,2,1,0],[1,9,2,2],[0,1,8,2]]),bg)
        q,offset,tail=score_distribution(weights,bg)
        scored=[(sum(q[i,b] for i,b in enumerate(seq)),np.prod(bg[list(seq)])) for seq in itertools.product(range(4),repeat=3)]
        for value in set(s for s,p in scored):
            self.assertAlmostEqual(tail[value-offset],sum(p for s,p in scored if s>=value),places=12)

    def test_reverse_hit_uses_complemented_reversed_weights(self):
        counts=np.ones((6,4));counts[np.arange(6),[0,0,1,1,2,3]]=100
        hits,_=scan('ACGGTT',[dict(database='test',matrix_id='test',name='test',counts=counts)],100,threshold=1.)
        hit=next(h for h in hits if h['strand']=='-')
        w=pwm(counts,np.full(4,.25))
        self.assertAlmostEqual(hit['score_bits'],sum(w[i,b] for i,b in enumerate([0,0,1,1,2,3])))
        self.assertAlmostEqual(variant_delta(hit,0,'A','C'),w[-1,2]-w[-1,3])
        self.assertIsNone(variant_delta(hit,6,'A','C'))

    def test_shared_at_excludes_changed_base_and_flanks(self):
        self.assertEqual(shared_positions('ATCGAT',0,6,1,0),[0,4,5])
        self.assertEqual(shared_positions('ATCGAT',0,6,1,2),[4,5])

    def test_standardization_uses_wt_weights_and_rejects_sparse_bins(self):
        y=np.r_[np.ones(10)*.2,np.ones(10)*.6];feature=np.r_[np.zeros(10),np.ones(10)]
        reference=(np.array([.5]),np.array([.75,.25]),[.1,.5])
        result,n,_=standardized(y,np.arange(20),feature,reference)
        self.assertAlmostEqual(result,10);self.assertEqual(n,10)
        self.assertIsNone(standardized(y,np.arange(19),feature,reference)[0])


if __name__=='__main__':unittest.main()

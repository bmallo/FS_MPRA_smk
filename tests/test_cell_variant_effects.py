import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from fsmpra.cell_concordance import classify
from fsmpra.io import write_json
from fsmpra.llr_coverage import coverage_from_calls
from fsmpra.site_effects import bh
from fsmpra.site_effects_clone import run as clone_run
from fsmpra.site_extract_hek import ELEM, oligo_id


def synthetic_extraction(path, rng, wt_clones=80, wt_depth=150, snvs=25, snv_depth=300, tau=0.05):
    """Clonal library: every clone gets its own accessibility offset (between-clone variance tau);
    every SNV is one clone with no true effect."""
    wt = ''.join(rng.choice(list('ACGT'), 211))
    variants = [dict(variant_id='lib:wt', source_variant='WT', oligo_pos1=None, ref=None, alt=None, sequence=wt, **{'class': 'WT'})]
    plan = [(0, wt_depth)] * wt_clones
    for p0 in rng.choice(np.arange(20, 190), snvs, replace=False):
        alt = next(b for b in 'ACGT' if b != wt[p0])
        variants.append(dict(variant_id=f'lib:{p0}{alt}', source_variant=f'{p0 + 21}{wt[p0]}>{alt}', oligo_pos1=int(p0 + 21),
                             ref=wt[p0], alt=alt, sequence=wt[:p0] + alt + wt[p0 + 1:], **{'class': 'SNV'}))
        plan.append((len(variants) - 1, snv_depth))
    at = np.array([b in 'AT' for b in wt])
    scores, known, cov, vi, cluster, distal = [], [], [], [], [], []
    for c, (v, n) in enumerate(plan):
        offset = rng.normal(0, tau)
        d = rng.uniform(0, 1, n)
        d[rng.random(n) < .12] = 0.
        p = np.where(d <= .01, .002, np.where(d >= .6, .8, np.clip(.15 + .5 * d + offset, .01, .95)))
        call = (rng.random((n, 211)) < p[:, None]) & at
        scores.append(np.where(call, 250, 0).astype(np.uint8))
        known.append(np.repeat(at[None], n, 0))
        msp = rng.random((n, 211)) < np.clip(.5 + offset, 0, 1)
        cov.append(np.where(msp, 2, 4).astype(np.uint8))
        vi += [v] * n
        cluster += [c] * n
        distal.append(d)
    n = len(vi)
    path.mkdir()
    np.save(path / 'scores.npy', np.concatenate(scores))
    np.save(path / 'known.npy', np.concatenate(known))
    np.save(path / 'qualities.npy', np.full((n, 211), 40, np.uint8))
    np.save(path / 'coverage.npy', np.concatenate(cov))
    np.savez(path / 'meta.npz', variant_index=np.array(vi, np.int16), cluster=np.array(cluster), distal_m6a=np.concatenate(distal),
             eligible=np.ones(n, bool), evaluated=np.ones(n, bool))
    write_json(path / 'variants.json', variants)
    write_json(path / 'audit.json', dict(schema_version='fsmpra.site_extract.v1', sources=[]))


def config(path):
    analysis = dict(cohort='eligible', baseq_min=20, ml_min=128, unmethylated_max_distal_m6a=.01, methylated_min_distal_m6a=.6,
                    labeling_strata=3, min_class_molecules=10, min_stratum_molecules=5, seed=7, mask_radius=4,
                    windows={'local': [5, 10], 'flank': [11, 30]}, regions0={},
                    primary_endpoints=['local:A_norm', 'flank:A_norm', 'msp_C_std', 'nuc_85_250_C_std'])
    write_json(path, dict(published_activity={}, analysis=analysis,
                          clone_null=dict(min_wt_clone_molecules=100, min_variant_molecules=200, within_bootstrap=60,
                                          calibration_positions=3)))


def table(path):
    with gzip.open(path, 'rt') as h:
        return list(csv.DictReader(h, delimiter='\t'))


class LlrCoverageTests(unittest.TestCase):
    def test_bits_msp_and_reversal(self):
        cov = coverage_from_calls(np.array([0, 0, 1]), np.array([0, 15, 5]), np.array([12, 40, 205]), 2, (5, 25))
        # Row 0: an 12-bp TF call (clipped to 5-12) and a 25-bp TF call; MSP everywhere (no call > 60 bp).
        self.assertTrue(((cov[0] & 2) > 0).all())
        self.assertEqual(int(((cov[0] & 1) > 0).sum()), 7 + 10)
        # Row 1: one 200-bp call covers the window: nucleosome-sized, never MSP.
        self.assertTrue(((cov[1] & 4) > 0).all() and not (cov[1] & 2).any())
        rev = coverage_from_calls(np.array([0]), np.array([5]), np.array([17]), 1, (5, 25), reverse=True)
        self.assertTrue(((rev[0, -12:] & 1) > 0).all() and not (rev[0, :-12] & 1).any())


class HekMappingTests(unittest.TestCase):
    def test_oligo_ids_complement_reverse_strand(self):
        last = ELEM[1] - 1
        self.assertEqual(oligo_id(last, 'G', 'A'), '21C>T')
        self.assertEqual(oligo_id(ELEM[0], 'A', 'G'), '231T>C')


class CloneNullTests(unittest.TestCase):
    def test_clone_null_is_calibrated_where_within_clone_se_is_not(self):
        rng = np.random.default_rng(11)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            synthetic_extraction(tmp / 'extract', rng)
            config(tmp / 'config.json')
            clone_run(tmp / 'config.json', tmp / 'extract', tmp / 'out')
            with open(tmp / 'out/calibration_summary.tsv') as h:
                gate = {r['endpoint']: r for r in csv.DictReader(h, delimiter='\t')}
            for ep in ('msp_C_std', 'element:C_std'):
                self.assertTrue(.75 <= float(gate[ep]['sd_z']) <= 1.3, (ep, gate[ep]['sd_z']))
            eff = [r for r in table(tmp / 'out/effects.tsv.gz') if r['endpoint'] == 'msp_C_std']
            self.assertEqual(len(eff), 25)
            z = np.array([float(r['z']) for r in eff])
            naive = np.array([float(r['delta']) / float(r['se_within']) for r in eff])
            self.assertLess(np.mean(np.abs(z) > 1.96), .2)
            self.assertGreater(np.std(naive), 2.)


class ConcordanceTests(unittest.TestCase):
    def test_partial_conjunction_controls_false_replication(self):
        rng = np.random.default_rng(5)
        hits = []
        for _ in range(200):
            # Real effects in HepG2 only; none in HEK293T.
            z1 = np.r_[rng.normal(4, 1, 30), rng.normal(0, 1, 70)]
            z2 = rng.normal(0, 1, 100)
            p1, p2 = 2 * (1 - _phi(np.abs(z1))), 2 * (1 - _phi(np.abs(z2)))
            p_rep = np.where(np.sign(z1) == np.sign(z2), np.maximum(p1, p2), 1.)
            hits.append((bh(p_rep) < .05).any())
        self.assertLess(np.mean(hits), .1)

    def test_classes(self):
        self.assertEqual(classify(.01, .02, 1., 2.), 'replicated')
        self.assertEqual(classify(.01, .02, 1., -2.), 'opposite')
        self.assertEqual(classify(.01, .5, 1., 2.), 'HepG2 only')
        self.assertEqual(classify(.5, .5, 1., 2.), 'neither')


def _phi(x):
    from scipy.stats import norm
    return norm.cdf(x)


if __name__ == '__main__':
    unittest.main()


class FootprintEffectsTests(unittest.TestCase):
    def test_masks_exclude_core_and_null_is_calibrated(self):
        from fsmpra.footprint_effects import footprint_masks, run as fp_run
        names, masks = footprint_masks(100, {'fp_local': (5, 20)}, 4)
        self.assertEqual(int(masks[0].sum()), 32)
        self.assertEqual(masks[0, 96:105].sum(), 0)
        rng = np.random.default_rng(3)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            synthetic_extraction(tmp / 'extract', rng, wt_clones=300, wt_depth=20, snvs=20, snv_depth=300, tau=0.)
            cov = np.load(tmp / 'extract/coverage.npy')
            tf = rng.random(cov.shape) < .05
            np.save(tmp / 'extract/coverage.npy', cov | tf.astype(np.uint8))
            write_json(tmp / 'extract/audit.json', dict(schema_version='fsmpra.site_extract.v1', sources=[], llr_coverage={}))
            analysis = dict(cohort='eligible', unmethylated_max_distal_m6a=.01, methylated_min_distal_m6a=.6, labeling_strata=3,
                            min_class_molecules=10, min_stratum_molecules=5, mask_radius=4, bootstrap=60)
            write_json(tmp / 'c.json', dict(published_activity={}, analysis=analysis, footprints=dict(
                tracks=['tf', 'nuc'], windows={'fp_local': [5, 20]}, primary_endpoints=['tf:fp_local:net'],
                min_molecules=200, aa_draws=120, seed=1)))
            fp_run(tmp / 'c.json', tmp / 'extract', tmp / 'out')
            with open(tmp / 'out/calibration_summary.tsv') as h:
                gate = {r['endpoint']: r for r in csv.DictReader(h, delimiter='\t')}
            self.assertTrue(.75 <= float(gate['tf:fp_local:net']['sd_z']) <= 1.3, gate['tf:fp_local:net']['sd_z'])

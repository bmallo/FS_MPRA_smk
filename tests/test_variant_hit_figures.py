import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from fsmpra.io import write_json
from tools import variant_hit_figures as vhf


def extraction(tmp, rng, per_allele=24):
    """Small synthetic extraction on the real WT sequence with every designed SNV."""
    wt = vhf.wt_sequence()
    variants = [dict(variant_id='lib:wt', source_variant='WT', oligo_pos1=None, ref=None, alt=None, sequence=wt, **{'class': 'WT'})]
    for p0, ref in enumerate(wt):
        for alt in 'ACGT':
            if alt != ref:
                variants.append(dict(variant_id=f'lib:{p0}{alt}', source_variant=f'{p0 + 21}{ref}>{alt}', oligo_pos1=p0 + 21, ref=ref,
                                     alt=alt, sequence=wt[:p0] + alt + wt[p0 + 1:], **{'class': 'SNV'}))
    deep = set(vhf.TOP_HITS) | {f'{p}{wt[p - 21]}>{b}' for p in vhf.COMPARE_POSITIONS for b in 'ACGT' if b != wt[p - 21]}
    n_per = np.array([40 * per_allele if v['class'] == 'WT' else 300 if v['source_variant'] in deep else per_allele for v in variants])
    vi = np.repeat(np.arange(len(variants)), n_per)
    n = len(vi)
    u = rng.random(n)
    distal = np.where(u < .1, rng.uniform(0, .008, n), np.where(u < .25, rng.uniform(.62, .8, n), rng.uniform(.05, .55, n)))
    known = np.array([[b in 'AT' for b in variants[i]['sequence']] for i in range(len(variants))])[vi]
    open_ = np.full((n, 211), .4)
    target = vi == next(i for i, v in enumerate(variants) if v['source_variant'] == '146C>G')
    open_[target] = .1                                             # one hit closes the element
    msp = rng.random((n, 211)) < open_
    cov = np.where(msp, 2, np.where(rng.random((n, 211)) < .5, 4, 1)).astype(np.uint8)
    called = known & msp & (rng.random((n, 211)) < .7)
    ext = tmp / 'extract'
    ext.mkdir()
    np.save(ext / 'scores.npy', np.where(called, 250, 0).astype(np.uint8))
    np.save(ext / 'known.npy', known)
    np.save(ext / 'qualities.npy', np.full((n, 211), 40, np.uint8))
    np.save(ext / 'coverage.npy', cov)
    np.savez(ext / 'meta.npz', variant_index=vi.astype(np.int16), cluster=np.arange(n), distal_m6a=distal.astype(np.float32),
             eligible=np.ones(n, bool), evaluated=np.ones(n, bool))
    write_json(ext / 'variants.json', variants)
    return ext, variants, target


class VariantHitFigureTests(unittest.TestCase):
    def test_motif_changes_score_both_alleles_and_flag_direction(self):
        rows = vhf.motif_changes(['56G>C', '55C>A'])
        self.assertTrue(rows)
        zn = [r for r in rows if r['variant'] == '56G>C' and r['name'] == 'ZN652']
        self.assertTrue(zn and all(r['change'] == 'loss' for r in zn))       # P1-core substitution weakens the ZNF652-like match
        for r in rows:
            self.assertLessEqual(min(r['wt_p'], r['alt_p']), 1e-3)
            self.assertLessEqual(r['oligo_start'], vhf.parse(r['variant'])[0])
            self.assertGreaterEqual(r['oligo_end'], vhf.parse(r['variant'])[0])

    def test_profiles_and_figures_run_on_real_names(self):
        rng = np.random.default_rng(3)
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            ext, variants, target = extraction(tmp, rng)
            vhf.profiles(ext, tmp / 'profiles', raster_n=20)
            prof = np.load(tmp / 'profiles' / 'profiles.npz')
            self.assertEqual(prof['msp_C_std'].shape, (len(variants), 211))
            names = list(prof['source_variant'])
            self.assertEqual(names[0], 'WT')
            hit = names.index('146C>G')
            self.assertLess(np.nanmean(prof['msp_C_std'][hit]), np.nanmean(prof['msp_C_std'][0]) - .2)
            self.assertTrue(np.allclose(np.nanmean(prof['msp_ALL'][0]), .4, atol=.02))
            shallow = names.index('22' + vhf.wt_sequence()[1] + '>' + next(b for b in 'ACGT' if b != vhf.wt_sequence()[1]))
            self.assertTrue(np.isnan(prof['msp_C_std'][shallow]).all())    # too few molecules per labeling stratum
            self.assertTrue(np.isfinite(prof['msp_ALL'][shallow]).all())
            raster = np.load(tmp / 'profiles' / 'rasters.npz')
            self.assertIn('146C>G', list(raster['names']))
            self.assertEqual(raster['state_0'].shape[1], 211)
            self.assertLessEqual(len(raster['state_0']), 20)
            audit = json.loads((tmp / 'profiles' / 'audit.json').read_text())
            self.assertEqual(audit['schema_version'], 'fsmpra.hit_profiles.v1')
            with self.assertRaises(FileExistsError):
                vhf.profiles(ext, tmp / 'profiles')                      # never overwrites
            motifs = tmp / 'motifs.tsv'
            from fsmpra.io import write_table
            rows = vhf.motif_changes(['146C>G', '56G>C', '125A>T'])
            write_table(motifs, rows, list(rows[0]))
            vhf.figures(tmp / 'profiles', motifs, tmp / 'figures', variants=['146C>G', '56G>C', '125A>T'])
            out = {p.name for p in (tmp / 'figures').iterdir()}
            for name in ('top_hits_overview.png', 'motif_sites.png', 'motif_sites.tsv', 'footprint_association.tsv',
                         '146C_to_G.png', '56G_to_C.png', '125A_to_T.png',
                         'position_146.png', 'position_56.png'):
                self.assertIn(name, out)
            with self.assertRaises(FileExistsError):
                vhf.figures(None, motifs, tmp / 'figures')


if __name__ == '__main__':
    unittest.main()

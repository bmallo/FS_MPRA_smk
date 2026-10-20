import csv
import gzip
import json
import tempfile
import unittest
from array import array
from pathlib import Path
import numpy as np
import pysam
from fsmpra.io import source_info, write_json
from fsmpra.site_effects import assign_groups, measures, run as effects_run, window_masks
from fsmpra.site_extract import run as extract_run

COMP = str.maketrans('ACGT', 'TGCA')


def library(rng):
    wt = ''.join(rng.choice(list('ACGT'), 211))
    variants = [dict(variant_id='lib:wt', source_variant='WT', **{'class': 'WT'}, oligo_pos1=None, ref=None, alt=None, sequence=wt)]
    for p0, ref in enumerate(wt):
        for alt in 'ACGT':
            if alt != ref:
                variants.append(dict(variant_id=f'lib:{p0}{alt}', source_variant=f'{p0 + 21}{ref}>{alt}', **{'class': 'SNV'},
                                     oligo_pos1=p0 + 21, ref=ref, alt=alt, sequence=wt[:p0] + alt + wt[p0 + 1:]))
    return variants


def mm_tags(stored, calls, reverse):
    """MM/ML for retained calls at stored-orientation positions."""
    forward = stored.translate(COMP)[::-1] if reverse else stored
    n = len(stored)
    groups = {'A': [], 'T': []}
    for q in calls:
        f = n - 1 - q if reverse else q
        groups[forward[f]].append(f)
    mm, ml = [], []
    for base, channel in (('A', 'A+a.'), ('T', 'T-a.')):
        opportunities = [i for i, b in enumerate(forward) if b == base]
        rank = {p: k for k, p in enumerate(opportunities)}
        chosen = sorted(groups[base])
        skips, last = [], -1
        for p in chosen:
            skips.append(rank[p] - last - 1)
            last = rank[p]
        mm.append(channel + ''.join(f',{s}' for s in skips) + ';')
        ml += [250] * len(chosen)
    return 'C+m.;' + ''.join(mm), ml


class SiteExtractTests(unittest.TestCase):
    def test_library_extraction_projects_calls_flags_and_skips_unselected(self):
        rng = np.random.default_rng(3)
        variants = library(rng)
        backbone = ''.join(rng.choice(list('ACGT'), 1600))
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            bam_path = tmp / 'fp.bam'
            plan = [('r1', 0, 'BC1', False, True, [1240, 1300, 1500]), ('r2', 5, '', True, False, [1250]),
                    ('skip', 0, 'BC1', False, True, []), ('r3', 0, 'BC1', False, False, [])]
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'plasmid', 'LN': 1600}]}
            expected = {}
            with pysam.AlignmentFile(str(bam_path), 'wb', header=header) as out:
                for name, vi, bc, reverse, msp, wanted in plan:
                    seq = backbone[:1234] + variants[vi]['sequence'] + backbone[1445:]
                    calls = [q for q in wanted if seq[q] in 'AT']
                    expected[name] = calls
                    r = pysam.AlignedSegment(out.header)
                    r.query_name, r.query_sequence, r.flag = name, seq, 16 if reverse else 0
                    r.reference_id, r.reference_start, r.mapping_quality, r.cigarstring = 0, 0, 60, '1600M'
                    r.query_qualities = array('B', [40] * 1600)
                    mm, ml = mm_tags(seq, calls, reverse)
                    r.set_tag('MM', mm)
                    r.set_tag('ML', array('B', ml))
                    r.set_tag('MA', '1600;msp+:1235-211' if msp else '1600;msp+:1-50')
                    out.write(r)
            pysam.index(str(bam_path))
            fields = ['read_id', 'observed_variant_id', 'corrected_barcode', 'barcode_assignment_status', 'passes', 'read_quality',
                      'mapq', 'read_length', 'strand', 'circular', 'footprint_status', 'footprint_analysis_eligible',
                      'footprint_analysis_eligible_with_barcode', 'whole_m6a_fraction', 'distal_m6a_fraction',
                      'distal_nucleosome_fraction', 'whole_nucleosome_segments', 'element_msp_fraction',
                      'element_tf_small_fraction', 'region_complete', 'm6a_status']
            with open(tmp / 'features.tsv', 'w', newline='') as h:
                w = csv.DictWriter(h, fieldnames=fields, delimiter='\t')
                w.writeheader()
                for name, vi, bc, reverse, msp, _ in plan:
                    w.writerow(dict(read_id=name, observed_variant_id=variants[vi]['variant_id'], corrected_barcode=bc,
                                    barcode_assignment_status='supported' if bc else '', passes=12, read_quality=.999, mapq=60,
                                    read_length=1600, strand='-' if reverse else '+', circular=1, footprint_status='evaluated',
                                    footprint_analysis_eligible=1, footprint_analysis_eligible_with_barcode=int(bool(bc)),
                                    whole_m6a_fraction=.3, distal_m6a_fraction='' if name == 'r3' else .3, distal_nucleosome_fraction=.5,
                                    whole_nucleosome_segments=4, element_msp_fraction=1. if msp else 0., element_tf_small_fraction=0.,
                                    region_complete=0 if name == 'skip' else 1, m6a_status='ok'))
            write_json(tmp / 'feature_audit.json', dict(sources=[source_info(bam_path)]))
            write_json(tmp / 'catalog.json', dict(library=dict(contig='plasmid', variants=variants)))
            write_json(tmp / 'contract.json', dict(bam=str(bam_path)))
            write_json(tmp / 'cfg.json', dict(schema_version='fsmpra.site_methylation_config.v1', features=str(tmp / 'features.tsv'),
                                              feature_audit=str(tmp / 'feature_audit.json'), catalog=str(tmp / 'catalog.json'),
                                              footprint_contract=str(tmp / 'contract.json')))
            extract_run(tmp / 'cfg.json', tmp / 'out')
            out = tmp / 'out'
            meta = np.load(out / 'meta.npz')
            scores, known = np.load(out / 'scores.npy'), np.load(out / 'known.npy')
            coverage = np.load(out / 'coverage.npy')
            self.assertEqual(scores.shape, (3, 211))
            np.testing.assert_array_equal(meta['variant_index'], [0, 5, 0])
            self.assertEqual(meta['cluster'][0], meta['cluster'][2])          # shared barcode
            self.assertNotEqual(meta['cluster'][1], meta['cluster'][0])       # unbarcoded singleton
            self.assertTrue(np.isnan(meta['distal_m6a'][2]))                   # missing stays missing
            for row, name, vi in ((0, 'r1', 0), (1, 'r2', 5)):
                element_calls = [q - 1234 for q in expected[name] if 1234 <= q < 1445]
                np.testing.assert_array_equal(np.flatnonzero(scores[row]), element_calls)
                np.testing.assert_array_equal(known[row], [b in 'AT' for b in variants[vi]['sequence']])
            self.assertTrue(((coverage[0] & 2) > 0).all())
            self.assertFalse(((coverage[2] & 2) > 0).any())
            audit = json.loads((out / 'audit.json').read_text())
            self.assertTrue(audit['complete'])
            self.assertEqual(audit['bam_primary_rows_scanned'], 4)
            extract_run(tmp / 'cfg.json', tmp / 'smoke', limit=2)               # smoke mode stops early
            self.assertEqual(np.load(tmp / 'smoke' / 'scores.npy').shape, (2, 211))
            np.testing.assert_array_equal(np.load(tmp / 'smoke' / 'scores.npy'), scores[:2])
            self.assertFalse(json.loads((tmp / 'smoke' / 'audit.json').read_text())['complete'])
            with self.assertRaises(FileExistsError):
                extract_run(tmp / 'cfg.json', tmp / 'out')                       # never overwrites


def simulate(tmp, rng, n_wt=20000, n_snv=1500, artifact=None, biology=None):
    """Synthetic extraction: per-site P(call) = fp + eff * access * labeling, by distal class."""
    wt = ''.join(rng.choice(list('ACGT'), 211))
    positions = list(range(25, 185, 3))[:52]
    variants = [dict(variant_id='lib:wt', source_variant='WT', oligo_pos1=None, ref=None, alt=None, sequence=wt, **{'class': 'WT'})]
    for p0 in positions:
        alt = next(b for b in 'GCAT' if b != wt[p0])
        variants.append(dict(variant_id=f'lib:{p0}', source_variant=f'{p0 + 21}{wt[p0]}>{alt}', oligo_pos1=p0 + 21, ref=wt[p0],
                             alt=alt, sequence=wt[:p0] + alt + wt[p0 + 1:], **{'class': 'SNV'}))
    at = np.array([b in 'AT' for b in wt])
    eff0 = np.where(at, .55 + .3 * rng.random(211), 0)
    blocks = {k: [] for k in ('variant_index', 'cluster', 'distal_m6a', 'scores', 'known', 'coverage')}
    cluster0 = 0
    for vi, v in enumerate(variants):
        n = n_wt if vi == 0 else n_snv
        p0 = None if vi == 0 else v['oligo_pos1'] - 21
        d = np.abs(np.arange(211) - p0) if p0 is not None else np.full(211, 999)
        eff, access = eff0.copy(), np.full(211, .5)
        if artifact and p0 == artifact:
            eff[(d >= 3) & (d <= 6)] *= 1.25
        if biology and p0 == biology:
            access[(d >= 3) & (d <= 10)] += .2
        known = np.array([b in 'AT' for b in v['sequence']])
        u = rng.random(n)
        cls = np.where(u < .1, 0, np.where(u < .25, 1, 2))
        distal = np.where(cls == 0, rng.uniform(0, .008, n), np.where(cls == 1, rng.uniform(.62, .8, n), rng.uniform(.05, .55, n)))
        labeling = np.where(cls == 0, 0., np.where(cls == 1, 1., distal / .7))
        acc = np.where(cls[:, None] == 1, 1., access[None, :])
        p = .004 + np.minimum(eff[None, :] * acc * labeling[:, None], .99)
        called = (rng.random((n, 211)) < p) & known[None, :]
        cov = np.where(rng.random((n, 211)) < .5 + .1 * (biology is not None and p0 == biology), 2, 0).astype(np.uint8)
        size = rng.integers(1, 6, n)
        clusters = cluster0 + np.repeat(np.arange(n), size)[:n]
        cluster0 = clusters.max() + 1
        for key, value in (('variant_index', np.full(n, vi)), ('cluster', clusters), ('distal_m6a', distal),
                           ('scores', np.where(called, 250, 0).astype(np.uint8)), ('known', np.repeat(known[None], n, 0)),
                           ('coverage', cov)):
            blocks[key].append(value)
    data = {k: np.concatenate(v) for k, v in blocks.items()}
    n = len(data['variant_index'])
    ext = tmp / 'extract'
    ext.mkdir()
    for name in ('scores', 'known', 'coverage'):
        np.save(ext / f'{name}.npy', data[name])
    np.save(ext / 'qualities.npy', np.full((n, 211), 40, np.uint8))
    np.savez(ext / 'meta.npz', variant_index=data['variant_index'].astype(np.int16), cluster=data['cluster'],
             distal_m6a=data['distal_m6a'].astype(np.float32), eligible=np.ones(n, bool), eligible_with_barcode=np.ones(n, bool))
    write_json(ext / 'variants.json', variants)
    write_json(ext / 'audit.json', dict(schema_version='fsmpra.site_extract.v1', sources=[]))
    act = tmp / 'activity.tsv'
    with open(act, 'w') as h:
        h.write('Chromosome\tPosition\tRef\tAlt\tValue\n')
        for v in variants[1:]:
            h.write(f"19\t{11089264 + v['oligo_pos1']}\t{v['ref']}\t{v['alt']}\t{rng.normal():.3f}\n")
    cfg = json.loads(Path('config/ldlr_site_methylation.json').read_text())
    cfg['published_activity'] = {'synthetic': str(act)}
    cfg['analysis'].update(bootstrap=60, aa_draws=40, figures=True)
    write_json(tmp / 'cfg.json', cfg)
    return variants, positions


class SiteEffectsTests(unittest.TestCase):
    def test_groups_use_distal_only_and_wt_quantile_strata(self):
        distal = np.array([0., .005, .02, .3, .59, .6, .9, np.nan])
        groups, weights, edges = assign_groups(distal, np.ones(8, bool), .01, .6, 2)
        np.testing.assert_array_equal(groups[[0, 1, 5, 6, 7]], [0, 0, 1, 1, -1])
        self.assertTrue((groups[2:5] >= 2).all())
        self.assertAlmostEqual(weights.sum(), 1)

    def test_windows_exclude_substituted_base_and_masked_flanks(self):
        names, masks = window_masks('A' * 211, 100, {'near': (1, 2), 'local': (3, 10)}, {'r': (95, 106)}, 2)
        rows = dict(zip(names, masks.astype(bool)))
        np.testing.assert_array_equal(np.flatnonzero(rows['near']), [98, 99, 101, 102])
        self.assertEqual(rows['local'].sum(), 16)
        self.assertFalse(rows['element'][98:103].any())
        np.testing.assert_array_equal(np.flatnonzero(rows['r']), [95, 96, 97, 103, 104, 105])

    def test_normalization_cancels_multiplicative_observability(self):
        masks = np.ones((1, 211))
        def sums(eff):
            g = 2 + 2
            cs, vs = np.zeros((g, 211)), np.full((g, 211), 1000.)
            cs[0], cs[1], cs[2], cs[3] = 5, 1000 * eff, 1000 * (.005 + .3 * eff), 1000 * (.005 + .5 * eff)
            return cs, vs, np.zeros((g, 3)), np.full(g, 1000.)
        a, _ = measures(*sums(.6), masks, np.array([.5, .5]))
        b, _ = measures(*sums(.75), masks, np.array([.5, .5]))
        self.assertGreater(b['M'][0] - a['M'][0], .1)
        self.assertGreater(b['C_std'][0] - a['C_std'][0], .05)
        self.assertAlmostEqual(b['A_norm'][0], a['A_norm'][0], places=2)

    def test_effects_separate_injected_artifact_from_chromatin_change(self):
        rng = np.random.default_rng(11)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            artifact, biology = 61, 121
            variants, positions = simulate(tmp, rng, artifact=artifact, biology=biology)
            effects_run(tmp / 'cfg.json', tmp / 'extract', tmp / 'out')
            out = tmp / 'out'
            with gzip.open(out / 'effects.tsv.gz', 'rt') as h:
                rows = list(csv.DictReader(h, delimiter='\t'))
            get = lambda p0, ep: next(r for r in rows if r['oligo_pos1'] == str(p0 + 21) and r['endpoint'] == ep)
            # Sequence-dependent observability: raised in M and C, cancelled by normalization.
            self.assertGreater(float(get(artifact, 'local:M')['z']), 4)
            self.assertGreater(float(get(artifact, 'local:C_std')['delta']), 0)
            self.assertLess(abs(float(get(artifact, 'local:A_norm')['z'])), 3.5)
            # Chromatin change: absent in M, retained after normalization.
            self.assertLess(abs(float(get(biology, 'local:M')['z'])), 3.5)
            self.assertGreater(float(get(biology, 'local:A_norm')['z']), 5)
            self.assertGreater(float(get(biology, 'msp_C_std')['delta']), .05)
            summary = json.loads((out / 'summary.json').read_text())
            self.assertAlmostEqual(summary['molecules']['wt']['U'] / 20000, .1, delta=.01)
            with open(out / 'calibration_summary.tsv') as h:
                calib = {r['endpoint']: r for r in csv.DictReader(h, delimiter='\t')}
            for ep in ('local:C_std', 'local:A_norm', 'element:M'):
                self.assertLess(abs(float(calib[ep]['sd_z']) - 1), .45, ep)
            with open(out / 'kernel.tsv') as h:
                kernel = list(csv.DictReader(h, delimiter='\t'))
            self.assertFalse(any(r['offset'] == '0' for r in kernel))
            self.assertTrue((out / 'kernel.png').exists() and (out / 'activity_association.tsv').exists())


if __name__ == '__main__':
    unittest.main()

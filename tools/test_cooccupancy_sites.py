"""P3.2 regression test — co-occupancy site-pair construction.

Deterministic synthetic checks of the geometry/logic (overlap, canonical
merge, FDR gate, site1 assignment, distal pairing) plus a real-rd check
of read_site_occupied + the abs->column coordinate mapping. No cluster,
no network. Run: pixi run python tools/test_cooccupancy_sites.py
"""
import sys

sys.path.insert(
    0, "/mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/"
       "FS_MPRA_smk/workflow/scripts")
import common as C  # noqa
import numpy as np  # noqa

# --- interval geometry ---
assert C._recip_overlap(0, 9, 0, 9) == 1.0
assert C._recip_overlap(0, 9, 5, 14) == 0.5
assert C._recip_overlap(0, 9, 20, 29) == 0.0
assert C._merge_bin_intervals(
    [(100, 119), (110, 129), (200, 219)], 0.5) == [(100, 129), (200, 219)]

# --- canonical sites: motif call + sig clusters, FDR gate ---
motif_result = {'motifs': [{'bin': 'TF', 'abs_start': 3402,
                            'abs_end': 3421}]}
avr = [
    {'variant_id': 'A', 'variant_fdr_q': 0.01, 'variant_pos0': 3405,
     'TF': {'significant_clusters': [
         {'abs_start': 3403, 'abs_end': 3420,
          'sum_abs_delta': 9.0, 'is_motif': True}]}},
    {'variant_id': 'B', 'variant_fdr_q': 0.02, 'variant_pos0': 3475,
     'sub_TF': {'significant_clusters': [
         {'abs_start': 3470, 'abs_end': 3481,
          'sum_abs_delta': 5.0, 'is_motif': True}]}},
    {'variant_id': 'C', 'variant_fdr_q': 0.9, 'variant_pos0': 3406,
     'TF': {'significant_clusters': [
         {'abs_start': 3403, 'abs_end': 3420,
          'sum_abs_delta': 99.0, 'is_motif': True}]}},  # FDR fail -> drop
]
cfg = C.build_cooccupancy_cfg(min_variant_instruments=1,
                              min_site_separation=20)
sites = C.build_canonical_sites(avr, motif_result, ['sub_TF', 'TF'], cfg)
assert len(sites) == 2, sites
assert any(s['bin'] == 'TF' and s['is_motif_site'] for s in sites)
assert any(s['bin'] == 'sub_TF' for s in sites)

pairs = C.build_site_pairs(avr, motif_result, None, ['sub_TF', 'TF'], cfg)
# A's SNV (3405) in TF site -> instrument; B's SNV (3475) in sub_TF site
assert set(pairs['site1_instruments']) == {'TF:3402-3421',
                                            'sub_TF:3470-3481'}, \
    pairs['site1_instruments']
# distal pair both directions (gap ~49 >= 20)
assert len(pairs['pairs']) == 2, pairs['pairs']
print("synthetic logic: PASS")

# --- real rd: coordinate mapping + read_site_occupied ---
BAM = ("/mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/"
       "FS_MPRA_smk/results/phase0_hmm/variants/LDLR_hmm.tagged.bam")
rd, *_ = C.parse_bam(BAM, C.build_analysis_bins(), target_chrom=None,
                     analysis_region=(3184, 3502))
site = {'site_id': 'TF:3300-3319', 'bin': 'TF',
        'abs_start': 3300, 'abs_end': 3319}
cols = C._site_columns(site, rd.analysis_start, rd.analysis_length)
assert cols[0] == 116 and cols[-1] == 135 and cols.size == 20, cols
occ = C.read_site_occupied(rd, site, rd.wt_indices, 0.5)
assert occ.dtype == bool and occ.shape[0] == len(rd.wt_indices)
out = {'site_id': 'x', 'bin': 'TF', 'abs_start': 9000, 'abs_end': 9010}
assert C.read_site_occupied(rd, out, rd.wt_indices[:50], 0.5).sum() == 0
print("real-rd coord/state: PASS")
print("ALL P3.2 CHECKS PASSED")

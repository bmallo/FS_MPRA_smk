"""P3.6 self-test — WT-vs-WT co-occupancy calibration harness.

Synthetic stub rd with correlated footprint structure (so WT itself has
a real site1<->site2 dependency — the harness must STILL be calibrated,
since pseudo-instruments and the conditional pool are both WT). Checks:
  (1) H0 (no injected dependency): pooled p ~ Uniform, BH FP@call_q
      ~ nominal, ~0 validated calls -> the gate does not manufacture
      false dependencies on structured real-like data.
  (2) TEETH (positive control): inject a fake cooperative loss ->
      the harness DOES produce significant calls, so a clean H0 result
      is meaningful.

Run: pixi run python tools/test_cooccupancy_wtvswt.py
"""
import sys
import numpy as np

sys.path.insert(
    0, "/mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/"
       "FS_MPRA_smk/workflow/scripts")
import common as C  # noqa


class StubRD:
    def __init__(self, cov, a_start):
        self.coverage_matrices = {'TF': cov}
        self.analysis_start = a_start
        self.analysis_length = cov.shape[1]
        self.wt_indices = np.arange(cov.shape[0])


A0, L, NW = 3000, 300, 9000
rng = np.random.default_rng(1)
# 4 "elements"; element occupancy is correlated through a per-read
# latent (some reads broadly more footprinted) -> non-trivial WT
# P(O2|O1) the harness must remain calibrated against.
elem = [(20, 39), (90, 109), (170, 189), (240, 259)]
latent = rng.normal(size=NW)
cov = np.zeros((NW, L), dtype=np.uint8)
for (s, e) in elem:
    p = 1.0 / (1.0 + np.exp(-(0.3 * latent + rng.normal(0, 0.4,
                                                        size=NW))))
    on = rng.random(NW) < p
    cov[on, s:e + 1] = 1
rd = StubRD(cov, A0)
sites = [{'site_id': f'TF:{A0+s}-{A0+e}', 'bin': 'TF',
          'abs_start': A0 + s, 'abs_end': A0 + e} for (s, e) in elem]
cfg = C.build_cooccupancy_cfg(min_site_separation=20, min_stratum=25,
                              call_q=0.10, min_consistency=0.70)

# (1) H0 — accumulate records over seeds, BH, count calls
recs = []
for sd in range(25):
    recs += C.cooccupancy_wt_vs_wt(rd, sites, cfg, n_instruments=6,
                                   instr_n=300, n_iter=1500, seed=sd)
p = np.array([r['p_two_sided'] for r in recs])
q = C.benjamini_hochberg(p)
mc = cfg['min_consistency']
calls = sum(1 for r, qq in zip(recs, q)
            if qq < cfg['call_q']
            and r['frac_instruments_consistent'] >= mc)
print(f"H0: {len(recs)} WT-vs-WT records  mean(p)={p.mean():.3f} "
      f"(~0.5)  FP@0.05={np.mean(p<0.05):.3f}  FP@0.10="
      f"{np.mean(p<0.10):.3f}  BH calls@q<{cfg['call_q']} & "
      f"consistency: {calls}/{len(recs)}")
assert 0.40 <= p.mean() <= 0.60, p.mean()
assert np.mean(p < 0.05) <= 0.12, np.mean(p < 0.05)
assert calls <= max(1, int(0.03 * len(recs))), calls
print("H0 calibration on structured real-like WT: PASS")

# (2) TEETH — inject a fake cooperative loss
recs2 = []
for sd in range(8):
    recs2 += C.cooccupancy_wt_vs_wt(rd, sites, cfg, n_instruments=6,
                                    instr_n=300, n_iter=1500, seed=sd,
                                    inject_dependency=0.7)
p2 = np.array([r['p_two_sided'] for r in recs2])
print(f"TEETH: injected dependency -> FP-style sig fraction "
      f"p<0.05 = {np.mean(p2 < 0.05):.2f} (want high)")
assert np.mean(p2 < 0.05) >= 0.5, np.mean(p2 < 0.05)
print("positive-control teeth: PASS")
print("ALL P3.6 SELF-TEST CHECKS PASSED")

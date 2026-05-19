"""P3.3 regression test — co-occupancy conditional-resampling test.

Synthetic ground truth, no rd/network. Validates:
  (1) CALIBRATION under H0: a variant that disrupts site 1 but whose
      site 2 follows the WT P(O2|O1) channel must NOT look like a
      dependency -> p-values ~Uniform, false-positive rate ~ nominal,
      mean excess ~ 0. (mini P3.6 preview)
  (2) POWER: a variant whose site 2 collapses BEYOND the WT channel
      when site 1 is lost -> significant, excess negative (site2_loss).
  (3) Degenerate WT stratum (< min) -> pair untestable (None).

Run: pixi run python tools/test_cooccupancy_stats.py
"""
import sys
import numpy as np

sys.path.insert(
    0, "/mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/"
       "FS_MPRA_smk/workflow/scripts")
import common as C  # noqa

P1, P0 = 0.60, 0.20          # WT P(O2=1|O1=1), P(O2=1|O1=0)
N_WT, B = 6000, 3000


def gen_states(rng, n, p_o1, p2_1, p2_0):
    o1 = rng.random(n) < p_o1
    p = np.where(o1, p2_1, p2_0)
    o2 = rng.random(n) < p
    return o1, o2


# (3) degenerate stratum -> None
rng = np.random.default_rng(0)
o1w = np.zeros(500, dtype=bool)          # no O1=1 reads
o2w = rng.random(500) < 0.3
assert C._wt_conditional_from_states(o1w, o2w, 100, rng, 25) is None
print("degenerate-stratum guard: PASS")

# (1) CALIBRATION under H0
rng = np.random.default_rng(42)
n_trials = 300
pvals, excesses = [], []
for _ in range(n_trials):
    o1w, o2w = gen_states(rng, N_WT, 0.5, P1, P0)
    wt = C._wt_conditional_from_states(o1w, o2w, B, rng, 25)
    assert wt is not None
    # variant disrupts site 1 (O1 mostly lost) but site 2 still follows
    # the SAME WT conditional -> H0 true (no excess dependency).
    o1v, o2v = gen_states(rng, 400, 0.12, P1, P0)
    r = C._cooccupancy_test_from_states(o1v, o2v, wt, rng)
    pvals.append(r['p_two_sided'])
    excesses.append(r['excess'])
pvals = np.array(pvals)
fp05 = np.mean(pvals < 0.05)
fp10 = np.mean(pvals < 0.10)
print(f"H0: n={n_trials}  FP@0.05={fp05:.3f} (~0.05)  "
      f"FP@0.10={fp10:.3f} (~0.10)  mean|excess|={np.mean(np.abs(excesses)):.4f}  "
      f"mean(excess)={np.mean(excesses):+.4f}")
assert fp05 <= 0.10, f"anti-conservative: FP@0.05={fp05:.3f}"
assert abs(np.mean(excesses)) < 0.02, np.mean(excesses)
print("calibration under H0: PASS (not anti-conservative, ~zero excess)")

# (2) POWER: cooperative — when site1 lost, site2 collapses below the
# WT channel (P(O2|O1=0) = 0.03 << WT's 0.20).
rng = np.random.default_rng(7)
hits = 0
NP = 60
for _ in range(NP):
    o1w, o2w = gen_states(rng, N_WT, 0.5, P1, P0)
    wt = C._wt_conditional_from_states(o1w, o2w, B, rng, 25)
    o1v, o2v = gen_states(rng, 400, 0.12, P1, 0.03)   # extra suppression
    r = C._cooccupancy_test_from_states(o1v, o2v, wt, rng)
    if r['p_two_sided'] < 0.05 and r['excess'] < 0 \
            and r['direction'] == 'site2_loss':
        hits += 1
power = hits / NP
print(f"POWER: cooperative detected {hits}/{NP} (power={power:.2f}, "
      f"p<0.05 & excess<0 & site2_loss)")
assert power >= 0.8, f"underpowered: {power:.2f}"
print("power on true cooperativity: PASS")
print("ALL P3.3 CHECKS PASSED")

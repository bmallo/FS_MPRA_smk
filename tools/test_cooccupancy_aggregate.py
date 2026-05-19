"""P3.4 regression test — cross-variant co-occupancy aggregation.

Synthetic ground truth, no rd/network. Validates the JOINT null
(instruments share the per-pair WT-conditional bootstrap):
  (1) CALIBRATION: K instruments all following the WT P(O2|O1) channel
      -> pair p ~Uniform, FP ~ nominal (the shared-WT correlation must
      not make it anti-conservative).
  (2) POWER: K instruments all showing extra site-2 suppression beyond
      the channel -> significant, site2_loss, high consistency.
  (3) SPECIFICITY TO CONSISTENCY: instruments that scatter (half down,
      half up) must NOT be CALLED — a call requires the directional-
      consistency gate (Phase-2 sign-consistency analogue), not just a
      significant mean-excess p. (The p may still be significant — that
      only says "site 2 deviates on average"; the gate enforces the
      biological multiple-instrument requirement.)

Run: pixi run python tools/test_cooccupancy_aggregate.py
"""
import sys
import numpy as np

sys.path.insert(
    0, "/mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/"
       "FS_MPRA_smk/workflow/scripts")
import common as C  # noqa

P1, P0 = 0.60, 0.20
N_WT, B, K, NK = 5000, 2000, 6, 300


def states(rng, n, p_o1, p2_1, p2_0):
    o1 = rng.random(n) < p_o1
    o2 = rng.random(n) < np.where(o1, p2_1, p2_0)
    return o1, o2


def wt(rng):
    o1w, o2w = states(rng, N_WT, 0.5, P1, P0)
    return C._wt_conditional_from_states(o1w, o2w, B, rng, 25)


# (1) CALIBRATION under H0
rng = np.random.default_rng(11)
pv = []
for _ in range(250):
    wc = wt(rng)
    inst = [states(rng, NK, 0.12, P1, P0) for _ in range(K)]
    r = C._cooccupancy_aggregate_from_states(inst, wc, rng)
    pv.append(r['p_two_sided'])
pv = np.array(pv)
fp05, fp10 = np.mean(pv < 0.05), np.mean(pv < 0.10)
print(f"H0 (K={K} consistent-null instruments): FP@0.05={fp05:.3f} "
      f"(~0.05) FP@0.10={fp10:.3f} (~0.10)")
assert fp05 <= 0.10, f"anti-conservative joint null: {fp05:.3f}"
print("aggregate calibration: PASS")

# (2) POWER: all instruments extra-suppress site 2 (0.20 -> 0.03)
rng = np.random.default_rng(3)
hits = 0
for _ in range(40):
    wc = wt(rng)
    inst = [states(rng, NK, 0.12, P1, 0.03) for _ in range(K)]
    r = C._cooccupancy_aggregate_from_states(inst, wc, rng)
    if (r['fdr_q'] if 'fdr_q' in r else r['p_two_sided']) < 0.05 \
            and r['direction'] == 'site2_loss' \
            and r['frac_instruments_consistent'] >= 0.8:
        hits += 1
print(f"POWER: {hits}/40 significant & site2_loss & consistent")
assert hits >= 32, f"underpowered: {hits}/40"
print("aggregate power: PASS")

# (3) SPECIFICITY: half instruments down, half up -> inconsistent ->
# the directional-consistency CALL gate must reject it (consistency
# ~0.5 << 0.70), even if the mean-excess p is significant.
rng = np.random.default_rng(99)
MIN_CONSIST = C.DEFAULT_COOCCUPANCY_CFG['min_consistency']
false_calls = 0
sig_p = 0
for _ in range(60):
    wc = wt(rng)
    inst = []
    for j in range(K):
        p2_0 = 0.03 if j % 2 == 0 else 0.45      # alternating sign
        inst.append(states(rng, NK, 0.12, P1, p2_0))
    r = C._cooccupancy_aggregate_from_states(inst, wc, rng)
    if r['p_two_sided'] < 0.05:
        sig_p += 1
    # the CALL = significant AND consistency gate (what run_cooccupancy
    # marks as is_call)
    if (r['p_two_sided'] < 0.05
            and r['frac_instruments_consistent'] >= MIN_CONSIST):
        false_calls += 1
print(f"SPECIFICITY: inconsistent scatter — p<0.05 in {sig_p}/60 "
      f"(mean-excess only), but CALLED (p<0.05 & consistency>="
      f"{MIN_CONSIST}) in {false_calls}/60 (want ~<=3)")
assert false_calls <= 4, f"call gate failed: {false_calls}/60"
print("aggregate specificity-to-consistency: PASS")
print("ALL P3.4 CHECKS PASSED")

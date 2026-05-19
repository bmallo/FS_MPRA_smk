"""P3.5 regression test — §2.5 secondary-mutation control.

Synthetic ground truth, no rd/network. Validates:
  (A) _secondary_free_mask: a read is dropped iff it carries a raw SNV
      inside the site-2 window (excluding the instrument's own id);
      out-of-window / excluded / unparseable handled.
  (B) ARTIFACT: the distal site-2 loss is carried ENTIRELY by reads
      bearing a recurrent site-2-local secondary SNV -> the
      secondary-free re-test collapses (not significant) -> the call
      does NOT survive -> flagged secondary_artifact.
  (C) GENUINE: site-2 loss spread across all reads with no site-2-local
      secondary mutation -> secondary-free re-test still significant &
      same direction -> validated_call.

Run: pixi run python tools/test_cooccupancy_secondary.py
"""
import sys
import numpy as np

sys.path.insert(
    0, "/mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/"
       "FS_MPRA_smk/workflow/scripts")
import common as C  # noqa

# (A) _secondary_free_mask --------------------------------------------
raw = [
    (),                                  # keep
    ("3405:t>C",),                       # site1 id (excluded) -> keep
    ("110:a>G",),                        # in window [90,130] -> drop
    ("200:c>T",),                        # out of window -> keep
    ("95:+T", "300:g>A"),                # 95 in window -> drop
    ("bogus",),                          # unparseable -> keep
]
mask = C._secondary_free_mask(raw, 90, 130, exclude_ids=("3405:t>C",))
assert mask.tolist() == [True, True, False, True, False, True], \
    mask.tolist()
print("_secondary_free_mask: PASS")

P1, P0 = 0.60, 0.20
B, K, NK = 2500, 6, 400
CALL_Q = C.DEFAULT_COOCCUPANCY_CFG['call_q']
MC = C.DEFAULT_COOCCUPANCY_CFG['min_consistency']


def states(rng, n, p_o1, p2_1, p2_0):
    o1 = rng.random(n) < p_o1
    o2 = rng.random(n) < np.where(o1, p2_1, p2_0)
    return o1, o2


def wt(rng):
    o1w, o2w = states(rng, 5000, 0.5, P1, P0)
    return C._wt_conditional_from_states(o1w, o2w, B, rng, 25)


def survives(full, sf):
    """Replicates the driver's §2.5 survival rule."""
    if sf is None:
        return False
    return (sf['p_two_sided'] < CALL_Q
            and sf['frac_instruments_consistent'] >= MC
            and np.sign(sf['weighted_mean_excess'])
            == np.sign(full['weighted_mean_excess']))


# (B) ARTIFACT: site-2 loss only on secondary-SNV-bearing reads -------
rng = np.random.default_rng(5)
artifact_flagged = 0
for _ in range(40):
    wc = wt(rng)
    full_inst, sf_inst = [], []
    for _k in range(K):
        # baseline: follows WT channel (no excess)
        o1, o2 = states(rng, NK, 0.12, P1, P0)
        # 25% of reads carry a recurrent site-2-local secondary SNV
        # and have site2 forced OFF -> drives a spurious distal loss
        sec = rng.random(NK) < 0.25
        o2 = o2 & ~sec
        full_inst.append((o1, o2))
        keep = ~sec                       # secondary-free subset
        if keep.sum() >= 25:
            sf_inst.append((o1[keep], o2[keep]))
    full = C._cooccupancy_aggregate_from_states(full_inst, wc, rng)
    sf = (C._cooccupancy_aggregate_from_states(sf_inst, wc, rng)
          if len(sf_inst) >= 3 else None)
    is_call = (full['p_two_sided'] < CALL_Q
               and full['frac_instruments_consistent'] >= MC)
    if is_call and not survives(full, sf):
        artifact_flagged += 1
print(f"ARTIFACT: flagged {artifact_flagged}/40 "
      f"(call collapses on secondary-free re-test)")
assert artifact_flagged >= 32, artifact_flagged
print("artifact detection: PASS")

# (C) GENUINE: cooperative loss across ALL reads, no secondary --------
rng = np.random.default_rng(8)
validated = 0
for _ in range(40):
    wc = wt(rng)
    full_inst, sf_inst = [], []
    for _k in range(K):
        o1, o2 = states(rng, NK, 0.12, P1, 0.03)   # real coop, all reads
        full_inst.append((o1, o2))
        sf_inst.append((o1, o2))                    # no secondary -> all kept
    full = C._cooccupancy_aggregate_from_states(full_inst, wc, rng)
    sf = C._cooccupancy_aggregate_from_states(sf_inst, wc, rng)
    is_call = (full['p_two_sided'] < CALL_Q
               and full['frac_instruments_consistent'] >= MC)
    if is_call and survives(full, sf):
        validated += 1
print(f"GENUINE: validated {validated}/40 (survives secondary screen)")
assert validated >= 32, validated
print("genuine survival: PASS")
print("ALL P3.5 CHECKS PASSED")

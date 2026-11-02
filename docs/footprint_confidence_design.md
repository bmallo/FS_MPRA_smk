# Footprint confidence calls: design note (2026-11-01)

Status: idea for the MPRA Fiber-seq pipeline revamp, to be built once Phase B data (naked-DNA
dose series B1, WT time course B2) exist. A display-only approximation is in the promoter explorer
(`tools/analysis_v2/explorer.html`, `explorer_data.py`).

## Idea

A footprint call is worth trusting only if the pattern of unmarked adenines it rests on is unlikely
when nothing is bound. Two levels:

1. **Per call (one molecule).** P_chance = probability that the call's unmarked adenines are all
   unmarked by chance, given that molecule's labeling. With first-order kinetics,
   p_j(d) = 1 − (1 − p_j,ref)^(d/d_ref), a molecule's dose d estimated from its backbone labeling L
   as d ∝ −ln(1 − L). Then log P_chance = (d/d_ref) · Σ_j log(1 − p_j,ref) over the unmarked adenines
   in the call. Report −log10 P_chance as a call quality.
2. **Per position and size (population).** Expected chance coverage E (the caller run on molecules
   simulated with no protein, at the real molecules' labeling), observed coverage O,
   local FDR = min(1, E / O), occupancy = (O − E) / spike-in detection (A8). Positions where E stays
   high at every labeling are "not resolvable", not "no footprint". Variant tests use confident
   positions only, or weight by 1 − local FDR.

The same E map computed from sequence alone is the sequence-information QC tool (deferred earlier).

## Evidence it separates real from chance (HepG2 run 3, exploratory)

| Labeling | 0.4–0.5 | 0.6–0.7 | 0.7–0.8 |
|---|---|---|---|
| 60–86 core, P(all 8 supporting adenines unmarked) | 0.10% | 0.0005% | 0.00001% |
| 114–138 call, P(all 4 supporting adenines unmarked) | 7.1% | 0.96% | 0.22% |
| local FDR at 73 (60–86) | 0.02 | 1e-4 | < 1e-4 |
| local FDR at 123 (114–138) | 1 | 0.43 | 0.35 |

(results/analysis_v2_20261028/fp60_86, fp114_138.)

## What it needs

- Null marking rates p_j,ref: B1 naked DNA at several Hia5 doses, per variant sequence where possible.
  Interim: in-cell M-molecule rates (conservative: protein lowers p_a at footprints, inflating chance)
  or naked-derived context tables (FiberHMM hexamer, A4 run-class lookup).
- The dose model and within-molecule independence: test with B2 (time course) and B1.
- Exact chance coverage: simulate the caller (it tolerates occasional marked adenines and scans many
  windows, so the product above is a lower bound on chance).
- Validation: naked controls give local FDR ≈ 1 everywhere; A8 spike-ins come back confident;
  60–86 confident, 114–138 and 182 not.

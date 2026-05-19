#!/usr/bin/env python
"""P2.9 — prove the WT-hoist + persistent-pool refactor is numerically
identical to the pre-refactor code (git HEAD b9fe79e).

Loads the OLD common.py (from git, /tmp/common_old.py) and the NEW
working-tree common.py as separate modules, parses a BAM once, and for
several variants/strata compares run_null_calibration output:

  OLD: common_old.run_null_calibration(... pool-per-call, internal build)
  NEW: common.run_null_calibration(... wt_ctx + persistent executor)

Asserts null_delta, wt_occ, pos_thresh, null_familywise_max,
null_max_cluster_sums, null_cluster_sums are bitwise-equal. The
refactor is purely build-once/reuse vs rebuild-each, so equivalence is
dataset-agnostic; the fast dev BAM is sufficient and deterministic.
"""
import importlib.util
import sys
import numpy as np

REPO = "/mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/FS_MPRA_smk"
sys.path.insert(0, f"{REPO}/workflow/scripts")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


NEW = _load("common", f"{REPO}/workflow/scripts/common.py")
OLD = _load("common_old", "/tmp/common_old.py")

BAM = f"{REPO}/results/phase0_hmm/variants/LDLR_hmm.tagged.bam"
PROM = (3184, 3502)               # internal half-open (== --promoter-region 3184-3501)
N_ITERS = 800
WORKERS = 4
MIN_READS = 20
SEED = 42

bins = NEW.build_analysis_bins()
rd, ref_len, ref_name, _ = NEW.parse_bam(
    BAM, bins, target_chrom=None, analysis_region=PROM)
wt_idx = rd.wt_indices
groups = NEW.group_variants(rd, include_multi=False, min_reads=MIN_READS)
print(f"WT={len(wt_idx)}  variants={len(groups)}  iters={N_ITERS} "
      f"workers={WORKERS}")

# Mirror run_variant_testing_parallel's stratum loop (stratify OFF =
# one stratum per variant); test the first few.
strata = NEW.build_strata(rd, groups, stratify=False)
test_strata = strata[:6]

# NEW path: build context + persistent pool ONCE (as the refactor does).
from concurrent.futures import ProcessPoolExecutor
wt_ctx = NEW.build_wt_null_context(rd, wt_idx, None, None, PROM,
                                   make_shared=True)
ex = ProcessPoolExecutor(max_workers=WORKERS)

ok = True
try:
    for si, st in enumerate(test_strata):
        rs = NEW.stable_variant_seed(SEED, f"stratum_{si}")
        common_kw = dict(
            min_nuc=None, max_nuc=None, n_iterations=N_ITERS,
            random_seed=rs, analysis_region=PROM,
            cluster_threshold_quantile=0.95,
            absolute_delta_threshold=None,
            gap_tolerance=2, merge_distance=5, n_workers=WORKERS)
        old_res = OLD.run_null_calibration(
            rd, wt_idx, st['rep_nc'], st['rep_n'], **common_kw)
        new_res = NEW.run_null_calibration(
            rd, wt_idx, st['rep_nc'], st['rep_n'],
            wt_ctx=wt_ctx, executor=ex, **common_kw)

        checks = []
        checks.append(("null_delta",
                       np.array_equal(old_res['null_delta'],
                                      new_res['null_delta'])))
        checks.append(("null_familywise_max",
                       np.array_equal(old_res['null_familywise_max'],
                                      new_res['null_familywise_max'])))
        for lbl in rd.bin_labels:
            checks.append((f"wt_occ[{lbl}]",
                           np.array_equal(old_res['wt_occ'][lbl],
                                          new_res['wt_occ'][lbl])))
            checks.append((f"pos_thresh[{lbl}]",
                           np.array_equal(old_res['pos_thresh'][lbl],
                                          new_res['pos_thresh'][lbl])))
            checks.append((f"null_max_cluster_sums[{lbl}]",
                           np.array_equal(
                               old_res['null_max_cluster_sums'][lbl],
                               new_res['null_max_cluster_sums'][lbl])))
            checks.append((f"null_cluster_sums[{lbl}]",
                           np.array_equal(
                               old_res['null_cluster_sums'][lbl],
                               new_res['null_cluster_sums'][lbl])))
        bad = [n for n, good in checks if not good]
        status = "OK" if not bad else f"MISMATCH: {bad}"
        if bad:
            ok = False
        print(f"stratum {si} N={st['rep_n']:>4}  {len(checks)} arrays  "
              f"{status}")
finally:
    ex.shutdown()
    NEW.cleanup_shared_memory(wt_ctx['shm_objects'])

print("\n" + ("ALL EQUAL — refactor is numerically identical"
               if ok else "EQUIVALENCE FAILED"))
sys.exit(0 if ok else 1)

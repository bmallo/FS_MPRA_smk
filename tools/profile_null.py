"""Decisive profile: where does run_null_calibration's per-variant time
go? The post-pool family-wise cluster loop is dataset-size-independent
(depends on n_iterations x bins x analysis_length, not WT size), so the
fast dev BAM reproduces that cost. If one null call at 10000 iters is
~hundreds of s even with tiny dev WT, the serial post-pool loop is the
bottleneck (not the WT setup my hoist addressed)."""
import cProfile
import io
import pstats
import sys
import time

REPO = "/mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/FS_MPRA_smk"
sys.path.insert(0, f"{REPO}/workflow/scripts")
import common as C  # noqa

BAM = f"{REPO}/results/phase0_hmm/variants/LDLR_hmm.tagged.bam"
PROM = (3184, 3502)
bins = C.build_analysis_bins()
rd, *_ = C.parse_bam(BAM, bins, target_chrom=None, analysis_region=PROM)
wt_idx = rd.wt_indices
groups = C.group_variants(rd, include_multi=False, min_reads=20)
strata = C.build_strata(rd, groups, stratify=False)
st = strata[0]
print(f"WT={len(wt_idx)} dev; one stratum N={st['rep_n']}, "
      f"iters=10000, workers=8")

from concurrent.futures import ProcessPoolExecutor
ctx = C.build_wt_null_context(rd, wt_idx, None, None, PROM, make_shared=True)
ex = ProcessPoolExecutor(max_workers=8)

pr = cProfile.Profile()
t0 = time.time()
pr.enable()
C.run_null_calibration(
    rd, wt_idx, st['rep_nc'], st['rep_n'],
    n_iterations=10000, random_seed=42, analysis_region=PROM,
    n_workers=8, wt_ctx=ctx, executor=ex)
pr.disable()
dt = time.time() - t0
ex.shutdown()
C.cleanup_shared_memory(ctx['shm_objects'])

print(f"\n=== one run_null_calibration: {dt:.1f} s "
      f"(WT tiny, so this ~= the dataset-independent serial cost) ===")
s = io.StringIO()
pstats.Stats(pr, stream=s).sort_stats("cumulative").print_stats(18)
print("\n".join(s.getvalue().splitlines()[:30]))

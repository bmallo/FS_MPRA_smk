"""Production-scale section timing. The dev profile was blind to
WT-scale costs. This parses the PRODUCTION BAM, builds the hoisted
context + persistent pool once, then times two consecutive
run_null_calibration calls (mimics variant 1 vs variant 2..N). If
call #2 is fast, the persistent-pool optimization works and variant 1
was just cold-pool + one-time build. If call #2 is also slow, the
per-variant production bottleneck is elsewhere and is shown here.
"""
import sys
import time
import cProfile
import io
import pstats
from concurrent.futures import ProcessPoolExecutor

REPO = "/mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/FS_MPRA_smk"
sys.path.insert(0, f"{REPO}/workflow/scripts")
import common as C  # noqa

BAM = f"{REPO}/results/phase0_hmm_full/variants_snv/LDLR_hmm_full.tagged.bam"
PROM = (3184, 3502)
WORKERS = 64

t = time.time()
bins = C.build_analysis_bins()
rd, *_ = C.parse_bam(BAM, bins, target_chrom=None, analysis_region=PROM)
print(f"parse_bam: {time.time()-t:.1f}s  WT={len(rd.wt_indices)}")
wt_idx = rd.wt_indices
groups = C.group_variants(rd, include_multi=False, min_reads=50)
strata = C.build_strata(rd, groups, stratify=False)
print(f"variants={len(groups)} strata={len(strata)}")

t = time.time()
ctx = C.build_wt_null_context(rd, wt_idx, None, None, PROM, make_shared=True)
print(f"build_wt_null_context (one-time): {time.time()-t:.1f}s")
t = time.time()
ex = ProcessPoolExecutor(max_workers=WORKERS)
# warm the pool (first .submit triggers worker spawn)
list(ex.map(int, range(WORKERS)))
print(f"pool create+warm (one-time): {time.time()-t:.1f}s")


def one(si, profile=False):
    st = strata[si]
    rs = C.stable_variant_seed(42, f"stratum_{si}")
    kw = dict(min_nuc=None, max_nuc=None, n_iterations=10000,
              random_seed=rs, analysis_region=PROM,
              cluster_threshold_quantile=0.95,
              absolute_delta_threshold=None, gap_tolerance=2,
              merge_distance=5, n_workers=WORKERS,
              wt_ctx=ctx, executor=ex)
    if profile:
        pr = cProfile.Profile()
        pr.enable()
    t0 = time.time()
    C.run_null_calibration(rd, wt_idx, st['rep_nc'], st['rep_n'], **kw)
    dt = time.time() - t0
    if profile:
        pr.disable()
        s = io.StringIO()
        pstats.Stats(pr, stream=s).sort_stats("cumulative").print_stats(12)
        print("\n".join(s.getvalue().splitlines()[4:22]))
    return dt, st['rep_n']


for si in range(4):
    dt, n = one(si, profile=(si == 1))
    print(f">>> call {si} (warm pool, N={n}): {dt:.1f}s")

ex.shutdown()
C.cleanup_shared_memory(ctx['shm_objects'])

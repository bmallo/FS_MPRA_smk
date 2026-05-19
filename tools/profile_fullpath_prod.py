"""Decisive end-to-end check: time the FULL run_variant_testing_parallel
per-variant path (null + nc_shift_null + _compute_variant_result +
nc_shift_stats) at PRODUCTION scale on the first N variants. Catches any
remaining WT-scale serial cost before committing the cluster to a run.
"""
import sys
import time

REPO = "/mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/FS_MPRA_smk"
sys.path.insert(0, f"{REPO}/workflow/scripts")
import common as C  # noqa

BAM = f"{REPO}/results/phase0_hmm_full/variants_snv/LDLR_hmm_full.tagged.bam"
PROM = (3184, 3502)
WORKERS = 64
N_VARIANTS = 5

t = time.time()
bins = C.build_analysis_bins()
rd, *_ = C.parse_bam(BAM, bins, target_chrom=None, analysis_region=PROM)
print(f"parse_bam: {time.time()-t:.1f}s  WT={len(rd.wt_indices)}")
groups = C.group_variants(rd, include_multi=False, min_reads=50)
print(f"variants={len(groups)}; timing first {N_VARIANTS} via "
      f"run_variant_testing_parallel")

# Trim to the first N variants (stable order) — same code path,
# just fewer strata so this finishes fast.
keys = list(groups.keys())[:N_VARIANTS]
sub = {k: groups[k] for k in keys}

t = time.time()
res = C.run_variant_testing_parallel(
    rd, rd.wt_indices, sub, PROM, PROM[0], PROM[1],
    random_seed=42, n_null_iterations=10000,
    n_workers=WORKERS, stratify=False)
dt = time.time() - t
print(f"\n{N_VARIANTS} variants full path: {dt:.1f}s  "
      f"=> {dt/N_VARIANTS:.1f}s/variant")
print(f"projected 149 variants: {dt/N_VARIANTS*149/60:.1f} min "
      f"(+ ~3 min parse)")
print(f"tested={len(res)}  sample keys={keys[:3]}")

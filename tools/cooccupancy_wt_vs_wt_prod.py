"""P3.6 production blocking gate — WT-vs-WT co-occupancy calibration
on REAL LDLR data + REAL footprint structure.

Parses the production tagged BAM (-> rd) and rebuilds the canonical
co-occupancy sites from the production Phase-2 HDF5 (FDR-significant
clusters + motif calls, merged per bin — same pooling as
build_canonical_sites). Then runs cooccupancy_wt_vs_wt over the real
site geometry. PASS = pooled p ~Uniform, BH calls @ call_q &
consistency ~ 0 (no false dependencies under H0). A positive-control
arm (inject_dependency) confirms teeth on real data.
"""
import argparse
import sys
import numpy as np
import h5py

sys.path.insert(
    0, "/mmfs1/gscratch/stergachislab/bmallo/large_home/git_repos/"
       "FS_MPRA_smk/workflow/scripts")
import common as C  # noqa


def canonical_sites_from_h5(h5_path, bin_labels, motif_fdr):
    """Pool FDR-significant clusters + motif calls per bin and merge —
    the build_canonical_sites pooling, reconstructed from the HDF5 so
    we don't re-run the ~11 min Stage 3."""
    with h5py.File(h5_path, 'r') as f:
        s = f['summary']
        svid = [v.decode() for v in s['variant_ids'][:]]
        sq = s['variant_fdr_q'][:]
        fdr_ok = {vid for vid, q in zip(svid, sq) if q < motif_fdr}
        by_bin = {b: [] for b in bin_labels}
        if 'clusters' in f:
            c = f['clusters']
            cvid = [v.decode() for v in c['variant_ids'][:]]
            cbin = [v.decode() for v in c['bin_labels'][:]]
            cs = c['abs_start'][:]
            ce = c['abs_end'][:]
            for vid, b, a0, a1 in zip(cvid, cbin, cs, ce):
                if vid in fdr_ok:
                    by_bin.setdefault(b, []).append((int(a0), int(a1)))
        motif_iv = set()
        if 'motifs' in f and 'calls' in f['motifs']:
            for k in f['motifs']['calls']:
                g = f['motifs']['calls'][k]
                b = g.attrs['bin']
                b = b.decode() if isinstance(b, bytes) else b
                a0, a1 = int(g.attrs['abs_start']), int(g.attrs['abs_end'])
                by_bin.setdefault(b, []).append((a0, a1))
                motif_iv.add((b, a0, a1))
    sites = []
    for b, ivs in by_bin.items():
        for (s0, s1) in C._merge_bin_intervals(ivs, 0.5):
            sites.append({'site_id': f"{b}:{s0}-{s1}", 'bin': b,
                          'abs_start': int(s0), 'abs_end': int(s1),
                          'is_motif_site': any(
                              mb == b and C._recip_overlap(
                                  s0, s1, ms, me) > 0
                              for (mb, ms, me) in motif_iv)})
    return sites


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam', required=True)
    ap.add_argument('--h5', required=True)
    ap.add_argument('--promoter-region', default='3184-3501')
    ap.add_argument('--n-instruments', type=int, default=6)
    ap.add_argument('--instr-n', type=int, default=60)
    ap.add_argument('--n-iter', type=int, default=4000)
    ap.add_argument('--seeds', type=int, default=15)
    ap.add_argument('--max-pairs', type=int, default=150)
    ap.add_argument('--motif-fdr', type=float, default=0.10)
    a = ap.parse_args()

    ps = a.promoter_region.split('-')
    region = (int(ps[0]), int(ps[1]) + 1)
    bins = C.build_analysis_bins()
    bin_labels = [b[0] for b in bins]
    rd, *_ = C.parse_bam(a.bam, bins, target_chrom=None,
                         analysis_region=region)
    print(f"WT reads: {len(rd.wt_indices):,}")
    sites = canonical_sites_from_h5(a.h5, bin_labels, a.motif_fdr)
    print(f"canonical sites from HDF5: {len(sites)}  "
          f"bins={sorted(set(s['bin'] for s in sites))}")
    cfg = C.build_cooccupancy_cfg()

    def sweep(inject):
        recs = []
        for sd in range(a.seeds):
            recs += C.cooccupancy_wt_vs_wt(
                rd, sites, cfg, a.n_instruments, a.instr_n,
                a.n_iter, seed=sd, max_pairs=a.max_pairs,
                inject_dependency=inject)
        return recs

    recs = sweep(0.0)
    if not recs:
        print("NO testable WT-vs-WT records — check sites/WT depth")
        sys.exit(2)
    p = np.array([r['p_two_sided'] for r in recs])
    q = C.benjamini_hochberg(p)
    mc = cfg['min_consistency']
    calls = sum(1 for r, qq in zip(recs, q)
                if qq < cfg['call_q']
                and r['frac_instruments_consistent'] >= mc)
    print("=" * 64)
    print(f"H0 WT-vs-WT (REAL data): n={len(recs)}  "
          f"mean(p)={p.mean():.3f} (~0.5)  "
          f"FP@0.05={np.mean(p < 0.05):.3f}  "
          f"FP@0.10={np.mean(p < 0.10):.3f}")
    print(f"  BH calls @ q<{cfg['call_q']} & consistency>={mc}: "
          f"{calls}/{len(recs)}")
    verdict = (0.40 <= p.mean() <= 0.60
               and np.mean(p < 0.05) <= 0.12
               and calls <= max(1, int(0.03 * len(recs))))
    print(f"  CALIBRATION VERDICT: "
          f"{'PASS' if verdict else 'FAIL'}")

    pc = sweep(0.7)
    if pc:
        ppc = np.array([r['p_two_sided'] for r in pc])
        print(f"positive control (inject=0.7): p<0.05 = "
              f"{np.mean(ppc < 0.05):.2f} (teeth — want high)")
    print("=" * 64)
    sys.exit(0 if verdict else 1)


if __name__ == '__main__':
    main()

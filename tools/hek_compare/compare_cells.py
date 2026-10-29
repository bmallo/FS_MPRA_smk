"""WT footprints in HEK293T vs HepG2 (replicate 3) over the shared 211-bp LDLR element.

Both cells were run through the same steps: WT analysis-eligible molecules, U/M/C proxy classes
on distal m6A, FiberHMM 2.12.0 linear re-call, plasmid-LLR (global, drawdown 25, per-position
table from each cell's own U/M classes, fold A). Comparison is on C-class reads.

Element frames (oligo o = 21..231):
  HepG2 rep3  LDLR_MPRA_NotI  ref0 = 1234 + (o - 21)            (element on + strand)
  HEK293T     LDLR (older)    ref0 = 3446 - (o - 21)            (element on - strand)

Threshold rule, the same as tools/llr_compare/compare_callers.py and fixed before looking at
chromatin reads: per cell, the primary min_llr is the smallest grid value whose fold-B M-class
false-positive bp rate over element +/-100 is <= that cell's FiberHMM 2.12.0. min_llr 5 is also
reported for both, so one shared threshold is always available.

Labeling-level control: HepG2 C reads are also reweighted to the HEK293T C-class distal-m6A
histogram (0.05 bins), since occupancy calls depend on how heavily a molecule was labeled.

Usage: python tools/hek_compare/compare_cells.py --out results/hek293t_vs_hepg2_20261028
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'llr_compare'))
import compare_callers as cc  # noqa: E402  (fiberhmm_calls, llr_calls, llrlib)

llrlib = cc.llrlib
OLIGO = np.arange(21, 232)
GRID = cc.GRID
CLASSES = cc.CLASSES
CELLS = {
    'HepG2': dict(workdir=Path('results/llr_vs_fiberhmm_20261026'), L=4647, ref0=1234 + (OLIGO - 21), elem100=(1134, 1545),
                  molecules=Path('results/rep3_site_methylation_20261021/extract/molecules.tsv.gz'), id_col='read_id'),
    'HEK293T': dict(workdir=Path('results/hek293t_llr_20261028'), L=4718, ref0=3446 - (OLIGO - 21), elem100=(3136, 3547),
                    molecules=Path('results/hek293t_llr_20261028/data/molecules.tsv.gz'), id_col='read_id'),
}
COLOR = {'HepG2': '#2B6CB0', 'HEK293T': '#C8641E'}
REGIONS = {'P1 core': (54, 62), 'R1 Sp1': (89, 104), 'R2 SRE1': (124, 139), 'R3 Sp1': (140, 155)}  # oligo, inclusive
PQ = 'llr_global_foldA_pp0_dd25'
DISTAL_BINS = np.arange(0.01, 0.6001, 0.05)
N_BOOT = 200
SEED = 20261028


def protected_matrix(calls, reads, ref0, L):
    """reads x 211 bool: oligo position protected; reads are full length so all positions are covered."""
    out = np.zeros((len(reads), len(ref0)), bool)
    for i, r in enumerate(reads):
        m = np.zeros(L, bool)
        for lo, hi in calls[r][2]:
            m[max(lo, 0):min(hi, L)] = True
        out[i] = m[ref0]
    return out


def fp_rate(calls, reads, region, L):
    a, b = region
    prot = cov = 0
    for r in reads:
        lo, hi, ivs = calls[r]
        m = np.zeros(L, bool)
        for s, e in ivs:
            m[max(s, 0):min(e, L)] = True
        c = np.zeros(L, bool); c[max(lo, 0):min(hi, L)] = True
        prot += int((m & c)[a:b].sum()); cov += int(c[a:b].sum())
    return prot / max(cov, 1)


def element_sizes(calls, reads, ref0):
    lo_e, hi_e = int(ref0.min()), int(ref0.max()) + 1
    counts = np.zeros(len(CLASSES))
    for r in reads:
        for s, e in calls[r][2]:
            if s < hi_e and e > lo_e:
                ln = e - s
                k = next((i for i, (_, x, y) in enumerate(CLASSES) if x <= ln <= y), None)
                if k is not None:
                    counts[k] += 1
    return counts / max(len(reads), 1)


def weighted_occ(mat, w):
    return (w[:, None] * mat).sum(0) / w.sum()


def boot_diff(m1, w1, m2, w2, rng):
    """Read-level Poisson bootstrap of occupancy(cell2) - occupancy(cell1)."""
    d = np.empty((N_BOOT, m1.shape[1]))
    f1, f2 = m1.astype(np.float32), m2.astype(np.float32)
    for b in range(N_BOOT):
        p1 = w1 * rng.poisson(1.0, len(w1)); p2 = w2 * rng.poisson(1.0, len(w2))
        d[b] = (p2 @ f2) / p2.sum() - (p1 @ f1) / p1.sum()
    return np.percentile(d, [2.5, 97.5], axis=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', required=True)
    out = Path(ap.parse_args().out)
    (out / 'figures').mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    res, summary = {}, {}
    for cell, c in CELLS.items():
        W, L, ref0 = c['workdir'], c['L'], c['ref0']
        data, calls_dir = W / 'data', W / 'results/calls'
        fh = {s: cc.fiberhmm_calls(data / f'{s}.fiberhmm.bam') for s in ('wt_meth', 'wt_unmeth', 'wt_chrom')}
        spans = {s: {r: (lo, hi) for r, (lo, hi, _) in d.items()} for s, d in fh.items()}
        foldB = {s: [r for r in spans[s] if llrlib.read_fold(r) == 1] for s in ('wt_meth', 'wt_unmeth')}
        bench = [dict(cell=cell, caller='fh2', min_llr=None, fp_bp=fp_rate(fh['wt_meth'], foldB['wt_meth'], c['elem100'], L),
                      tp_bp=fp_rate(fh['wt_unmeth'], foldB['wt_unmeth'], c['elem100'], L))]
        for t in GRID:
            bench.append(dict(cell=cell, caller='llr_measured', min_llr=t,
                              fp_bp=fp_rate(cc.llr_calls(calls_dir / f'wt_meth.{PQ}.parquet', t, spans['wt_meth']), foldB['wt_meth'], c['elem100'], L),
                              tp_bp=fp_rate(cc.llr_calls(calls_dir / f'wt_unmeth.{PQ}.parquet', t, spans['wt_unmeth']), foldB['wt_unmeth'], c['elem100'], L)))
        bench = pd.DataFrame(bench)
        fh2_fp = float(bench.fp_bp.iloc[0])
        ok = bench[(bench.caller == 'llr_measured') & (bench.fp_bp <= fh2_fp)]
        primary = float(ok.min_llr.min()) if len(ok) else float(max(GRID))
        chrom = sorted(spans['wt_chrom'])
        mol = pd.read_csv(c['molecules'], sep='\t', usecols=[c['id_col'], 'distal_m6a_fraction'])
        distal = mol.set_index(c['id_col']).distal_m6a_fraction.reindex(chrom).values
        callers = {f'llr@{primary:g}': cc.llr_calls(calls_dir / f'wt_chrom.{PQ}.parquet', primary, spans['wt_chrom']),
                   'llr@5': cc.llr_calls(calls_dir / f'wt_chrom.{PQ}.parquet', 5.0, spans['wt_chrom']),
                   'fh2': fh['wt_chrom']}
        mats = {k: protected_matrix(v, chrom, ref0, L) for k, v in callers.items()}
        sizes = {k: element_sizes(v, chrom, ref0) for k, v in callers.items()}
        em = pd.read_csv(W / 'results/emissions.tsv.gz', sep='\t')
        em = em[em.fold == 'A'].set_index('ref_pos')
        pa = em.p_a.reindex(ref0 + 1).values; pp = em.p_p.reindex(ref0 + 1).values
        floor = pd.read_csv(W / 'results/detection_floor.tsv', sep='\t').iloc[0].to_dict()
        res[cell] = dict(mats=mats, sizes=sizes, distal=distal, pa=pa, pp=pp, primary=primary, bench=bench)
        summary[cell] = dict(primary_min_llr=primary, fh2_fold_b_fp=fh2_fp, chromatin_reads=len(chrom),
                             fold_b_reads={k: len(v) for k, v in foldB.items()}, detection_floor=floor,
                             distal_m6a_median_C=float(np.nanmedian(distal)),
                             protected_frac_element={k: float(m.mean()) for k, m in mats.items()})
        print(cell, json.dumps(summary[cell], default=str), flush=True)
    pd.concat([r['bench'] for r in res.values()]).to_csv(out / 'control_benchmark.tsv', sep='\t', index=False)

    # Weights: HepG2 C reads reweighted to the HEK293T distal histogram.
    h_hek = np.histogram(res['HEK293T']['distal'], DISTAL_BINS)[0] / len(res['HEK293T']['distal'])
    b_hep = np.clip(np.digitize(res['HepG2']['distal'], DISTAL_BINS) - 1, 0, len(DISTAL_BINS) - 2)
    h_hep = np.bincount(b_hep, minlength=len(DISTAL_BINS) - 1) / len(b_hep)
    w_matched = (h_hek / np.where(h_hep > 0, h_hep, np.nan))[b_hep]
    w_matched = np.nan_to_num(w_matched)
    ones = {cell: np.ones(len(r['distal'])) for cell, r in res.items()}
    p_hep, p_hek = res['HepG2']['primary'], res['HEK293T']['primary']
    compare = {'llr_primary': (f'llr@{p_hep:g}', f'llr@{p_hek:g}'), 'llr@5': ('llr@5', 'llr@5'), 'fh2': ('fh2', 'fh2')}
    tab = pd.DataFrame({'oligo': OLIGO, 'hepg2_ref0': CELLS['HepG2']['ref0'], 'hek293t_ref0': CELLS['HEK293T']['ref0'],
                        'hepg2_p_a': res['HepG2']['pa'], 'hek293t_p_a': res['HEK293T']['pa'],
                        'hepg2_p_p': res['HepG2']['pp'], 'hek293t_p_p': res['HEK293T']['pp']})
    for name, (k1, k2) in compare.items():
        m1, m2 = res['HepG2']['mats'][k1], res['HEK293T']['mats'][k2]
        tab[f'{name}_hepg2'] = m1.mean(0)
        tab[f'{name}_hepg2_matched'] = weighted_occ(m1, w_matched)
        tab[f'{name}_hek293t'] = m2.mean(0)
        tab[f'{name}_diff'] = tab[f'{name}_hek293t'] - tab[f'{name}_hepg2']
        tab[f'{name}_diff_matched'] = tab[f'{name}_hek293t'] - tab[f'{name}_hepg2_matched']
        lo, hi = boot_diff(m1, w_matched, m2, ones['HEK293T'], rng)
        tab[f'{name}_diff_matched_lo'], tab[f'{name}_diff_matched_hi'] = lo, hi
    tab.to_csv(out / 'occupancy_oligo.tsv', sep='\t', index=False)
    sizes = pd.DataFrame([dict(cell=cell, caller=k, **{cl[0]: v[i] for i, cl in enumerate(CLASSES)})
                          for cell, r in res.items() for k, v in r['sizes'].items()])
    sizes.to_csv(out / 'element_size_classes.tsv', sep='\t', index=False)
    reg = []
    for rname, (a, b) in REGIONS.items():
        sel = (OLIGO >= a) & (OLIGO <= b)
        row = dict(region=rname, oligo=f'{a}-{b}')
        for name in compare:
            row[f'{name}_hepg2_matched'] = float(tab[f'{name}_hepg2_matched'][sel].mean())
            row[f'{name}_hek293t'] = float(tab[f'{name}_hek293t'][sel].mean())
            row[f'{name}_diff_matched'] = float(tab[f'{name}_diff_matched'][sel].mean())
        reg.append(row)
    pd.DataFrame(reg).to_csv(out / 'region_occupancy.tsv', sep='\t', index=False)
    ok = np.isfinite(tab.hepg2_p_a) & np.isfinite(tab.hek293t_p_a)
    summary['emissions'] = dict(at_sites=int(ok.sum()), pearson_p_a=float(np.corrcoef(tab.hepg2_p_a[ok], tab.hek293t_p_a[ok])[0, 1]),
                                median_p_a=dict(HepG2=float(tab.hepg2_p_a[ok].median()), HEK293T=float(tab.hek293t_p_a[ok].median())),
                                dead_sites_below_0_10=dict(HepG2=int((tab.hepg2_p_a[ok] < 0.10).sum()), HEK293T=int((tab.hek293t_p_a[ok] < 0.10).sum())))
    for name in compare:
        d = tab[f'{name}_diff_matched']
        summary[f'diff_{name}'] = dict(pearson_occupancy=float(np.corrcoef(tab[f'{name}_hepg2_matched'], tab[f'{name}_hek293t'])[0, 1]),
                                       mean_diff=float(d.mean()), max_abs_diff=float(d.abs().max()),
                                       max_abs_diff_oligo=int(tab.oligo[d.abs().idxmax()]),
                                       positions_ci_excludes_0=int(((tab[f'{name}_diff_matched_lo'] > 0) | (tab[f'{name}_diff_matched_hi'] < 0)).sum()),
                                       positions_abs_diff_ge_0_10=int((d.abs() >= 0.10).sum()))
    summary['weights'] = dict(bins=DISTAL_BINS.tolist(), hek293t_hist=h_hek.tolist(), hepg2_hist=h_hep.tolist(),
                              hepg2_effective_n=float(w_matched.sum() ** 2 / (w_matched ** 2).sum()))
    json.dump(summary, open(out / 'summary.json', 'w'), indent=1, default=str)

    # Figure: occupancy and difference, primary LLR and FiberHMM 2.12.0.
    plt.rcParams.update({'figure.facecolor': cc.BG, 'axes.facecolor': cc.BG, 'savefig.facecolor': cc.BG})
    fig, axes = plt.subplots(3, 1, figsize=(13, 10), sharex=True, gridspec_kw=dict(height_ratios=[3, 2, 2]))
    ax = axes[0]
    for cell, key in (('HepG2', 'llr_primary_hepg2_matched'), ('HEK293T', 'llr_primary_hek293t')):
        ax.plot(OLIGO, tab[key], color=COLOR[cell], lw=2, label=f'{cell} (LLR, measured table)')
    for cell, key in (('HepG2', 'fh2_hepg2_matched'), ('HEK293T', 'fh2_hek293t')):
        ax.plot(OLIGO, tab[key], color=COLOR[cell], lw=1.2, ls='--', label=f'{cell} (FiberHMM 2.12.0)')
    for rname, (a, b) in REGIONS.items():
        ax.axvspan(a, b + 1, color=cc.GRIDC, alpha=0.5, lw=0)
        ax.text((a + b + 1) / 2, 1.0, rname, ha='center', va='bottom', fontsize=10, color=cc.MUTED)
    ax.set(ylim=(0, 1.05), ylabel='fraction protected', title='WT occupancy, C-class reads (HepG2 matched to HEK293T labeling)')
    ax.legend(loc='lower left', fontsize=10, ncol=2, frameon=False)
    for ax, name, title in ((axes[1], 'llr_primary', 'HEK293T − HepG2, LLR measured table'),
                            (axes[2], 'fh2', 'HEK293T − HepG2, FiberHMM 2.12.0')):
        ax.fill_between(OLIGO, tab[f'{name}_diff_matched_lo'], tab[f'{name}_diff_matched_hi'], color=cc.GRIDC, lw=0)
        ax.plot(OLIGO, tab[f'{name}_diff_matched'], color=cc.INK, lw=1.5)
        ax.axhline(0, color=cc.MUTED, lw=0.8)
        ax.set(ylabel='Δ occupancy', title=title)
    axes[2].set_xlabel('oligo position (shared 211-bp element, 21–231)')
    fig.tight_layout()
    fig.savefig(out / 'figures/occupancy_hek_vs_hepg2.png', dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(tab.hepg2_p_a, tab.hek293t_p_a, s=14, color=cc.INK, alpha=0.7)
    ax.plot([0, 1], [0, 1], color=cc.MUTED, lw=0.8)
    ax.set(xlabel='HepG2 P(m6A | accessible), M class', ylabel='HEK293T P(m6A | accessible), M class',
           title='Per-site labeling, shared element', xlim=(0, 1), ylim=(0, 1))
    fig.tight_layout(); fig.savefig(out / 'figures/emissions_hek_vs_hepg2.png', dpi=200); plt.close(fig)
    print(json.dumps({k: v for k, v in summary.items() if k != 'weights'}, indent=1, default=str))


if __name__ == '__main__':
    main()

"""Slide-sized figures for the FiberHMM vs plasmid-LLR comparison.

Reads the outputs of compare_callers.py and the measured emission table; computes FiberHMM's
hexamer emission probabilities at each plasmid A/T with FiberHMM's own encoder.

Usage (FS_MPRA_smk root, plasmid-llr python env):
  python tools/llr_compare/deck_figures.py --workdir results/llr_vs_fiberhmm_20261026
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

sys.path.insert(0, '/mmfs1/gscratch/stergachislab/bmallo/large_home/code/mine/plasmid-llr/llr')
import llrlib  # noqa: E402

ELEM0 = 1234                      # element start, 0-based plasmid coordinate; oligo = ref0 - 1234 + 21
OLIGO = lambda ref0: ref0 - ELEM0 + 21
REGIONS = {'P1 core': (54, 62), 'R1 Sp1': (89, 104), 'R2 SRE1': (124, 139), 'R3 Sp1': (140, 155)}  # oligo, inclusive
CALLERS = [('fh2', 'FiberHMM 2.12.0', '#2B6CB0', '-'), ('fh3', 'FiberHMM 3.0.2', '#1B2A41', '--'),
           ('llr_kmer@5', 'LLR, hexamer table', '#E9A36B', '-'), ('llr_measured@5', 'LLR, measured table', '#C8641E', '-')]
BG, INK, MUTED, GRID = '#F7F6F1', '#1B2A41', '#4A5568', '#DCDDD5'
BLUE, ORANGE, GRAY = '#2B6CB0', '#C8641E', '#8A94A3'
plt.rcParams.update({'figure.facecolor': BG, 'axes.facecolor': BG, 'savefig.facecolor': BG, 'font.size': 14,
                     'axes.titlesize': 15, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'axes.labelcolor': MUTED,
                     'axes.edgecolor': GRID, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'text.color': INK, 'axes.titlecolor': INK, 'axes.grid': True, 'grid.color': GRID,
                     'grid.linewidth': 0.7, 'axes.axisbelow': True})


def hexamer_rates(seq):
    """FiberHMM hia5_pacbio P(m6A | accessible), P(m6A | protected) at each reference position."""
    import fiberhmm
    from fiberhmm.core.bam_reader import encode_from_query_sequence
    from fiberhmm.core.model_io import load_model
    from fiberhmm.inference.tf_recaller import N_CTX
    model = load_model(str(Path(fiberhmm.__file__).parent / 'models/hia5_pacbio.json'))
    ep = np.asarray(model.emissionprob_, float)  # state 0 protected, 1 accessible; columns = observation codes
    codes = encode_from_query_sequence(seq, set(range(len(seq))), edge_trim=0, mode='pacbio-fiber', context_size=3)
    ok = (codes >= 0) & (codes < N_CTX)
    pa, pp = np.full(len(seq), np.nan), np.full(len(seq), np.nan)
    pa[ok], pp[ok] = ep[1, codes[ok]], ep[0, codes[ok]]
    return pa, pp


def shade(ax, ymax=1.0):
    for name, (a, b) in REGIONS.items():
        ax.axvspan(a - 0.5, b + 0.5, color='#ECEAE2', zorder=0)
        ax.text((a + b) / 2, 0.965, name, ha='center', va='top', fontsize=11, color=MUTED, transform=ax.get_xaxis_transform())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workdir', required=True)
    W = Path(ap.parse_args().workdir)
    out = W / 'results/deck'
    out.mkdir(parents=True, exist_ok=True)
    seq, at_mask = llrlib.load_reference(str(W / 'ref/LDLR_MPRA_NotI.fa'))
    em = pd.read_csv(W / 'results/emissions.tsv.gz', sep='\t')
    em = em[em.fold == 'A']
    ref0 = em.ref_pos.values - 1                       # table is 1-based
    hk_a, hk_p = hexamer_rates(seq)
    pa_meas, pp_meas = em.p_a.values, em.p_p.values
    occ = pd.read_csv(W / 'results/comparison/chromatin_occupancy.tsv.gz', sep='\t', index_col=0)
    summ = pd.read_csv(W / 'results/comparison/chromatin_summary.tsv', sep='\t').set_index('caller')
    bench = pd.read_csv(W / 'results/comparison/control_benchmark.tsv', sep='\t')
    stats = {}

    # 1. Measured vs hexamer emission rates.
    fig, (a, b) = plt.subplots(1, 2, figsize=(16, 6.2), gridspec_kw={'wspace': 0.25})
    dead = pa_meas < 0.10
    a.scatter(hk_a[ref0][~dead], pa_meas[~dead], s=7, color=BLUE, alpha=0.5, linewidths=0, label='measured ≥ 0.10')
    a.scatter(hk_a[ref0][dead], pa_meas[dead], s=12, color=ORANGE, alpha=0.8, linewidths=0, label=f'measured < 0.10 (n = {dead.sum()})')
    a.plot([0, 1], [0, 1], color=GRAY, ls=':')
    a.set(xlim=(0, 1), ylim=(-0.02, 1.02), xlabel='FiberHMM hexamer P(m6A | accessible)', ylabel='measured on M-class molecules',
          title='Labeling rate when accessible, per A/T site')
    a.legend(frameon=False, fontsize=12, loc='upper left')
    ok = np.isfinite(hk_p[ref0])
    b.scatter(hk_p[ref0][ok], np.maximum(pp_meas[ok], 1e-5), s=7, color=GRAY, alpha=0.5, linewidths=0)
    b.plot([1e-5, 1], [1e-5, 1], color=GRAY, ls=':')
    b.set(xscale='log', yscale='log', xlim=(1e-4, 0.3), ylim=(8e-6, 0.5), xlabel='FiberHMM hexamer P(m6A | protected)',
          ylabel='measured on U-class molecules', title='False-call rate when protected, per A/T site')
    stats['emissions'] = dict(sites=int(len(pa_meas)), dead_measured=int(dead.sum()), dead_hexamer=int((hk_a[ref0] < 0.10).sum()),
                              median_pa_measured=float(np.median(pa_meas)), median_pa_hexamer=float(np.nanmedian(hk_a[ref0])),
                              median_pp_measured=float(np.median(pp_meas)), median_pp_hexamer=float(np.nanmedian(hk_p[ref0])),
                              r_pa=float(np.corrcoef(hk_a[ref0][ok], pa_meas[ok])[0, 1]))
    fig.savefig(out / 'deck_emissions.png', dpi=200, bbox_inches='tight'); plt.close(fig)

    # 2. Control benchmark: false-positive bp on M-class (held out) by caller.
    fig, ax = plt.subplots(figsize=(13, 5.2))
    for which, label, col in (('llr_measured', 'LLR, measured table', ORANGE), ('llr_kmer', 'LLR, hexamer table', '#E9A36B')):
        d = bench[bench.caller == which].sort_values('min_llr')
        ax.plot(d.min_llr, 100 * d.fp_bp, '-o', color=col, lw=2.2, label=label)
    for which, label, col, ls in (('fh2', 'FiberHMM 2.12.0 (TF min_llr 5)', BLUE, '-'), ('fh3', 'FiberHMM 3.0.2 (TF min_llr 5)', INK, '--')):
        v = 100 * float(bench[bench.caller == which].fp_bp.iloc[0])
        ax.axhline(v, color=col, ls=ls, lw=2, label=f'{label}: {v:.1f}%')
    ax.set_xscale('log'); ax.set_xticks([3, 5, 8, 12, 18, 25, 35, 50, 75], ['3', '5', '8', '12', '18', '25', '35', '50', '75'])
    ax.set(xlabel='LLR threshold (min_llr, nats)', ylabel='% of bp called protected', ylim=(0, None),
           title='False footprints on held-out M-class (near-saturated) reads, element ±100 bp')
    ax.legend(frameon=False, fontsize=12, loc='upper right')
    fig.savefig(out / 'deck_control_fp.png', dpi=200, bbox_inches='tight'); plt.close(fig)

    # 3. Size classes per read.
    classes = ['sub-TF 1-19', 'TF 20-39', 'PIC 40-69', 'nucleosome 70-1000', '>1000']
    nice = ['1–19 bp', '20–39 bp (TF)', '40–69 bp', '70–1,000 bp (nucleosome)', '> 1,000 bp']
    fig, ax = plt.subplots(figsize=(16, 5.6))
    x = np.arange(len(classes)); w = 0.2
    for k, (key, label, col, _) in enumerate(CALLERS):
        v = [summ.loc[key, f'per_read {c}'] for c in classes]
        bars = ax.bar(x + (k - 1.5) * w, v, width=w, color=col, edgecolor=BG, linewidth=1.5, label=label)
        for bx, vv in zip(bars, v):
            ax.text(bx.get_x() + bx.get_width() / 2, vv + 0.15, f'{vv:.1f}', ha='center', fontsize=10, color=INK)
    ax.set_xticks(x, nice); ax.grid(axis='x', visible=False)
    ax.set(ylabel='footprints per read', title='Footprints per read by length, 60,000 C-class WT reads (LLR at min_llr 5)')
    ax.set_ylim(0, 13.5)
    ax.legend(frameon=False, fontsize=12, ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.1))
    fig.savefig(out / 'deck_sizes.png', dpi=200, bbox_inches='tight'); plt.close(fig)

    # 4. Element occupancy in oligo coordinates, with dead sites marked.
    lo, hi = ELEM0 - 30, ELEM0 + 211 + 30
    xs = OLIGO(np.arange(lo, hi))
    dead_ref = ref0[dead]
    fig, ax = plt.subplots(figsize=(16, 5.8))
    shade(ax)
    for key, label, col, ls in CALLERS:
        ax.plot(xs, occ[key].values[lo:hi], color=col, ls=ls, lw=2.2, label=label)
    dr = dead_ref[(dead_ref >= lo) & (dead_ref < hi)]
    ax.scatter(OLIGO(dr), np.full(len(dr), 0.03), marker='|', s=160, color=ORANGE, linewidths=2, label='A/T never labeled (measured < 0.10)')
    ax.set(xlim=(xs[0], xs[-1]), ylim=(0, 1), xlabel='oligo position (element = 21–231)', ylabel='fraction of reads protected',
           title='Footprint occupancy across the element, C-class WT reads')
    ax.legend(frameon=False, fontsize=12, ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.16))
    fig.savefig(out / 'deck_occupancy.png', dpi=200, bbox_inches='tight'); plt.close(fig)

    # 5. Repeat 1 zoom: occupancy and the per-site labeling ceiling underneath.
    zlo, zhi = ELEM0 + (70 - 21), ELEM0 + (112 - 21)
    zx = OLIGO(np.arange(zlo, zhi))
    fig, (a, b) = plt.subplots(2, 1, figsize=(16, 7.4), sharex=True, gridspec_kw={'height_ratios': [1.5, 1], 'hspace': 0.15})
    a.axvspan(89 - 0.5, 104 + 0.5, color='#ECEAE2', zorder=0)
    a.text(96.5, 0.95, 'Repeat 1 (Sp1)', ha='center', va='top', fontsize=12, color=MUTED, transform=a.get_xaxis_transform())
    for key, label, col, ls in CALLERS:
        a.plot(zx, occ[key].values[zlo:zhi], color=col, ls=ls, lw=2.4, label=label)
    a.set(ylim=(0.4, 0.85), ylabel='fraction protected', title='Repeat 1 edge: extra FiberHMM protection over oligo 86–91')
    a.legend(frameon=False, fontsize=11, ncol=2, loc='upper left')
    sel = (ref0 >= zlo) & (ref0 < zhi)
    zs = OLIGO(ref0[sel]); zp = pa_meas[sel]; zh = hk_a[ref0[sel]]
    b.bar(zs - 0.2, zp, width=0.4, color=np.where(zp < 0.10, ORANGE, BLUE), label='measured (M-class)')
    b.bar(zs + 0.2, zh, width=0.4, color=GRAY, label='FiberHMM hexamer')
    for p in (86, 89, 90, 91):
        b.text(p, 1.02, seq[ELEM0 + p - 21], ha='center', fontsize=12, color=INK, fontweight='bold')
    b.set(ylim=(0, 1.12), ylabel='P(m6A | accessible)', xlabel='oligo position')
    b.legend(frameon=False, fontsize=11, loc='upper right', ncol=2)
    fig.savefig(out / 'deck_repeat1.png', dpi=200, bbox_inches='tight'); plt.close(fig)
    i90 = ELEM0 + 90 - 21
    stats['repeat1'] = dict(oligo_86_89_91_measured_pa=[float(pa_meas[ref0 == ELEM0 + p - 21][0]) for p in (86, 89, 90, 91)],
                            oligo_86_89_91_hexamer_pa=[float(hk_a[ELEM0 + p - 21]) for p in (86, 89, 90, 91)],
                            occupancy_at_90={k: float(occ[k].values[i90]) for k, *_ in CALLERS})
    json.dump(stats, open(out / 'deck_stats.json', 'w'), indent=1)
    print(json.dumps(stats, indent=1))


if __name__ == '__main__':
    main()

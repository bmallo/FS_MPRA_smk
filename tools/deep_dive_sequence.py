"""Where is the 211-bp element information-rich or -poor for plasmid Fiber-seq?

Combines the WT sequence (A/T opportunity, homopolymers, GC) with what the
replicate-3 extraction actually measured on WT molecules (per-site retained m6A
in the U / C / M distal-labeling classes) and with the per-SNV readout quality
of the mask-4 run (usable sites and SE of local:A_norm).

Usage: python tools/deep_dive_sequence.py --output results/deep_dive_20261026/sequence
"""
import argparse
import csv
import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

EXTRACT = Path('results/rep3_site_methylation_20261021/extract')
MASK4 = Path('docs/evidence/site_methylation_mask4_20261021/effects/effects.tsv.gz')
CONFIG = Path('config/ldlr_site_methylation_mask4.json')
U_MAX, M_MIN = 0.01, 0.60  # prespecified primary thresholds
OLIGO_OFFSET = 21  # oligo_pos1 = element index (0-based) + 21
WINDOW = 10  # +/- bp for local A/T and GC density

BG, INK, MUTED, GRID = '#F7F6F1', '#1B2A41', '#4A5568', '#DCDDD5'
BLUE, ORANGE, GRAY = '#2B6CB0', '#C8641E', '#8A94A3'
BASE_COLORS = {'A': '#2B6CB0', 'T': '#6FA0D8', 'G': '#C8641E', 'C': '#E9A36B'}
plt.rcParams.update({
    'figure.facecolor': BG, 'axes.facecolor': BG, 'savefig.facecolor': BG, 'font.size': 13,
    'axes.titlesize': 14, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'axes.labelcolor': MUTED,
    'axes.edgecolor': GRID, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.color': MUTED,
    'ytick.color': MUTED, 'text.color': INK, 'axes.titlecolor': INK, 'axes.grid': True, 'grid.color': GRID,
    'grid.linewidth': 0.7, 'axes.axisbelow': True,
})


def runs(seq):
    """Homopolymer run length at each position."""
    out, i = np.zeros(len(seq), int), 0
    while i < len(seq):
        j = i
        while j + 1 < len(seq) and seq[j + 1] == seq[i]:
            j += 1
        out[i:j + 1] = j - i + 1
        i = j + 1
    return out


def window_mean(x, w):
    k = np.ones(2 * w + 1)
    return np.convolve(x, k, 'same') / np.convolve(np.ones_like(x, float), k, 'same')


def wt_site_rates():
    """Per-position retained-call fraction among A/T opportunities, WT molecules, by distal class."""
    meta = np.load(EXTRACT / 'meta.npz')
    variants = json.load(open(EXTRACT / 'variants.json'))
    wt_index = next(i for i, v in enumerate(variants) if v['class'] == 'WT')
    rows = np.flatnonzero((meta['variant_index'] == wt_index) & meta['eligible'])
    distal = meta['distal_m6a'][rows]
    scores = np.load(EXTRACT / 'scores.npy', mmap_mode='r')
    known = np.load(EXTRACT / 'known.npy', mmap_mode='r')
    out = {}
    for name, sel in (('U', distal <= U_MAX), ('M', distal >= M_MIN), ('C', (distal > U_MAX) & (distal < M_MIN))):
        r = rows[sel]
        called = np.zeros(211); opp = np.zeros(211)
        for start in range(0, len(r), 200_000):  # bounded memory
            idx = r[start:start + 200_000]
            k = known[idx] > 0
            called += ((scores[idx] > 0) & k).sum(0)
            opp += k.sum(0)
        out[name] = dict(molecules=int(len(r)), rate=np.where(opp > 0, called / np.maximum(opp, 1), np.nan), opp=opp)
    return variants[wt_index]['sequence'], out


def snv_readout():
    by = {}
    with gzip.open(MASK4, 'rt') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            if r['endpoint'] in ('local:A_norm', 'local:M', 'flank:A_norm'):
                by.setdefault(int(r['oligo_pos1']), []).append(r)
    return by


def deck_figures(seq, rates, at, run, run_list, at_density, gc_density, local_sites, regions, out):
    """Two slide-sized figures: sequence context, and A-run blindness in saturated molecules."""
    n = len(seq)
    x = np.arange(n) + OLIGO_OFFSET
    short = {'P1_core': 'P1 core', 'Repeat1_Sp1': 'R1 Sp1', 'Repeat2_SRE1': 'R2 SRE1', 'Repeat3_Sp1': 'R3 Sp1'}
    fig, ax = plt.subplots(3, 1, figsize=(16, 7.2), sharex=True, gridspec_kw={'height_ratios': [0.55, 1, 1], 'hspace': 0.45})
    for i, b in enumerate(seq):
        ax[0].bar(x[i], 1, width=1, color=BLUE if b in 'AT' else ORANGE, edgecolor=BG, linewidth=0.3)
    for r in run_list:
        ax[0].plot([r['oligo_start'] - 0.5, r['oligo_end'] + 0.5], [1.18, 1.18], color=INK, lw=3.5, solid_capstyle='butt')
    ax[0].set(ylim=(0, 1.35), yticks=[], title='WT sequence: A/T (blue) can carry m6A, G/C (orange) cannot; black bars = runs ≥ 3')
    ax[0].grid(False)
    for k, (a, b) in regions.items():
        if k in short:
            for axx in ax[1:]:
                axx.axvspan(a + OLIGO_OFFSET - 0.5, b + OLIGO_OFFSET - 0.5, color='#ECEAE2', zorder=0)
            ax[1].text((a + b) / 2 + OLIGO_OFFSET - 0.5, 1.03, short[k], ha='center', fontsize=11, color=MUTED, transform=ax[1].get_xaxis_transform())
    ax[1].plot(x, at_density, color=BLUE, lw=2.2)
    ax[1].plot(x, gc_density, color=ORANGE, lw=2.2)
    ax[1].text(x[-1] + 2, at_density[-1], 'A/T', color=BLUE, fontsize=12, va='center')
    ax[1].text(x[-1] + 2, gc_density[-1] - 0.08, 'G/C', color=ORANGE, fontsize=12, va='center')
    ax[1].set(ylabel='fraction ±10 bp', ylim=(0, 1))
    ax[2].bar(x, local_sites, width=0.9, color=np.where(local_sites < 3, ORANGE, BLUE))
    ax[2].set(ylabel='A/T sites', xlabel='oligo position (element = 21–231)',
              title='A/T sites available to a local readout (5–10 bp either side of each SNV); orange = fewer than 3')
    fig.savefig(out / 'deck_sequence_context.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    m = rates['M']['rate']
    ctx = np.array(['G/C' if not at[i] else 'run ≥ 3' if run[i] >= 3 else 'run 2' if run[i] == 2 else 'isolated' for i in range(n)])
    colors = {'isolated': BLUE, 'run 2': '#6FA0D8', 'run ≥ 3': ORANGE}
    fig, (a, b) = plt.subplots(1, 2, figsize=(16, 5.6), gridspec_kw={'width_ratios': [2.6, 1], 'wspace': 0.18})
    for k, c in colors.items():
        sel = ctx == k
        a.scatter(x[sel], m[sel], s=34 if k == 'run ≥ 3' else 22, color=c, label=f'A/T {k}', zorder=3)
    for r in run_list:
        if r['base'] in 'AT':
            a.axvspan(r['oligo_start'] - 0.5, r['oligo_end'] + 0.5, color='#F6E3D3', zorder=0)
    i125 = 125 - OLIGO_OFFSET
    a.annotate('125', (125, m[i125]), xytext=(125, 0.42), ha='center', fontsize=12, color=INK,
               arrowprops=dict(arrowstyle='-', color=INK, lw=0.8))
    a.set(title='Retained m6A per A/T site in near-saturated (M) WT molecules', xlabel='oligo position', ylabel='fraction of molecules called', ylim=(-0.03, 1.03))
    a.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), frameon=False, fontsize=12, ncol=3)
    groups = [(k, m[ctx == k]) for k in colors]
    rng = np.random.default_rng(1)
    for j, (k, v) in enumerate(groups):
        b.scatter(j + rng.uniform(-0.18, 0.18, len(v)), v, s=22, color=colors[k], alpha=0.85)
        b.plot([j - 0.3, j + 0.3], [np.median(v)] * 2, color=INK, lw=2.5)
        b.text(j, 1.06, f'median {np.median(v):.2f}\nn = {len(v)}', ha='center', fontsize=11)
    b.set_xticks(range(len(groups)), [k for k, _ in groups])
    b.set(ylim=(-0.03, 1.2), title='By A/T context', ylabel='')
    b.grid(axis='x', visible=False)
    fig.savefig(out / 'deck_arun_blind.png', dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', required=True)
    out = Path(ap.parse_args().output)
    out.mkdir(parents=True, exist_ok=True)
    seq, rates = wt_site_rates()
    regions = json.load(open(CONFIG))['analysis']['regions0']
    n = len(seq)
    at = np.array([b in 'AT' for b in seq], float)
    gc = 1 - at
    run = runs(seq)
    at_run = (run >= 3) & (at > 0)
    gc_run = (run >= 3) & (at == 0)
    at_density = window_mean(at, WINDOW)
    gc_density = window_mean(gc, WINDOW)
    cpg = np.array([seq[i:i + 2] == 'CG' or seq[i - 1:i + 1] == 'CG' for i in range(n)])
    # A/T sites usable by a mask-4 local readout (5-10 bp either side, unchanged A/T).
    local_sites = np.array([sum(at[j] for j in range(n) if 5 <= abs(j - i) <= 10) for i in range(n)])
    readout = snv_readout()

    rows = []
    for i in range(n):
        p = i + OLIGO_OFFSET
        se = [float(r['se']) for r in readout.get(p, []) if r['endpoint'] == 'local:A_norm' and r['se']]
        rows.append(dict(element0=i, oligo_pos1=p, genomic_pos1=11089264 + p, base=seq[i], is_at=int(at[i]),
                         run_length=int(run[i]), at_density_pm10=round(at_density[i], 3), gc_density_pm10=round(gc_density[i], 3),
                         cpg=int(cpg[i]), local_at_sites_mask4=int(local_sites[i]),
                         wt_rate_U=rates['U']['rate'][i], wt_rate_C=rates['C']['rate'][i], wt_rate_M=rates['M']['rate'][i],
                         snv_local_Anorm_se_median=float(np.median(se)) if se else None,
                         regions=';'.join(k for k, (a, b) in regions.items() if a <= i < b)))
    with open(out / 'per_position.tsv', 'w', newline='') as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow({k: ('' if isinstance(v, float) and not np.isfinite(v) else v) for k, v in r.items()})

    # Hia5 / caller context in near-saturated (M) WT molecules: base and run context.
    m = rates['M']['rate']
    ctx = {}
    for i in range(n):
        if not at[i] or not np.isfinite(m[i]):
            continue
        key = f"{seq[i]} in {'A/T run ≥3' if at_run[i] else 'run 2' if run[i] == 2 else 'isolated'}"
        ctx.setdefault(key, []).append(m[i])
    context = {k: dict(sites=len(v), median_M_rate=float(np.median(v)), min=float(np.min(v)), max=float(np.max(v))) for k, v in sorted(ctx.items())}

    # Information-poor stretches: windows with < 3 A/T sites in +/-10 bp.
    poor = np.flatnonzero(local_sites < 3) + OLIGO_OFFSET
    run_list = []
    i = 0
    while i < n:
        if run[i] >= 3:
            run_list.append(dict(oligo_start=i + OLIGO_OFFSET, oligo_end=i + run[i] - 1 + OLIGO_OFFSET, base=seq[i], length=int(run[i])))
            i += run[i]
        else:
            i += 1
    summary = dict(
        length=n, at_sites=int(at.sum()), gc_fraction=float(gc.mean()), cpg_dinucleotides=int(seq.count('CG')),
        homopolymers_ge3=run_list, at_runs_ge3=sum(1 for r in run_list if r['base'] in 'AT'), gc_runs_ge3=sum(1 for r in run_list if r['base'] in 'GC'),
        wt_molecules=rates['U']['molecules'] + rates['C']['molecules'] + rates['M']['molecules'],
        wt_class_molecules={k: v['molecules'] for k, v in rates.items()},
        median_site_rate={k: float(np.nanmedian(v['rate'])) for k, v in rates.items()},
        M_rate_range=[float(np.nanmin(m)), float(np.nanmax(m))], M_rate_by_context=context,
        snv_positions_with_lt3_local_sites=[int(p) for p in poor],
        corr_local_sites_vs_se=float(np.corrcoef(*np.array([(r['local_at_sites_mask4'], r['snv_local_Anorm_se_median']) for r in rows if r['snv_local_Anorm_se_median']]).T)[0, 1]),
        regions={k: dict(oligo=[a + OLIGO_OFFSET, b - 1 + OLIGO_OFFSET], at_sites=int(at[a:b].sum()), length=b - a,
                         gc_fraction=float(gc[a:b].mean())) for k, (a, b) in regions.items()},
    )
    json.dump(summary, open(out / 'summary.json', 'w'), indent=1, default=lambda o: o.item())

    # Figure: tracks along the element.
    x = np.arange(n) + OLIGO_OFFSET
    fig, ax = plt.subplots(5, 1, figsize=(16, 11), sharex=True, gridspec_kw={'height_ratios': [0.5, 1, 1.3, 1.3, 1], 'hspace': 0.35})
    for i, b in enumerate(seq):
        ax[0].bar(x[i], 1, width=1, color=BASE_COLORS[b], edgecolor=BG, linewidth=0.3)
    for r in run_list:
        ax[0].plot([r['oligo_start'] - 0.5, r['oligo_end'] + 0.5], [1.15, 1.15], color=INK, lw=3, solid_capstyle='butt')
    ax[0].set(ylim=(0, 1.3), yticks=[], title='WT sequence (blue A/T = m6A-observable, orange G/C = blind); black bars = homopolymer runs ≥ 3')
    ax[0].grid(False)
    for k, (a, b) in regions.items():
        if k.startswith(('P1_core', 'Repeat')):
            for axx in ax[1:]:
                axx.axvspan(a + OLIGO_OFFSET - 0.5, b + OLIGO_OFFSET - 0.5, color='#ECEAE2', zorder=0)
            ax[1].text((a + b) / 2 + OLIGO_OFFSET, 1.02, {'P1_core': 'P1 core', 'Repeat1_Sp1': 'R1 Sp1', 'Repeat2_SRE1': 'R2 SRE1', 'Repeat3_Sp1': 'R3 Sp1'}[k], ha='center', fontsize=10, color=MUTED, transform=ax[1].get_xaxis_transform())
    ax[1].plot(x, at_density, color=BLUE, lw=2, label='A/T fraction ±10 bp')
    ax[1].plot(x, gc_density, color=ORANGE, lw=2, label='G/C fraction ±10 bp')
    ax[1].set(ylabel='fraction', ylim=(0, 1))
    ax[1].legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=10, frameon=False)
    for name, col in (('M', ORANGE), ('C', BLUE), ('U', GRAY)):
        ax[2].scatter(x, rates[name]['rate'], s=12, color=col, label=f'{name} class (n = {rates[name]["molecules"]:,})')
    ax[2].set(ylabel='retained m6A\nper A/T site', title='WT per-site m6A by distal-labeling class (M ≈ observability ceiling, U ≈ false-positive floor)')
    ax[2].legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=10, frameon=False)
    ax[3].bar(x, local_sites, width=0.9, color=np.where(local_sites < 3, ORANGE, BLUE))
    ax[3].set(ylabel='A/T sites', title='Usable A/T sites in a mask-4 local window (5–10 bp either side); orange < 3')
    se = np.array([r['snv_local_Anorm_se_median'] if r['snv_local_Anorm_se_median'] else np.nan for r in rows])
    ax[4].scatter(x, 100 * se, s=12, color=INK)
    ax[4].set(ylabel='SE, pp', xlabel='oligo position (element = 21–231)', title='Per-SNV standard error of local:A_norm (median of 3 substitutions, mask-4 run)')
    fig.savefig(out / 'sequence_tracks.png', dpi=170, bbox_inches='tight')
    fig.savefig(out / 'sequence_tracks.svg', bbox_inches='tight')
    plt.close(fig)
    deck_figures(seq, rates, at, run, run_list, at_density, gc_density, local_sites, regions, out)
    print(json.dumps(summary, indent=1, default=lambda o: o.item()))


if __name__ == '__main__':
    main()

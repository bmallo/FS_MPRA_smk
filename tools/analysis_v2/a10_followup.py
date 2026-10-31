"""Step A10 follow-up figures (exploratory, after A10; HepG2 run 3 only), from existing outputs only.

  a10f_positions.png  where the SNVs with and without a significant E1 change sit (A7 and +/-5 bp masks)
  a10f_shift.png      each SNV's mean z more than 50 bp away vs its change in 81+ bp call coverage there; and the
                      far-field z distribution before and after removing each SNV's own shift
  a10f_evidence.png   A8 spike detection (M molecules, perfect spikes) vs the summed evidence of the spiked
                      adenines, sum of ln((1 - p_p) / (1 - p_a)) from the run-3 fold-A emissions, and vs their number

Usage: python tools/analysis_v2/a10_followup.py --a10 <a10 dir> --a8 <a8 dir> --extract <A6 HepG2_run3 dir> \
           --emissions <emissions.tsv.gz> --output <new dir>
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('a7f', Path(__file__).with_name('a7_figures.py'))
a7f = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a7f)
_spec = importlib.util.spec_from_file_location('a7s', Path(__file__).with_name('a7_snv_edits.py'))
a7s = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a7s)
INK, BODY, BLUE, ORANGE, GREY = a7f.INK, a7f.BODY, a7f.BLUE, a7f.ORANGE, a7f.GREY
MOTIFS = [('ZNF652-like', 52, 62), ('SP4', 82, 102), ('SREBF1/2', 127, 136), ('SP1/SP3/KLF', 140, 156),
          ('CREB1/CREM', 182, 192)]
FAR = 50
EL0, OLIGO0 = 1234, 21


def robust_sd(x):
    x = np.asarray(x, float)
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def positions(a10, out):
    fig, ax = plt.subplots(2, 1, figsize=(16, 7.6), sharex=True, gridspec_kw=dict(hspace=0.45))
    for a, (run, rl) in zip(ax, [('HepG2_run3', 'A7 mask'), ('HepG2_run3_pm5', '±5 bp mask')]):
        s = pd.read_csv(a10 / run / 'snvs.tsv', sep='\t')
        s = s[s.endpoint == 'E1']
        y = -np.log10(s.simes_p.clip(lower=1e-30))
        for sel, c, lab in ((True, '#A7AFBA', 'changes protection somewhere'), (False, ORANGE, 'no change found')):
            g = s.selected == sel
            a.scatter(s.snv_oligo[g], y[g], s=16 if sel else 34, color=c, zorder=3 if not sel else 2,
                      label=f'{lab} ({int(g.sum())})')
        for name, lo, hi in MOTIFS:
            a.axvspan(lo - 0.5, hi + 0.5, color='#DCE8F5', lw=0, zorder=0)
            a.text((lo + hi) / 2, 1.02, name, transform=a.get_xaxis_transform(), ha='center', va='bottom',
                   color=BODY, fontsize=10)
        a.set_ylim(0, 31)
        a.set_ylabel('-log10 Simes p')
        a.set_title(f'HepG2 run 3 · adenine unmarked · {rl}', loc='left', color=INK, fontsize=13, pad=22)
        a.spines[['top', 'right']].set_visible(False)
        a.legend(frameon=False, fontsize=11, loc='upper right', bbox_to_anchor=(1, 1.22), ncol=2)
    ax[-1].set_xlabel('SNV position (oligo) · values above 30 drawn at 30 · blue bands: motif windows')
    fig.savefig(out / 'a10f_positions.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def shift(a10, out):
    e = pd.read_parquet(a10 / 'HepG2_run3_pm5' / 'effects.parquet')
    e['dist'] = (e.oligo - e.snv_oligo).abs()
    f = e[(e.endpoint == 'E1') & e.testable & (e.dist > FAR)].copy()
    sh = f.groupby('snv').z.mean()
    nuc = e[(e.endpoint == 'E2') & (e['bin'] == 'NUC') & e.testable & (e.dist > FAR)].groupby('snv').D.mean()
    j = pd.concat([sh.rename('z'), nuc.rename('nuc')], axis=1).dropna()
    r = j.corr().iloc[0, 1]
    f['zc'] = f.z - f.snv.map(sh)
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.2), gridspec_kw=dict(wspace=0.25))
    a = ax[0]
    a.scatter(100 * j.nuc, j.z, s=12, color=BLUE, alpha=0.7)
    a.axhline(0, color=GREY, lw=0.8)
    a.axvline(0, color=GREY, lw=0.8)
    a.set_xlabel('change in 81+ bp call coverage > 50 bp away (points)')
    a.set_ylabel('mean z of adenines > 50 bp away')
    a.set_title(f'One point per SNV · r = {r:.2f} (n = {len(j)})', loc='left', color=INK, fontsize=13)
    a = ax[1]
    bins = np.linspace(-6, 6, 61)
    a.hist(f.z.clip(-6, 6), bins=bins, density=True, histtype='step', lw=2, color=ORANGE,
           label=f'as tested: robust SD {robust_sd(f.z):.2f}, {100 * (f.z.abs() > 1.96).mean():.1f}% |z|>1.96')
    a.hist(f.zc.clip(-6, 6), bins=bins, density=True, histtype='step', lw=2, color=BLUE,
           label=f'each SNV\'s shift removed: {robust_sd(f.zc):.2f}, {100 * (f.zc.abs() > 1.96).mean():.1f}%')
    x = np.linspace(-6, 6, 300)
    a.plot(x, np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi), color=GREY, lw=1.2, ls='--', label='chance: 1.00, 5.0%')
    a.set_xlabel('z of adenines > 50 bp from their SNV')
    a.set_yticks([])
    a.legend(frameon=False, fontsize=10.5, loc='upper left')
    a.set_title('HepG2 run 3 · ±5 bp mask', loc='left', color=INK, fontsize=13)
    for a in ax:
        a.spines[['top', 'right']].set_visible(False)
    fig.savefig(out / 'a10f_shift.png', dpi=100, bbox_inches='tight')
    plt.close(fig)
    return dict(r=r, sd_before=robust_sd(f.z), sd_after=robust_sd(f.zc))


def evidence(a8, extract, emissions, out):
    pos = np.load(Path(extract) / 'molecules.npz')['pos']
    _hit, miss = a7s.emission_steps(emissions, pos)
    sp = pd.read_csv(a8 / 'spikes' / 'HepG2_run3' / 'spikes.tsv', sep='\t')
    sp = sp[(sp.emissions == 'position') & (sp['mode'] == 'perfect') & (sp.group == 'M')]
    sp = sp[(sp.centre_oligo >= OLIGO0) & (sp.centre_oligo <= OLIGO0 + 210) & (sp['size'] <= 40)].copy()
    s0 = sp.centre_oligo - OLIGO0 + EL0 - (sp['size'] - 1) // 2
    w, k = [], []
    for a_, b_ in zip(s0, s0 + sp['size']):
        m = (pos >= a_) & (pos < b_)
        w.append(miss[m].sum())
        k.append(int(m.sum()))
    sp['W'], sp['K'] = w, k
    fig, ax = plt.subplots(1, 2, figsize=(15, 5), gridspec_kw=dict(wspace=0.25))
    a = ax[0]
    for L, c in ((10, ORANGE), (15, '#9A4A12'), (20, BLUE), (30, INK), (40, GREY)):
        g = sp[sp['size'] == L]
        a.scatter(g.W, 100 * g.net, s=12, color=c, alpha=0.7, label=f'{L} bp')
    a.axvline(a7s.MIN_LLR, color=GREY, lw=1, ls='--')
    a.text(a7s.MIN_LLR, 50, ' caller threshold (5)', color=BODY, fontsize=11)
    a.set_xscale('symlog', linthresh=10)
    a.set_xlabel('summed evidence of the protected adenines')
    a.set_ylabel('detected in M molecules (%)')
    a.legend(frameon=False, fontsize=10.5, title='footprint', title_fontsize=10.5, loc='lower right')
    a = ax[1]
    g = sp[sp['size'] <= 20]
    kk = sorted(g.K.unique())
    data = [100 * g[g.K == q].net.to_numpy() for q in kk]
    a.boxplot(data, positions=kk, widths=0.6, showfliers=True, flierprops=dict(markersize=3, markeredgecolor=GREY),
              medianprops=dict(color=ORANGE, lw=2))
    a.set_xlabel('A/T bases inside the footprint (10-20 bp footprints)')
    a.set_ylabel('detected in M molecules (%)')
    for a in ax:
        a.spines[['top', 'right']].set_visible(False)
    ax[0].set_title('HepG2 run 3 · A8 perfect spikes', loc='left', color=INK, fontsize=13)
    ax[1].set_title('Same spikes, by count of A/T', loc='left', color=INK, fontsize=13)
    fig.savefig(out / 'a10f_evidence.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--a10', required=True)
    ap.add_argument('--a8', required=True)
    ap.add_argument('--extract', required=True)
    ap.add_argument('--emissions', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    a7f.a4f.style()
    positions(Path(a.a10), out)
    print(shift(Path(a.a10), out))
    evidence(Path(a.a8), a.extract, a.emissions, out)


if __name__ == '__main__':
    main()

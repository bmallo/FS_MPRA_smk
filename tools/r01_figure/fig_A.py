"""R01 figure option A: multi-modal evidence for LDLR promoter variants (display rules in figA_rules.json).

a  HepG2 MPRA expression change along the promoter, ClinVar class by marker, regulatory regions shaded.
b  Region x readout: share of each region's SNVs in the top 20% of each readout (MPRA, plasmid Fiber-seq,
   AlphaGenome; PromoterAI pending). Shows where readouts agree (R3 Sp1) and where each is blind.
c  Sp1 site (143-151): nucleosome-size protection change vs MPRA, ClinVar class by marker.
d  HepG2 vs HEK293T local inaccessibility change (SNVs with >= 3 HEK clones); separate cells, separate axes values.
Writes fig_A.pdf/.png and panel_*.tsv to --output.
Usage: python tools/r01_figure/fig_A.py --root results/analysis_v2_20261028 --work results/r01_figure_20261102 --output <dir>
"""
import argparse
import importlib.util
import json
import subprocess
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

HERE = Path(__file__).parent
_s = importlib.util.spec_from_file_location('fig_repeat3', HERE / 'fig_repeat3.py')
fr = importlib.util.module_from_spec(_s)
_s.loader.exec_module(fr)
INK, INK2, MUTED, GRID = fr.INK, fr.INK2, fr.MUTED, fr.GRID
ORANGE, BLUE, AQUA = '#eb6834', '#2a78d6', '#1baf7a'
RULES = json.load(open(HERE / 'figA_rules.json'))
PLP = ['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic']
BLB = ['Benign', 'Likely benign']
# panel a annotations (panel_a_rules.json; sites from results/.../panel_a_sites): (label, oligo lo, oligo hi, kind)
SITES = [('SP-like\n(weak match)', 89, 100, 'motif', 0),
         ('SREBP\n(SRE)', 127, 136, 'motif', 0), ('SP1 / KLF\n(GC box)', 143, 151, 'motif', 0)]
FP_MOTIF = {'F04': 'no known motif', 'F07': 'at SP-like site', 'F09': 'spans SP-like site and SRE', 'F18': ''}


def footprint_sites(root):
    """Confident TF-size footprints (panel_a_rules.json footprint_version) as SITES entries."""
    fp = pd.read_csv(Path(root) / 'fp_features/footprints.tsv', sep='\t')
    fp = fp[(fp.len_lo <= 45) & (fp.len_hi >= 17) & (fp.confident_share >= 0.8)].sort_values('span_lo')
    out, last_hi = [], -1
    for r in fp.itertuples():
        lo, hi = (66, 81) if r.feature == 'F04' else (int(r.span_lo), int(r.span_hi))
        row = 1 if lo <= last_hi else 0
        last_hi = max(last_hi, hi) if row == 0 else last_hi
        lab = f'~{int(r.len_peak)} bp\n{100 * r.wt_occupancy_confident:.2g}%'
        out.append((lab, lo, hi, 'fp' if r.feature == 'F04' else 'fp_motif', row))
    return out
C0, TSS_C = 285, -86                     # c. = oligo - 285 (c < 0); NM_000527.5 TSS at c.-86 (chr19:11,089,463)
VIOLET, VIOLET_BG, AQUA_BG = '#6250d6', '#ebe8fb', '#ddf3ea'
BLUES = LinearSegmentedColormap.from_list('blues', ['#f4f7fc', '#a9c6ef', '#2a78d6', '#103e7a'])


def clin_class(c):
    if c in PLP:
        return 'plp'
    if c in BLB:
        return 'blb'
    if isinstance(c, str) and c.startswith('Conflicting'):
        return 'conf'
    if c == 'Uncertain significance':
        return 'vus'
    return 'none'


def scatter_clin(ax, x, y, cls, small=4):
    st = {'none': dict(s=small, color='#c9c7c1', lw=0, zorder=2),
          'vus': dict(s=small + 6, color=BLUE, lw=0, zorder=3),
          'conf': dict(s=18, marker='D', facecolor='white', edgecolor=ORANGE, lw=0.9, zorder=4),
          'plp': dict(s=20, marker='D', color=ORANGE, edgecolor='white', lw=0.4, zorder=5),
          'blb': dict(s=16, marker='o', facecolor='white', edgecolor=AQUA, lw=1.0, zorder=4)}
    lab = {'none': 'Not in ClinVar', 'vus': 'ClinVar VUS', 'conf': 'Conflicting', 'plp': 'Pathogenic / likely path.',
           'blb': 'Benign / likely benign'}
    for k in ('none', 'vus', 'conf', 'blb', 'plp'):
        m = cls == k
        ax.scatter(x[m], y[m], label=f'{lab[k]} ({int(m.sum())})', **st[k])


def draw_panel_a(fig, spec, t, sites=SITES, gnomad=False, track_h=0.34):
    sa = spec.subgridspec(3, 1, height_ratios=[track_h, 1, 0.13 if gnomad else 0.001], hspace=0.08)
    axt = fig.add_subplot(sa[0])
    ax = fig.add_subplot(sa[1], sharex=axt)
    axr = fig.add_subplot(sa[2], sharex=axt)
    lo_c, hi_c = 21 - C0 - 1, 231 - C0 + 1
    for name, lo, hi, kind, row in sites:
        c_lo, c_hi = lo - C0, hi - C0
        fc, ec = (AQUA_BG, AQUA) if kind.startswith('fp') else (VIOLET_BG, VIOLET)
        if row == 0:
            ax.axvspan(c_lo - 0.5, c_hi + 0.5, color=fc, lw=0, zorder=0)
            y0, h = 0.05, 0.32
        else:
            y0, h = 0.84, 0.13
        axt.add_patch(plt.Rectangle((c_lo - 0.5, y0), c_hi - c_lo + 1, h, facecolor=fc, edgecolor=ec, lw=0.8,
                                    hatch='////' if kind == 'fp' else None))
        if row == 0:
            axt.text((c_lo + c_hi) / 2, 0.42, name, ha='center', va='bottom', fontsize=5.6, color=INK, linespacing=1.0)
        else:
            axt.text(c_hi + 2, y0 + h / 2, name.replace('\n', ', '), ha='left', va='center', fontsize=5.6, color=INK)
    ax.axvline(TSS_C, color=INK2, lw=0.6, ls=(0, (2, 2)), zorder=1)
    axt.annotate('', xy=(TSS_C + 14, 0.22), xytext=(TSS_C, 0.22),
                 arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.9, shrinkA=0, shrinkB=0))
    axt.plot([TSS_C, TSS_C], [0.0, 0.22], color=INK, lw=0.9)
    axt.text(TSS_C, 0.47, 'TSS', ha='center', va='bottom', fontsize=5.8, color=INK, fontweight='bold')
    axt.set_ylim(0, 1)
    axt.axis('off')
    fr.tidy(ax)
    scatter_clin(ax, t.c_dot.values, t.act.values, t.cls.values)
    ax.axhline(0, color=INK2, lw=0.5)
    ax.set_ylabel('Expression change\n(log2, HepG2 MPRA)')

    g = t[t.gnomad_ac > 0] if gnomad else t.iloc[:0]
    axr.vlines(g.c_dot, 0, 1, color=np.where(g.gnomad_ac >= 2, BLUE, '#a9c6ef'), lw=0.9)
    axr.set_ylim(0, 1)
    axr.set_yticks([0.5])
    axr.set_yticklabels([f'gnomAD\n({len(g)})' if gnomad else ''], fontsize=5.6, linespacing=0.9)
    xa = axr if gnomad else ax
    axr.tick_params(axis='y', length=0)
    for sp in ('top', 'right', 'left'):
        axr.spines[sp].set_visible(False)
    axr.set_xlim(lo_c, hi_c)
    ticks = list(range(-250, -54, 25))
    if gnomad:
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(axis='x', length=0)
    else:
        axr.axis('off')
    xa.set_xticks(ticks)
    xa.set_xticklabels([f'c.{k}' for k in ticks])
    xa.tick_params(axis='x', labelbottom=True)
    xa.set_xlabel('Position relative to the start codon (HGVS c., NM_000527.5)')
    from matplotlib.lines import Line2D
    h, l = ax.get_legend_handles_labels()
    if gnomad:
        h += [Line2D([], [], marker='|', ls='none', color=BLUE, ms=7, mew=1.2),
              Line2D([], [], marker='|', ls='none', color='#a9c6ef', ms=7, mew=1.2), Line2D([], [], ls='none')]
        l += ['gnomAD v4.1, seen ≥ 2 times', 'gnomAD v4.1, seen once', 'GTEx: no eQTL in this region']
    if any(k.startswith('fp') for _, _, _, k, _ in sites) and not any(k == 'motif' for _, _, _, k, _ in sites):
        from matplotlib.patches import Patch
        h += [Patch(facecolor=AQUA_BG, edgecolor=AQUA, hatch='////'), Patch(facecolor=AQUA_BG, edgecolor=AQUA)]
        l += ['Fiber-seq footprint, no known motif', 'Fiber-seq footprint, motif match nearby']
    y = xa.get_position().y0 - 0.40 / fig.get_figheight()
    fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.5 * (axr.get_position().x0 + axr.get_position().x1), y),
               ncol=5 if len(h) <= 5 else 4, frameon=False, fontsize=5.8, handletextpad=0.3, columnspacing=1.0)
    fr.panel_letter(axt, 'a')
    return ax


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--work', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, w, out = Path(a.root), Path(a.work), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()
    stamp = f'draft {date.today()} · git {commit}+uncommitted'

    t = pd.read_csv(w / 'evidence_v3/evidence.tsv', sep='\t')
    t['cls'] = t.clinvar.map(clin_class).fillna('none')
    corr = pd.read_csv(w / 'correlations_v2/pairs.tsv', sep='\t')

    # readout scores for panel b
    fe = pd.read_parquet(root / 'fp_features/effects.parquet').query('mask == "A7" and calls == "confident"')
    fe = fe[fe.p.notna()].assign(z=lambda d: (d.D / d.se).abs())
    fpz = fe.groupby('snv').z.max()
    t['score_mpra'] = t.act.abs()
    t['score_fs'] = np.fmax(t.nuc_z.abs(), t.snv.map(fpz))
    t['score_pai'] = t.pai_score.abs()
    t['score_ag'] = np.fmax((t.ag_dnase_hepg2 / t.ag_dnase_hepg2.std()).abs(),
                            (t.ag_rna_seq_hepg2 / t.ag_rna_seq_hepg2.std()).abs())
    readouts = [('MPRA\n(HepG2)', 'score_mpra'), ('Plasmid\nFiber-seq\n(HepG2)', 'score_fs'),
                ('AlphaGenome\n(HepG2)', 'score_ag'), ('PromoterAI', 'score_pai')]
    for _, col in readouts:
        if col:
            t[col + '_top'] = t[col] >= t[col].quantile(0.8)
    rows = RULES['panel_b']['rows']
    inreg = np.zeros(len(t), bool)
    bt = []
    for name, rng in rows.items():
        if isinstance(rng, list):
            m = t.oligo.between(*rng).values
            inreg |= m
        else:
            m = ~inreg
        r = dict(region=name, n=int(m.sum()))
        for lab, col in readouts:
            r[lab.replace('\n', ' ')] = t.loc[m, col + '_top'].mean() if col else np.nan
        bt.append(r)
    bt = pd.DataFrame(bt)
    bt.to_csv(out / 'panel_b_regions.tsv', sep='\t', index=False, float_format='%.4g')

    fig = plt.figure(figsize=(7.2, 7.6))
    gs_a = fig.add_gridspec(1, 1, left=0.085, right=0.985, top=0.965, bottom=0.715)
    gs_l = fig.add_gridspec(2, 1, left=0.165, right=0.53, top=0.56, bottom=0.075, height_ratios=[1, 1.05],
                            hspace=0.62)
    gs_r = fig.add_gridspec(1, 1, left=0.63, right=0.985, top=0.56, bottom=0.075)

    # a
    draw_panel_a(fig, gs_a[0, 0], t)
    t[['snv', 'oligo', 'c_dot', 'act', 'clinvar', 'cls', 'gnomad_ac']].to_csv(out / 'panel_a.tsv', sep='\t', index=False)
    fa = plt.figure(figsize=(7.2, 3.0))
    draw_panel_a(fa, fa.add_gridspec(1, 1, left=0.085, right=0.985, top=0.92, bottom=0.3)[0, 0], t)
    fa.savefig(out / 'panel_a.pdf')
    fa.savefig(out / 'panel_a.png', dpi=300)
    ff = plt.figure(figsize=(7.2, 3.0))
    draw_panel_a(ff, ff.add_gridspec(1, 1, left=0.085, right=0.985, top=0.92, bottom=0.3)[0, 0], t,
                 sites=footprint_sites(root), track_h=0.42)
    ff.savefig(out / 'panel_a_footprints.pdf')
    ff.savefig(out / 'panel_a_footprints.png', dpi=300)

    # b (inset grid inside the left half of row 2, upper part)
    axb = fig.add_subplot(gs_l[0, 0])
    M = bt[[c for c in bt.columns if c not in ('region', 'n')]].values.astype(float)
    axb.imshow(np.ma.masked_invalid(M), cmap=BLUES, vmin=0, vmax=1, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            txt, col = ('pending', MUTED) if np.isnan(v) else (f'{100 * v:.0f}%', 'white' if v > 0.55 else INK)
            axb.text(j, i, txt, ha='center', va='center', fontsize=6 if not np.isnan(v) else 5, color=col)
    axb.set_xticks(range(M.shape[1]))
    axb.set_xticklabels([r[0] for r in readouts], fontsize=5.8)
    axb.xaxis.tick_top()
    axb.set_yticks(range(M.shape[0]))
    axb.set_yticklabels([f'{r} ({n})' for r, n in zip(bt.region, bt.n)], fontsize=6)
    axb.tick_params(length=0)
    for s in axb.spines.values():
        s.set_visible(False)
    axb.set_xticks(np.arange(-0.5, M.shape[1]), minor=True)
    axb.set_yticks(np.arange(-0.5, M.shape[0]), minor=True)
    axb.grid(which='minor', color='white', lw=1.5)
    axb.text(0.5, -0.08, 'share of SNVs among each readout\'s top 20% (expected 20%)', transform=axb.transAxes,
             ha='center', va='top', fontsize=5.8, color=INK2)
    fr.panel_letter(axb, 'b')

    # d: HepG2 vs HEK293T
    axd = fig.add_subplot(gs_l[1, 0])
    fr.tidy(axd, 'both')
    h = t[(t.hek_clones >= 3)].dropna(subset=['local_inacc_pts', 'hek_local_pts'])
    site = h.oligo.between(143, 151)
    axd.scatter(h.local_inacc_pts[~site], h.hek_local_pts[~site], s=7, color='#9a9893', lw=0, alpha=0.8,
                label='other SNVs')
    axd.scatter(h.local_inacc_pts[site], h.hek_local_pts[site], s=14, color=ORANGE, edgecolor='white', lw=0.3,
                label='Sp1 site')
    axd.axhline(0, color=INK2, lw=0.5)
    axd.axvline(0, color=INK2, lw=0.5)
    r = corr[(corr.x == 'local_inacc_pts') & (corr.y == 'hek_local_pts') & (corr.subset == 'all')].iloc[0]
    axd.text(0.03, 0.97, f'Spearman ρ = {r.rho:.2f} (95% CI {r.lo:.2f}–{r.hi:.2f})\n{int(r.n)} SNVs with ≥3 HEK293T clones',
             transform=axd.transAxes, ha='left', va='top', fontsize=5.8, color=INK)
    axd.set_xlabel('HepG2: inaccessibility change ≤10 bp from SNV (pts)')
    axd.set_ylabel('HEK293T: same measure (pts)')
    axd.legend(loc='lower right', frameon=False, fontsize=5.6, handletextpad=0.2)
    fr.panel_letter(axd, 'd')
    h[['snv', 'oligo', 'local_inacc_pts', 'hek_local_pts', 'hek_clones', 'hek_molecules']].to_csv(
        out / 'panel_d_hek.tsv', sep='\t', index=False)

    # c: Sp1 site
    axc = fig.add_subplot(gs_r[0, 0])
    fr.tidy(axc, 'both')
    s = t[t.oligo.between(143, 151)].copy()
    nse = pd.read_csv(w / 'gnomad_check/nucleosome_all.tsv', sep='\t').query(
        'endpoint == "nucleosome" and calls == "confident" and mask == "A7"').set_index('snv')
    s['se'] = 100 * s.snv.map(nse.se)
    axc.errorbar(s.nuc_D_pts, s.act, xerr=1.96 * s.se, fmt='none', ecolor='#c9c7c1', lw=0.6, zorder=1)
    scatter_clin(axc, s.nuc_D_pts.values, s.act.values, s.cls.values, small=10)
    for v, dx, dy in (('149C>T', -4.6, 0.35), ('150C>G', -1.0, 0.75), ('145C>T', -1.2, 0.55), ('146C>G', -6.2, -0.05),
                      ('147T>G', 0.8, -0.12), ('147T>A', 0.8, 0.15)):
        rr = s[s.snv == v].iloc[0]
        axc.annotate(f'{v} (c.{rr.c_dot})', (rr.nuc_D_pts, rr.act), xytext=(rr.nuc_D_pts + dx, rr.act + dy),
                     fontsize=5.5, color=INK2, arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.4))
    axc.axhline(0, color=INK2, lw=0.5)
    axc.axvline(0, color=INK2, lw=0.5)
    rc = corr[(corr.x == 'act') & (corr.y == 'nuc_D_pts') & (corr.subset == 'sp1_site_143_151')].iloc[0]
    nv = int(((s.cls == 'vus') & (s.nuc_D_pts / s.se > 1.96) & (s.act < -1)).sum())
    axc.text(0.97, 0.97, f'Sp1 site, {len(s)} SNVs\nSpearman ρ = {rc.rho:.2f}\n{nv} of {int((s.cls == "vus").sum())} VUS '
             f'match the\npathogenic pattern', transform=axc.transAxes, ha='right', va='top', fontsize=5.8, color=INK)
    axc.set_xlim(-5, 12.5)
    axc.set_xlabel('Plasmid Fiber-seq: change in molecules with a\n≥81-bp protected stretch (pts, 95% CI)')
    axc.set_ylabel('Expression change (log2, HepG2 MPRA)')
    axc.legend(loc='lower left', frameon=False, fontsize=5.5, handletextpad=0.2, bbox_to_anchor=(0, 0.3))
    fr.panel_letter(axc, 'c')
    s[['snv', 'c_dot', 'cls', 'nuc_D_pts', 'se', 'act', 'ag_dnase_hepg2']].to_csv(out / 'panel_c_site.tsv', sep='\t',
                                                                                  index=False)

    fig.text(0.985, 0.004, stamp, ha='right', va='bottom', fontsize=5, color=MUTED)
    fig.savefig(out / 'fig_A.pdf', metadata={'Title': 'LDLR promoter multi-modal variant evidence', 'Subject': stamp})
    fig.savefig(out / 'fig_A.png', dpi=300)


if __name__ == '__main__':
    main()

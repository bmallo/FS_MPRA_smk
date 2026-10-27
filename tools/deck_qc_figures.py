"""QC figures for the replicate-3 walkthrough deck, drawn from the pipeline's own tables.

Usage: python tools/deck_qc_figures.py --output docs/figures/deck_qc_20261026
Every number plotted is read from a file listed in SOURCES; nothing is typed in.
"""
import argparse
import csv
import gzip
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

REP3 = Path('/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/fiberseq/HepG2/replicate_3')
SOURCES = {
    'random50k': REP3 / 'qc/random50k.reads.tsv.gz',
    'read_classes': REP3 / 'report/filters/read_class_summary.tsv',
    'promoter_classes': REP3 / 'report/class_summary.tsv',
    'per_read_m6a': REP3 / 'qc/m6a_nuc/plasmid_per_read.tsv.gz',
    'ledger_summary': Path('results/rep3_ledger_20261016/summary.json'),
    'variant_coverage': Path('results/rep3_ledger_20261016/barcodes/variant_barcode_coverage.tsv'),
    'barcode_map': Path('results/rep3_ledger_20261016/barcodes/barcode_map.tsv.gz'),
    'distal_hist': Path('results/rep3_site_methylation_20261021/effects/distal_m6a_histogram.tsv'),
    'site_summary': Path('results/rep3_site_methylation_20261021/effects/summary.json'),
}
PLASMID_LEN = 4647

# Deck palette (validated: tools/dataviz validator, light surface #F7F6F1).
BG, INK, MUTED, GRID = '#F7F6F1', '#1B2A41', '#4A5568', '#DCDDD5'
BLUE, ORANGE, GRAY, LIGHT = '#2B6CB0', '#C8641E', '#8A94A3', '#B9BFC8'

plt.rcParams.update({
    'figure.facecolor': BG, 'axes.facecolor': BG, 'savefig.facecolor': BG,
    'font.family': 'DejaVu Sans', 'font.size': 15, 'axes.titlesize': 17, 'axes.titleweight': 'bold',
    'axes.titlelocation': 'left', 'axes.labelsize': 15, 'axes.labelcolor': MUTED, 'axes.edgecolor': GRID,
    'axes.spines.top': False, 'axes.spines.right': False, 'xtick.color': MUTED, 'ytick.color': MUTED,
    'xtick.labelsize': 13, 'ytick.labelsize': 13, 'text.color': INK, 'axes.titlecolor': INK,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8, 'axes.axisbelow': True,
})


def tsv(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as h:
        return list(csv.DictReader(h, delimiter='\t'))


def fmt(n):
    return f'{n / 1e6:.2f}M' if n >= 1e6 else f'{n / 1e3:.0f}k' if n >= 1e4 else f'{n:,.0f}'


def save(fig, out, name):
    fig.savefig(out / f'{name}.png', dpi=200, bbox_inches='tight')
    fig.savefig(out / f'{name}.svg', bbox_inches='tight')
    plt.close(fig)
    print('wrote', name)


def hbars(ax, labels, values, colors, total=None):
    y = np.arange(len(labels))[::-1]
    ax.barh(y, values, color=colors, height=0.7, edgecolor=BG, linewidth=2)
    ax.set_yticks(y, labels)
    ax.grid(axis='y', visible=False)
    ax.tick_params(axis='y', length=0, labelsize=15, labelcolor=INK)
    xmax = max(values)
    for yi, v in zip(y, values):
        text = fmt(v) + ((f'  ·  {100 * v / total:.1f}%' if v / total >= 0.001 else '  ·  <0.1%') if total else '')
        ax.text(v + xmax * 0.01, yi, text, va='center', fontsize=14, color=INK)
    ax.set_xlim(0, xmax * 1.32)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: fmt(v) if v else '0'))


def delivery(out):
    rows = tsv(SOURCES['random50k'])
    length = np.array([int(r['len']) for r in rows])
    passes = np.array([int(r['np']) for r in rows if r['np']])
    flanked = np.array([int(r['np']) for r in rows if r['np'] and r['insert_len'] != '-1'])
    fig, (a, b) = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw={'wspace': 0.25})
    a.hist(np.clip(length, 0, 12000), bins=np.arange(0, 12001, 100), color=BLUE, edgecolor=BG, linewidth=0.3)
    a.set_yscale('log')
    a.axvline(PLASMID_LEN, color=ORANGE, lw=1.5, ls='--')
    a.text(PLASMID_LEN + 250, a.get_ylim()[1] * 0.35, f'plasmid\n{PLASMID_LEN:,} bp', color=ORANGE, fontsize=14)
    a.set(title=f'Read length (random 50k sample)', xlabel='read length, bp (clipped at 12 kb)', ylabel='reads (log)')
    b.hist(np.clip(flanked, 0, 60), bins=np.arange(0, 62, 2), color=BLUE, edgecolor=BG, linewidth=0.3)
    top = b.get_ylim()[1]
    for t, lab, yy in ((10, 'np ≥ 10 (cohort)', 1.02), (20, 'np ≥ 20', 0.94)):
        b.axvline(t, color=ORANGE, lw=1.5, ls='--')
        b.text(t + 0.8, top * yy, lab, color=ORANGE, fontsize=13)
    b.set_ylim(0, top * 1.1)
    b.set(title=f'Passes per flanked read (n = {len(flanked):,})', xlabel='CCS passes (clipped at 60)', ylabel='reads')
    save(fig, out, 'delivery_length_passes')
    return dict(n=len(length), median_len=float(np.median(length)), flanked=len(flanked),
                np_ge10=float((flanked >= 10).mean()), np_ge20=float((flanked >= 20).mean()))


def alignment(out):
    rc = {r['pc']: int(r['reads']) for r in tsv(SOURCES['read_classes'])}
    order = [('plasmid_full', 'Plasmid, full length', BLUE), ('plasmid_partial', 'Plasmid, partial', BLUE),
             ('plasmid_concatemer', 'Plasmid, concatemer', BLUE), ('chrM', 'chrM', ORANGE),
             ('nuclear', 'Nuclear hg38', GRAY), ('other_genomic', 'Other genomic', GRAY), ('unmapped', 'Unmapped', GRAY)]
    total = sum(rc.values())
    fig, ax = plt.subplots(figsize=(13, 5))
    hbars(ax, [o[1] for o in order], [rc[o[0]] for o in order], [o[2] for o in order], total)
    ax.set_title(f'Primary reads by alignment target (n = {total:,})')
    ax.set_xlabel('reads')
    save(fig, out, 'alignment_targets')
    return dict(total=total, **rc)


def genotypes(out):
    all_ = {r['']: int(r['all']) for r in tsv(SOURCES['promoter_classes'])}
    order = [('WT', 'WT', BLUE), ('SNV', 'SNV', BLUE), ('INS', 'Insertion', GRAY), ('DEL', 'Deletion', GRAY),
             ('MULTI', 'Two or more edits', ORANGE), ('LARGE_INDEL', 'Large indel', GRAY),
             ('PARTIAL', 'Partial insert', LIGHT), ('NO_INSERT', 'No insert', LIGHT)]
    total = sum(all_.values())
    fig, ax = plt.subplots(figsize=(13, 5.4))
    hbars(ax, [o[1] for o in order], [all_[o[0]] for o in order], [o[2] for o in order], total)
    ax.set_title(f'Plasmid reads by promoter genotype (n = {total:,})')
    ax.set_xlabel('reads')
    save(fig, out, 'genotype_classes')
    return dict(total=total, **all_)


def funnel(out):
    s = json.load(open(SOURCES['ledger_summary']))
    q, v = s['extraction_qc'], s['validated_counts']
    steps = [('Primary plasmid reads', q['plasmid_primary'], GRAY),
             ('Exact design, element bases ≥ Q20', q['genotype_exact_design'], GRAY),
             ('Analysis-eligible WT + SNV', v['analysis_eligible'], BLUE),
             ('  … with supported barcode', v['analysis_eligible_with_barcode'], ORANGE)]
    fig, ax = plt.subplots(figsize=(13, 3.8))
    hbars(ax, [x[0] for x in steps], [x[1] for x in steps], [x[2] for x in steps], q['plasmid_primary'])
    ax.set_title('From plasmid reads to the analysis cohort')
    ax.set_xlabel('molecules')
    save(fig, out, 'cohort_funnel')
    return {x[0]: x[1] for x in steps}


def snv_rows(cohort='analysis_eligible'):
    rows = [r for r in tsv(SOURCES['variant_coverage']) if r['cohort'] == cohort and re.fullmatch(r'\d+[ACGT]>[ACGT]', r['source_variant'])]
    pos = np.array([int(re.match(r'\d+', r['source_variant']).group()) for r in rows])
    return rows, pos


def coverage(out):
    rows, pos = snv_rows()
    mol = np.array([float(r['molecules']) for r in rows])
    alt = [r['source_variant'][-1] for r in rows]
    med = np.median(mol)
    fig, (a, b) = plt.subplots(1, 2, figsize=(15, 5), gridspec_kw={'width_ratios': [1, 1.6], 'wspace': 0.2})
    a.hist(mol, bins=np.arange(0, 8001, 200), color=BLUE, edgecolor=BG, linewidth=0.5)
    a.axvline(med, color=ORANGE, lw=2, ls='--')
    a.text(med + 120, a.get_ylim()[1] * 0.9, f'median {med:,.0f}', color=ORANGE, fontsize=14)
    a.set(title=f'Molecules per SNV (n = {len(mol)})', xlabel='analysis-eligible molecules', ylabel='SNVs')
    b.scatter(pos, mol, s=14, color=BLUE, alpha=0.75, linewidths=0)
    low = mol < 100
    b.scatter(pos[low], mol[low], s=40, color=ORANGE, zorder=3)
    for k, (p, m, r) in enumerate(sorted(zip(pos[low], mol[low], np.array(rows, dtype=object)[low]), key=lambda t: t[0])):
        b.annotate(f"{r['source_variant']} ({m:.0f})", (p, m), xytext=(-40, 18 + 26 * (k % 2)), textcoords='offset points',
                   fontsize=12, color=ORANGE, arrowprops=dict(arrowstyle='-', color=ORANGE, lw=0.8),
                   bbox=dict(fc=BG, ec='none', pad=1.5))
    b.axhline(med, color=GRAY, lw=1.2, ls=':')
    b.set(title='Molecules per SNV along the element', xlabel='oligo position (element = 21–231)', ylabel='molecules', ylim=(0, None))
    save(fig, out, 'snv_coverage')
    q = np.percentile(mol, [10, 25, 50, 75, 90])
    return dict(n=len(mol), p10_25_50_75_90=q.tolist(), lt100=int(low.sum()), fold80=float(mol.mean() / np.percentile(mol, 20)))


def barcodes(out):
    rows, pos = snv_rows()
    bc = np.array([float(r['unique_barcodes']) for r in rows])
    mol = np.array([float(r['molecules']) for r in rows])
    reads = np.array([int(r['observed_reads']) for r in tsv(SOURCES['barcode_map']) if r['status'] == 'supported'])
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(18, 5.2), gridspec_kw={'wspace': 0.3})
    a.hist(bc, bins=np.arange(0, 1701, 50), color=ORANGE, edgecolor=BG, linewidth=0.5)
    a.axvline(np.median(bc), color=INK, lw=1.5, ls='--')
    a.text(np.median(bc) + 30, a.get_ylim()[1] * 0.9, f'median {np.median(bc):,.0f}', fontsize=14)
    a.set(title='Supported barcodes per SNV', xlabel='supported barcodes', ylabel='SNVs')
    b.scatter(bc, mol, s=14, color=BLUE, alpha=0.7, linewidths=0)
    b.set(title='Barcodes vs molecules', xlabel='supported barcodes', ylabel='eligible molecules')
    edges = np.arange(1.5, 41.5, 1)
    c.hist(np.clip(reads, 2, 40), bins=np.concatenate([[1.5], edges[1:]]), color=BLUE, edgecolor=BG, linewidth=0.5)
    c.set(title=f'Reads per barcode (n = {len(reads):,})', xlabel='reads (all classes; clipped at 40)', ylabel='barcodes')
    c.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: fmt(v) if v else '0'))
    save(fig, out, 'barcode_coverage')
    return dict(median_bc=float(np.median(bc)), p10_bc=float(np.percentile(bc, 10)), ge50=int((bc >= 50).sum()),
                supported_barcodes=len(reads), median_reads_per_bc=float(np.median(reads)), r_bc_mol=float(np.corrcoef(bc, mol)[0, 1]))


def methylation(out):
    h = tsv(SOURCES['distal_hist'])
    lo = np.array([float(r['bin_low']) for r in h]); hi = np.array([float(r['bin_high']) for r in h])
    wt = np.array([float(r['wt_molecules']) for r in h])
    th = json.load(open(SOURCES['site_summary']))['thresholds']
    u, m = th['unmethylated_max_distal_m6a'], th['methylated_min_distal_m6a']
    per = [r for r in tsv(SOURCES['per_read_m6a']) if r['structure'] == 'full_length']
    x = np.array([float(r['m6a_at']) for r in per]); y = np.array([float(r['n_nuc']) for r in per])
    fig, (a, b) = plt.subplots(1, 2, figsize=(16, 5.4), gridspec_kw={'wspace': 0.25})
    frac = wt / wt.sum()
    cls = np.where(hi <= u + 1e-9, GRAY, np.where(lo >= m - 1e-9, ORANGE, BLUE))
    a.bar(lo, frac, width=hi - lo, align='edge', color=cls, edgecolor=BG, linewidth=0.4)
    a.set_yscale('log')
    for t in (u, m):
        a.axvline(t, color=INK, lw=1.2, ls='--')
    shares = [frac[hi <= u + 1e-9].sum(), frac[(lo >= u - 1e-9) & (hi <= m + 1e-9)].sum(), frac[lo >= m - 1e-9].sum()]
    a.set_ylim(3e-7, 3)
    for xx, lab, col in ((0.02, f'U ≤ {u}\n{100 * shares[0]:.1f}%', MUTED), (0.3, f'C\n{100 * shares[1]:.1f}%', BLUE),
                         (0.62, f'M ≥ {m}\n{100 * shares[2]:.1f}%', ORANGE)):
        a.text(xx, 2.2, lab, color=col, fontsize=14, fontweight='bold', va='top')
    a.set(title=f'Distal m6A, WT molecules (n = {wt.sum():,.0f})', xlabel='distal retained m6A fraction', ylabel='fraction of molecules (log)', xlim=(0, 0.85))
    hb = b.hexbin(x, np.clip(y, 0, 35), gridsize=(45, 18), extent=(0, 0.85, 0, 35), bins='log', cmap='Blues', mincnt=1, linewidths=0)
    b.set(title=f'm6A vs nucleosomes (1% sample, n = {len(x):,})', xlabel='m6A / A+T per read', ylabel='nucleosome calls per read')
    b.grid(False)
    cb = fig.colorbar(hb, ax=b, pad=0.02)
    cb.set_label('reads (log)', color=MUTED)
    save(fig, out, 'm6a_labeling')
    return dict(wt=float(wt.sum()), U=shares[0], C=shares[1], M=shares[2], per_read_n=len(x), per_read_median_m6a=float(np.median(x)))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', required=True)
    out = Path(ap.parse_args().output)
    out.mkdir(parents=True, exist_ok=True)
    stats = {f.__name__: f(out) for f in (delivery, alignment, genotypes, funnel, coverage, barcodes, methylation)}
    stats['sources'] = {k: str(v) for k, v in SOURCES.items()}
    json.dump(stats, open(out / 'figure_stats.json', 'w'), indent=1)
    print(json.dumps(stats, indent=1))


if __name__ == '__main__':
    main()

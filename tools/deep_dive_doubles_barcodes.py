"""Two-SNV molecules and barcode collisions in the replicate-3 ledger.

doubles:  reads whose 211-bp element differs from WT at exactly two positions.
          Are the pairs random, recombination-like, clonal? How much SNV
          coverage could they add to local readouts?
barcodes: are collisions above the chance rate for N15, what are they made of,
          and why are most molecules without a barcode assignment?

Usage: python tools/deep_dive_doubles_barcodes.py --output results/deep_dive_20261026
"""
import argparse
import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LEDGER = Path('results/rep3_ledger_20261016/ledger/molecule_ledger.tsv.gz')
SUPPORT = Path('results/rep3_ledger_20261016/barcodes/barcode_variant_support.tsv.gz')
CATALOG = Path('results/rep3_ledger_20261016/catalog/catalog.json')
OLIGO_OFFSET = 21
# Partner distance beyond which a second SNV is outside a mask-4 local (5-10 bp) or
# flank (11-30 bp) window plus the ~±8 bp observability kernel.
RESCUE_LOCAL, RESCUE_FLANK = 20, 40

BG, INK, MUTED, GRID = '#F7F6F1', '#1B2A41', '#4A5568', '#DCDDD5'
BLUE, ORANGE, GRAY, LIGHT = '#2B6CB0', '#C8641E', '#8A94A3', '#B9BFC8'
plt.rcParams.update({
    'figure.facecolor': BG, 'axes.facecolor': BG, 'savefig.facecolor': BG, 'font.size': 13,
    'axes.titlesize': 14, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'axes.labelcolor': MUTED,
    'axes.edgecolor': GRID, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.color': MUTED,
    'ytick.color': MUTED, 'text.color': INK, 'axes.titlecolor': INK, 'axes.grid': True, 'grid.color': GRID,
    'grid.linewidth': 0.7, 'axes.axisbelow': True,
})


def jdump(obj, path):
    json.dump(obj, open(path, 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))


def load():
    cols = ['observed_sequence', 'observed_variant_id', 'genotype_class', 'genotype_status', 'barcode', 'barcode_status',
            'passes', 'read_quality', 'mapq', 'circular', 'sequence_qc', 'methylation_fraction', 'analysis_eligible',
            'barcode_assignment_status', 'assigned_variant_id']
    df = pd.read_csv(LEDGER, sep='\t', usecols=cols, dtype=str, keep_default_na=False)
    for c in ('passes', 'read_quality', 'mapq', 'methylation_fraction'):
        df[c] = pd.to_numeric(df[c], errors='coerce')
    for c in ('circular', 'sequence_qc', 'analysis_eligible'):
        df[c] = df[c] == '1'
    print('ledger rows', len(df), flush=True)
    return df


def catalog():
    lib = json.load(open(CATALOG))['library']
    variants = {v['variant_id']: v for v in lib['variants']}
    wt = next(v['sequence'] for v in variants.values() if v['class'] == 'WT')
    return variants, wt


def save_panels(fig, axes, names, out):
    """Also save each panel on its own, for slides."""
    fig.canvas.draw()
    inv = fig.dpi_scale_trans.inverted()
    axes = list(axes)
    boxes = [a.get_tightbbox(fig.canvas.get_renderer()).transformed(inv).padded(0.12) for a in axes]
    for a, name, box in zip(axes, names, boxes):
        for other in axes:
            other.set_visible(other is a)
        fig.savefig(out / f'panel_{name}.png', dpi=200, bbox_inches=box)
    for other in axes:
        other.set_visible(True)


# ---------------------------------------------------------------- doubles
def doubles(df, wt, out):
    out.mkdir(parents=True, exist_ok=True)
    qc = (df.passes >= 10) & (df.read_quality >= 0.99) & (df.mapq >= 20) & df.circular
    L = len(wt)
    same_len = df.observed_sequence.str.len() == L
    sub = df[qc & same_len]
    arr = np.frombuffer(''.join(sub.observed_sequence).encode(), dtype=np.uint8).reshape(-1, L)
    wt_arr = np.frombuffer(wt.encode(), dtype=np.uint8)
    diff = arr != wt_arr
    nd = diff.sum(1)
    print('hamming counts', Counter(nd.tolist()).most_common(6), flush=True)
    idx1, idx2 = np.flatnonzero(nd == 1), np.flatnonzero(nd == 2)
    # Singles: abundance per (pos, alt).
    p1 = diff[idx1].argmax(1)
    alt1 = arr[idx1, p1]
    single = Counter(zip(p1.tolist(), alt1.tolist()))
    pos_single = np.bincount(p1, minlength=L).astype(float)
    # Doubles: the two positions and alts.
    d2 = diff[idx2]
    pa = d2.argmax(1)
    pb = L - 1 - d2[:, ::-1].argmax(1)
    alta, altb = arr[idx2, pa], arr[idx2, pb]
    dist = pb - pa
    passes2 = sub.passes.values[idx2]
    # Null distance distributions from single-SNV position abundances.
    fpos = pos_single / pos_single.sum()
    null_ind = np.zeros(L); null_rec = np.zeros(L)
    for d in range(1, L):
        w = (fpos[:-d] * fpos[d:]).sum()
        null_ind[d] = w
        null_rec[d] = w * d
    null_ind /= null_ind.sum(); null_rec /= null_rec.sum()
    obs = np.bincount(dist, minlength=L).astype(float); obs /= obs.sum()
    cdf = lambda x: np.cumsum(x)
    ks_ind, ks_rec = float(np.abs(cdf(obs) - cdf(null_ind)).max()), float(np.abs(cdf(obs) - cdf(null_rec)).max())
    # Pair-level: observed pair counts vs product of single abundances (x distance).
    pairs = Counter(zip(pa.tolist(), alta.tolist(), pb.tolist(), altb.tolist()))
    xs, ys, ds = [], [], []
    for (a, aa, b, bb), n in pairs.items():
        fa, fb = single.get((a, aa), 0), single.get((b, bb), 0)
        if fa and fb:
            xs.append(fa * fb); ys.append(n); ds.append(b - a)
    xs, ys, ds = map(np.array, (xs, ys, ds))
    r_ind = float(np.corrcoef(np.log(xs), np.log(ys))[0, 1])
    r_rec = float(np.corrcoef(np.log(xs * ds), np.log(ys))[0, 1])
    # Per-position participation in doubles relative to singles.
    pos_double = np.bincount(np.concatenate([pa, pb]), minlength=L).astype(float)
    ratio = pos_double / np.maximum(pos_single, 1)
    # Clonality via barcodes: what else carries a double read's barcode?
    seqs = sub.observed_sequence.values
    bc = sub.barcode.values
    bstat = sub.barcode_status.values
    dkey = {i: seqs[i] for i in idx2}
    by_bc = defaultdict(list)
    ok = bstat == 'ok'
    for i in np.flatnonzero(ok):
        by_bc[bc[i]].append(i)
    clon = Counter()
    for k, i in enumerate(idx2):
        if not ok[i]:
            clon['barcode unreadable'] += 1
            continue
        mates = [j for j in by_bc[bc[i]] if j != i]
        if not mates:
            clon['barcode seen once'] += 1
            continue
        mseq = {seqs[j] for j in mates}
        parts = {wt[:pa[k]] + chr(alta[k]) + wt[pa[k] + 1:], wt[:pb[k]] + chr(altb[k]) + wt[pb[k] + 1:]}
        if mseq == {seqs[i]}:
            clon['same double on every read'] += 1
        elif mseq & parts and not (mseq - parts - {seqs[i]}):
            clon['shared with a constituent single'] += 1
        elif wt in mseq and not (mseq - {wt, seqs[i]} - parts):
            clon['shared with WT'] += 1
        else:
            clon['shared with unrelated genotype'] += 1
    # Rescue potential per SNV.
    gain_local, gain_flank = Counter(), Counter()
    for a, aa, b, bb in zip(pa, alta, pb, altb):
        d = b - a
        if d > RESCUE_LOCAL:
            gain_local[(a, aa)] += 1; gain_local[(b, bb)] += 1
        if d > RESCUE_FLANK:
            gain_flank[(a, aa)] += 1; gain_flank[(b, bb)] += 1
    keys = sorted(single)
    s_cnt = np.array([single[k] for k in keys], float)
    gl = np.array([gain_local.get(k, 0) for k in keys], float)
    gf = np.array([gain_flank.get(k, 0) for k in keys], float)
    meth = dict(single=float(np.nanmedian(sub.methylation_fraction.values[idx1])), double=float(np.nanmedian(sub.methylation_fraction.values[idx2])),
                wt=float(np.nanmedian(sub.methylation_fraction.values[nd == 0])))
    summary = dict(
        qc='passes >= 10, read quality >= 0.99, MAPQ >= 20, full-length plasmid, element length 211 (no indel)',
        qc_reads=int(len(sub)), hamming_counts={str(k): int(v) for k, v in sorted(Counter(nd.tolist()).items()) if k <= 6},
        hamming_ge7=int((nd >= 7).sum()), singles=int(len(idx1)), doubles=int(len(idx2)), doubles_np_ge20=int((passes2 >= 20).sum()),
        doubles_per_single=float(len(idx2) / len(idx1)), distinct_pairs=len(pairs),
        distance=dict(median=float(np.median(dist)), frac_lt10=float((dist < 10).mean()), frac_le20=float((dist <= RESCUE_LOCAL).mean()),
                      null_ind_median=float(np.searchsorted(np.cumsum(null_ind), 0.5)), null_rec_median=float(np.searchsorted(np.cumsum(null_rec), 0.5)),
                      ks_vs_independent=ks_ind, ks_vs_recombination=ks_rec,
                      frac_lt10_null_ind=float(null_ind[:10].sum()), frac_lt10_null_rec=float(null_rec[:10].sum())),
        pair_count_correlation=dict(log_product=r_ind, log_product_times_distance=r_rec, pairs_tested=int(len(xs))),
        position_ratio=dict(median=float(np.median(ratio)), p10=float(np.percentile(ratio, 10)), p90=float(np.percentile(ratio, 90)),
                            top=[dict(oligo_pos=int(p + OLIGO_OFFSET), ratio=float(ratio[p])) for p in np.argsort(-ratio)[:8]]),
        barcode_context=dict(clon),
        methylation_fraction_median=meth,
        rescue=dict(rule_local=f'partner > {RESCUE_LOCAL} bp away', rule_flank=f'partner > {RESCUE_FLANK} bp away',
                    snvs=len(keys), median_single_reads=float(np.median(s_cnt)), median_gain_local=float(np.median(gl)),
                    median_gain_flank=float(np.median(gf)), median_pct_gain_local=float(np.median(100 * gl / s_cnt)),
                    median_pct_gain_flank=float(np.median(100 * gf / s_cnt)), total_rescuable_reads_local=int((dist > RESCUE_LOCAL).sum())),
    )
    jdump(summary, out / 'summary.json')
    with open(out / 'pair_counts.tsv', 'w') as h:
        h.write('oligo_pos_a\talt_a\toligo_pos_b\talt_b\tdistance\treads\tsingle_reads_a\tsingle_reads_b\n')
        for (a, aa, b, bb), n in sorted(pairs.items(), key=lambda t: -t[1]):
            h.write(f'{a + OLIGO_OFFSET}\t{chr(aa)}\t{b + OLIGO_OFFSET}\t{chr(bb)}\t{b - a}\t{n}\t{single.get((a, aa), 0)}\t{single.get((b, bb), 0)}\n')

    fig, ax = plt.subplots(2, 2, figsize=(16, 10), gridspec_kw={'hspace': 0.4, 'wspace': 0.22})
    a = ax[0, 0]
    x = np.arange(1, L)
    a.bar(x, obs[1:], width=1, color=LIGHT, label=f'observed (n = {len(idx2):,})')
    a.plot(x, null_ind[1:], color=BLUE, lw=2, label='random pairing of singles')
    a.plot(x, null_rec[1:], color=ORANGE, lw=2, label='pairing × distance (recombination)')
    a.set(title='Distance between the two SNVs', xlabel='bp apart', ylabel='fraction of double reads')
    a.legend(frameon=False, fontsize=11)
    a = ax[0, 1]
    xp = np.arange(L) + OLIGO_OFFSET
    a.scatter(xp, ratio, s=12, color=INK)
    a.axhline(np.median(ratio), color=GRAY, ls=':')
    a.set(title='Double-read participation per position ÷ single reads', xlabel='oligo position', ylabel='ratio')
    a = ax[1, 0]
    labels = ['same double on every read', 'barcode seen once', 'shared with a constituent single', 'shared with WT', 'shared with unrelated genotype', 'barcode unreadable']
    vals = [clon.get(k, 0) for k in labels]
    y = np.arange(len(labels))[::-1]
    a.barh(y, vals, color=[BLUE, GRAY, ORANGE, ORANGE, ORANGE, LIGHT], edgecolor=BG, linewidth=2)
    a.set_yticks(y, labels); a.grid(axis='y', visible=False)
    for yi, v in zip(y, vals):
        a.text(v, yi, f'  {v:,} ({100 * v / len(idx2):.1f}%)', va='center', fontsize=11)
    a.set_xlim(0, max(vals) * 1.45)
    a.set(title="What else carries a double read's barcode", xlabel='double reads')
    a = ax[1, 1]
    a.hist(100 * gl / s_cnt, bins=40, color=BLUE, edgecolor=BG, label=f'local readout (partner > {RESCUE_LOCAL} bp)')
    a.hist(100 * gf / s_cnt, bins=40, color=ORANGE, edgecolor=BG, alpha=0.8, label=f'flank readout (partner > {RESCUE_FLANK} bp)')
    a.set(title='Extra molecules per SNV if far-apart doubles were rescued', xlabel='% of single-SNV reads added', ylabel='SNVs')
    a.legend(frameon=False, fontsize=11)
    save_panels(fig, ax.flat, ['distance', 'position_ratio', 'barcode_context', 'rescue'], out)
    fig.savefig(out / 'doubles.png', dpi=170, bbox_inches='tight')
    fig.savefig(out / 'doubles.svg', bbox_inches='tight')
    plt.close(fig)
    return summary


# ---------------------------------------------------------------- barcodes
def edit_relation(s1, s2):
    if len(s1) == len(s2):
        h = sum(x != y for x, y in zip(s1, s2))
        return '1 substitution' if h == 1 else '2 substitutions' if h == 2 else '≥ 3 substitutions (shifted indels)'
    if abs(len(s1) - len(s2)) == 1:
        long, short = (s1, s2) if len(s1) > len(s2) else (s2, s1)
        if any(long[:i] + long[i + 1:] == short for i in range(len(long))):
            return '1-bp indel'
        return 'indel + other'
    return 'larger difference'


def barcodes(df, variants, wt, out):
    out.mkdir(parents=True, exist_ok=True)
    usable = (df.genotype_status == 'exact_design') & df.sequence_qc & (df.barcode_status == 'ok')
    okb = df.barcode[(df.barcode_status == 'ok')]
    # Composition and effective barcode space.
    distinct = pd.unique(okb)
    distinct = [b for b in distinct if len(b) == 15]
    comp = np.zeros((15, 4))
    lut = {c: i for i, c in enumerate('ACGT')}
    for b in distinct:
        for i, c in enumerate(b):
            comp[i, lut[c]] += 1
    comp /= comp.sum(1, keepdims=True)
    q = float(np.prod((comp ** 2).sum(1)))
    # Richness from usable reads per barcode.
    ub = df.barcode[usable].value_counts()
    f = Counter(ub.values.tolist())
    s_obs = int(len(ub))
    chao1 = s_obs + f[1] ** 2 / (2 * f[2]) if f[2] else float('nan')
    exp_collide_obs = s_obs ** 2 * q / 2
    exp_collide_lib = chao1 ** 2 * q / 2
    # Collision anatomy.
    sup = pd.read_csv(SUPPORT, sep='\t', dtype={'barcode': str, 'variant_id': str}, keep_default_na=False)
    nvar = sup.groupby('barcode').variant_id.transform('size')
    col = sup[nvar > 1].sort_values(['barcode', 'supporting_molecules'], ascending=[True, False])
    rel, minority, dominance, classes, balanced = Counter(), Counter(), [], Counter(), 0
    resolvable_bc, resolvable_mol = 0, 0
    for b, g in col.groupby('barcode', sort=False):
        n = g.supporting_molecules.values
        v = g.variant_id.values
        top, second = variants[v[0]], variants[v[1]]
        rel[edit_relation(top['sequence'], second['sequence'])] += 1
        classes[' vs '.join(sorted([top['class'], second['class']]))] += 1
        minority[int(min(n[1:].sum(), 5))] += 1
        dominance.append(n[0] / n.sum())
        if n[1] >= 2:
            balanced += 1
        if n[0] >= 2 and n[0] >= 5 * n[1:].sum():
            resolvable_bc += 1; resolvable_mol += int(n[0])
    # Why molecules are not assigned.
    st = df.barcode_assignment_status
    reason = np.select(
        [st == 'supported', st == 'collision', st == 'genotype_disagreement', st == 'genotype_unresolved',
         (st == 'insufficient_support') & (df.barcode_status != 'ok'),
         (st == 'insufficient_support') & (df.genotype_status == 'off_design'),
         (st == 'insufficient_support') & (df.genotype_status == 'low_base_quality'),
         (st == 'insufficient_support') & (df.genotype_status != 'exact_design'),
         (st == 'insufficient_support') & ~df.sequence_qc,
         (st == 'insufficient_support') & usable],
        ['supported', 'barcode collision', 'genotype disagrees with barcode', 'off-design read on an assigned barcode',
         'barcode unreadable', 'element is not a design (doubles, large indels)', 'element base below Q20', 'element flanks incomplete',
         'low passes or read quality', 'only usable read of its barcode'],
        default='other')
    reasons = Counter(reason.tolist())
    singleton_eligible = int(((reason == 'only usable read of its barcode') & df.analysis_eligible).sum())
    singleton_seed = None
    summary = dict(
        barcode_length=15, theoretical_space=4 ** 15, effective_space=1 / q, collision_prob_per_pair=q,
        composition=[{b: round(float(comp[i, j]), 3) for j, b in enumerate('ACGT')} for i in range(15)],
        distinct_readable_barcodes=len(distinct), usable_barcodes=s_obs, usable_reads_per_barcode={str(k): int(v) for k, v in sorted(f.items())[:10]},
        singleton_fraction_of_usable_barcodes=f[1] / s_obs, chao1_library_barcodes=chao1, sampled_fraction=s_obs / chao1,
        expected_random_collisions=dict(among_observed=exp_collide_obs, among_library_estimate=exp_collide_lib),
        observed_collision_barcodes=int(col.barcode.nunique()), collisions_with_minority_ge2_reads=balanced,
        collision_minority_reads={('5+' if k == 5 else str(k)): int(v) for k, v in sorted(minority.items())},
        collision_median_majority_share=float(np.median(dominance)),
        collision_pair_relation=dict(rel), collision_pair_classes=dict(classes.most_common(10)),
        resolvable_by_majority=dict(rule='majority >= 2 reads and >= 5x all others', barcodes=resolvable_bc, molecules=resolvable_mol),
        molecule_status=dict(reasons), singleton_usable_reads_analysis_eligible=singleton_eligible,
    )
    jdump(summary, out / 'summary.json')

    fig, ax = plt.subplots(2, 2, figsize=(19, 10), gridspec_kw={'hspace': 0.42, 'wspace': 0.85})
    a = ax[0, 0]
    k = np.array(sorted(f)); v = np.array([f[i] for i in k])
    a.bar(k[k <= 30], v[k <= 30], color=BLUE, edgecolor=BG)
    a.set_yscale('log')
    a.set(title=f'Usable reads per barcode ({s_obs:,} barcodes)', xlabel='usable reads carrying the barcode', ylabel='barcodes (log)')
    a.text(0.98, 0.72, f'singletons {100 * f[1] / s_obs:.0f}%\nChao1 library ≈ {chao1 / 1e6:.1f}M barcodes\nsampled ≈ {100 * s_obs / chao1:.0f}%', transform=a.transAxes, ha='right', va='top', fontsize=12)
    a = ax[0, 1]
    order = ['supported', 'element is not a design (doubles, large indels)', 'barcode unreadable', 'only usable read of its barcode',
             'off-design read on an assigned barcode', 'element base below Q20', 'barcode collision', 'element flanks incomplete',
             'genotype disagrees with barcode', 'low passes or read quality']
    vals = [reasons.get(o, 0) for o in order]
    y = np.arange(len(order))[::-1]
    a.barh(y, vals, color=[BLUE, ORANGE, ORANGE, ORANGE] + [GRAY] * 6, edgecolor=BG, linewidth=2)
    a.set_yticks(y, order); a.grid(axis='y', visible=False)
    tot = sum(reasons.values())
    for yi, vv in zip(y, vals):
        a.text(vv, yi, f'  {vv / 1e6:.2f}M ({100 * vv / tot:.1f}%)', va='center', fontsize=11)
    a.set_xlim(0, max(vals) * 1.5)
    a.set(title='All 8.72M plasmid molecules by barcode outcome', xlabel='molecules')
    a.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v / 1e6:.0f}M' if v else '0'))
    a = ax[1, 0]
    labs = list(rel.keys()); vv = [rel[l] for l in labs]
    o = np.argsort(vv)
    a.barh(np.arange(len(labs)), np.array(vv)[o], color=ORANGE, edgecolor=BG)
    a.set_yticks(np.arange(len(labs)), np.array(labs)[o]); a.grid(axis='y', visible=False)
    a.set(title='How the top two genotypes on a colliding barcode differ', xlabel='collision barcodes')
    a = ax[1, 1]
    kk = sorted(minority); a.bar([str(x) if x < 5 else '5+' for x in kk], [minority[x] for x in kk], color=GRAY, edgecolor=BG)
    a.set(title='Reads supporting the minority genotype(s)', xlabel='minority reads on the colliding barcode', ylabel='collision barcodes')
    a.text(0.98, 0.95, f'expected by chance ≈ {exp_collide_obs:,.0f}\nobserved {col.barcode.nunique():,}', transform=a.transAxes, ha='right', va='top', fontsize=12)
    save_panels(fig, ax.flat, ['reads_per_barcode', 'outcomes', 'collision_relation', 'minority'], out)
    fig.savefig(out / 'barcodes.png', dpi=170, bbox_inches='tight')
    fig.savefig(out / 'barcodes.svg', bbox_inches='tight')
    plt.close(fig)
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', required=True)
    ap.add_argument('--part', choices=('doubles', 'barcodes', 'both'), default='both')
    args = ap.parse_args()
    out = Path(args.output)
    variants, wt = catalog()
    df = load()
    if args.part in ('doubles', 'both'):
        print(json.dumps(doubles(df, wt, out / 'doubles'), indent=1, default=str))
    if args.part in ('barcodes', 'both'):
        print(json.dumps(barcodes(df, variants, wt, out / 'barcodes'), indent=1, default=str))


if __name__ == '__main__':
    main()

"""Compare plasmid-LLR with FiberHMM 2.12.0 / 3.0.2 on the replicate-3 WT read sets.

Callers (all on the same reads, reference coordinates):
  llr_measured   plasmid-LLR, global mode, drawdown 25, per-position table from the proxy controls
  llr_kmer       the same scan with FiberHMM's hexamer table (isolates the table's effect)
  fh2            FiberHMM 2.12.0, linear re-call (MA nuc + tf)
  fh3            FiberHMM 3.0.2,  linear re-call (MA nuc + tf)
  fh2_orig       FiberHMM 2.12.0 as originally run (circular for full-length reads)

Threshold rule, fixed before looking at chromatin-read results: the primary LLR min_llr is
the smallest grid value whose false-positive bp rate on held-out (fold B) M-control reads,
over the element +/-100 bp, is <= FiberHMM 2.12.0's on the same reads. min_llr 5 and 12 are
also reported.

Usage (from the FS_MPRA_smk root, plasmid-llr's python env):
  python tools/llr_compare/compare_callers.py --workdir results/llr_vs_fiberhmm_20261026
"""
import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, '/mmfs1/gscratch/stergachislab/bmallo/large_home/code/mine/plasmid-llr/llr')
import llrlib  # noqa: E402

L = 4647
ELEM = (1234, 1445)            # 0-based half-open, the 211-bp mutable element
ELEM100 = (1134, 1545)         # element +/-100 bp: not used to define the U/M/C classes
CLASSES = [('sub-TF 1-19', 1, 19), ('TF 20-39', 20, 39), ('PIC 40-69', 40, 69), ('nucleosome 70-1000', 70, 1000), ('>1000', 1001, 10**6)]
GRID = [3, 5, 8, 12, 18, 25, 35, 50, 75]
LABEL = {'llr_measured': 'LLR, measured table', 'llr_kmer': 'LLR, FiberHMM hexamer table', 'fh2': 'FiberHMM 2.12.0',
         'fh3': 'FiberHMM 3.0.2', 'fh2_orig': 'FiberHMM 2.12.0 (original circular run)'}
COLOR = {'llr_measured': '#C8641E', 'llr_kmer': '#E9A36B', 'fh2': '#2B6CB0', 'fh3': '#1B2A41', 'fh2_orig': '#8A94A3'}
BG, INK, MUTED, GRIDC = '#F7F6F1', '#1B2A41', '#4A5568', '#DCDDD5'
plt.rcParams.update({'figure.facecolor': BG, 'axes.facecolor': BG, 'savefig.facecolor': BG, 'font.size': 13,
                     'axes.titlesize': 14, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'axes.labelcolor': MUTED,
                     'axes.edgecolor': GRIDC, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'text.color': INK, 'axes.titlecolor': INK, 'axes.grid': True, 'grid.color': GRIDC,
                     'grid.linewidth': 0.7, 'axes.axisbelow': True})


def fiberhmm_calls(bam):
    """{read: (ref_lo, ref_hi, [(lo, hi), ...])} from MA nuc + tf, reference coords, half-open.

    FiberHMM 3.x declares `@CO fiberhmm:coord=molecular`: ns/nl/as/al/MA are in original-
    fiber coordinates, which run opposite to stored SEQ on reverse-aligned reads. 2.x writes
    stored-SEQ coordinates. Molecular intervals are flipped to stored-SEQ before projection.
    """
    af = pysam.AlignmentFile(str(bam), 'rb')
    molecular = any(c.strip() == 'fiberhmm:coord=molecular' for c in af.header.to_dict().get('CO', []))
    print(f'{bam.name}: coordinates {"molecular" if molecular else "stored-SEQ"}', flush=True)
    out = {}
    for rec in af:
        if not llrlib.usable(rec):
            continue
        ma = llrlib.ma_intervals(rec)
        q = [(s, s + l) for s, l in ma.get('nuc', []) + ma.get('tf', [])]
        if molecular and rec.is_reverse:
            ql = rec.query_length
            q = [(ql - e, ql - s) for s, e in q]
        out[rec.query_name] = (rec.reference_start, rec.reference_end, llrlib.query_intervals_to_ref(rec, q))
    return out


def llr_calls(parquet, min_llr, spans):
    df = pd.read_parquet(parquet, columns=['min_llr', 'read_name', 'start', 'end'])
    df = df[np.isclose(df.min_llr, min_llr)]
    by = defaultdict(list)
    for r, s, e in zip(df.read_name.values, df.start.values, df.end.values):
        by[r].append((int(s) - 1, int(e)))  # 1-based inclusive -> 0-based half-open
    return {r: (lo, hi, by.get(r, [])) for r, (lo, hi) in spans.items()}


def mask(ivs):
    m = np.zeros(L, bool)
    for lo, hi in ivs:
        m[max(lo, 0):min(hi, L)] = True
    return m


def summarize(calls, reads, region=ELEM100):
    a, b = region
    prot, cov, n = 0, 0, 0
    occ, depth = np.zeros(L), np.zeros(L)
    cls = np.zeros(len(CLASSES)); cls_elem = np.zeros(len(CLASSES))
    lengths = []
    for r in reads:
        lo, hi, ivs = calls[r]
        m = mask(ivs)
        c = np.zeros(L, bool); c[max(lo, 0):min(hi, L)] = True
        occ += m & c; depth += c
        rc = c[a:b]
        prot += int((m[a:b] & rc).sum()); cov += int(rc.sum()); n += 1
        for s, e in ivs:
            ln = e - s
            lengths.append(ln)
            k = next(i for i, (_, x, y) in enumerate(CLASSES) if x <= ln <= y) if ln >= 1 else None
            if k is None:
                continue
            cls[k] += 1
            if s < ELEM[1] and e > ELEM[0]:
                cls_elem[k] += 1
    return dict(reads=n, protected_bp_fraction=prot / max(cov, 1), occupancy=occ / np.maximum(depth, 1),
                per_read={c[0]: cls[i] / max(n, 1) for i, c in enumerate(CLASSES)},
                per_read_overlapping_element={c[0]: cls_elem[i] / max(n, 1) for i, c in enumerate(CLASSES)},
                lengths=np.array(lengths))


def jaccard(c1, c2, reads, region):
    a, b = region
    out = []
    for r in reads:
        m1, m2 = mask(c1[r][2])[a:b], mask(c2[r][2])[a:b]
        u = (m1 | m2).sum()
        if u:
            out.append((m1 & m2).sum() / u)
    return float(np.mean(out)) if out else float('nan'), len(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workdir', required=True)
    W = Path(ap.parse_args().workdir)
    out = W / 'results/comparison'
    out.mkdir(parents=True, exist_ok=True)
    data, calls_dir = W / 'data', W / 'results/calls'
    samples = ('wt_meth', 'wt_unmeth', 'wt_chrom')

    fh = {s: {'fh2': fiberhmm_calls(data / f'{s}.fiberhmm.bam'), 'fh3': fiberhmm_calls(data / f'{s}.fiberhmm3.bam')} for s in samples}
    fh['wt_chrom']['fh2_orig'] = fiberhmm_calls(data / 'wt_chrom.orig.bam')
    print('loaded FiberHMM calls', {s: {k: len(v) for k, v in d.items()} for s, d in fh.items()}, flush=True)
    spans = {s: {r: (lo, hi) for r, (lo, hi, _) in fh[s]['fh2'].items()} for s in samples}
    pq = {'llr_measured': 'llr_global_foldA_pp0_dd25', 'llr_kmer': 'llr_global_kmer_dd25'}

    def llr(s, which, t):
        return llr_calls(calls_dir / f'{s}.{pq[which]}.parquet', t, spans[s])

    foldB = {s: [r for r in spans[s] if llrlib.read_fold(r) == 1] for s in ('wt_meth', 'wt_unmeth')}
    # Control benchmark (held-out fold B), element +/-100.
    bench = []
    for which in ('fh2', 'fh3'):
        fp = summarize(fh['wt_meth'][which], foldB['wt_meth'])['protected_bp_fraction']
        tp = summarize(fh['wt_unmeth'][which], foldB['wt_unmeth'])['protected_bp_fraction']
        bench.append(dict(caller=which, min_llr=None, fp_bp=fp, tp_bp=tp))
    for which in pq:
        for t in GRID:
            fp = summarize(llr('wt_meth', which, t), foldB['wt_meth'])['protected_bp_fraction']
            tp = summarize(llr('wt_unmeth', which, t), foldB['wt_unmeth'])['protected_bp_fraction']
            bench.append(dict(caller=which, min_llr=t, fp_bp=fp, tp_bp=tp))
    bench = pd.DataFrame(bench)
    bench.to_csv(out / 'control_benchmark.tsv', sep='\t', index=False)
    fh2_fp = float(bench[bench.caller == 'fh2'].fp_bp.iloc[0])
    ok = bench[(bench.caller == 'llr_measured') & (bench.fp_bp <= fh2_fp)]
    primary = float(ok.min_llr.min()) if len(ok) else float(max(GRID))
    print(f'primary min_llr = {primary} (FiberHMM 2.12.0 fold-B FP = {fh2_fp:.4f})', flush=True)

    # Chromatin reads.
    chrom = sorted(spans['wt_chrom'])
    rows, occ, lens = [], {}, {}
    callers = {}
    for t in sorted({primary, 5.0, 12.0}):
        for which in pq:
            callers[f'{which}@{t:g}'] = llr('wt_chrom', which, t)
    for which in ('fh2', 'fh3', 'fh2_orig'):
        callers[which] = fh['wt_chrom'][which]
    for name, c in callers.items():
        s = summarize(c, chrom)
        se = summarize(c, chrom, ELEM)
        occ[name], lens[name] = s['occupancy'], s['lengths']
        rows.append(dict(caller=name, reads=s['reads'], protected_frac_element_pm100=s['protected_bp_fraction'],
                         protected_frac_element=se['protected_bp_fraction'],
                         **{f'per_read {k}': v for k, v in s['per_read'].items()},
                         **{f'per_read_in_element {k}': v for k, v in s['per_read_overlapping_element'].items()}))
    table = pd.DataFrame(rows)
    table.to_csv(out / 'chromatin_summary.tsv', sep='\t', index=False)
    pd.DataFrame({k: v for k, v in occ.items()}).rename_axis('ref0').to_csv(out / 'chromatin_occupancy.tsv.gz', sep='\t')
    main_set = [f'llr_measured@{primary:g}', f'llr_kmer@{primary:g}', 'fh2', 'fh3', 'fh2_orig']
    jac = []
    for i, a in enumerate(main_set):
        for b in main_set[i + 1:]:
            for rname, reg in (('element_pm100', ELEM100), ('whole_read', (0, L))):
                j, n = jaccard(callers[a], callers[b], chrom, reg)
                jac.append(dict(a=a, b=b, region=rname, mean_jaccard=j, reads=n))
    pd.DataFrame(jac).to_csv(out / 'chromatin_agreement.tsv', sep='\t', index=False)
    json.dump(dict(primary_min_llr=primary, fh2_fp_bp_element_pm100=fh2_fp, chromatin_reads=len(chrom),
                   fold_b_reads={k: len(v) for k, v in foldB.items()}), open(out / 'summary.json', 'w'), indent=1)

    # Figures.
    base = lambda n: n.split('@')[0]
    fig, ax = plt.subplots(figsize=(8, 6))
    for which in pq:
        d = bench[bench.caller == which].sort_values('min_llr')
        ax.plot(d.fp_bp, d.tp_bp, '-o', color=COLOR[which], label=LABEL[which], ms=5)
        for _, r in d.iterrows():
            ax.annotate(f'{r.min_llr:g}', (r.fp_bp, r.tp_bp), xytext=(4, -10), textcoords='offset points', fontsize=9, color=COLOR[which])
    for which in ('fh2', 'fh3'):
        d = bench[bench.caller == which]
        ax.scatter(d.fp_bp, d.tp_bp, s=110, marker='D', color=COLOR[which], label=LABEL[which], zorder=4)
    ax.set_xscale('log')
    ax.set(title='Held-out proxy controls, element ±100 bp', xlabel='protected bp called on M-class reads (false positive)',
           ylabel='protected bp called on U-class reads')
    ax.legend(frameon=False, fontsize=11)
    fig.savefig(out / 'control_roc.png', dpi=200, bbox_inches='tight'); plt.close(fig)

    fig, ax = plt.subplots(2, 1, figsize=(16, 8), gridspec_kw={'hspace': 0.35})
    for name in main_set:
        for a_, (lo, hi) in ((ax[0], (0, L)), (ax[1], (ELEM100[0], ELEM100[1]))):
            x = np.arange(lo, hi)
            a_.plot(x + 1, occ[name][lo:hi], color=COLOR[base(name)], lw=1.6, label=LABEL[base(name)] + (f' (min_llr {name.split("@")[1]})' if '@' in name else ''))
    for a_ in ax:
        a_.axvspan(ELEM[0] + 1, ELEM[1], color='#ECEAE2', zorder=0)
        a_.set(ylabel='fraction of reads protected', ylim=(0, 1))
    ax[0].set(title=f'Footprint occupancy, C-class WT reads (n = {len(chrom):,}); shaded = 211-bp element', xlabel='plasmid position')
    ax[1].set(title='Element ±100 bp', xlabel='plasmid position (element = 1235–1445)')
    ax[0].legend(frameon=False, fontsize=10, ncol=3, loc='upper left', bbox_to_anchor=(0, -0.15))
    fig.savefig(out / 'chromatin_occupancy.png', dpi=200, bbox_inches='tight'); plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 5))
    bins = np.arange(0, 301, 5)
    for name in main_set:
        ax.hist(np.clip(lens[name], 0, 300), bins=bins, histtype='step', lw=2, color=COLOR[base(name)],
                weights=np.full(len(lens[name]), 1 / len(chrom)), label=LABEL[base(name)])
    ax.set(title='Footprint lengths, C-class WT reads', xlabel='footprint length, bp (clipped at 300)', ylabel='footprints per read per 5-bp bin')
    ax.legend(frameon=False, fontsize=11)
    fig.savefig(out / 'chromatin_lengths.png', dpi=200, bbox_inches='tight'); plt.close(fig)
    print(table.to_string()); print(pd.DataFrame(jac).to_string())


if __name__ == '__main__':
    main()

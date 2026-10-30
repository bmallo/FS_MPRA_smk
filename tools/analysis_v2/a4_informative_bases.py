"""Step A4: informative bases. What can the sequence (and the caller) let us see, and how does each variant change it?

Hia5 methylates adenine (N6-methyladenine) and a PacBio m6A caller reads it on both strands, so every A/T base
pair carries exactly one methylatable adenine (top strand for A, bottom strand for T) and every G/C pair carries
none. A footprint is seen only as a run of unmethylated A/T, so A/T placement sets the resolution.

Per element position (211 bp, oligo_pos1 = index + 21):
  sequence view: A/T or not, adenine strand, A/T counts in centred 11- and 21-bp windows, length of the A/T-free
      gap the position sits in, and the shortest window covering the position that holds >= k A/T (k = the
      number of unmethylated A/T the caller needs: ceil(min_llr / median llr_miss))
  caller view (plasmid-LLR as run on HepG2): llr_miss = ln((1 - p_p) / (1 - p_a)) per A/T from the fold-A
      emission table, and the shortest covering window whose summed llr_miss reaches min_llr. The emission table
      comes from in-sample proxy classes (M/U), not naked controls (step B1).
Per designed variant (SNV, 1-bp DEL, 1-bp INS):
  sequence class: A/T loss, A/T gain, A<->T swap (adenine moves strand), G<->C swap, DEL/INS of A/T or G/C
  caller view: plasmid-LLR keeps only reference A/T positions where the read base equals the reference base
      (llrlib.read_observations), so any substitution at an A/T (loss or swap) and any deleted A/T drops that
      position's evidence, and no gained or inserted A/T is ever used
  local effect: A/T count and summed llr_miss within +-10 bp, WT vs variant
Optionally, analysis-ready molecules per variant from a run's fsmpra ledger (rules from A3).

Usage: python tools/analysis_v2/a4_informative_bases.py --catalog <catalog.json> \
           --emissions <emissions.tsv.gz> --ledger <ledger dir> --output <new dir>
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

AT = set('AT')
OLIGO_OFFSET = 21
# LDLR promoter repeats 1-3 (Sudhof et al. 1987), oligo positions, inclusive
REGIONS = {'Repeat 1 (Sp1)': (89, 104), 'Repeat 2 (SRE-1)': (124, 139), 'Repeat 3 (Sp1)': (140, 155)}
LOCAL = 10


def shortest_cover(weights, need):
    """For each position i: length of the shortest window [s, e] with s <= i <= e and sum(weights[s:e+1]) >= need
    (np.nan when no window reaches it)."""
    w = np.asarray(weights, float)
    n = len(w)
    c = np.concatenate([[0.0], np.cumsum(w)])
    # emin[s]: last index of the shortest window starting at s that reaches `need` (n when none does)
    emin = np.searchsorted(c, c[:n] + need - 1e-12, side='left') - 1
    best = np.full(n, np.nan)
    for i in range(n):
        s = np.arange(i + 1)
        ok = emin[:i + 1] < n
        if ok.any():
            best[i] = (np.maximum(emin[:i + 1][ok], i) - s[ok] + 1).min()
    return best


def gap_length(is_at):
    """Length of the A/T-free run containing each position (0 at A/T)."""
    out = np.zeros(len(is_at), int)
    i = 0
    while i < len(is_at):
        if is_at[i]:
            i += 1
            continue
        j = i
        while j < len(is_at) and not is_at[j]:
            j += 1
        out[i:j] = j - i
        i = j
    return out


def window_sum(x, half):
    x = np.asarray(x, float)
    c = np.concatenate([[0.0], np.cumsum(x)])
    i = np.arange(len(x))
    lo, hi = np.maximum(i - half, 0), np.minimum(i + half + 1, len(x))
    return c[hi] - c[lo]


def seq_class(cls, ref, alt):
    if cls == 'SNV':
        if ref in AT and alt in AT:
            return 'A<->T swap'
        if ref in AT:
            return 'A/T loss'
        if alt in AT:
            return 'A/T gain'
        return 'G<->C swap'
    if cls == 'DEL':
        return 'DEL of A/T' if ref in AT else 'DEL of G/C'
    if cls == 'INS':
        return 'INS of A/T' if alt in AT else 'INS of G/C'
    raise ValueError(cls)


def region_of(pos1):
    for name, (a, b) in REGIONS.items():
        if a <= pos1 <= b:
            return name
    return ''


def emissions(path, start0, wt, pp_floor=0.0):
    t = pd.read_csv(path, sep='\t')
    t = t[(t.fold == 'A') & (t.ref_pos > start0) & (t.ref_pos <= start0 + len(wt))].copy()
    t['i'] = t.ref_pos - 1 - start0
    for _, r in t.iterrows():
        if wt[r.i] != r.ref_base:
            raise ValueError(f'emission base {r.ref_base} != WT {wt[r.i]} at element {r.i}')
    p_p = np.maximum(t.p_p.to_numpy(), pp_floor)
    t['llr_miss'] = np.log((1 - p_p) / (1 - t.p_a.to_numpy()))
    return t.set_index('i')


def positions_table(wt, em, k, min_llr):
    n = len(wt)
    is_at = np.array([b in AT for b in wt])
    llr = np.zeros(n)
    p_a = np.full(n, np.nan)
    low = np.zeros(n, bool)
    llr[em.index.to_numpy()] = em.llr_miss.to_numpy()
    p_a[em.index.to_numpy()] = em.p_a.to_numpy()
    if 'low_cov' in em:
        low[em.index.to_numpy()] = em.low_cov.astype(str).str.lower().eq('true').to_numpy()
    missing = is_at & ~np.isin(np.arange(n), em.index.to_numpy())
    return pd.DataFrame(dict(
        index0=np.arange(n), oligo_pos1=np.arange(n) + OLIGO_OFFSET, base=list(wt), is_at=is_at,
        adenine_strand=['top' if b == 'A' else 'bottom' if b == 'T' else '' for b in wt],
        at_11bp=window_sum(is_at, 5).astype(int), at_21bp=window_sum(is_at, 10).astype(int),
        gap_bp=gap_length(is_at), min_fp_seq_bp=shortest_cover(is_at, k),
        p_a=p_a, llr_miss=llr, emission_low_cov=low, emission_missing=missing,
        min_fp_llr_bp=shortest_cover(llr, min_llr),
        region=[region_of(i + OLIGO_OFFSET) for i in range(n)]))


def variants_table(variants, wt, pos):
    is_at = pos.is_at.to_numpy()
    llr = pos.llr_miss.to_numpy()
    rows = []
    for v in variants:
        if v['class'] == 'WT':
            continue
        cls, ref, alt, p1 = v['class'], v['ref'] or '', v['alt'] or '', v['oligo_pos1']
        i = p1 - OLIGO_OFFSET
        if cls in ('SNV', 'DEL') and wt[i] != ref:
            raise ValueError(f"{v['source_variant']}: WT base {wt[i]} at {i} != ref {ref}")
        if cls == 'SNV':
            changed = wt[:i] + alt + wt[i + 1:]
        elif cls == 'DEL':
            changed = wt[:i] + wt[i + 1:]
        else:
            changed = wt[:i] + alt + wt[i:]
        if changed != v['sequence']:
            raise ValueError(f"{v['source_variant']}: rebuilt sequence differs from catalog")
        sc = seq_class(cls, ref, alt)
        d_seq = (alt in AT) - (ref in AT)
        # caller view: a substituted or deleted reference A/T is dropped; inserted or gained A/T never counted
        dropped = cls in ('SNV', 'DEL') and is_at[i]
        lo, hi = max(i - LOCAL, 0), min(i + LOCAL + 1, len(wt))
        at_local = int(is_at[lo:hi].sum())
        llr_local = float(llr[lo:hi].sum())
        rows.append(dict(
            source_variant=v['source_variant'], variant_id=v['variant_id'], cls=cls, oligo_pos1=p1, index0=i,
            ref=ref, alt=alt, seq_class=sc, at_count_change=d_seq, caller_drops_at=dropped,
            caller_llr_removed=float(llr[i]) if dropped else 0.0,
            at_local_wt=at_local, at_local_variant=at_local + d_seq,
            llr_local_wt=llr_local, llr_local_frac_removed=(float(llr[i]) / llr_local if dropped and llr_local else 0.0),
            min_fp_seq_bp_wt=pos.min_fp_seq_bp.iat[min(i, len(wt) - 1)],
            min_fp_llr_bp_wt=pos.min_fp_llr_bp.iat[min(i, len(wt) - 1)],
            region=region_of(p1)))
    return pd.DataFrame(rows)


def ledger_counts(ledger_dir):
    """Analysis-ready molecules per variant id (A3 rules: exact design, sequence QC, full-length, m6A ok)."""
    cols = ['observed_variant_id', 'genotype_status', 'sequence_qc', 'length_qc', 'methylation_qc']
    n = pd.Series(dtype=float)
    for ch in pd.read_csv(Path(ledger_dir) / 'ledger/molecule_ledger.tsv.gz', sep='\t', usecols=cols,
                          keep_default_na=False, dtype=str, chunksize=2_000_000):
        ok = (ch.genotype_status == 'exact_design') & (ch.sequence_qc == '1') & (ch.length_qc == '1') & \
             (ch.methylation_qc == '1')
        n = n.add(ch[ok].observed_variant_id.value_counts(), fill_value=0)
    return n.astype(int)


def summarize(wt, pos, var, k, min_llr, med):
    gaps = []
    is_at = pos.is_at.to_numpy()
    i = 0
    while i < len(wt):
        if is_at[i]:
            i += 1
            continue
        j = i
        while j < len(wt) and not is_at[j]:
            j += 1
        gaps.append(dict(oligo_start=i + OLIGO_OFFSET, oligo_end=j - 1 + OLIGO_OFFSET, bp=j - i, seq=wt[i:j]))
        i = j
    gaps = sorted(gaps, key=lambda g: -g['bp'])
    fs, fl = pos.min_fp_seq_bp, pos.min_fp_llr_bp
    out = dict(
        element_bp=len(wt), at_bp=int(is_at.sum()), at_fraction=float(is_at.mean()),
        adenine_top=int((pos.base == 'A').sum()), adenine_bottom=int((pos.base == 'T').sum()),
        k_at_needed=k, min_llr=min_llr, median_llr_miss=med,
        gaps_ge_5bp=sum(g['bp'] >= 5 for g in gaps), longest_gaps=gaps[:6],
        min_fp_seq={'median': float(fs.median()), 'max': float(fs.max()), 'positions_gt_20bp': int((fs > 20).sum())},
        min_fp_llr={'median': float(fl.median()), 'max': float(fl.max()), 'positions_gt_20bp': int((fl > 20).sum()),
                    'positions_gt_40bp': int((fl > 40).sum()), 'positions_none': int(fl.isna().sum())},
        emission_low_cov=int(pos.emission_low_cov.sum()), emission_missing=int(pos.emission_missing.sum()),
        llr_miss_quantiles={q: float(np.quantile(pos.llr_miss[pos.is_at], q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9)},
        variant_classes=var.groupby(['cls', 'seq_class']).size().rename('n').reset_index().to_dict('records'),
        caller_drops={c: int(n) for c, n in var[var.caller_drops_at].groupby('cls').size().items()},
        variants_by_cls={c: int(n) for c, n in var.groupby('cls').size().items()},
        regions={})
    for name in REGIONS:
        p = pos[pos.region == name]
        out['regions'][name] = dict(at_bp=int(p.is_at.sum()), bp=len(p), min_fp_llr_median=float(p.min_fp_llr_bp.median()),
                                    llr_sum=float(p.llr_miss.sum()))
    if 'molecules' in var:
        m = var.groupby('cls').molecules.sum()
        md = var[var.caller_drops_at].groupby('cls').molecules.sum()
        out['molecules'] = {c: dict(total=int(m[c]), caller_drops_at=int(md.get(c, 0))) for c in m.index}
    return out


def figure(pos, var, out, k, min_llr):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ink, body, blue, orange, grey = '#1B2430', '#4A5568', '#2F6DB5', '#D9782D', '#A7AFB9'
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 13, 'axes.edgecolor': grey,
                         'axes.labelcolor': body, 'xtick.color': body, 'ytick.color': body})
    x = pos.oligo_pos1.to_numpy()
    fig, ax = plt.subplots(4, 1, figsize=(18, 9.2), sharex=True,
                           gridspec_kw=dict(height_ratios=[0.8, 1.6, 1.6, 1.4], hspace=0.18))
    for a in ax:
        for name, (s, e) in REGIONS.items():
            a.axvspan(s - 0.5, e + 0.5, color='#E3EDF8', zorder=0, lw=0)
        a.spines[['top', 'right']].set_visible(False)
    # 1. A/T map by adenine strand
    a = ax[0]
    top, bot = pos.base == 'A', pos.base == 'T'
    a.vlines(x[top], 0.05, 0.95, color=blue, lw=2.2)
    a.vlines(x[bot], -0.95, -0.05, color=blue, lw=2.2, alpha=0.55)
    a.set_ylim(-1.1, 1.1)
    a.set_yticks([0.5, -0.5], ['A (top)', 'T (bottom A)'])
    for name, (s, e) in REGIONS.items():
        a.text((s + e) / 2, 1.25, name, ha='center', va='bottom', color=body, fontsize=12)
    # 2. evidence per A/T (caller view)
    a = ax[1]
    a.bar(x, pos.llr_miss, width=0.8, color=np.where(pos.is_at, blue, grey))
    a.set_ylabel('evidence if\nunmethylated\n(llr_miss, nats)')
    # 3. smallest callable footprint
    a = ax[2]
    a.step(x, pos.min_fp_seq_bp, where='mid', color=grey, lw=2, label=f'sequence only (>= {k} A/T)')
    a.step(x, pos.min_fp_llr_bp, where='mid', color=ink, lw=2, label=f'LLR caller (sum >= {min_llr} nats)')
    a.axhline(20, color=orange, lw=1, ls='--')
    a.text(x[-1] + 1, 20, 'TF\n20 bp', color=orange, va='center', fontsize=11)
    a.set_ylabel('smallest callable\nfootprint (bp)')
    a.set_ylim(0, max(60, np.nanmax(pos.min_fp_llr_bp) + 5))
    a.legend(loc='upper left', frameon=False, ncol=2, fontsize=12)
    # 4. variant classes by position
    a = ax[3]
    snv = var[var.cls == 'SNV']
    colors = {'A/T loss': orange, 'A<->T swap': ink, 'A/T gain': blue, 'G<->C swap': grey}
    ylev = {'A/T loss': 3, 'A<->T swap': 2, 'A/T gain': 1, 'G<->C swap': 0}
    for c, col in colors.items():
        s = snv[snv.seq_class == c]
        a.scatter(s.oligo_pos1, np.full(len(s), ylev[c]), s=22, color=col, label=c)
    a.set_yticks(list(ylev.values()), [f'{c} ({(snv.seq_class == c).sum()})' for c in ylev])
    a.set_ylim(-0.6, 3.6)
    a.set_xlabel('oligo position')
    a.text(x[0], 3.75, 'SNV class (A/T loss and A<->T swap both remove caller evidence)', color=body, fontsize=12)
    fig.savefig(out / 'a4_tracks.png', dpi=110, bbox_inches='tight', facecolor='#F7F6F2')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--emissions', required=True, help='plasmid-LLR emission table (fold A is used)')
    ap.add_argument('--min-llr', type=float, default=5.0)
    ap.add_argument('--pp-floor', type=float, default=0.0)
    ap.add_argument('--ledger', help='fsmpra ledger dir: add analysis-ready molecules per variant')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    lib = json.load(open(a.catalog))['library']
    wt = next(v['sequence'] for v in lib['variants'] if v['class'] == 'WT')
    if wt != lib['mutable_sequence']:
        raise ValueError('WT design differs from catalog mutable_sequence')
    em = emissions(a.emissions, lib['mutable_ref_start0'], wt, a.pp_floor)
    med = float(em.llr_miss.median())
    k = math.ceil(a.min_llr / med)
    pos = positions_table(wt, em, k, a.min_llr)
    var = variants_table(lib['variants'], wt, pos)
    if a.ledger:
        n = ledger_counts(a.ledger)
        var['molecules'] = var.variant_id.map(n).fillna(0).astype(int)
        wt_id = next(v['variant_id'] for v in lib['variants'] if v['class'] == 'WT')
        (out / 'wt_molecules.txt').write_text(f"{int(n.get(wt_id, 0))}\n")
    pos.to_csv(out / 'positions.tsv', sep='\t', index=False)
    var.to_csv(out / 'variants.tsv', sep='\t', index=False)
    s = summarize(wt, pos, var, k, a.min_llr, med)
    s['inputs'] = dict(catalog=a.catalog, emissions=a.emissions, ledger=a.ledger, pp_floor=a.pp_floor)
    json.dump(s, open(out / 'summary.json', 'w'), indent=1, default=str)
    figure(pos, var, out, k, a.min_llr)
    print(json.dumps(s, indent=1, default=str))


if __name__ == '__main__':
    main()

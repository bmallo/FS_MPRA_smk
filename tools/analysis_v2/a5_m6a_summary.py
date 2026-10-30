"""Step A5 (part 2): methylation-call QC tables from the a5_m6a_extract tallies (one directory per run).

Checks, each per run (cells kept separate):
  ml          distribution of ML over counted adenines by class (unlisted = no call)
  element     per-position m6A rate in U (false-positive proxy), C and M over the element, WT reads, >= 10 passes
  passes      U-element false-positive rate, M-element rate and class fractions by pass bin
  gatc        dam check: m6A at GATC adenines vs adenines in the same 3-mer (GAT top / ATC bottom, not in GATC)
              and vs all isolated adenines; distal positions >= 60 bp from the reference ends only
  orientation per-position M rate from forward- vs reverse-aligned reads (both strands of every A/T are called on
              every read, so the two should agree)
  strands     M rate of top- vs bottom-strand adenines by A-run class (informs the swap sensitivity rule)

Usage: python tools/analysis_v2/a5_m6a_summary.py --run HepG2_run3=<extract dir> [--run ...] --output <new dir>
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

_spec = importlib.util.spec_from_file_location('a4v', Path(__file__).with_name('a4_variant_visibility.py'))
a4v = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a4v)
_spec = importlib.util.spec_from_file_location('a5x', Path(__file__).with_name('a5_m6a_extract.py'))
a5x = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a5x)

ELIGIBLE_BINS = slice(1, None)           # pass bins >= 10 (rules.json min_passes)
END_EXCLUDE = 60                         # ft no-call zone near read ends is ~45 bp
COMP = str.maketrans('ACGT', 'TGCA')


def revcomp(s):
    return s.translate(COMP)[::-1]


def gatc_sites(ref):
    """Adenine positions of every GATC: top-strand A at i+1 (ref A), bottom-strand A at i+2 (ref T)."""
    top, bottom = [], []
    i = ref.find('GATC')
    while i >= 0:
        top.append(i + 1)
        bottom.append(i + 2)
        i = ref.find('GATC', i + 1)
    return top, bottom


def same_context_controls(ref):
    """Adenines in the GATC 3-mer context but not in a GATC: top A of GAT (not GATC), bottom A of ATC (not GATC)."""
    top = [i + 1 for i in range(len(ref) - 3) if ref[i:i + 3] == 'GAT' and ref[i:i + 4] != 'GATC']
    bottom = [i for i in range(1, len(ref) - 1) if ref[i - 1:i + 2] == 'ATC' and ref[max(0, i - 2):i + 2] != 'GATC']
    return top, bottom


def rate(counts, idx, k=1):
    """Pooled rate and per-position rates for positions idx; counts is (3 [obs, m128, m204], L)."""
    idx = np.asarray(idx, int)
    obs, m = counts[0, idx], counts[k, idx]
    per = np.where(obs > 0, m / np.maximum(obs, 1), np.nan)
    return (m.sum() / obs.sum() if obs.sum() else np.nan), per, int(obs.sum())


def load(d):
    z = np.load(Path(d) / 'counts.npz')
    meta = json.load(open(Path(d) / 'extract.json'))
    ref = a5x.ref_seq(meta['reference'])
    reads = pd.read_csv(Path(d) / 'reads.tsv.gz', sep='\t')
    return z, meta, ref, reads


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', action='append', required=True, help='LABEL=extract dir')
    ap.add_argument('--hepg2-element', default='HepG2_run3', help='run whose element orientation defines oligo coords')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    runs = dict(x.split('=', 1) for x in a.run)
    data = {k: load(v) for k, v in runs.items()}
    hep_el = None
    if a.hepg2_element in data:
        _, meta, ref, _ = data[a.hepg2_element]
        hep_el = ref[meta['element'][0]:meta['element'][1]]
    run_rows, el_rows, pass_rows, ml_rows, g_rows, gsite_rows, or_rows, st_rows = [], [], [], [], [], [], [], []
    for run, (z, meta, ref, reads) in data.items():
        C = z['counts']                                   # class, strand, pass bin, obs/m128/m204, pos
        L = len(ref)
        s, e = meta['element']
        el = ref[s:e]
        reverse = hep_el is not None and el != hep_el
        if reverse and revcomp(el) != hep_el:
            raise ValueError(f'{run}: element is neither the HepG2 element nor its reverse complement')
        oligo = (lambda p: e - 1 - p + 21) if reverse else (lambda p: p - s + 21)
        at = np.array([b in 'AT' for b in ref])
        elig = C[:, :, ELIGIBLE_BINS].sum(axis=(1, 2))      # class, 3, L
        rd = reads[reads.passes >= 10]
        # --- element per-position rates (WT reads only by construction of the extract)
        epos = [p for p in range(s, e) if at[p]]
        for p in epos:
            row = dict(run=run, ref_pos0=p, oligo_pos1=oligo(p), base=ref[p],
                       adenine_strand='top' if ref[p] == 'A' else 'bottom')
            for ci, cl in enumerate(a5x.CLASSES):
                row[f'n_{cl}'] = int(elig[ci, 0, p])
                row[f'rate_{cl}'] = elig[ci, 1, p] / elig[ci, 0, p] if elig[ci, 0, p] else np.nan
                row[f'rate_{cl}_204'] = elig[ci, 2, p] / elig[ci, 0, p] if elig[ci, 0, p] else np.nan
            el_rows.append(row)
        ed = pd.DataFrame([r for r in el_rows if r['run'] == run])
        u_pool, _, u_n = rate(elig[0], epos)
        u_pool204, _, _ = rate(elig[0], epos, 2)
        m_pool, _, _ = rate(elig[2], epos)
        m_pool204, _, _ = rate(elig[2], epos, 2)
        wt = rd[rd.wt == 1]
        run_rows.append(dict(
            run=run, sampled_primary=meta['sampled_primary'], fraction=meta['fraction'], kept=meta['kept'],
            **{f'status_{k}': v for k, v in meta['status'].items()},
            kept_ge10_passes=len(rd), wt_ge10_passes=len(wt),
            frac_U=float((rd['class'] == 'U').mean()), frac_C=float((rd['class'] == 'C').mean()),
            frac_M=float((rd['class'] == 'M').mean()), frac_zero_distal_calls=float((rd.distal_m6a == 0).mean()),
            wt_U_molecules=int((wt['class'] == 'U').sum()), wt_M_molecules=int((wt['class'] == 'M').sum()),
            element_U_rate=u_pool, element_U_rate_204=u_pool204, element_M_rate=m_pool, element_M_rate_204=m_pool204,
            element_U_obs=u_n, U_rate_median_position=float(ed.rate_U.median()),
            U_rate_p95_position=float(ed.rate_U.quantile(0.95)), U_rate_max_position=float(ed.rate_U.max()),
            U_positions_ge_0p02=int((ed.rate_U >= 0.02).sum()), element_at=len(epos),
            M_over_U_median_position=float((ed.rate_M / ed.rate_U.where(ed.rate_U > 0)).median())))
        # --- pass bins
        for b, lab in enumerate(a5x.PASS_LABELS):
            cb = C[:, :, b].sum(axis=1)
            rb = reads[reads.passes.apply(a5x.pass_bin) == b]
            wb = rb[rb.wt == 1]
            pass_rows.append(dict(run=run, passes=lab, reads=len(rb), wt_reads=len(wb),
                                  frac_U=float((rb['class'] == 'U').mean()) if len(rb) else np.nan,
                                  frac_M=float((rb['class'] == 'M').mean()) if len(rb) else np.nan,
                                  distal_m6a_median=float(rb.distal_m6a.median()) if len(rb) else np.nan,
                                  wt_U_molecules=int((wb['class'] == 'U').sum()),
                                  element_U_rate=rate(cb[0], epos)[0], element_U_obs=rate(cb[0], epos)[2],
                                  element_M_rate=rate(cb[2], epos)[0], element_C_rate=rate(cb[1], epos)[0]))
        # --- ML distribution (>= 10 passes)
        h = z['ml_hist'][:, ELIGIBLE_BINS].sum(axis=1)
        n = z['n_at'][:, ELIGIBLE_BINS].sum(axis=1)
        for ci, cl in enumerate(a5x.CLASSES):
            ml_rows.append(dict(run=run, cls=cl, adenines=int(n[ci]), frac_no_call=1 - h[ci].sum() / n[ci],
                                **{f'ml_{lo}_{hi - 1}': h[ci, lo:hi].sum() / n[ci]
                                   for lo, hi in [(0, 64), (64, 128), (128, 204), (204, 250), (250, 256)]}))
        # --- GATC (dam) check on distal positions away from the ends
        dmask = a5x.distal_mask(L, meta['element'], meta['barcode'])
        ok = dmask.copy()
        ok[:END_EXCLUDE] = ok[L - END_EXCLUDE:] = False
        gt, gb = gatc_sites(ref)
        ct, cb_ = same_context_controls(ref)
        cls_iso = a4v.adenine_classes(ref)
        iso = [k for (strand, k), c in cls_iso.items() if c[1] == 'isolated']
        gset = set(gt) | set(gb)
        groups = {'GATC': [p for p in gt + gb if ok[p]],
                  'same_3mer': [p for p in ct + cb_ if ok[p] and p not in gset],
                  'isolated_other': [p for p in iso if ok[p] and p not in gset]}
        for ci, cl in enumerate(a5x.CLASSES):
            pooled = {g: rate(elig[ci], idx) for g, idx in groups.items()}
            row = dict(run=run, cls=cl, gatc_sites=len(gt), gatc_adenines_used=len(groups['GATC']),
                       gatc_in_distal_total=sum(ok[p] for p in gt))
            for g, (pr, per, nobs) in pooled.items():
                row[f'{g}_n'] = len(groups[g])
                row[f'{g}_rate'] = pr
                row[f'{g}_median_site'] = float(np.nanmedian(per))
            row['gatc_over_same_3mer'] = row['GATC_rate'] / row['same_3mer_rate'] if row['same_3mer_rate'] else np.nan
            row['mannwhitney_p_gatc_vs_3mer'] = float(mannwhitneyu(pooled['GATC'][1], pooled['same_3mer'][1]).pvalue)
            g_rows.append(row)
            for g, idx in groups.items():
                for p, v in zip(idx, pooled[g][1]):
                    gsite_rows.append(dict(run=run, cls=cl, group=g, ref_pos0=p, rate=v, n=int(elig[ci, 0, p])))
        # --- orientation: forward- vs reverse-aligned reads, M class, every counted A/T
        fw, rv = C[2, 0, ELIGIBLE_BINS].sum(axis=0), C[2, 1, ELIGIBLE_BINS].sum(axis=0)
        use = at & (fw[0] >= 200) & (rv[0] >= 200)
        pf, pr_ = fw[1, use] / fw[0, use], rv[1, use] / rv[0, use]
        pool = (fw[1, use] + rv[1, use]) / (fw[0, use] + rv[0, use])
        se = np.sqrt(pool * (1 - pool) * (1 / fw[0, use] + 1 / rv[0, use]))
        zz = np.where(se > 0, (pf - pr_) / np.where(se > 0, se, 1), 0)
        or_rows.append(dict(run=run, positions=int(use.sum()), pearson=float(np.corrcoef(pf, pr_)[0, 1]),
                            median_abs_diff=float(np.median(np.abs(pf - pr_))), mean_diff_fw_minus_rv=float(np.mean(pf - pr_)),
                            positions_abs_z_gt_5=int((np.abs(zz) > 5).sum()),
                            reads_fw=int((rd.strand == '+').sum()), reads_rv=int((rd.strand == '-').sum())))
        # --- top vs bottom strand adenines by run class, M class, distal away from ends
        mrate = np.where(elig[2, 0] > 0, elig[2, 1] / np.maximum(elig[2, 0], 1), np.nan)
        iso_med = np.nanmedian([mrate[k] for (st, k), c in cls_iso.items() if c[1] == 'isolated' and ok[k]])
        tab = pd.DataFrame([dict(strand=st, run_bin=c[0], pos_in_run=c[1], rate=mrate[k])
                            for (st, k), c in cls_iso.items() if ok[k]])
        g = tab.groupby(['run_bin', 'strand']).rate.agg(['size', 'median']).reset_index()
        g['rel_to_isolated_median'] = g['median'] / iso_med
        st_rows.append(g.assign(run=run))
    pd.DataFrame(run_rows).to_csv(out / 'runs.tsv', sep='\t', index=False)
    pd.DataFrame(el_rows).to_csv(out / 'element_positions.tsv', sep='\t', index=False)
    pd.DataFrame(pass_rows).to_csv(out / 'passes.tsv', sep='\t', index=False)
    pd.DataFrame(ml_rows).to_csv(out / 'ml_distribution.tsv', sep='\t', index=False)
    pd.DataFrame(g_rows).to_csv(out / 'gatc.tsv', sep='\t', index=False)
    pd.DataFrame(gsite_rows).to_csv(out / 'gatc_sites.tsv', sep='\t', index=False)
    pd.DataFrame(or_rows).to_csv(out / 'orientation.tsv', sep='\t', index=False)
    pd.concat(st_rows).to_csv(out / 'strands.tsv', sep='\t', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        for f in ['runs', 'passes', 'ml_distribution', 'gatc', 'orientation', 'strands']:
            t = pd.read_csv(out / f'{f}.tsv', sep='\t')
            print(f'== {f}\n{(t.T if f == "runs" else t).round(5).to_string()}\n')


if __name__ == '__main__':
    main()

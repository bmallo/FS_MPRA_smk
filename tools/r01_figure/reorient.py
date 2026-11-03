"""R01 figure, step 1 (re-orient): collect what analysis v2 already says about each panel. No new statistics.

Reads stored analysis v2 outputs only (fp_features, A10, A11, fp60_86, HEK293T extract) and the published
HepG2 MPRA (Kircher 2019). Writes, to --output:
  regions.tsv      one row per HepG2 footprint F01-F24: box, WT occupancy, confident share, flag, SNVs whose
                   position lies in its span with their published activity, significant footprint changes
                   (A7 mask, confident calls) at that footprint by direction, and HEK293T coverage of SNVs in
                   its span (molecules, clones; per-SNV HEK statistics are not valid).
  snv_footprint.tsv every SNV x footprint pair significant with the A7 mask and confident calls (fp_features
                   'sig': BH across footprints within an SNV), with the +/-5 bp mask and all-call results for
                   the same pair, activity in both elements, distance from the SNV to the footprint span, the
                   SNV's caller class, and HEK293T molecules/clones. 'pattern' = activator-like when the
                   footprint change and the activity change have the same sign, repressor-like when opposite.
  motifs.tsv       A11 motif dose-response rows plus HEK293T clone-mean motif results (descriptive copy).
  leads.tsv        the SNVs named in the brief (60-86 alleles, A10 protection-lowering SNVs, top-activity SNVs)
                   with activity, A11 shifts and their footprint-change rows.
  summary.json     counts quoted in the step-1 report.

Usage: python tools/r01_figure/reorient.py --root results/analysis_v2_20261028 --output <new dir>
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

MPRA_DIR = Path('/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/mpra_published')
OLIGO_OFFSET = 11089264                     # oligo k = chr19:11,089,264 + k (hg38)
HEK_DEEP = 300                              # molecules: the brief's "about 105 SNVs with >= 300 molecules"
BRIEF_A10 = ['80G>C', '110G>C', '115A>C', '117A>C', '117A>T', '126A>T', '157A>G', '209G>C']
BRIEF_F04 = ['71T>A', '71T>C', '71T>G', '77A>T', '77A>G', '77A>C']


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def activity():
    """Published log2 effect and P per SNV for both elements (SNVs only; 1-bp deletions kept apart)."""
    out = []
    for el, f, sep in (('', 'GRCh38_LDLR.tsv', '\t'), ('2', 'GRCh38_LDLR.2.csv', ',')):
        t = pd.read_csv(MPRA_DIR / f, sep=sep)
        t['oligo'] = t['Position'] - OLIGO_OFFSET
        t = t[(t.oligo >= 21) & (t.oligo <= 231) & (t.Alt != '-')]
        t['snv'] = t.oligo.astype(str) + t.Ref + '>' + t.Alt
        out.append(t.set_index('snv')[['Value', 'P-Value', 'Tags']].rename(
            columns={'Value': f'act{el}', 'P-Value': f'act{el}_p', 'Tags': f'act{el}_tags'}))
    return out[0].join(out[1], how='outer')


def hek_coverage(root):
    c = pd.read_csv(root / 'a9_statistics/extract/HEK293T_run1/clones.tsv', sep='\t')
    c = c[c.label != 'WT']
    return c.groupby('label').agg(hek_molecules=('molecules', 'sum'), hek_clones=('clone', 'size'),
                                  hek_top_clone_share=('molecules', lambda m: m.max() / m.sum()))


def size_class(r):
    if r.len_hi <= 19:
        return 'sub-TF (10-19)'
    if r.len_lo <= 45 and r.len_hi >= 17:
        return 'TF (~20-40)'
    if r.len_hi <= 80:
        return 'PIC-size (41-80)'
    return 'long (81+)'


def span_distance(pos, lo, hi):
    return np.where(pos < lo, lo - pos, np.where(pos > hi, pos - hi, 0))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)

    inputs = {k: root / v for k, v in dict(
        footprints='fp_features/footprints.tsv', fp_effects='fp_features/effects.parquet',
        function='a11_motifs_function/function_snvs.tsv', dose='a11_motifs_function/dose.tsv',
        hek_motifs='a10_effects/HEK293T_run1/motifs.tsv', hek_clones='a9_statistics/extract/HEK293T_run1/clones.tsv',
        visibility='a4_informative/visibility/variants_visibility.tsv', fp60_86='fp60_86/variants.tsv').items()}
    inputs['mpra_LDLR'] = MPRA_DIR / 'GRCh38_LDLR.tsv'
    inputs['mpra_LDLR.2'] = MPRA_DIR / 'GRCh38_LDLR.2.csv'
    checksums = {k: dict(path=str(p), md5=md5(p)) for k, p in inputs.items()}

    fp = pd.read_csv(inputs['footprints'], sep='\t').set_index('feature')
    fp['size_class'] = fp.apply(size_class, axis=1)
    act = activity()
    fn = pd.read_csv(inputs['function'], sep='\t').set_index('snv')
    hek = hek_coverage(root)
    vis = pd.read_csv(inputs['visibility'], sep='\t')
    vis = vis[vis.cls == 'SNV'].set_index('source_variant')
    snv = fn[['pos', 'region', 'C1', 'C2', 'far_hits']].join(act).join(hek)
    snv['drops_its_base'] = vis['caller_drops_at'].reindex(snv.index).fillna(False).astype(bool)
    snv[['hek_molecules', 'hek_clones']] = snv[['hek_molecules', 'hek_clones']].fillna(0).astype(int)

    fe = pd.read_parquet(inputs['fp_effects'])
    key = ['snv', 'footprint']
    wide = None
    for (mask, calls), g in fe.groupby(['mask', 'calls']):
        g = g.set_index(key)[['occ', 'n', 'C', 'D', 'se', 'p', 'sig']].add_suffix(f'_{mask}_{calls}')
        wide = g if wide is None else wide.join(g, how='outer')
    wide = wide.reset_index()

    # regions
    rows = []
    for f, r in fp.iterrows():
        inspan = snv[(snv.pos >= r.span_lo) & (snv.pos <= r.span_hi)]
        w = wide[(wide.footprint == f) & wide['sig_A7_confident'].astype(bool)]
        dec, inc = w[w.D_A7_confident < 0], w[w.D_A7_confident > 0]
        rows.append(dict(
            footprint=f, size_class=r.size_class, flag=r.flag, span=f'{r.span_lo}-{r.span_hi}',
            mid_box=f'{r.oligo_lo}-{r.oligo_hi}', len_box=f'{r.len_lo}-{r.len_hi}',
            wt_occ_pct=100 * r.wt_occupancy, wt_occ_conf_pct=100 * r.wt_occupancy_confident,
            confident_share=r.confident_share,
            snvs_in_span=len(inspan), act_median_in_span=inspan.act.median(),
            act_min_in_span=inspan.act.min(), act_max_in_span=inspan.act.max(),
            n_act_up_p05=int(((inspan.act > 0) & (inspan.act_p < 0.05)).sum()),
            n_act_down_p05=int(((inspan.act < 0) & (inspan.act_p < 0.05)).sum()),
            sig_loss=len(dec), sig_gain=len(inc),
            sig_loss_median_act=dec.snv.map(snv.act).median(), sig_gain_median_act=inc.snv.map(snv.act).median(),
            sig_loss_pm5_also=int(dec['sig_pm5_confident'].fillna(False).astype(bool).sum()),
            sig_gain_pm5_also=int(inc['sig_pm5_confident'].fillna(False).astype(bool).sum()),
            hek_snvs_in_span_observed=int((inspan.hek_molecules > 0).sum()),
            hek_snvs_in_span_ge300=int((inspan.hek_molecules >= HEK_DEEP).sum()),
            hek_molecules_in_span=int(inspan.hek_molecules.sum()),
            hek_clones_median=float(inspan.loc[inspan.hek_molecules > 0, 'hek_clones'].median())
            if (inspan.hek_molecules > 0).any() else np.nan))
    regions = pd.DataFrame(rows)
    regions.to_csv(out / 'regions.tsv', sep='\t', index=False, float_format='%.4g')

    # SNV x footprint pairs significant with the A7 mask and confident calls
    sig = wide[wide['sig_A7_confident'].astype(bool)].copy()
    sig = sig.join(snv, on='snv')
    sig = sig.join(fp[['size_class', 'flag', 'span_lo', 'span_hi', 'confident_share']], on='footprint')
    sig['distance_to_span'] = span_distance(sig.pos.values, sig.span_lo.values, sig.span_hi.values)
    sig['pattern'] = np.where(np.sign(sig.D_A7_confident) == np.sign(sig.act), 'activator-like', 'repressor-like')
    sig.loc[sig.act.isna() | (sig.act == 0), 'pattern'] = 'no activity change'
    cols = ['snv', 'pos', 'footprint', 'size_class', 'flag', 'span_lo', 'span_hi', 'confident_share',
            'distance_to_span', 'drops_its_base', 'region',
            'occ_A7_confident', 'D_A7_confident', 'se_A7_confident', 'p_A7_confident', 'n_A7_confident',
            'C_A7_confident', 'D_pm5_confident', 'p_pm5_confident', 'sig_pm5_confident',
            'D_A7_all', 'p_A7_all', 'sig_A7_all', 'act', 'act_p', 'act2', 'act2_p', 'pattern', 'C1', 'C2',
            'hek_molecules', 'hek_clones', 'hek_top_clone_share']
    sig = sig[cols].sort_values(['footprint', 'pos', 'snv'])
    sig.to_csv(out / 'snv_footprint.tsv', sep='\t', index=False, float_format='%.4g')

    # motifs (descriptive copies)
    dose = pd.read_csv(inputs['dose'], sep='\t')
    hm = pd.read_csv(inputs['hek_motifs'], sep='\t')
    pd.concat([dose.assign(source='A11 dose-response (HepG2)'),
               hm.assign(source='A10 clone-mean (HEK293T)')], ignore_index=True).to_csv(
        out / 'motifs.tsv', sep='\t', index=False, float_format='%.4g')

    # leads named in the brief
    top_act = snv.sort_values('act', ascending=False).head(15).index.tolist()
    lead_sets = [('F04 alleles', BRIEF_F04), ('A10 protection-lowering', BRIEF_A10), ('top activity', top_act),
                 ('Repeat 3 (143-151)', snv.index[snv.region.astype(bool)].tolist())]
    lrows = []
    for name, lst in lead_sets:
        for v in lst:
            s = snv.loc[v] if v in snv.index else None
            w = wide[(wide.snv == v)]
            best = w.loc[w['p_A7_confident'].idxmin()] if w['p_A7_confident'].notna().any() else None
            sigs = w[w['sig_A7_confident'].astype(bool)]
            lrows.append(dict(
                set=name, snv=v, pos=s.pos if s is not None else np.nan,
                act=s.act if s is not None else np.nan, act_p=s.act_p if s is not None else np.nan,
                act2=s.act2 if s is not None else np.nan, C1=s.C1 if s is not None else np.nan,
                C2=s.C2 if s is not None else np.nan, drops_its_base=bool(s.drops_its_base) if s is not None else None,
                sig_footprints_A7_conf=';'.join(f'{x.footprint}:{100 * x.D_A7_confident:+.2f}'
                                                f'{"*" if x.sig_pm5_confident else ""}' for x in sigs.itertuples()),
                min_p_footprint=best.footprint if best is not None else '',
                min_p_D_pts=100 * best.D_A7_confident if best is not None else np.nan,
                min_p=best.p_A7_confident if best is not None else np.nan,
                hek_molecules=s.hek_molecules if s is not None else 0,
                hek_clones=s.hek_clones if s is not None else 0))
    pd.DataFrame(lrows).to_csv(out / 'leads.tsv', sep='\t', index=False, float_format='%.4g')

    # repressor scan: SNVs with raised activity (log2 > 0.5, P < 0.05) x footprints within 20 bp of the SNV,
    # significant or not (direction only; most per-SNV footprint tests are underpowered)
    up = snv[(snv.act > 0.5) & (snv.act_p < 0.05)]
    rs = wide[wide.snv.isin(up.index)].join(snv[['pos', 'act', 'act2', 'C1', 'C2', 'drops_its_base', 'hek_molecules',
                                                'hek_clones']], on='snv')
    rs = rs.join(fp[['size_class', 'span_lo', 'span_hi', 'confident_share']], on='footprint')
    rs['distance_to_span'] = span_distance(rs.pos.values, rs.span_lo.values, rs.span_hi.values)
    rs = rs[(rs.distance_to_span <= 20) & rs.size_class.isin(['TF (~20-40)', 'sub-TF (10-19)'])]
    rs[['snv', 'pos', 'act', 'act2', 'C1', 'C2', 'drops_its_base', 'footprint', 'size_class', 'confident_share',
        'distance_to_span', 'D_A7_confident', 'p_A7_confident', 'sig_A7_confident', 'D_pm5_confident',
        'p_pm5_confident', 'D_A7_all', 'p_A7_all', 'hek_molecules', 'hek_clones']].sort_values(
        ['pos', 'snv', 'footprint']).to_csv(out / 'repressor_scan.tsv', sep='\t', index=False, float_format='%.4g')

    summ = dict(
        inputs=checksums,
        snvs=len(snv), snvs_with_activity=int(snv.act.notna().sum()),
        act_up_p05=int(((snv.act > 0) & (snv.act_p < 0.05)).sum()),
        act_down_p05=int(((snv.act < 0) & (snv.act_p < 0.05)).sum()),
        act_up_gt0p5=int((snv.act > 0.5).sum()), act_down_lt_m0p5=int((snv.act < -0.5).sum()),
        sig_pairs_A7_conf=len(sig), sig_snvs_A7_conf=int(sig.snv.nunique()),
        sig_pairs_by_footprint=sig.footprint.value_counts().to_dict(),
        sig_pairs_by_pattern=sig.pattern.value_counts().to_dict(),
        sig_pairs_pm5_also=int(sig.sig_pm5_confident.fillna(False).astype(bool).sum()),
        sig_pairs_beyond_10bp=int((sig.distance_to_span > 10).sum()),
        hek_snvs_observed=int((snv.hek_molecules > 0).sum()),
        hek_snvs_ge300=int((snv.hek_molecules >= HEK_DEEP).sum()),
        hek_snvs_ge3clones=int((snv.hek_clones >= 3).sum()))
    json.dump(summ, open(out / 'summary.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()

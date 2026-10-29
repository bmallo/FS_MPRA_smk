"""Join function with chromatin architecture per SNV, label candidate activators/repressors, score motifs.

Inputs (config/ldlr_footprints.json plus the per-cell result directories):
- published activity, two sources (r = 0.97, so one measurement): effect = mean log2 value, hit =
  BH q < 0.05 in both sources with the same sign;
- HepG2 site effects with LLR coverage (accessibility msp_C_std, nucleosome-sized nuc_85_250_C_std;
  primary-family q) and HepG2 footprint effects (tf:fp_local:net etc.; footprint-family q);
- HEK293T effects and footprints for the shared SNVs (replication columns, descriptive for footprints).

Labels (config 'classification'): a footprint label when the local TF-footprint change is
significant, otherwise an accessibility label. Loss of a footprint with lower activity marks a
candidate activator site; loss with higher activity a candidate repressor site; a gained footprint
marks a candidate created site. These are hypotheses about the bound factor's role, not
identifications.

Motifs: JASPAR 2024 + HOCOMOCO v12 PWMs (extras/data) scanned on the WT element with the local
exact-tail scanner (fsmpra.biology_motifs; per-window p <= 1e-4). For each strongest variant (a
functional hit with any significant chromatin change) every overlapping WT site gets a REF->ALT
score change, and the ALT sequence is scanned for sites created at the substitution. For each WT
site overlapped by a strongest variant, Spearman of score change vs activity and vs chromatin
across all SNVs inside the site ranks which motif best explains the region.

Usage: python -m fsmpra.variant_classes --config config/ldlr_footprints.json --dir results/variant_effects_20261028
"""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy import stats

from . import biology_motifs as bm
from .io import atomic_output, read_json, write_json, write_table
from .site_effects import GENOMIC_OFFSET, bh

MOTIF_FILES = [('extras/data/JASPAR2024_CORE_vertebrates_nr.jaspar', 'JASPAR2024'),
               ('extras/data/HOCOMOCOv12_H12CORE_pcm.txt', 'HOCOMOCOv12')]


def read_effects(path, endpoints, q_field):
    out = {}
    with gzip.open(path, 'rt') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            if r['endpoint'] in endpoints:
                f = lambda k: float(r[k]) if r.get(k) not in ('', None) else np.nan
                out[(r['source_variant'], r['endpoint'])] = dict(delta=f('delta'), z=f('z'), q=f(q_field),
                                                                 molecules=f('molecules'), context=r['context'],
                                                                 at_change=r['at_change'])
    return out


def observability_residual(path, track, window):
    """Footprint change not predicted by the near-naked (M) change: r = dC - b * dM, b = library slope of dC on dM.

    Post hoc (2026-10-28): the prespecified net = dC - dM assumed slope 1, but footprint changes on M
    molecules run 2-3x those on chromatin molecules (the observability kernel scales with labeling),
    so C - M turned shared losses into apparent gains. SE ignores the C/M covariance from shared
    barcode clusters (positive), which makes it conservative."""
    rows = {}
    with gzip.open(path, 'rt') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            if r['track'] == track and r['window'] == window and r['measure'] in ('C_std', 'M') and r['delta'] not in ('', None) \
                    and r['se'] not in ('', None):
                rows.setdefault(r['source_variant'], {})[r['measure']] = (float(r['delta']), float(r['se']))
    keys = [k for k, v in rows.items() if len(v) == 2]
    C = np.array([rows[k]['C_std'][0] for k in keys]); M = np.array([rows[k]['M'][0] for k in keys])
    sC = np.array([rows[k]['C_std'][1] for k in keys]); sM = np.array([rows[k]['M'][1] for k in keys])
    b = float(np.polyfit(M, C, 1)[0])
    r = C - b * M
    z = r / np.sqrt(sC ** 2 + (b * sM) ** 2)
    p = 2 * stats.norm.sf(np.abs(z))
    q = bh(p)
    return {k: dict(delta=r[i], z=z[i], q=q[i], c=C[i], m=M[i]) for i, k in enumerate(keys)}, b


def activity(paths):
    act = {}
    for src, path in paths.items():
        delim = ',' if str(path).endswith('.csv') else '\t'
        rows = [r for r in csv.DictReader(open(path), delimiter=delim) if r['Alt'] in ('A', 'C', 'G', 'T')]
        q = bh([float(r['P-Value']) for r in rows])
        for r, qq in zip(rows, q):
            act.setdefault((int(r['Position']), r['Ref'], r['Alt']), {})[src] = (float(r['Value']), qq)
    return act


def label(fp_q, fp_d, msp_q, msp_d, func_hit, a, alpha):
    if not func_hit:
        return ''
    down = a < 0
    if np.isfinite(fp_q) and fp_q < alpha:
        if fp_d < 0:
            return 'activator_footprint' if down else 'repressor_footprint'
        return 'created_repressor_footprint' if down else 'created_activator_footprint'
    if np.isfinite(msp_q) and msp_q < alpha:
        if (msp_d < 0) == down:
            return 'activator_accessibility' if down else 'repressor_accessibility'
        return 'discordant_accessibility'
    return ''


class Scorer:
    """Cached PWM tables so ALT windows can be scanned without recomputing score distributions."""

    def __init__(self, matrices):
        uniform = np.full(4, .25)
        self.m = []
        for rec in matrices:
            w = bm.pwm(rec['counts'], uniform)
            q, off, tail = bm.score_distribution(w, uniform)
            self.m.append((rec, w, q, off, tail))

    def best_over(self, seq, p0, ref_seq):
        """Sites overlapping p0 in seq that pass p <= 1e-4 while the same matrix/start/strand does not in ref_seq."""
        s = np.array([bm.BASES.index(b) for b in seq])
        r = np.array([bm.BASES.index(b) for b in ref_seq])
        out = []
        for rec, w, q, off, tail in self.m:
            L = len(w)
            for start in range(max(0, p0 - L + 1), min(p0, len(s) - L) + 1):
                for strand in '+-':
                    win, rwin = s[start:start + L], r[start:start + L]
                    if strand == '-':
                        win, rwin = 3 - win[::-1], 3 - rwin[::-1]
                    p = tail[q[np.arange(L), win].sum() - off]
                    if p <= 1e-4 and tail[q[np.arange(L), rwin].sum() - off] > 1e-4:
                        out.append(dict(database=rec['database'], matrix_id=rec['matrix_id'], name=rec['name'], start0=start, length=L,
                                        strand=strand, alt_p=float(p), ref_p=float(tail[q[np.arange(L), rwin].sum() - off]),
                                        alt_bits=float(w[np.arange(L), win].sum()), ref_bits=float(w[np.arange(L), rwin].sum())))
        return out


def run(config, root, output):
    cfg = read_json(config)
    alpha = cfg['classification']['chromatin_q']
    V = Path(root)
    variants = read_json(V / 'hepg2/extract_llr/variants.json')
    wt = next(v['sequence'] for v in variants if v['class'] == 'WT')
    snvs = [v for v in variants if v['class'] == 'SNV']
    arch = ['msp_C_std', 'nuc_85_250_C_std']
    fpe = ['tf:fp_local:net', 'tf:fp_flank:net', 'nuc:fp_local:net', 'tf:fp_local:C_std', 'tf:fp_local:M']
    hep_site = read_effects(V / 'hepg2/effects/effects.tsv.gz', arch, 'q_bh_primary_family')
    hep_fp = read_effects(V / 'hepg2/footprints/effects.tsv.gz', fpe, 'q_bh_footprint_family')
    resid, slopes = {}, {}
    for tr in ('tf', 'nuc'):
        resid[tr], slopes[tr] = observability_residual(V / 'hepg2/footprints/effects.tsv.gz', tr, 'fp_local')
    hek_site = read_effects(V / 'hek293t/effects/effects.tsv.gz', arch, 'q_bh_primary_family')
    hek_fp = read_effects(V / 'hek293t/footprints/effects.tsv.gz', fpe, 'q_bh_footprint_family')
    gates = {}
    for cell, sub in (('hepg2', 'effects'), ('hepg2', 'footprints'), ('hek293t', 'effects')):
        with open(V / cell / sub / 'calibration_summary.tsv') as h:
            for r in csv.DictReader(h, delimiter='\t'):
                gates[(cell, r['endpoint'])] = r['gate_pass'] == 'True'
    act = activity({k: Path(v) for k, v in cfg['published_activity'].items()})
    rows = []
    for v in snvs:
        sv = v['source_variant']
        a = act.get((GENOMIC_OFFSET + v['oligo_pos1'], v['ref'], v['alt']), {})
        vals = [a[k][0] for k in a]
        a_mean = float(np.mean(vals)) if vals else np.nan
        func_hit = len(a) == 2 and all(x[1] < .05 for x in a.values()) and len({np.sign(x[0]) for x in a.values()}) == 1
        g = lambda d, ep, k: d.get((sv, ep), {}).get(k, np.nan)
        row = dict(source_variant=sv, oligo_pos1=v['oligo_pos1'], genomic_pos1=GENOMIC_OFFSET + v['oligo_pos1'], ref=v['ref'], alt=v['alt'],
                   context=hep_site.get((sv, 'msp_C_std'), {}).get('context', ''),
                   at_change=hep_site.get((sv, 'msp_C_std'), {}).get('at_change', ''),
                   activity=a_mean, activity_q_max=max([x[1] for x in a.values()], default=np.nan), functional_hit=func_hit,
                   hepg2_molecules=g(hep_site, 'msp_C_std', 'molecules'))
        for name, d, ep in (('msp', hep_site, 'msp_C_std'), ('nuc', hep_site, 'nuc_85_250_C_std'),
                            ('tf_local_net_prespecified', hep_fp, 'tf:fp_local:net'), ('nuc_local_net_prespecified', hep_fp, 'nuc:fp_local:net')):
            row[f'hepg2_{name}_delta'], row[f'hepg2_{name}_z'], row[f'hepg2_{name}_q'] = g(d, ep, 'delta'), g(d, ep, 'z'), g(d, ep, 'q')
        for tr, name in (('tf', 'tf_local'), ('nuc', 'nuc_local')):
            x = resid[tr].get(sv, {})
            row[f'hepg2_{name}_delta'], row[f'hepg2_{name}_z'], row[f'hepg2_{name}_q'] = (x.get('delta', np.nan), x.get('z', np.nan),
                                                                                       x.get('q', np.nan))
        row['hepg2_tf_local_C_std_delta'] = g(hep_fp, 'tf:fp_local:C_std', 'delta')
        row['hepg2_tf_local_M_delta'] = g(hep_fp, 'tf:fp_local:M', 'delta')
        row['hek293t_molecules'] = g(hek_site, 'msp_C_std', 'molecules')
        for name, d, ep in (('msp', hek_site, 'msp_C_std'), ('nuc', hek_site, 'nuc_85_250_C_std'), ('tf_local', hek_fp, 'tf:fp_local:net')):
            row[f'hek293t_{name}_delta'] = g(d, ep, 'delta')
        row['hek293t_msp_z'] = g(hek_site, 'msp_C_std', 'z')
        chrom_hit = any(np.isfinite(row[f'hepg2_{k}_q']) and row[f'hepg2_{k}_q'] < alpha for k in ('msp', 'nuc', 'tf_local', 'nuc_local'))
        row['chromatin_hit'] = chrom_hit
        row['class'] = label(row['hepg2_tf_local_q'], row['hepg2_tf_local_delta'], row['hepg2_msp_q'], row['hepg2_msp_delta'],
                             func_hit, a_mean, alpha)
        row['strongest'] = func_hit and chrom_hit
        rows.append(row)

    # Position summary: classes and consistency across the three substitutions.
    positions = []
    for p in sorted({r['oligo_pos1'] for r in rows}):
        rs = [r for r in rows if r['oligo_pos1'] == p]
        cls = [r['class'] for r in rs if r['class']]
        positions.append(dict(oligo_pos1=p, genomic_pos1=GENOMIC_OFFSET + p, wt_base=wt[p - 21], snvs=len(rs),
                              functional_hits=sum(r['functional_hit'] for r in rs), classed=len(cls),
                              classes=','.join(sorted(set(cls))), consensus=cls[0] if cls and len(set(cls)) == 1 and len(cls) >= 2 else '',
                              activity_mean=float(np.nanmean([r['activity'] for r in rs])),
                              msp_delta_mean=float(np.nanmean([r['hepg2_msp_delta'] for r in rs])),
                              tf_local_delta_mean=float(np.nanmean([r['hepg2_tf_local_delta'] for r in rs]))))

    # Motifs.
    matrices = []
    for path, db in MOTIF_FILES:
        matrices += bm.read_matrices(path, db)
    hits, settings = bm.scan(wt, matrices, GENOMIC_OFFSET + 21)
    by_sv = {r['source_variant']: r for r in rows}
    strongest = [r for r in rows if r['strongest']]
    broken = []
    for r in strongest:
        p0 = r['oligo_pos1'] - 21
        for h in hits:
            d = bm.variant_delta(h, p0, r['ref'], r['alt'])
            if d is not None:
                broken.append(dict(source_variant=r['source_variant'], oligo_pos1=r['oligo_pos1'], class_=r['class'],
                                   activity=r['activity'], database=h['database'], matrix_id=h['matrix_id'], name=h['name'],
                                   site_oligo=f"{h['start0'] + 21}-{h['end0'] + 20}", strand=h['strand'],
                                   wt_site_p=h['sequence_tail_p'], gc_matched_p=h['gc_matched_tail_p'], score_change_bits=d))
    scorer = Scorer(matrices)
    created = []
    for r in strongest:
        p0 = r['oligo_pos1'] - 21
        lo, hi = max(0, p0 - 30), min(211, p0 + 31)
        ref_w = wt[lo:hi]
        alt_w = ref_w[:p0 - lo] + r['alt'] + ref_w[p0 - lo + 1:]
        for c in scorer.best_over(alt_w, p0 - lo, ref_w):
            created.append(dict(source_variant=r['source_variant'], oligo_pos1=r['oligo_pos1'], class_=r['class'], activity=r['activity'],
                                database=c['database'], matrix_id=c['matrix_id'], name=c['name'],
                                site_oligo=f"{lo + c['start0'] + 21}-{lo + c['start0'] + c['length'] + 20}",
                                strand=c['strand'], alt_p=c['alt_p'], ref_p=c['ref_p'], score_change_bits=c['alt_bits'] - c['ref_bits']))
    # Rank WT sites hit by strongest variants: does disrupting the site track activity / chromatin across all its SNVs?
    site_rank = []
    touched = {(b['database'], b['matrix_id'], b['site_oligo'], b['strand']) for b in broken}
    for h in hits:
        key = (h['database'], h['matrix_id'], f"{h['start0'] + 21}-{h['end0'] + 20}", h['strand'])
        if key not in touched:
            continue
        xs, ya, ym, yt = [], [], [], []
        for r in rows:
            d = bm.variant_delta(h, r['oligo_pos1'] - 21, r['ref'], r['alt'])
            if d is None or not np.isfinite(r['activity']):
                continue
            xs.append(d); ya.append(r['activity']); ym.append(r['hepg2_msp_delta']); yt.append(r['hepg2_tf_local_delta'])
        if len(xs) < 6:
            continue
        sp = lambda y: stats.spearmanr(xs, y, nan_policy='omit')
        sa, sm, st = sp(ya), sp(ym), sp(yt)
        site_rank.append(dict(database=h['database'], matrix_id=h['matrix_id'], name=h['name'], site_oligo=key[2], strand=h['strand'],
                              site_genomic=f"chr19:{h['genomic_start1']}-{h['genomic_end1']}", wt_site_p=h['sequence_tail_p'],
                              gc_matched_p=h['gc_matched_tail_p'], snvs=len(xs),
                              rho_activity=float(sa.correlation), p_activity=float(sa.pvalue),
                              rho_accessibility=float(sm.correlation), p_accessibility=float(sm.pvalue),
                              rho_tf_footprint=float(st.correlation), p_tf_footprint=float(st.pvalue),
                              strongest_snvs=','.join(sorted({b['source_variant'] for b in broken if (b['database'], b['matrix_id'], b['site_oligo'], b['strand']) == key}))))
    site_rank.sort(key=lambda s: s['p_activity'])

    def clean(rs):
        return [{k: (None if isinstance(x, (float, np.floating)) and not np.isfinite(x) else
                     bool(x) if isinstance(x, np.bool_) else x) for k, x in r.items()} for r in rs]
    with atomic_output(output) as out:
        write_table(out / 'variants.tsv', clean(rows), list(rows[0]))
        write_table(out / 'positions.tsv', clean(positions), list(positions[0]))
        if broken:
            write_table(out / 'motifs_broken.tsv', clean(broken), list(broken[0]))
        if created:
            write_table(out / 'motifs_created.tsv', clean(created), list(created[0]))
        if site_rank:
            write_table(out / 'motif_sites_ranked.tsv', clean(site_rank), list(site_rank[0]))
        counts = {}
        for r in rows:
            if r['class']:
                counts[r['class']] = counts.get(r['class'], 0) + 1
        write_json(out / 'summary.json', dict(
            schema_version='fsmpra.variant_classes.v1', snvs=len(rows), functional_hits=sum(r['functional_hit'] for r in rows),
            chromatin_hits=sum(r['chromatin_hit'] for r in rows), strongest=len(strongest), classes=counts,
            gates={f'{c}:{e}': p for (c, e), p in gates.items() if e in arch + fpe}, motif_scan=settings,
            observability_slopes=slopes, footprint_measure='slope-scaled residual dC - b*dM (post hoc; see observability_residual)',
            wt_motif_sites=len(hits), rules=cfg['classification'], functional=cfg['functional'],
            interpretation='Labels are hypotheses from sign logic (loss of footprint or accessibility with activity change); '
                           'they do not identify the bound protein or its mechanism. Motif tails are unadjusted per-window '
                           'sequence matches.'))
    print(json.dumps(dict(functional=sum(r['functional_hit'] for r in rows), chromatin=sum(r['chromatin_hit'] for r in rows),
                          strongest=len(strongest), classes=counts), indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', required=True)
    ap.add_argument('--dir', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    run(a.config, a.dir, a.output)


if __name__ == '__main__':
    main()

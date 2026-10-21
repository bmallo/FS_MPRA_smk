"""Evidence figures for the top replicate-3 variant hits.

Three steps, so the heavy part runs once on Klone and plotting can be redone anywhere:

    # here or on Klone: motif scores for WT and ALT at each hit (both strands, both databases)
    python -m tools.variant_hit_figures motifs --output docs/evidence/top_hits_20261021/motif_changes.tsv
    # Klone only (needs the site extraction): per-base profiles for all 634 alleles + molecule rasters
    python -m tools.variant_hit_figures profiles --extract results/rep3_site_methylation_20261021/extract \\
        --output results/rep3_site_methylation_20261021/hit_profiles
    # anywhere: one figure per hit plus an overview
    python -m tools.variant_hit_figures figures --profiles <hit_profiles dir> --output <figures dir>

Descriptive figures from one preparation. Motif scores are sequence matches, not binding.
"""
import argparse
import csv
import gzip
import textwrap
from pathlib import Path
import numpy as np

GENOMIC_OFFSET = 11089264
EVIDENCE = Path('docs/evidence')
EFFECTS = EVIDENCE / 'site_methylation_mask4_20261021' / 'effects' / 'effects.tsv.gz'
TOP_HITS = ['146C>G', '146C>T', '149C>G', '98C>G', '129C>G', '133C>G', '56G>C', '56G>T', '58G>T', '117A>T', '117A>C', '115A>C',
            '191G>T', '73T>C', '125A>T']
COMPARE_POSITIONS = [56, 58, 73, 117, 146]   # all three substitutions side by side
HIT_POSITIONS = [54, 56, 58, 59, 73, 74, 98, 115, 117, 125, 129, 133, 146, 149, 191]   # motif table covers every substitution here
ROLES = {'146C>G': 'Repeat 3 (Sp1): accessibility loss', '146C>T': 'Repeat 3 (Sp1): accessibility loss, adds A/T',
         '149C>G': 'Repeat 3 (Sp1): accessibility loss', '98C>G': 'Repeat 1 (Sp1): accessibility loss',
         '129C>G': 'SRE-1: activity loss without accessibility change; TF-sized call loss, proxy also moves',
         '133C>G': 'SRE-1: activity loss without accessibility loss',
         '56G>C': 'P1 core: left-half opening with activity loss; TF-sized calls unchanged',
         '56G>T': 'P1 core: left-half opening with activity loss; adds an A/T inside the TF-sized footprint',
         '117A>T': 'Oligo 117: opening with activity gain; creates AP-1/CRE-like match',
         '117A>C': 'Oligo 117: opening with smaller activity gain; no motif created',
         '115A>C': 'Oligo 115: opening with activity gain; creates AP-1/CRE-like match',
         '191G>T': 'Oligo 191: local opening without function',
         '58G>T': 'P1 core: activity loss with small chromatin change',
         '73T>C': 'Artifact example: methylated proxy drops more than chromatin at the A run 8-10 bp away',
         '125A>T': 'Oligo 125 (A run): observability control; largest TF-sized call loss, no function'}
# Named elements and the motif model that best matches each in WT (chosen after seeing the data); oligo search range.
MOTIF_SITES = [('P1 core', 'HOCOMOCOv12', 'ZN652', (45, 70)), ('Repeat 1 (Sp1)', 'HOCOMOCOv12', 'SP4', (82, 108)),
               ('SRE-1', 'JASPAR2024', 'SREBF1', (120, 140)), ('Repeat 3 (Sp1)', 'HOCOMOCOv12', 'WT1', (138, 158))]
STATES = ['unannotated', 'accessible (MSP)', 'nucleosome-sized 85-250 bp', 'long protection >250 bp', 'TF-sized footprint 10-60 bp']


def wt_sequence():
    with open(EVIDENCE / '20261020' / 'validation_panel.tsv') as h:
        return next(r for r in csv.DictReader(h, delimiter='\t') if r['source_variant'] == 'WT')['mutable_sequence']


def parse(name):
    return int(name[:-3]), name[-3], name[-1]


# ---------------------------------------------------------------- motifs
def motif_changes(variants, threshold=1e-3):
    """Best WT and ALT match of every matrix over windows that contain the substituted base."""
    from fsmpra.biology_motifs import BASES, read_matrices, pwm, score_distribution
    seq = wt_sequence()
    mats = read_matrices('extras/data/JASPAR2024_CORE_vertebrates_nr.jaspar', 'JASPAR2024') + \
        read_matrices('extras/data/HOCOMOCOv12_H12CORE_pcm.txt', 'HOCOMOCOv12')
    uniform = np.full(4, .25)
    prepared = []
    for m in mats:
        w = pwm(m['counts'], uniform)
        q, offset, tail = score_distribution(w, uniform)
        prepared.append((m, w, q, offset, tail))
    rows = []
    for v in variants:
        pos, ref, alt = parse(v)
        p0 = pos - 21
        if seq[p0] != ref:
            raise ValueError('Reference mismatch for ' + v)
        alleles = {'wt': np.array([BASES.index(b) for b in seq]),
                   'alt': np.array([BASES.index(b) for b in seq[:p0] + alt + seq[p0 + 1:]])}
        for m, w, q, offset, tail in prepared:
            L = len(w)
            starts = range(max(0, p0 - L + 1), min(p0, 211 - L) + 1)
            best = {}
            for allele, s in alleles.items():
                top = None
                for start in starts:
                    window = s[start:start + L]
                    for strand, frag in (('+', window), ('-', 3 - window[::-1])):
                        total = int(q[np.arange(L), frag].sum())
                        p = float(tail[total - offset])
                        if top is None or p < top[0]:
                            top = (p, float(w[np.arange(L), frag].sum()), start, strand)
                best[allele] = top
            if best and min(best['wt'][0], best['alt'][0]) <= threshold:
                lw, la = np.log10(best['wt'][0]), np.log10(best['alt'][0])
                start, strand = best['wt'][2:] if best['wt'][0] <= best['alt'][0] else best['alt'][2:]
                rows.append(dict(variant=v, database=m['database'], matrix_id=m['matrix_id'], name=m['name'],
                                 wt_p=best['wt'][0], alt_p=best['alt'][0], wt_bits=round(best['wt'][1], 2),
                                 alt_bits=round(best['alt'][1], 2), log10_p_change=round(la - lw, 2),
                                 change='gain' if la < lw else 'loss' if la > lw else 'none',
                                 oligo_start=start + 21, oligo_end=start + L + 20, strand=strand,
                                 wt_site=seq[start:start + L]))
    rows.sort(key=lambda r: (r['variant'], -abs(r['log10_p_change'])))
    return rows


# ---------------------------------------------------------------- profiles (Klone)
def profiles(extract, output, config='config/ldlr_site_methylation_mask4.json', raster_n=150, seed=20261021):
    from fsmpra.io import atomic_output, read_json, write_json
    from fsmpra.site_effects import COHORTS, assign_groups, grouped_sum, ratio
    extract = Path(extract)
    spec = read_json(config)['analysis']
    variants = read_json(extract / 'variants.json')
    wt_i = next(i for i, v in enumerate(variants) if v['class'] == 'WT')
    meta = dict(np.load(extract / 'meta.npz'))
    key = COHORTS[spec['cohort']]                         # same cohort as the effects run; unevaluated tracks are missing
    cohort = meta[key].astype(bool) if key else np.ones(len(meta['variant_index']), bool)
    cohort &= meta.get('evaluated', np.ones(len(cohort), bool)).astype(bool)
    vi = meta['variant_index'].astype(int)
    groups, weights, edges = assign_groups(meta['distal_m6a'].astype(float), (vi == wt_i) & cohort,
                                           spec['unmethylated_max_distal_m6a'], spec['methylated_min_distal_m6a'], spec['labeling_strata'])
    G, V = 2 + len(weights), len(variants)
    arr = {n: np.load(extract / f'{n}.npy', mmap_mode='r') for n in ('scores', 'known', 'qualities', 'coverage')}
    keep = np.flatnonzero(cohort & (groups >= 0))
    acc = {k: np.zeros((V * G, 211)) for k in ('called', 'valid', 'msp', 'nuc', 'long', 'tf')}
    cnt = np.zeros(V * G)
    for a in range(0, len(keep), 200000):
        rows = keep[a:a + 200000]
        key = vi[rows] * G + groups[rows]
        valid = arr['known'][rows] & (arr['qualities'][rows] >= spec['baseq_min'])
        cov = arr['coverage'][rows]
        parts = dict(called=valid & (arr['scores'][rows] >= spec['ml_min']), valid=valid, msp=(cov & 2) > 0,
                     nuc=(cov & 4) > 0, long=(cov & 8) > 0, tf=(cov & 1) > 0)
        for k, x in parts.items():
            acc[k] += grouped_sum(key, x, V * G)
        cnt += np.bincount(key, minlength=V * G)
    acc = {k: x.reshape(V, G, 211) for k, x in acc.items()}
    cnt = cnt.reshape(V, G)
    out = {}
    # Same minimum depths as the calibrated effects: every labeling stratum and each proxy class.
    strata_ok = (cnt[:, 2:] >= spec['min_stratum_molecules']).all(1)[:, None]
    enough_u = (cnt[:, 0] >= spec['min_class_molecules'])[:, None]
    enough_m = (cnt[:, 1] >= spec['min_class_molecules'])[:, None]
    with np.errstate(invalid='ignore', divide='ignore'):
        for k in ('msp', 'nuc', 'long', 'tf'):
            per = ratio(acc[k][:, 2:], cnt[:, 2:, None])
            out[f'{k}_C_std'] = np.where(strata_ok, np.einsum('vsj,s->vj', per, weights), np.nan)
            out[f'{k}_ALL'] = ratio(acc[k].sum(1), cnt.sum(1)[:, None])
        per = ratio(acc['called'][:, 2:], acc['valid'][:, 2:])
        out['m6a_C_std'] = np.where(strata_ok, np.einsum('vsj,s->vj', per, weights), np.nan)
        out['m6a_U'] = np.where(enough_u, ratio(acc['called'][:, 0], acc['valid'][:, 0]), np.nan)
        out['m6a_M'] = np.where(enough_m, ratio(acc['called'][:, 1], acc['valid'][:, 1]), np.nan)
        out['m6a_ALL'] = ratio(acc['called'].sum(1), acc['valid'].sum(1))
    # Molecule rasters: random chromatinized (C-class) molecules for WT and each top hit.
    rng = np.random.default_rng(seed)
    names = {v['source_variant']: i for i, v in enumerate(variants)}
    raster, raster_names = {}, []
    seq = next(v['sequence'] for v in variants if v['class'] == 'WT')
    extra = [f'{p}{seq[p - 21]}>{b}' for p in COMPARE_POSITIONS for b in 'ACGT' if b != seq[p - 21]]
    for name in ['WT'] + list(dict.fromkeys(TOP_HITS + extra)):
        if name not in names:
            continue
        idx = np.flatnonzero(cohort & (vi == names[name]) & (groups >= 2))
        idx = np.sort(rng.choice(idx, min(raster_n, len(idx)), replace=False))
        cov = np.asarray(arr['coverage'][idx])
        state = np.zeros(cov.shape, np.int8)
        state[(cov & 2) > 0] = 1
        state[(cov & 8) > 0] = 3
        state[(cov & 4) > 0] = 2
        state[(cov & 1) > 0] = 4
        called = (np.asarray(arr['scores'][idx]) >= spec['ml_min']) & np.asarray(arr['known'][idx])
        order = np.argsort(-(state == 1).mean(1), kind='stable')
        raster[f'state_{len(raster_names)}'], raster[f'm6a_{len(raster_names)}'] = state[order], called[order]
        raster_names.append(name)
    with atomic_output(output) as o:
        np.savez_compressed(o / 'profiles.npz', **{k: v.astype(np.float32) for k, v in out.items()},
                            molecules=cnt.astype(np.int64), source_variant=np.array([v['source_variant'] for v in variants]))
        np.savez_compressed(o / 'rasters.npz', names=np.array(raster_names), **raster)
        write_json(o / 'audit.json', dict(schema_version='fsmpra.hit_profiles.v1', extract=str(extract), config=config,
                   strata_edges=edges.tolist(), strata_weights=weights.tolist(), groups='0 U, 1 M, 2.. C strata',
                   cohort=spec['cohort'] + ' and footprint evaluated', raster_per_allele=raster_n, seed=seed,
                   states=STATES, interpretation='Per-base fractions; C_std is standardized over WT distal-labeling strata. '
                   'Rasters are random C-class molecules ordered by element MSP fraction. Descriptive; one preparation.'))
    print(output)


# ---------------------------------------------------------------- figures
def load_effects():
    rows = {}
    with gzip.open(EFFECTS, 'rt') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            rows.setdefault(r['source_variant'], {})[r['endpoint']] = r
    return rows


def load_activity():
    act = {}
    for name, path, delim in (('LDLR', '20261020/published_GRCh38_LDLR.tsv', '\t'), ('LDLR.2', '20261020/published_GRCh38_LDLR.2.csv', ',')):
        with open(EVIDENCE / path) as h:
            for r in csv.DictReader(h, delimiter=delim):
                if r['Alt'] in 'ACGT' and len(r['Alt']) == 1:
                    act.setdefault(name, {})[(int(r['Position']), r['Ref'], r['Alt'])] = float(r['Value'])
    return act


def load_motifs(path):
    with open(path) as h:
        return list(csv.DictReader(h, delimiter='\t'))


def proxy_dominated(rows):
    """Local change explained by observability: the methylated proxy moves (|z| > 3) at least as much as chromatin."""
    m, c = rows.get('local:M'), rows.get('local:C_std')
    if not m or not c or not m['delta'] or not c['delta'] or not m['z']:
        return False
    return abs(float(m['z'])) > 3 and abs(float(m['delta'])) >= abs(float(c['delta']))


def strong_motifs(motifs, variant, n=2, best=1e-4, min_change=1.):
    """Largest changes among matches reaching p <= best in WT or variant, at least 10-fold."""
    rows = [r for r in motifs if r['variant'] == variant and min(float(r['wt_p']), float(r['alt_p'])) <= best
            and abs(float(r['log10_p_change'])) >= min_change]
    return sorted(rows, key=lambda r: -abs(float(r['log10_p_change'])))[:n]


def element_motif(motifs, variant):
    """Change in the named element motif (MOTIF_SITES) when the substitution falls inside that match, whatever its size."""
    for label, database, name, (a, b) in MOTIF_SITES:
        for r in motifs:
            if r['variant'] == variant and r['database'] == database and r['name'] == name and a <= int(r['oligo_start']) and int(r['oligo_end']) <= b:
                fold = 10 ** abs(float(r['log10_p_change']))
                return f"{label} {name}: p {float(r['wt_p']):.0e} → {float(r['alt_p']):.0e} ({fold:.{0 if fold >= 10 else 1}f}× {'weaker' if r['change'] == 'loss' else 'stronger'})"
    return None


def footprint_association(effects, min_molecules=300):
    """Does local chromatin change beyond broad accessibility/nucleosome change track activity?

    Partial correlation of the local 5-10 bp change with published activity after regressing both on the MSP and
    nucleosome-sized changes, for SNVs outside A runs (common set, >= min_molecules), with the circular position-shift null."""
    from fsmpra.library_stats import circular_shift
    out = []
    base = []
    for v, r in effects.items():
        if v == 'WT' or 'msp_C_std' not in r or int(r['msp_C_std']['molecules']) < min_molecules or r['msp_C_std']['context'] != 'other':
            continue
        vals = {ep: float(r[ep]['delta']) if ep in r and r[ep]['delta'] else np.nan
                for ep in ('local:A_norm', 'local:C_std', 'msp_C_std', 'nuc_85_250_C_std')}
        vals.update(variant=v, pos=int(r['msp_C_std']['oligo_pos1']), proxy_dominated=proxy_dominated(r),
                    activity=float(r['msp_C_std']['activity_LDLR']) if r['msp_C_std']['activity_LDLR'] else np.nan)
        base.append(vals)
    for subset in ('outside A runs', 'outside A runs, proxy-dominated removed'):
        rows = [b for b in base if not (subset.endswith('removed') and b['proxy_dominated'])]
        for ep in ('local:A_norm', 'local:C_std'):
            d = [b for b in rows if all(np.isfinite(b[k]) for k in (ep, 'msp_C_std', 'nuc_85_250_C_std', 'activity'))]
            X = np.c_[np.ones(len(d)), [b['msp_C_std'] for b in d], [b['nuc_85_250_C_std'] for b in d]]
            x, y = np.array([b[ep] for b in d]), np.array([b['activity'] for b in d])
            rx = x - X @ np.linalg.lstsq(X, x, rcond=None)[0]
            ry = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
            cs = circular_shift(np.array([b['pos'] for b in d]), rx, ry)
            out.append(dict(subset=subset, endpoint=ep, snvs=len(d), partial_r=round(float(np.corrcoef(rx, ry)[0, 1]), 4),
                            raw_r=round(float(np.corrcoef(x, y)[0, 1]), 4), positions=cs['positions'],
                            position_partial_r=round(cs['position_pearson'], 4), rotations=cs['rotations'],
                            null_abs_max=round(cs['null_abs_max'], 4), p_circular_shift=round(cs['p_circular_shift'], 4)))
    return out


def profile_changes(prof, effects, min_molecules=300):
    """Per-SNV changes from the per-base profiles (pp), with observability flags for footprint readouts.

    site = +/-5 bp, region = +/-30 bp, element = all 211 bases. tf_site_null_sd is the spread of the same window's change
    over SNVs > 40 bp away (sampling noise plus distal effects). m6A columns average unchanged A/T sites within 10 bp."""
    seq = wt_sequence()
    names = list(prof['source_variant'])
    wt = names.index('WT')
    rows = []
    for v, r in effects.items():
        if v == 'WT' or v not in names or 'msp_C_std' not in r or int(r['msp_C_std']['molecules']) < min_molecules:
            continue
        i = names.index(v)
        if not np.isfinite(prof['tf_C_std'][i]).all():
            continue
        pos, ref, alt = parse(v)
        d = {k: 100 * (prof[k][i] - prof[k][wt]) for k in ('tf_C_std', 'msp_C_std', 'nuc_C_std', 'm6a_C_std', 'm6a_M')}
        site, region = slice(max(pos - 26, 0), pos - 15), slice(max(pos - 51, 0), pos + 10)
        at = [k for k in range(pos - 31, pos - 10) if 0 <= k < 211 and seq[k] in 'AT' and k != pos - 21]
        with np.errstate(invalid='ignore'):
            dc, dm = float(np.nanmean(d['m6a_C_std'][at])), float(np.nanmean(d['m6a_M'][at]))
        rows.append(dict(variant=v, oligo_pos=pos, at_change='gain' if ref in 'CG' and alt in 'AT' else 'loss' if ref in 'AT' and alt in 'CG'
                         else 'none', context=r['msp_C_std']['context'], activity_log2=float(r['msp_C_std']['activity_LDLR']),
                         wt_tf_site_pct=100 * float(prof['tf_C_std'][wt][site].mean()), tf_site_pp=float(d['tf_C_std'][site].mean()),
                         msp_site_pp=float(d['msp_C_std'][site].mean()), msp_region_pp=float(d['msp_C_std'][region].mean()),
                         nuc_region_pp=float(d['nuc_C_std'][region].mean()), msp_element_pp=float(d['msp_C_std'].mean()),
                         nuc_element_pp=float(d['nuc_C_std'].mean()), m6a_chromatin_10bp_pp=dc, m6a_proxy_10bp_pp=dm,
                         proxy_clean=bool(abs(dm) < max(2., abs(dc)))))
    pos = np.array([r['oligo_pos'] for r in rows])
    D = np.array([100 * (prof['tf_C_std'][names.index(r['variant'])] - prof['tf_C_std'][wt]) for r in rows])
    for r in rows:
        w = D[:, max(r['oligo_pos'] - 26, 0):r['oligo_pos'] - 15].mean(1)
        sd = float(np.std(w[np.abs(pos - r['oligo_pos']) > 40], ddof=1))
        r.update(tf_site_null_sd=sd, tf_site_z=r['tf_site_pp'] / sd,
                 tf_resolved=r['at_change'] == 'none' and r['context'] == 'other' and r['proxy_clean'])
    return rows


def profile_association(rows):
    """Partial correlation of site/region changes with activity after regressing both on element MSP and nucleosome change."""
    from fsmpra.library_stats import circular_shift
    subsets = {'all': lambda r: True, 'A/T-neutral': lambda r: r['at_change'] == 'none',
               'A/T-neutral, outside A runs': lambda r: r['at_change'] == 'none' and r['context'] == 'other',
               'A/T-neutral, outside A runs, proxy-clean': lambda r: r['tf_resolved'],
               'A/T-neutral, WT TF-sized occupancy >= 4%': lambda r: r['at_change'] == 'none' and r['wt_tf_site_pct'] >= 4}
    out = []
    for name, keep in subsets.items():
        d = [r for r in rows if keep(r)]
        if len({r['oligo_pos'] for r in d}) < 30:                   # too few positions for the rotation null
            out += [dict(subset=name, readout=key, snvs=len(d), partial_r=np.nan, raw_r=np.nan, positions=len({r['oligo_pos'] for r in d}),
                         position_partial_r=np.nan, rotations=0, null_abs_max=np.nan, p_circular_shift=np.nan)
                    for key in ('tf_site_pp', 'msp_site_pp', 'msp_region_pp')]
            continue
        X = np.c_[np.ones(len(d)), [r['msp_element_pp'] for r in d], [r['nuc_element_pp'] for r in d]]
        y = np.array([r['activity_log2'] for r in d])
        ry = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
        for key in ('tf_site_pp', 'msp_site_pp', 'msp_region_pp'):
            x = np.array([r[key] for r in d])
            rx = x - X @ np.linalg.lstsq(X, x, rcond=None)[0]
            cs = circular_shift(np.array([r['oligo_pos'] for r in d]), rx, ry)
            out.append(dict(subset=name, readout=key, snvs=len(d), partial_r=round(float(np.corrcoef(rx, ry)[0, 1]), 4),
                            raw_r=round(float(np.corrcoef(x, y)[0, 1]), 4), positions=cs['positions'],
                            position_partial_r=round(cs['position_pearson'], 4), rotations=cs['rotations'],
                            null_abs_max=round(cs['null_abs_max'], 4), p_circular_shift=round(cs['p_circular_shift'], 4)))
    return out


def _style(plt):
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 8.5, 'axes.edgecolor': '#c3c2b7', 'axes.labelcolor': '#52514e',
                         'xtick.color': '#898781', 'ytick.color': '#898781', 'text.color': '#0b0b0b', 'axes.titlesize': 9.5,
                         'axes.titlelocation': 'left'})


INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#898781', '#e1e0d9'
WT_C, VAR_C = '#52514e', '#2a78d6'
STATE_COLORS = ['#f2f1ec', '#cfe0f5', '#52514e', '#9b9a94', '#eb6834']


def overview(effects, activity, motifs, output, plt, prows=None):
    """One row per top hit: activity, calibrated chromatin changes with 95% intervals, proxy check, TF-sized footprint change
    (when per-base profiles are given), top motif change."""
    hits = [h for h in TOP_HITS if h in effects]
    endpoints = [('msp_C_std', 'Accessibility'), ('nuc_85_250_C_std', 'Nucleosome-sized'), ('local:A_norm', 'Local, normalized'),
                 ('local:C_std', 'Local, chromatin'), ('local:M', 'Local, proxy')]
    n_tf = 1 if prows else 0
    fig = plt.figure(figsize=(17 + 1.6 * n_tf, 0.6 * len(hits) + 1.8))
    gs = fig.add_gridspec(1, 8 + n_tf, width_ratios=[2.6] + [1.05] * (6 + n_tf) + [2.3], wspace=.14)
    fig.subplots_adjust(top=.9)
    y = np.arange(len(hits))[::-1]
    ax0 = fig.add_subplot(gs[0])
    ax0.set_xlim(0, 1)
    ax0.set_ylim(-.6, len(hits) - .4)
    ax0.axis('off')
    for yi, h in zip(y, hits):
        pos, ref, alt = parse(h)
        r = effects[h]['msp_C_std']
        ax0.text(0, yi + .12, f'{h}  ·  chr19:{GENOMIC_OFFSET + pos:,}', fontsize=9, fontweight='bold', va='center')
        role = textwrap.fill(ROLES[h] + ('; A run' if r['context'] == 'homopolymer' else ''), 52)
        ax0.text(0, yi - .08, role, fontsize=7.2, color=INK2, va='top', linespacing=1.1)
    ax0.set_title('Variant', loc='left')
    # activity
    ax = fig.add_subplot(gs[1], sharey=ax0)
    for yi, h in zip(y, hits):
        pos, ref, alt = parse(h)
        for k, (src, mk) in enumerate((('LDLR', 'o'), ('LDLR.2', 'D'))):  # noqa: B007
            v = activity[src].get((GENOMIC_OFFSET + pos, ref, alt))
            if v is not None:
                ax.plot([0, v], [yi + (k - .5) * .18] * 2, color=MUTED, lw=1.2)
                ax.plot(v, yi + (k - .5) * .18, mk, color=INK2, ms=4.5)
    ax.axvline(0, color='#c3c2b7', lw=1)
    ax.set_title('Activity\n(log2)')
    _tidy(ax)
    for k, (ep, label) in enumerate(endpoints):
        ax = fig.add_subplot(gs[2 + k], sharey=ax0)
        for yi, h in zip(y, hits):
            r = effects[h].get(ep)
            if not r or not r['delta']:
                continue
            d, se = 100 * float(r['delta']), 100 * float(r['se']) if r['se'] else np.nan
            q = float(r['q_bh_primary_family']) if r.get('q_bh_primary_family') else np.nan
            unresolved = ep in ('local:A_norm', 'local:C_std') and (r['context'] == 'homopolymer' or proxy_dominated(effects[h]))
            if ep == 'local:M':
                color = INK2
            elif ep == 'local:C_std':
                color = VAR_C if abs(float(r['z'])) > 3 and not unresolved else MUTED
            else:
                color = VAR_C if q < .05 and not unresolved else MUTED
            ax.plot([d - 1.96 * se, d + 1.96 * se], [yi, yi], color=color, lw=1.6, solid_capstyle='round')
            ax.plot(d, yi, 'o', color=color if not unresolved else 'white', mec=color, ms=5.5)
        ax.axvline(0, color='#c3c2b7', lw=1)
        ax.set_title(f'{label}\n(Δ pp, 95%)')
        _tidy(ax)
    if prows:
        ax = fig.add_subplot(gs[7], sharey=ax0)
        for yi, h in zip(y, hits):
            r = prows.get(h)
            if not r:
                continue
            band = 1.96 * r['tf_site_null_sd']
            ax.plot([-band, band], [yi, yi], color=GRID, lw=6, solid_capstyle='butt', zorder=1)
            color = VAR_C if abs(r['tf_site_z']) > 3 and r['tf_resolved'] else MUTED
            ax.plot(r['tf_site_pp'], yi, 'o', color=color if r['tf_resolved'] else 'white', mec=color, ms=5.5, zorder=3)
        ax.axvline(0, color='#c3c2b7', lw=1)
        ax.set_title('TF-sized footprint\n±5 bp (Δ pp)')
        _tidy(ax)
    ax = fig.add_subplot(gs[7 + n_tf], sharey=ax0)
    ax.set_xlim(0, 1)
    ax.axis('off')
    ax.set_title('Largest motif changes\n(best p ≤ 1e-4, ≥10-fold)', loc='left')
    for yi, h in zip(y, hits):
        lines = []
        for r in strong_motifs(motifs, h):
            arrow = '↓' if r['change'] == 'loss' else '↑'
            lines.append(f"{arrow} {r['name']} ({r['database'][:4]})  p {float(r['wt_p']):.0e} → {float(r['alt_p']):.0e}")
        em = element_motif(motifs, h)
        if em and not any(em.split(':')[0].split()[-1] in line for line in lines):
            lines.append(em)
        ax.text(0, yi, '\n'.join(lines) if lines else 'no strong motif change', fontsize=7.3, color=INK2, va='center')
    tf_note = ('TF-sized footprint: per-base FiberHMM 10–60 bp calls in chromatinized molecules, ±5 bp of the substitution; grey band, 95% range '
               'of the same window over SNVs > 40 bp away; blue, |z| > 3; hollow, not interpretable (adds or removes an A/T, A run, or proxy m6A moves '
               'at least as much as chromatin within 10 bp). ' if prows else '')
    fig.text(.01, .01, tf_note + 'Activity: ● LDLR, ◆ LDLR.2 published MPRA. Blue: called at q < 0.05 in the prespecified primary family (mask-4 run); for local chromatin, |z| > 3. '
             'Grey: not called. Hollow: local readout unresolved (substitution in or next to an A run, or the methylated proxy changes at least as much as chromatin '
             'with |z| > 3). Dark grey: heavily methylated proxy class, an observability check, not a result. Motifs: best match over windows containing the base, '
             'JASPAR 2024 / HOCOMOCO v12, uniform background; sequence matches, not binding.', fontsize=7, color=MUTED, wrap=True)
    fig.suptitle('Top replicate-3 variant hits: function, chromatin and motif evidence', x=.01, ha='left', fontsize=12, y=.97)
    for ext in ('png', 'svg'):
        fig.savefig(Path(output) / f'top_hits_overview.{ext}', dpi=170, bbox_inches='tight')
    plt.close(fig)


def _tidy(ax):
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.tick_params(left=False, labelleft=False)
    ax.grid(axis='x', color=GRID, lw=.8)
    ax.set_axisbelow(True)


def variant_figure(h, prof, raster, effects, activity, motifs, output, plt, prow=None):
    """WT vs variant: per-base occupancy, per-site m6A in chromatin and proxy classes, molecules, function and motifs."""
    seq = wt_sequence()
    names = list(prof['source_variant'])
    wi, vi = names.index('WT'), names.index(h)
    pos, ref, alt = parse(h)
    p0 = pos - 21
    x = np.arange(21, 232)
    fig = plt.figure(figsize=(13, 12))
    gs = fig.add_gridspec(7, 2, height_ratios=[.9, 1, 1, 1, 1.1, 1.6, 1.6], width_ratios=[3, 1], hspace=.55, wspace=.12)
    r = effects[h]
    # header
    ax = fig.add_subplot(gs[0, :])
    ax.axis('off')
    stat = lambda ep: (100 * float(r[ep]['delta']), float(r[ep]['q_bh_primary_family']) if r[ep].get('q_bh_primary_family') else np.nan)
    a, n = stat('msp_C_std'), stat('nuc_85_250_C_std')
    l, lm = stat('local:A_norm'), 100 * float(r['local:M']['delta'])
    act = [activity[s].get((GENOMIC_OFFSET + pos, ref, alt)) for s in ('LDLR', 'LDLR.2')]
    ax.text(0, .95, f'{h}  ·  chr19:{GENOMIC_OFFSET + pos:,} {ref}>{alt}  ·  {ROLES[h]}', fontsize=13, fontweight='bold', va='top')
    ax.text(0, .52, f"Published activity log2 {act[0]:+.2f} / {act[1]:+.2f}   ·   accessibility Δ {a[0]:+.2f} pp (q {a[1]:.2g})   ·   "
            f"nucleosome-sized Δ {n[0]:+.2f} pp (q {n[1]:.2g})   ·   local 5–10 bp Δ {l[0]:+.2f} pp (q {l[1]:.2g}; proxy {lm:+.2f})"
            + ('   ·   A run: local readout unresolved' if r['msp_C_std']['context'] == 'homopolymer' else
               '   ·   proxy-dominated: local readout unresolved' if proxy_dominated(r) else ''),
            fontsize=8.5, color=INK2, va='top')
    ax.text(0, .1, f"Molecules {r['msp_C_std']['molecules']} (U {r['msp_C_std']['molecules_U']}, M {r['msp_C_std']['molecules_M']}, "
            f"C {r['msp_C_std']['molecules_C']}). Profiles: chromatinized molecules standardized over WT labeling strata."
            + (f"   TF-sized footprint ±5 bp Δ {prow['tf_site_pp']:+.2f} pp (WT {prow['wt_tf_site_pct']:.1f}%, z {prow['tf_site_z']:+.1f}"
               + (')' if prow['tf_resolved'] else '; not interpretable: ' + ('adds/removes an A/T' if prow['at_change'] != 'none' else
                  'A run' if prow['context'] != 'other' else 'proxy m6A moves as much as chromatin') + ')') if prow else ''),
            fontsize=7.5, color=MUTED, va='bottom')
    tracks = [('msp_C_std', 'Accessible (MSP)'), ('nuc_C_std', 'Nucleosome-sized\n85–250 bp'), ('tf_C_std', 'TF-sized footprint\n10–60 bp')]
    axes = []
    for k, (key, label) in enumerate(tracks):
        ax = fig.add_subplot(gs[1 + k, 0], sharex=axes[0] if axes else None)
        axes.append(ax)
        ax.plot(x, 100 * prof[key][wi], color=WT_C, lw=1.6, label='WT')
        ax.plot(x, 100 * prof[key][vi], color=VAR_C, lw=1.6, label=h)
        ax.fill_between(x, 100 * prof[key][wi], 100 * prof[key][vi], color=VAR_C, alpha=.12, lw=0)
        ax.set_ylabel(label + '\n(% molecules)', rotation=0, ha='right', va='center', fontsize=8)
        if k == 0:
            ax.legend(frameon=False, fontsize=8, ncol=2, loc='upper right')
    ax = fig.add_subplot(gs[4, 0], sharex=axes[0])
    axes.append(ax)
    at = np.array([b in 'AT' for b in seq]) & (np.arange(211) != p0)
    for key, mk, color, label in (('m6a_C_std', 'o', None, 'chromatinized'), ('m6a_M', 's', INK2, 'methylated proxy')):
        d = 100 * (prof[key][vi] - prof[key][wi])
        ax.vlines(x[at], 0, d[at], color=color or VAR_C, lw=1, alpha=.6)
        ax.plot(x[at], d[at], mk, ms=3.2, color=color or VAR_C, label=f'{label}', alpha=.9)
    ax.axhline(0, color='#c3c2b7', lw=1)
    ax.set_ylabel('m6A at unchanged\nA/T sites, Δ pp', rotation=0, ha='right', va='center', fontsize=8)
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc='lower right', bbox_to_anchor=(1, .98), handletextpad=.2)
    for ax in axes:
        ax.axvspan(pos - 4.5, pos + 4.5, color=MUTED, alpha=.12, lw=0)
        ax.axvline(pos, color=VAR_C, lw=1, ls=(0, (2, 2)))
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        ax.grid(axis='y', color=GRID, lw=.7)
        ax.set_axisbelow(True)
    for ax in axes[:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)
    axes[-1].set_xlim(20.5, 231.5)
    axes[-1].set_xlabel('Oligo position (chr19 position = 11,089,264 + oligo)', fontsize=8)
    # molecule rasters
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(STATE_COLORS)
    for k, name in enumerate(['WT', h]):
        ax = fig.add_subplot(gs[5 + k, 0], sharex=axes[0])
        if name in raster['names']:
            j = list(raster['names']).index(name)
            st = raster[f'state_{j}']
            ax.imshow(st, aspect='auto', cmap=cmap, vmin=0, vmax=4, interpolation='nearest', extent=(20.5, 231.5, len(st), 0))
            ax.axvline(pos, color=INK, lw=.8, ls=(0, (2, 2)))
        ax.set_ylabel(f'{name}\n{len(st) if name in raster["names"] else 0} molecules', rotation=0, ha='right', va='center', fontsize=8)
        ax.set_yticks([])
        plt.setp(ax.get_xticklabels(), visible=k == 1)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in STATE_COLORS]
    fig.legend(handles, STATES, loc='upper left', bbox_to_anchor=(.12, .085), frameon=False, fontsize=7.5, ncol=5,
               title='Molecule state (per base; random chromatinized molecules, most accessible first)', title_fontsize=7.5)
    # sequence and motifs (right column top)
    ax = fig.add_subplot(gs[1:3, 1])
    ax.axis('off')
    window = range(max(0, p0 - 12), min(211, p0 + 13))
    ax.text(0, .98, 'Sequence (oligo %d–%d)' % (window[0] + 21, window[-1] + 21), fontsize=8.5, fontweight='bold', va='top')
    for k, j in enumerate(window):
        ax.text(k / 26, .84, seq[j], fontsize=9, family='monospace', color=VAR_C if j == p0 else INK, fontweight='bold' if j == p0 else 'normal', ha='left')
        if j == p0:
            ax.text(k / 26, .74, alt, fontsize=9, family='monospace', color=VAR_C, fontweight='bold', ha='left')
    m = strong_motifs(motifs, h, n=4)
    ax.text(0, .62, 'Largest motif changes (p ≤ 1e-4, ≥10-fold)', fontsize=8, fontweight='bold', va='top')
    for k, row in enumerate(m):
        arrow = '↓ loss' if row['change'] == 'loss' else '↑ gain'
        ax.text(0, .52 - k * .13, f"{arrow}  {row['name']} ({row['database'][:4]}, {row['strand']})\n      p {float(row['wt_p']):.1e} → {float(row['alt_p']):.1e}",
                fontsize=7.3, color=INK2, va='top')
    if not m:
        ax.text(0, .52, 'No strong motif change', fontsize=7.5, color=INK2, va='top')
    em = element_motif(motifs, h)
    if em:
        ax.text(0, .52 - max(len(m), 1) * .13, 'Element motif: ' + em, fontsize=7.3, color=VAR_C, va='top', wrap=True)
    # function at neighbouring positions (right column)
    ax = fig.add_subplot(gs[3:5, 1])
    for src, mk, off in (('LDLR', 'o', -.12), ('LDLR.2', 'D', .12)):
        for j in range(pos - 10, pos + 11):
            for b in 'ACGT':
                if 21 <= j <= 231 and b != seq[j - 21]:
                    v = activity[src].get((GENOMIC_OFFSET + j, seq[j - 21], b))
                    if v is not None:
                        hl = j == pos and b == alt
                        ax.plot(j + off, v, mk, ms=5 if hl else 3, color=VAR_C if hl else MUTED, alpha=1 if hl else .7, zorder=3 if hl else 2)
    ax.axhline(0, color='#c3c2b7', lw=1)
    ax.set_title('Published activity, ±10 bp (log2)', fontsize=8.5)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    # calibrated statistics (right column bottom)
    ax = fig.add_subplot(gs[5:7, 1])
    eps = [('msp_C_std', 'Accessibility'), ('nuc_85_250_C_std', 'Nucleosome-sized'), ('flank:A_norm', 'Methylation 11–30 bp'),
           ('local:A_norm', 'Local 5–10 bp'), ('local:M', 'Local, methylated proxy')]
    for k, (ep, label) in enumerate(eps[::-1]):
        rr = r.get(ep)
        if not rr or not rr['delta']:
            continue
        d, se = 100 * float(rr['delta']), 100 * float(rr['se'])
        q = float(rr['q_bh_primary_family']) if rr.get('q_bh_primary_family') else np.nan
        unresolved = ep == 'local:A_norm' and (rr['context'] == 'homopolymer' or proxy_dominated(r))
        color = INK2 if ep == 'local:M' else VAR_C if q < .05 and not unresolved else MUTED
        ax.plot([d - 1.96 * se, d + 1.96 * se], [k, k], color=color, lw=2, solid_capstyle='round')
        ax.plot(d, k, 'o', color='white' if unresolved else color, mec=color, ms=6)
        ax.text(1.02, k, f'{label}', transform=ax.get_yaxis_transform(), fontsize=7.5, va='center', color=INK2)
    ax.axvline(0, color='#c3c2b7', lw=1)
    ax.set_yticks([])
    ax.set_title('Calibrated changes vs WT (pp, 95%)', fontsize=8.5)
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    fig.text(.01, .02, 'Replicate 3 only; site-methylation mask-4 run. Blue: q < 0.05 (primary family); hollow: local readout unresolved. Grey band: ±4 bp around the substitution (adjacent-adenine artifact reach). '
             'Methylated proxy: heavily methylated plasmids, a same-reaction observability check, not a prepared naked control.', fontsize=7, color=MUTED)
    safe = h.replace('>', '_to_')
    for ext in ('png', 'svg'):
        fig.savefig(Path(output) / f'{safe}.{ext}', dpi=150, bbox_inches='tight')
    plt.close(fig)


def position_figure(pos, prof, effects, activity, motifs, output, plt):
    """All three substitutions at one position: per-site m6A change (chromatin vs proxy), occupancy change, function, motifs."""
    seq = wt_sequence()
    names = list(prof['source_variant'])
    wi = names.index('WT')
    subs = [f'{pos}{seq[pos - 21]}>{b}' for b in 'ACGT' if b != seq[pos - 21]]
    subs = [v for v in subs if v in names and v in effects]
    lo, hi = max(21, pos - 30), min(231, pos + 30)
    x = np.arange(lo, hi + 1)
    j = x - 21
    at = np.array([seq[k] in 'AT' for k in j]) & (x != pos)
    fig, axes = plt.subplots(len(subs), 3, figsize=(13, 2.4 * len(subs) + 1.2), gridspec_kw=dict(width_ratios=[2.2, 2.2, 1.3], wspace=.25, hspace=.6),
                             squeeze=False)
    for row, v in enumerate(subs):
        vi = names.index(v)
        ax = axes[row, 0]
        for key, mk, color, label in (('m6a_C_std', 'o', VAR_C, 'chromatinized'), ('m6a_M', 's', INK2, 'methylated proxy')):
            d = 100 * (prof[key][vi][j] - prof[key][wi][j])
            ax.vlines(x[at], 0, d[at], color=color, lw=1, alpha=.55)
            ax.plot(x[at], d[at], mk, ms=3.2, color=color, label=label)
        ax.set_ylabel(v, rotation=0, ha='right', va='center', fontsize=10, fontweight='bold')
        ax = axes[row, 1]
        for key, color, label in (('msp_C_std', VAR_C, 'accessible'), ('nuc_C_std', INK2, 'nucleosome-sized'), ('tf_C_std', '#eb6834', 'TF-sized footprint')):
            ax.plot(x, 100 * (prof[key][vi][j] - prof[key][wi][j]), color=color, lw=1.6, label=label)
        ax = axes[row, 2]
        pp, ref, alt = parse(v)
        vals = [activity[s].get((GENOMIC_OFFSET + pp, ref, alt), np.nan) for s in ('LDLR', 'LDLR.2')]
        ax.barh([1, 0], vals, color=[INK2, MUTED], height=.55)
        ax.set_yticks([1, 0])
        ax.set_yticklabels(['LDLR', 'LDLR.2'], fontsize=7.5)
        ax.axvline(0, color='#c3c2b7', lw=1)
        ax.set_xlim(-3, 2)
        gains = sorted([r for r in motifs if r['variant'] == v and float(r['log10_p_change']) <= -1 and min(float(r['wt_p']), float(r['alt_p'])) <= 1e-4],
                       key=lambda r: float(r['log10_p_change']))
        losses = sorted([r for r in motifs if r['variant'] == v and float(r['log10_p_change']) >= 1 and min(float(r['wt_p']), float(r['alt_p'])) <= 1e-4],
                        key=lambda r: -float(r['log10_p_change']))
        txt = []
        if gains:
            txt.append('↑ ' + ', '.join(dict.fromkeys(r['name'] for r in gains[:4])))
        if losses:
            txt.append('↓ ' + ', '.join(dict.fromkeys(r['name'] for r in losses[:4])))
        r = effects[v]
        flag = ' (A run)' if r['msp_C_std']['context'] == 'homopolymer' else ' (proxy-dominated)' if proxy_dominated(r) else ''
        em = element_motif(motifs, v)
        ax.set_title(f"access. {100 * float(r['msp_C_std']['delta']):+.1f} pp · local {100 * float(r['local:A_norm']['delta']):+.1f} pp{flag}\n"
                     + (textwrap.shorten('; '.join(txt), 60) if txt else 'no strong motif change') + (f'\n{em}' if em else ''), fontsize=7.3, color=INK2)
        for a in axes[row, :2]:
            a.axvspan(pos - 4.5, pos + 4.5, color=MUTED, alpha=.12, lw=0)
            a.axvline(pos, color=VAR_C, lw=1, ls=(0, (2, 2)))
            a.axhline(0, color='#c3c2b7', lw=1)
        for a in axes[row]:
            for side in ('top', 'right'):
                a.spines[side].set_visible(False)
    axes[0, 0].set_title('m6A at unchanged A/T sites, variant − WT (pp)')
    axes[0, 1].set_title('Per-base occupancy, variant − WT (pp)')
    handles = axes[0, 0].get_legend_handles_labels()[0] + axes[0, 1].get_legend_handles_labels()[0]
    labels = ['m6A, chromatinized', 'm6A, methylated proxy', 'accessible (MSP)', 'nucleosome-sized 85–250 bp', 'TF-sized footprint 10–60 bp']
    fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(.1, .02), ncol=5, frameon=False, fontsize=7.5)
    axes[-1, 0].set_xlabel('Oligo position')
    axes[-1, 1].set_xlabel('Oligo position')
    axes[-1, 2].set_xlabel('Published activity (log2)')
    fig.suptitle(f'Oligo {pos} (chr19:{GENOMIC_OFFSET + pos:,}): every substitution compared with WT', x=.01, ha='left', fontsize=12)
    for ext in ('png', 'svg'):
        fig.savefig(Path(output) / f'position_{pos}.{ext}', dpi=150, bbox_inches='tight')
    plt.close(fig)


def motif_site_rows(effects, activity, prof=None):
    """For each named site: every substitution inside the best WT match, with motif change and chromatin/function effects.

    With per-base profiles, also the TF-sized footprint and MSP change over the whole match +/-5 bp (blank otherwise)."""
    from fsmpra.biology_motifs import BASES, read_matrices, pwm, score_distribution
    seq = wt_sequence()
    db = {'HOCOMOCOv12': read_matrices('extras/data/HOCOMOCOv12_H12CORE_pcm.txt', 'HOCOMOCOv12'),
          'JASPAR2024': read_matrices('extras/data/JASPAR2024_CORE_vertebrates_nr.jaspar', 'JASPAR2024')}
    out = []
    for label, database, name, (a, b) in MOTIF_SITES:
        best = None
        for m in [m for m in db[database] if m['name'] == name]:
            w = pwm(m['counts'], np.full(4, .25))
            L = len(w)
            q, off, tail = score_distribution(w, np.full(4, .25))

            def score(s, lo, hi):
                s = np.array([BASES.index(c) for c in s])
                top = (2., None)
                for st in range(lo, hi - L + 2):
                    for strand, f in (('+', s[st:st + L]), ('-', 3 - s[st:st + L][::-1])):
                        pv = float(tail[int(q[np.arange(L), f].sum()) - off])
                        if pv < top[0]:
                            top = (pv, st)
                return top
            pw, st = score(seq, a - 21, b - 21)
            if best is None or pw < best[0]:
                best = (pw, st, L, m['matrix_id'], score)
        pw, st, L, mid, score = best
        span = slice(max(st - 5, 0), st + L + 5)
        names = list(prof['source_variant']) if prof is not None else []
        for j in range(st, st + L):
            for alt in 'ACGT':
                if alt == seq[j]:
                    continue
                v = f'{j + 21}{seq[j]}>{alt}'
                pa = score(seq[:j] + alt + seq[j + 1:], max(0, st - L + 1), min(211, st + 2 * L - 1) - 1)[0]
                r = effects.get(v, {})
                get = lambda ep: 100 * float(r[ep]['delta']) if ep in r and r[ep]['delta'] and int(r[ep]['molecules']) >= 300 else np.nan
                deep = v in names and np.isfinite(get('msp_C_std'))
                site = lambda k: 100 * float(prof[k][names.index(v)][span].mean() - prof[k][names.index('WT')][span].mean()) if deep else np.nan
                out.append(dict(site=label, database=database, motif=name, matrix_id=mid, site_start=st + 21, site_end=st + L + 20, wt_p=pw,
                                variant=v, oligo_pos=j + 21, ref=seq[j], alt=alt, motif_log10_p_change=np.log10(pa) - np.log10(pw),
                                activity_log2=activity['LDLR'].get((GENOMIC_OFFSET + j + 21, seq[j], alt), np.nan),
                                local_chromatin_pp=get('local:C_std'), local_proxy_pp=get('local:M'), accessibility_pp=get('msp_C_std'),
                                nucleosome_pp=get('nuc_85_250_C_std'), proxy_dominated=bool(r) and proxy_dominated(r),
                                a_run=bool(r) and r['msp_C_std']['context'] == 'homopolymer',
                                at_change=(seq[j] in 'AT') != (alt in 'AT'), tf_site_pp=site('tf_C_std'), msp_site_pp=site('msp_C_std')))
    return out


def motif_site_figure(rows, output, plt):
    """Columns: motif sites. Rows: motif disruption, activity, local chromatin, accessibility, per substitution in the site."""
    from scipy import stats
    sites = list(dict.fromkeys(r['site'] for r in rows))
    tracks = [('motif_log10_p_change', 'Motif disruption\n(log10 p, alt − WT)'), ('activity_log2', 'Published activity\n(log2)'),
              ('local_chromatin_pp', 'Local chromatin\n5–10 bp (Δ pp)'), ('accessibility_pp', 'Accessibility\nMSP (Δ pp)')]
    profiled = any(np.isfinite(r.get('tf_site_pp', np.nan)) for r in rows)
    if profiled:
        tracks.insert(3, ('tf_site_pp', 'TF-sized footprint\nover match ±5 bp (Δ pp)'))
    fig, axes = plt.subplots(len(tracks), len(sites), figsize=(4.1 * len(sites), 2.4 * len(tracks) - .1), squeeze=False,
                             gridspec_kw=dict(hspace=.35, wspace=.28))
    seq = wt_sequence()
    for c, site in enumerate(sites):
        rs = [r for r in rows if r['site'] == site]
        x = np.array([r['oligo_pos'] for r in rs]) + np.array([{'A': -.24, 'C': -.08, 'G': .08, 'T': .24}[r['alt']] for r in rs])
        motif = np.array([r['motif_log10_p_change'] for r in rs])
        for k, (key, label) in enumerate(tracks):
            ax = axes[k, c]
            y = np.array([r[key] for r in rs], float)
            unresolved = np.array([r['a_run'] or r['proxy_dominated'] for r in rs]) & (key == 'local_chromatin_pp')
            if key == 'tf_site_pp':                                      # A/T changes alter footprint-call observability
                unresolved = np.array([r['a_run'] or r['at_change'] for r in rs])
            ok = np.isfinite(y)
            ax.axhline(0, color='#c3c2b7', lw=1)
            ax.scatter(x[ok & ~unresolved], y[ok & ~unresolved], s=22, c=np.clip(motif[ok & ~unresolved], -1, 3), cmap='Blues', vmin=-1, vmax=3,
                       edgecolor=INK2, lw=.4, zorder=3)
            ax.scatter(x[ok & unresolved], y[ok & unresolved], s=22, facecolor='none', edgecolor=MUTED, lw=.8, zorder=2)
            if key != 'motif_log10_p_change':
                use = ok & ~unresolved
                if use.sum() > 5:
                    rho, pv = stats.spearmanr(motif[use], y[use])
                    ax.text(.02, .97, f'ρ with motif disruption {rho:+.2f} (p {pv:.0e}, n {use.sum()})', transform=ax.transAxes,
                            fontsize=7.3, va='top', color=INK2)
                    lo_, hi_ = ax.get_ylim()
                    ax.set_ylim(lo_, hi_ + .25 * (hi_ - lo_))                # room for the label above the points
            if c == 0:
                ax.set_ylabel(label, fontsize=8.5)
            for side in ('top', 'right'):
                ax.spines[side].set_visible(False)
            ax.grid(axis='y', color=GRID, lw=.7)
            ax.set_axisbelow(True)
            lo, hi = rs[0]['site_start'], rs[0]['site_end']
            ax.set_xlim(lo - .6, hi + .6)
            ax.set_xticks(range(lo, hi + 1))
            ax.set_xticklabels([seq[p - 21] for p in range(lo, hi + 1)] if k == len(tracks) - 1 else [], fontsize=7.5, family='monospace')
        r0 = rs[0]
        mid = r0['matrix_id'].split('.')[0]
        axes[0, c].set_title(f"{site}: {r0['motif']} ({r0['database']}{'' if mid == r0['motif'] else ' ' + mid})\noligo {r0['site_start']}–{r0['site_end']}, "
                             f"WT match p {r0['wt_p']:.0e}", fontsize=9)
        axes[-1, c].set_xlabel(f"WT base (oligo {r0['site_start']}–{r0['site_end']}); dots: A, C, G, T left to right", fontsize=7.5)
    fig.suptitle('Motif disruption, footprint, accessibility and function at four promoter elements (replicate 3)', x=.01, ha='left', fontsize=12)
    fig.text(.01, .005, 'Each dot is one substitution inside the best WT match, shaded by how much it weakens the motif (darker = more). Hollow: readout unresolved '
             '(local chromatin: A run or proxy-dominated' + ('; TF-sized footprint: A run, or the substitution adds or removes an A/T' if profiled else '') +
             '). Spearman ρ excludes unresolved points.\nMotif matches are sequence models, not identified proteins; sites were chosen after seeing the data.'
             + (' TF-sized footprint: per-base FiberHMM 10–60 bp calls in chromatinized molecules, averaged over the match ±5 bp.' if profiled else ''),
             fontsize=7, color=MUTED)
    for ext in ('png', 'svg'):
        fig.savefig(Path(output) / f'motif_sites.{ext}', dpi=160, bbox_inches='tight')
    plt.close(fig)


def figures(profiles_dir, motifs_path, output, variants=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    _style(plt)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    from fsmpra.io import write_table
    effects, activity, motifs = load_effects(), load_activity(), load_motifs(motifs_path)
    prof = dict(np.load(Path(profiles_dir) / 'profiles.npz')) if profiles_dir else None
    prows = None
    if prof is not None:
        changes = profile_changes(prof, effects)
        write_table(output / 'profile_changes.tsv', changes, list(changes[0]))
        assoc = profile_association(changes)
        write_table(output / 'profile_association.tsv', assoc, list(assoc[0]))
        prows = {r['variant']: r for r in changes}
    overview(effects, activity, motifs, output, plt, prows)
    assoc = footprint_association(effects)
    write_table(output / 'footprint_association.tsv', assoc, list(assoc[0]))
    rows = motif_site_rows(effects, activity, prof)
    write_table(output / 'motif_sites.tsv', rows, list(rows[0]))
    motif_site_figure(rows, output, plt)
    if prof is not None:
        raster = dict(np.load(Path(profiles_dir) / 'rasters.npz'))
        for h in variants or TOP_HITS:
            if h in effects and h in list(prof['source_variant']):
                variant_figure(h, prof, raster, effects, activity, motifs, output, plt, prows.get(h))
        for pos in COMPARE_POSITIONS:
            position_figure(pos, prof, effects, activity, motifs, output, plt)
    print(output)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = p.add_subparsers(dest='cmd', required=True)
    m = sub.add_parser('motifs')
    m.add_argument('--output', required=True)
    pr = sub.add_parser('profiles')
    pr.add_argument('--extract', required=True)
    pr.add_argument('--output', required=True)
    f = sub.add_parser('figures')
    f.add_argument('--profiles', help='hit_profiles directory from the profiles step (omit for the overview only)')
    f.add_argument('--motifs', default=str(EVIDENCE / 'top_hits_20261021' / 'motif_changes.tsv'))
    f.add_argument('--output', required=True)
    a = p.parse_args()
    if a.cmd == 'motifs':
        from fsmpra.io import write_table
        seq = wt_sequence()
        rows = motif_changes([f'{p}{seq[p - 21]}>{b}' for p in HIT_POSITIONS for b in 'ACGT' if b != seq[p - 21]])
        Path(a.output).parent.mkdir(parents=True, exist_ok=True)
        write_table(a.output, rows, list(rows[0]))
        print(a.output, len(rows))
    elif a.cmd == 'profiles':
        profiles(a.extract, a.output)
    else:
        figures(a.profiles, a.motifs, a.output)


if __name__ == '__main__':
    main()

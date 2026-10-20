"""Independent library-level checks computed only from the committed evidence snapshot.

Run from the repository root (numpy/scipy only; no BAMs or results/ needed):

    python -m tools.independent_review_20261020 --output results/independent_review_20261020

The checks ask what the existing replicate-3 summaries already support at the
library level: split-half reproducibility of SNV effects, how much reproducible
local signal tracks A/T-opportunity change, and whether chromatin changes track
published activity under a spatial (circular position-shift) null. Barcode
halves come from one preparation, so reproducibility here includes any
deterministic sequence-dependent measurement bias.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy import stats
from fsmpra.io import atomic_output, open_text, source_info, write_json, write_table

GENOMIC_OFFSET = 11089264  # genomic_pos1 = offset + oligo_pos1 (forward strand, GRCh38)
LITERATURE_ELEMENTS = {  # LDLR promoter repeats (Südhof et al. 1987); oligo 1-based inclusive
    'Repeat1_Sp1': 'AAACTCCTCCTCTTGC', 'Repeat2_SRE1': 'AAAATCACCCCACTGC', 'Repeat3_Sp1': 'AAACTCCTCCCCCTGC'}
DATA_ELEMENTS = {'P1_core': (54, 62)}  # data-selected; activity-sensitive core of WT_P1
PEAKS = ['WT_P1', 'WT_P2', 'WT_P3', 'WT_P4', 'adjacent_143_150']
MIN_BARCODES = 50


def number(value):
    return float(value) if value not in ('', None) else np.nan


def finite(value):
    return float(value) if np.isfinite(value) else None


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    keep = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[keep], y[keep])[0, 1]) if keep.sum() > 2 else np.nan


def split_half(a, b):
    """Split-half correlation, Spearman-Brown reliability and signal/noise SDs (pp)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    r = float(np.corrcoef(a, b)[0, 1])
    cov = float(np.cov(a, b)[0, 1])
    return dict(n=len(a), pearson_halves=r, spearman_halves=float(stats.spearmanr(a, b)[0]),
                reliability_full=2 * r / (1 + r), signal_sd_pp=float(np.sqrt(max(cov, 0.))),
                noise_sd_full_pp=float(np.sqrt(np.var(a - b, ddof=1) / 4)))


def position_means(pos, values):
    pos, values = np.asarray(pos), np.asarray(values, float)
    keep = np.isfinite(values)
    upos = np.unique(pos[keep])
    return upos, np.array([values[keep & (pos == p)].mean() for p in upos])


def circular_shift(pos, chromatin, activity, margin=10):
    """Exact circular position-shift test of position-level correlation.

    Position means keep the three substitutions at a site together; rotating one
    profile against the other preserves each profile's spatial autocorrelation.
    Every rotation with offset in [margin, L-margin] is enumerated, so the
    smallest attainable p is 1/(number of rotations + 1)."""
    pos, chromatin, activity = np.asarray(pos), np.asarray(chromatin, float), np.asarray(activity, float)
    keep = np.isfinite(chromatin) & np.isfinite(activity)
    pc, xc = position_means(pos[keep], chromatin[keep])
    pa, ya = position_means(pos[keep], activity[keep])
    if not np.array_equal(pc, pa):
        raise ValueError('Position sets differ')
    observed = float(np.corrcoef(xc, ya)[0, 1])
    shifts = range(margin, len(xc) - margin + 1)
    null = np.array([np.corrcoef(np.roll(xc, s), ya)[0, 1] for s in shifts])
    p = (1 + np.sum(np.abs(null) >= abs(observed))) / (len(null) + 1)
    return dict(positions=len(xc), position_pearson=observed, rotations=len(null),
                null_abs_max=float(np.abs(null).max()), p_circular_shift=float(p))


def load(evidence):
    evidence = Path(evidence)
    wt_seq = None
    with open(evidence / 'validation_panel.tsv') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            if r['source_variant'] == 'WT':
                wt_seq = r['mutable_sequence']
    if wt_seq is None or len(wt_seq) != 211:
        raise ValueError('WT mutable sequence unavailable')
    activity = {}
    for name, path, delim in (('act1', 'published_GRCh38_LDLR.tsv', '\t'), ('act2', 'published_GRCh38_LDLR.2.csv', ',')):
        with open(evidence / path) as h:
            for r in csv.DictReader(h, delimiter=delim):
                if r['Alt'] in ('A', 'C', 'G', 'T'):
                    activity.setdefault(f"{r['Position']}:{r['Ref']}>{r['Alt']}", {})[name] = float(r['Value'])
    rows = []
    with open_text(evidence / 'barcode_sensitivity.tsv.gz') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            if r['source_variant'] == 'WT':
                continue
            name = r['source_variant']
            pos, ref, alt = int(name[:-3]), name[-3], name[-1]
            if wt_seq[pos - 21] != ref:
                raise ValueError(f'Reference mismatch for {name}')
            act = activity.get(f'{GENOMIC_OFFSET + pos}:{ref}>{alt}', {})
            rows.append(dict(variant_id=r['variant_id'], source_variant=name, pos=pos, ref=ref, alt=alt,
                             region=r['region_id'], in_region=r['mutation_in_region'] == 'True',
                             n_barcodes=int(r['n_barcodes']), at_change=r['at_change'],
                             delta=number(r['delta_pp']), split0=number(r['split_0_delta_pp']),
                             split1=number(r['split_1_delta_pp']),
                             masked2=number(r['shared_at_radius2_delta_pp']),
                             distal_m6a_std=number(r['distal_m6a_standardized_delta_pp']),
                             distal_nuc_std=number(r['distal_nuc_standardized_delta_pp']),
                             act1=act.get('act1', np.nan), act2=act.get('act2', np.nan)))
    landscape, wt = {}, {}
    with open(evidence / 'landscape_categories.tsv') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            (wt if r['source_variant'] == 'WT' else landscape.setdefault(r['variant_id'], {}))[r['track']] = number(r['mean_coverage'])
    landscape = {v: {t: 100 * (c - wt[t]) for t, c in tracks.items()} for v, tracks in landscape.items()}
    return wt_seq, rows, landscape


def run(evidence, output):
    evidence = Path(evidence)
    inputs = ['barcode_sensitivity.tsv.gz', 'landscape_categories.tsv', 'published_GRCh38_LDLR.tsv',
              'published_GRCh38_LDLR.2.csv', 'validation_panel.tsv', 'manifest.json']
    sources = [source_info(evidence / f, True) for f in inputs] + [source_info(Path(__file__), True)]
    wt_seq, rows, landscape = load(evidence)
    rows = [r for r in rows if r['n_barcodes'] >= MIN_BARCODES]
    by_region = {}
    for r in rows:
        by_region.setdefault(r['region'], []).append(r)

    reliability = []
    for region, group in by_region.items():
        for subset, keep in (('all', lambda r: True), ('mutation_in_region', lambda r: r['in_region']),
                             ('mutation_outside', lambda r: not r['in_region'])):
            g = [r for r in group if keep(r) and np.isfinite(r['split0']) and np.isfinite(r['split1'])]
            if len(g) >= 10 and not (region.startswith('whole') and subset != 'all'):
                reliability.append(dict(region=region, subset=subset,
                                        **split_half([r['split0'] for r in g], [r['split1'] for r in g])))

    at_class = []
    for region in PEAKS:
        g = [r for r in by_region[region] if r['in_region']]
        classes = sorted({r['at_change'] for r in g})
        mean0 = {c: np.mean([r['split0'] for r in g if r['at_change'] == c]) for c in classes}
        mean1 = {c: np.mean([r['split1'] for r in g if r['at_change'] == c]) for c in classes}
        signal = np.cov([r['split0'] for r in g], [r['split1'] for r in g])[0, 1]
        between = np.cov([mean0[r['at_change']] for r in g], [mean1[r['at_change']] for r in g])[0, 1]
        for c in classes:
            h = [r for r in g if r['at_change'] == c]
            at_class.append(dict(region=region, at_change=c, n=len(h),
                                 median_delta_pp=float(np.median([r['delta'] for r in h])),
                                 mean_masked_pm2_delta_pp=float(np.nanmean([r['masked2'] for r in h])),
                                 fraction_negative=float(np.mean([r['delta'] < 0 for r in h])),
                                 at_class_share_of_reproducible_variance=float(between / signal)))

    msp = {r['variant_id']: r for r in by_region['whole_msp']}
    variants = list(msp.values())
    for v in variants:
        for track in ('nuc_85_250', 'tf_small', 'msp'):
            v[f'land_{track}'] = landscape[v['variant_id']][track]
    endpoints = [('msp_delta_pp', 'delta'), ('msp_distal_m6a_standardized_pp', 'distal_m6a_std'),
                 ('msp_distal_protection_standardized_pp', 'distal_nuc_std'), ('msp_barcode_half_0_pp', 'split0'),
                 ('msp_barcode_half_1_pp', 'split1'), ('nuc_85_250_delta_pp', 'land_nuc_85_250'),
                 ('tf_small_element_delta_pp', 'land_tf_small')]
    subsets = [('all', lambda v: True), ('at_gain', lambda v: v['at_change'] == 'gain'),
               ('at_loss', lambda v: v['at_change'] == 'loss'),
               ('at_retained', lambda v: v['at_change'] == 'retained_count'),
               ('excluding_oligo_133_160', lambda v: v['pos'] < 133 or v['pos'] > 160),
               ('excluding_abs_activity_gt1', lambda v: abs(v['act1']) <= 1)]
    association = []
    for label, key in endpoints:
        for subset, keep in subsets:
            if subset != 'all' and key in ('split0', 'split1'):
                continue
            g = [v for v in variants if keep(v)]
            x = [v[key] for v in g]
            for act in ('act1', 'act2'):
                y = [v[act] for v in g]
                association.append(dict(endpoint=label, activity_source=act, subset=subset, variants=len(g),
                                        variant_pearson=pearson(x, y), **circular_shift([v['pos'] for v in g], x, y)))

    elements = {}
    for name, motif in LITERATURE_ELEMENTS.items():
        start = wt_seq.find(motif)
        if start < 0 or wt_seq.find(motif, start + 1) >= 0:
            raise ValueError(f'{name} must occur once in the WT element')
        elements[name] = (start + 21, start + 20 + len(motif))
    elements.update(DATA_ELEMENTS)
    pos = np.array([v['pos'] for v in variants])
    element_rows = []
    for name, (a, b) in elements.items():
        inside = (pos >= a) & (pos <= b)
        width = b - a + 1
        null = [np.nanmean([v['distal_m6a_std'] for v, p in zip(variants, pos) if s <= p < s + width])
                for s in range(21, 232 - width + 1) if s + width <= a or s > b]
        g = [v for v, k in zip(variants, inside) if k]
        mean = lambda key: float(np.nanmean([v[key] for v in g]))
        element_rows.append(dict(element=name, oligo_start1=a, oligo_end1=b, genomic_start1=GENOMIC_OFFSET + a,
                                 genomic_end1=GENOMIC_OFFSET + b, sequence=wt_seq[a - 21:b - 20], variants=len(g),
                                 activity_log2=mean('act1'), msp_delta_pp=mean('delta'), msp_half_0_pp=mean('split0'),
                                 msp_half_1_pp=mean('split1'), msp_distal_m6a_standardized_pp=mean('distal_m6a_std'),
                                 standardized_msp_rank_among_other_windows=float((np.sum(np.array(null) <= mean('distal_m6a_std')) + 1) / (len(null) + 1)),
                                 other_windows=len(null), nuc_85_250_delta_pp=mean('land_nuc_85_250'),
                                 selection='data-selected' if name in DATA_ELEMENTS else 'literature-defined'))

    upos = sorted(set(pos.tolist()))
    position_map = []
    local = {(r['region'], r['variant_id']): r['delta'] for r in rows if r['region'] in PEAKS and r['in_region']}
    for p in upos:
        g = [v for v in variants if v['pos'] == p]
        loc = [local[k] for k in local if k[1] in {v['variant_id'] for v in g}]
        position_map.append(dict(oligo_pos1=p, genomic_pos1=GENOMIC_OFFSET + p, ref=wt_seq[p - 21], variants=len(g),
                                 activity_log2=float(np.mean([v['act1'] for v in g])),
                                 msp_delta_pp=float(np.mean([v['delta'] for v in g])),
                                 msp_distal_m6a_standardized_pp=float(np.nanmean([v['distal_m6a_std'] for v in g])),
                                 nuc_85_250_delta_pp=float(np.mean([v['land_nuc_85_250'] for v in g])),
                                 local_peak_tf_delta_pp=finite(np.mean(loc)) if loc else None))

    def clean(rows):
        return [{k: (finite(v) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()} for r in rows]
    with atomic_output(output) as out:
        for name, table in (('split_half_reliability.tsv', reliability), ('at_class_local_effects.tsv', at_class),
                            ('activity_association.tsv', association), ('promoter_elements.tsv', element_rows),
                            ('position_map.tsv', position_map)):
            table = clean(table)
            write_table(out / name, table, list(table[0]))
        write_json(out / 'summary.json', dict(
            schema_version='fsmpra.independent_review.v1', sources=sources, min_barcodes=MIN_BARCODES,
            snvs=len(variants), activity_sources_pearson=pearson([v['act1'] for v in variants], [v['act2'] for v in variants]),
            interpretation='Within-preparation reproducibility and post hoc library-level associations; '
                           'not biological replication, sequence-matched calibration or causal mediation.'))
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--evidence', default='docs/evidence/20261020')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(run(args.evidence, args.output))


if __name__ == '__main__':
    main()

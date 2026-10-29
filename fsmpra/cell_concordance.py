"""Do SNVs change chromatin in the same direction in HepG2 and HEK293T?

Joins the per-cell effect tables (site_effects for HepG2, site_effects_clone for HEK293T) on
the HepG2 SNV ID and endpoint. Only endpoints that pass the calibration gate in both cells and
SNVs with enough molecules in both enter the per-SNV tests. The library-wide correlation, whose
permutation and rotation nulls do not depend on z calibration, is also reported for endpoints
that fail a gate, flagged `usable = False` (supplementary, decided after the HEK293T gates were
seen: 124 clones cannot bound an exactly calibrated rate below 0.10).

Library-wide (per endpoint): Pearson/Spearman of HEK293T vs HepG2 effects in z units, with a
label-permutation p and the exact circular position-shift p (fsmpra.library_stats), plus the
sign agreement among SNVs significant in HepG2 (binomial vs 0.5).

Per SNV: replicability by partial conjunction, p_rep = max(p_HepG2, p_HEK293T) when the signs
agree and 1 otherwise, BH across SNVs within each endpoint and across the primary family; a
difference z (effects scaled by each cell's own SE) is descriptive only, because cell type,
library construct, clone, reaction and m6A caller are confounded between the two cells.

Usage: python -m fsmpra.cell_concordance --config config/ldlr_cell_concordance.json
"""
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from .io import atomic_output, open_text, read_json, source_info, write_json, write_table
from .library_stats import circular_shift
from .site_effects import bh

BG, INK, MUTED, GRID = '#F7F6F1', '#1B2A41', '#4A5568', '#DCDDD5'
BLUE, ORANGE, GRAY = '#2B6CB0', '#C8641E', '#8A94A3'
CLASS_STYLE = {'replicated': (ORANGE, 'same direction, replicated'), 'HepG2 only': (BLUE, 'significant in HepG2 only'),
               'HEK293T only': (INK, 'significant in HEK293T only'), 'opposite': ('#A8521A', 'significant, opposite signs'),
               'neither': (GRAY, 'neither')}


def num(x):
    return float(x) if x not in ('', None) else np.nan


def load(path):
    rows = {}
    with open_text(Path(path) / 'effects.tsv.gz') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            rows[(r['source_variant'], r['endpoint'])] = r
    gate = {}
    with open(Path(path) / 'calibration_summary.tsv') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            gate[r['endpoint']] = r['gate_pass'] == 'True'
    return rows, gate


def classify(q1, q2, d1, d2, alpha=.05):
    s1, s2 = q1 < alpha, q2 < alpha
    if s1 and s2:
        return 'replicated' if np.sign(d1) == np.sign(d2) else 'opposite'
    return 'HepG2 only' if s1 else 'HEK293T only' if s2 else 'neither'


def run(config):
    cfg = read_json(config)
    cc, spec = cfg['concordance'], cfg['analysis']
    primary = spec['primary_endpoints']
    paths = {c: cfg['cells'][c]['effects'] for c in ('HepG2', 'HEK293T')}
    (hep, gate_hep), (hek, gate_hek) = load(paths['HepG2']), load(paths['HEK293T'])
    endpoints = sorted({ep for _, ep in hek} & {ep for _, ep in hep})
    candidates = [ep for ep in endpoints if not ep.endswith((':U', ':M'))]
    usable = [ep for ep in candidates if gate_hep.get(ep) and gate_hek.get(ep)]
    rng = np.random.default_rng(cc['seed'])
    rows = []
    for (sv, ep), b in hek.items():
        a = hep.get((sv, ep))
        if a is None or ep not in candidates:
            continue
        if num(a['molecules']) < cc['min_molecules'] or num(b['molecules']) < cc['min_molecules']:
            continue
        z1, z2, p1, p2 = num(a['z']), num(b['z']), num(a['p']), num(b['p'])
        if not all(np.isfinite([z1, z2, p1, p2])):
            continue
        rows.append(dict(source_variant=sv, endpoint=ep, oligo_pos1=int(a['oligo_pos1']), genomic_pos1=int(a['genomic_pos1']),
                         at_change=a['at_change'], context=a['context'],
                         hepg2_molecules=int(num(a['molecules'])), hek293t_molecules=int(num(b['molecules'])),
                         hek293t_top_clone_fraction=num(b.get('top_clone_fraction')),
                         hepg2_delta=num(a['delta']), hepg2_se=num(a['se']), hepg2_z=z1, hepg2_p=p1,
                         hepg2_q_primary=num(a.get('q_bh_primary_family')),
                         hek293t_delta=num(b['delta']), hek293t_se=num(b['se']), hek293t_z=z2, hek293t_p=p2,
                         hek293t_q_primary=num(b.get('q_bh_primary_family')),
                         p_rep=max(p1, p2) if np.sign(z1) == np.sign(z2) else 1.,
                         difference_z=(z2 - z1) / np.sqrt(2),
                         activity_LDLR=num(a.get('activity_LDLR')), activity_LDLR_2=num(a.get('activity_LDLR_2')),
                         usable=ep in usable))
    for ep in usable:
        sel = [r for r in rows if r['endpoint'] == ep]
        for r, q in zip(sel, bh([r['p_rep'] for r in sel])):
            r['q_rep_within_endpoint'] = float(q)
    prim = [r for r in rows if r['endpoint'] in primary and r['usable']]
    for r, q in zip(prim, bh([r['p_rep'] for r in prim])):
        r['q_rep_primary_family'] = float(q)
    for r in rows:
        r['class'] = (classify(r['hepg2_q_primary'], r['hek293t_q_primary'], r['hepg2_z'], r['hek293t_z'])
                      if r['endpoint'] in primary and r['usable'] else '')
    library = []
    for ep in sorted({r['endpoint'] for r in rows}):
        sel = [r for r in rows if r['endpoint'] == ep]
        x, y = np.array([r['hepg2_z'] for r in sel]), np.array([r['hek293t_z'] for r in sel])
        if len(sel) < 20:
            continue
        r_obs = float(np.corrcoef(x, y)[0, 1])
        perm = np.array([np.corrcoef(x, rng.permutation(y))[0, 1] for _ in range(cc['permutations'])])
        rot = circular_shift([r['oligo_pos1'] for r in sel], x, y)
        sig = [r for r in sel if r['hepg2_q_primary'] < .05] if ep in primary else [r for r in sel if r['hepg2_p'] < .05]
        agree = sum(np.sign(r['hepg2_z']) == np.sign(r['hek293t_z']) for r in sig)
        library.append(dict(endpoint=ep, primary=ep in primary, usable=ep in usable, gate_hepg2=bool(gate_hep.get(ep)),
                            gate_hek293t=bool(gate_hek.get(ep)), snvs=len(sel), pearson_z=r_obs,
                            spearman_z=float(stats.spearmanr(x, y).correlation),
                            pearson_delta=float(np.corrcoef([r['hepg2_delta'] for r in sel], [r['hek293t_delta'] for r in sel])[0, 1]),
                            p_permutation=float((1 + np.sum(np.abs(perm) >= abs(r_obs))) / (len(perm) + 1)),
                            position_pearson=rot['position_pearson'], p_circular_shift=rot['p_circular_shift'],
                            rotations=rot['rotations'], hepg2_significant=len(sig), sign_agree=int(agree),
                            sign_agree_fraction=agree / len(sig) if sig else np.nan,
                            p_sign_binomial=float(stats.binomtest(int(agree), len(sig), .5).pvalue) if sig else np.nan,
                            replicated_q05=sum(r.get('q_rep_primary_family' if ep in primary else 'q_rep_within_endpoint', 1) < .05
                                               for r in sel) if ep in usable else None))
    out_dir = Path(cfg['output_dir'])
    fields = list(dict.fromkeys(k for r in rows for k in r))

    def clean(rs):
        return [{k: (None if isinstance(v, (float, np.floating)) and not np.isfinite(v) else v) for k, v in r.items()} for r in rs]
    with atomic_output(out_dir) as out:
        write_table(out / 'concordance.tsv', clean(rows), fields)
        write_table(out / 'library_concordance.tsv', clean(library), list(library[0]))
        write_json(out / 'summary.json', dict(
            schema_version='fsmpra.cell_concordance.v1',
            sources=[source_info(Path(config), True)] + [source_info(Path(p) / f, True) for p in paths.values()
                                                         for f in ('effects.tsv.gz', 'calibration_summary.tsv')],
            endpoints_both_cells=len(endpoints), endpoints_usable=usable,
            endpoints_failing_gate={'HepG2': sorted(e for e in endpoints if not gate_hep.get(e)),
                                    'HEK293T': sorted(e for e in endpoints if not gate_hek.get(e))},
            snvs=len({r['source_variant'] for r in rows}),
            primary_classes={ep: {k: sum(r['class'] == k for r in rows if r['endpoint'] == ep) for k in CLASS_STYLE}
                             for ep in primary},
            interpretation='Concordance (same direction in both cells) is robust to the confounds between cells; a '
                           'cell-specific or opposite call is not by itself evidence of cell-type biology.'))
        figures(out, rows, library, primary)
    print(f'{len({r["source_variant"] for r in rows})} SNVs, {len(usable)} usable endpoints', flush=True)


def figures(out, rows, library, primary):
    plt.rcParams.update({'figure.facecolor': BG, 'axes.facecolor': BG, 'savefig.facecolor': BG, 'font.size': 13,
                         'axes.titlesize': 14, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left',
                         'axes.labelcolor': MUTED, 'axes.edgecolor': GRID, 'axes.spines.top': False,
                         'axes.spines.right': False, 'xtick.color': MUTED, 'ytick.color': MUTED, 'text.color': INK,
                         'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': .7, 'axes.axisbelow': True})
    eps = [e for e in primary if any(r['endpoint'] == e for r in rows)]
    if not eps:
        return
    lib = {d['endpoint']: d for d in library}
    fig, axes = plt.subplots(1, len(eps), figsize=(5.2 * len(eps), 5.2), squeeze=False)
    for ax, ep in zip(axes[0], eps):
        sel = [r for r in rows if r['endpoint'] == ep]
        if not sel[0]['usable']:
            ax.scatter([r['hepg2_z'] for r in sel], [r['hek293t_z'] for r in sel], s=22, color=GRAY, edgecolor=BG, lw=.8,
                       label=f'per-SNV tests withheld: gate failed ({len(sel)})')
        for cls, (color, label) in CLASS_STYLE.items():
            pts = [r for r in sel if r['class'] == cls]
            if pts:
                ax.scatter([r['hepg2_z'] for r in pts], [r['hek293t_z'] for r in pts], s=34 if cls != 'neither' else 18,
                           color=color, edgecolor=BG, lw=.8, label=f'{label} ({len(pts)})', zorder=3 if cls != 'neither' else 2)
        ax.axhline(0, color=GRAY, lw=.8)
        ax.axvline(0, color=GRAY, lw=.8)
        d = lib.get(ep, {})
        ax.set(title=f"{ep}\nr = {d.get('pearson_z', np.nan):.2f}, rotation p = {d.get('p_circular_shift', np.nan):.3f}",
               xlabel='HepG2 z', ylabel='HEK293T z')
        ax.legend(fontsize=9, frameon=False, loc='best')
    fig.tight_layout()
    fig.savefig(out / 'concordance_scatter.png', dpi=200)
    plt.close(fig)
    fig, axes = plt.subplots(len(eps), 1, figsize=(14, 2.6 * len(eps)), sharex=True, squeeze=False)
    for ax, ep in zip(axes[:, 0], eps):
        sel = sorted((r for r in rows if r['endpoint'] == ep), key=lambda r: r['oligo_pos1'])
        x = np.array([r['oligo_pos1'] for r in sel])
        ax.scatter(x - .15, [r['hepg2_z'] for r in sel], s=14, color=BLUE, label='HepG2')
        ax.scatter(x + .15, [r['hek293t_z'] for r in sel], s=14, color=ORANGE, label='HEK293T')
        ax.axhline(0, color=GRAY, lw=.8)
        ax.set(ylabel='z', title=ep)
    axes[0, 0].legend(frameon=False, ncol=2, loc='upper right')
    axes[-1, 0].set_xlabel('oligo position')
    fig.tight_layout()
    fig.savefig(out / 'concordance_track.png', dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', required=True)
    run(ap.parse_args().config)


if __name__ == '__main__':
    main()

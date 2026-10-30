"""Step A4 (part 3): predicted change in visible adenines per variant, from naked-DNA sequence bias.

A variant changes evidence two ways: it adds or removes an A/T, and it changes the A-run context of nearby
adenines (breaking AAAAA into AAGAA makes the remaining adenines far more visible on naked DNA). Visibility of
an adenine = its naked-DNA methylation rate relative to isolated adenines, looked up by (run length bin, position
in run) from the naked Hia5-methylated pEF-GFP control (a4_context_bias.load). Per variant:
  visible_wt, visible_variant = sum of visibility over all adenines (both strands) of the element
  delta_visible = variant - WT (in "isolated-adenine equivalents")
  neighbours_changed = adenines other than the edited base whose (run, position) class changed
This is a prediction from a different plasmid and caller version; step A7 measures the effect on our data and
B1 (naked controls of this library) replaces it.

Usage: python tools/analysis_v2/a4_variant_visibility.py --catalog <catalog.json> --naked <emissions.tsv.gz> \
           --variants <a4 variants.tsv> --output <new dir>
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('a4c', Path(__file__).with_name('a4_context_bias.py'))
a4c = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a4c)


def adenine_classes(seq, max_bin=5):
    """{(strand, index): (run_bin, pos_in_run)} for every adenine of both strands."""
    out = {}
    n = len(seq)
    i = 0
    while i < n:
        b = seq[i]
        j = i
        while j < n and seq[j] == b:
            j += 1
        if b in 'AT':
            L = j - i
            for k in range(i, j):
                rank = (k - i + 1) if b == 'A' else (j - k)
                pos = 'isolated' if L == 1 else "5' end" if rank == 1 else "3' end" if rank == L else 'interior'
                out[('top' if b == 'A' else 'bottom', k)] = (min(L, max_bin), pos)
        i = j
    return out


def lookup_table(naked_path):
    t = a4c.load(naked_path)
    return t.groupby(['run_bin', 'pos_in_run']).p_rel.median().to_dict()


def visibility(seq, lut):
    return {key: lut[c] for key, c in adenine_classes(seq).items()}


def compare(wt, var_seq, cls, i, lut):
    """Sum of visibility in WT and variant; adenines (other than the edited base) whose class changed.
    Variant coordinates are mapped back to WT: a DEL at i shifts later bases left, an INS at i right."""
    vw, vv = visibility(wt, lut), visibility(var_seq, lut)
    cw, cv = adenine_classes(wt), adenine_classes(var_seq)

    def to_wt(k):
        if cls == 'DEL':
            return k + 1 if k >= i else k
        if cls == 'INS':
            return None if k == i else (k - 1 if k > i else k)
        return k
    changed = 0
    for (s, k), c in cv.items():
        kw = to_wt(k)
        if kw is None or (cls == 'SNV' and kw == i):
            continue
        if cw.get((s, kw)) != c:
            changed += 1
    return sum(vw.values()), sum(vv.values()), changed


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--naked', required=True)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    lib = json.load(open(a.catalog))['library']
    wt = lib['mutable_sequence']
    seqs = {v['variant_id']: v['sequence'] for v in lib['variants']}
    lut = lookup_table(a.naked)
    var = pd.read_csv(a.variants, sep='\t')
    r = [compare(wt, seqs[v.variant_id], v.cls, v.index0, lut) for v in var.itertuples()]
    var['visible_wt'] = [x[0] for x in r]
    var['visible_variant'] = [x[1] for x in r]
    var['delta_visible'] = var.visible_variant - var.visible_wt
    var['neighbours_changed'] = [x[2] for x in r]
    var.to_csv(out / 'variants_visibility.tsv', sep='\t', index=False)
    pd.Series({f'{k[0]}|{k[1]}': v for k, v in lut.items()}).rename('p_rel').to_csv(out / 'naked_lookup.tsv', sep='\t')
    g = var.groupby(['cls', 'seq_class']).agg(
        n=('delta_visible', 'size'), delta_visible_median=('delta_visible', 'median'),
        delta_visible_min=('delta_visible', 'min'), delta_visible_max=('delta_visible', 'max'),
        frac_neighbours_changed=('neighbours_changed', lambda x: float((x > 0).mean())),
        frac_gain_ge_0p5=('delta_visible', lambda x: float((x >= 0.5).mean())),
        frac_loss_le_m0p5=('delta_visible', lambda x: float((x <= -0.5).mean())))
    g.reset_index().to_csv(out / 'by_class.tsv', sep='\t', index=False)
    print(g.round(3).to_string())
    print('visible_wt', round(var.visible_wt.iat[0], 2), 'adenines', len(adenine_classes(wt)))
    top = var.reindex(var.delta_visible.abs().sort_values(ascending=False).index).head(15)
    print(top[['source_variant', 'cls', 'seq_class', 'delta_visible', 'neighbours_changed', 'region']].to_string(index=False))


if __name__ == '__main__':
    main()

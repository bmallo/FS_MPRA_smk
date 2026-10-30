"""Step A4 (part 2): is the per-adenine methylation rate set by local sequence (A-runs)?

Uses plasmid-LLR emission tables (fold A; p_a = methylation rate of each reference A/T in the "accessible"
class). Each table's 11-mer context gives, for the central adenine on its own strand (A for top, T for bottom),
the length of the same-base run it sits in (max 11) and its rank from the run's 5' end on that strand.

Comparing a truly naked Hia5-methylated plasmid (pEF-GFP benchmark control) with the in-sample M classes of
the LDLR runs separates sequence bias (present on naked DNA) from protection (absent on naked DNA). Different
plasmids and m6A-caller versions, so only the shape by context is compared, as p_a relative to each table's
median over isolated adenines (run length 1).

Usage: python tools/analysis_v2/a4_context_bias.py --table naked_pEFGFP=<emissions.tsv.gz> \
           --table HepG2_M=<emissions.tsv.gz> [--element HepG2_M=1234:1445] --output <new dir>
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

RUN_BINS = [1, 2, 3, 4, 5]


def run_context(context, strand):
    """(run length, rank from the run's 5' end on the adenine's strand) for the centre base of an odd k-mer."""
    c = len(context) // 2
    base = 'A' if strand == 'top' else 'T'
    if context[c] != base:
        raise ValueError(f'centre of {context} is not {base}')
    lo = c
    while lo > 0 and context[lo - 1] == base:
        lo -= 1
    hi = c
    while hi < len(context) - 1 and context[hi + 1] == base:
        hi += 1
    n = hi - lo + 1
    rank = (c - lo + 1) if strand == 'top' else (hi - c + 1)   # bottom strand reads right to left
    return n, rank


def load(path, element=None):
    t = pd.read_csv(path, sep='\t')
    t = t[(t.fold == 'A') & ~t.low_cov.astype(str).str.lower().eq('true')].copy()
    t = t[t.context.str.len() == t.context.str.len().max()]   # drop truncated contexts at the reference ends
    if element:
        a, b = element
        t = t[(t.ref_pos > a) & (t.ref_pos <= b)]
    rc = [run_context(c, s) for c, s in zip(t.context, t.strand)]
    t['run'] = [r[0] for r in rc]
    t['rank'] = [r[1] for r in rc]
    t['run_bin'] = np.minimum(t.run, RUN_BINS[-1])
    t['pos_in_run'] = np.where(t.run == 1, 'isolated', np.where(t['rank'] == 1, "5' end",
                               np.where(t['rank'] == t.run, "3' end", 'interior')))
    t['p_rel'] = t.p_a / t.loc[t.run == 1, 'p_a'].median()
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--table', action='append', required=True, help='LABEL=emissions.tsv.gz')
    ap.add_argument('--element', action='append', default=[], help='LABEL=start0:end0 (restrict to the element)')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    el = {k: tuple(map(int, v.split(':'))) for k, v in (e.split('=', 1) for e in a.element)}
    rows, pos_rows, summ = [], [], {}
    for spec in a.table:
        label, path = spec.split('=', 1)
        for scope in (['plasmid', 'element'] if label in el else ['plasmid']):
            t = load(path, el[label] if scope == 'element' else None)
            g = t.groupby('run_bin').agg(n=('p_a', 'size'), p_a_median=('p_a', 'median'), p_rel_median=('p_rel', 'median'),
                                         frac_below_0p2=('p_a', lambda x: float((x < 0.2).mean())))
            rows.append(g.reset_index().assign(table=label, scope=scope))
            h = t[t.run >= 3].groupby('pos_in_run').agg(n=('p_a', 'size'), p_rel_median=('p_rel', 'median'))
            pos_rows.append(h.reset_index().assign(table=label, scope=scope))
            summ[f'{label}:{scope}'] = dict(source=path, n=len(t), median_p_a_isolated=float(t.loc[t.run == 1, 'p_a'].median()),
                                            spearman_run_vs_p_rel=float(t[['run', 'p_rel']].corr('spearman').iat[0, 1]))
    by_run = pd.concat(rows)[['table', 'scope', 'run_bin', 'n', 'p_a_median', 'p_rel_median', 'frac_below_0p2']]
    by_pos = pd.concat(pos_rows)[['table', 'scope', 'pos_in_run', 'n', 'p_rel_median']]
    by_run.to_csv(out / 'p_a_by_run_length.tsv', sep='\t', index=False)
    by_pos.to_csv(out / 'p_a_by_position_in_run.tsv', sep='\t', index=False)
    json.dump(summ, open(out / 'summary.json', 'w'), indent=1)
    print(by_run.to_string(index=False))
    print(by_pos.to_string(index=False))
    print(json.dumps(summ, indent=1))


if __name__ == '__main__':
    main()

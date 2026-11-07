"""AlphaGenome-only proxy benign counts in a 350 bp window per element (rules.json "v2_alphagenome_350bp").

Uses snv_proxy.tsv from proxy.py (AlphaGenome classes with the 1 kb element's top-20% cutoffs; ClinVar >= 1-star
benchmark flags) and <root>/v2_350bp/windows_350bp.tsv (MANE Select TSS from UCSC). No AlphaGenome rerun.
Writes <root>/v2_350bp/windows_final.tsv, snv_350bp.tsv, counts.tsv
Usage: python tools/proxy_benign/proxy_350.py --root results/proxy_benign_20261106
"""
import argparse
from pathlib import Path

import pandas as pd

L, UP = 350, 175
ANCHOR = {'APC': ('promoter 1B TSS NM_001127511.3', 112707498, '+')}   # no MANE Select TSS in the element


def place(r):
    if r.gene == 'NIPBL':                                    # NIPBL TSS just past the 3' edge: TSS-side 350 bp
        return r.end - L + 1, r.end
    lo = r.tss - UP if r.strand == '+' else r.tss - (L - UP - 1)
    lo = min(max(lo, r.start), r.end - L + 1)                # slide inward to stay inside the element
    return lo, lo + L - 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    a = ap.parse_args()
    root = Path(a.root)
    d = root / 'v2_350bp'
    t = pd.read_csv(d / 'windows_350bp.tsv', sep='\t')
    e = pd.read_csv(Path(__file__).parent / 'elements.tsv', sep='\t')
    t = e[['gene', 'chrom', 'start', 'end']].merge(t.drop(columns=['chrom', 'start', 'end', 'w_start', 'w_end',
                                                                   'inside_element', 'tss_offset_in_element']),
                                                   on='gene', how='left')
    for g, (lab, tss, st) in ANCHOR.items():
        t.loc[t.gene == g, ['anchor', 'tss', 'strand']] = [lab, tss, st]
    t.loc[t.gene == 'NIPBL', 'anchor'] = "element 3' edge (NIPBL TSS just outside)"
    t[['w_start', 'w_end']] = [list(place(r)) for r in t.itertuples()]
    t['shift_bp'] = [(r.w_start - (r.tss - UP if r.strand == '+' else r.tss - (L - UP - 1))) if r.gene != 'NIPBL'
                     else None for r in t.itertuples()]
    t['tss_in_window'] = (t.tss >= t.w_start) & (t.tss <= t.w_end)
    t.to_csv(d / 'windows_final.tsv', sep='\t', index=False)
    s = pd.read_csv(root / 'snv_proxy.tsv', sep='\t', low_memory=False)
    s['clinvar_benchmark'] = s.clinvar_benchmark.fillna('')
    W = t.set_index('gene')
    s = s[[W.loc[g, 'w_start'] <= p <= W.loc[g, 'w_end'] for g, p in zip(s.window, s.pos)]]
    s.to_csv(d / 'snv_350bp.tsv', sep='\t', index=False, float_format='%.6g')
    rows = []
    for g, x in s.groupby('window', sort=False):
        P = x[x.clinvar_benchmark == '']
        rows.append(dict(element=g, possible_snvs=len(x), benchmark_snvs=int((x.clinvar_benchmark != '').sum()),
                         proxy_benign=int(P.ag_no_change.sum()),
                         proxy_benign_pct=round(100 * P.ag_no_change.sum() / len(x), 1),
                         ag_rna_large=int(P.ag_rna_large.sum()), ag_dnase_large=int(P.ag_dnase_large.sum())))
    c = pd.DataFrame(rows)
    c.to_csv(d / 'counts.tsv', sep='\t', index=False)
    print(t[['gene', 'anchor', 'strand', 'tss', 'start', 'end', 'w_start', 'w_end', 'shift_bp', 'tss_in_window']]
          .to_string(index=False))
    print(c.to_string(index=False))


if __name__ == '__main__':
    main()

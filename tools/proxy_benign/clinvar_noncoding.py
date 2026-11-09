"""ClinVar records in the proxy benign elements, flagged as protein-changing or not and labelled by MANE Select region
(rules.json "v3_noncoding_clinvar").

Inputs: <root>/clinvar/<element>.json (fetch.py clinvar esummary), elements.tsv. MANE Select transcripts come from the
UCSC hg38 mane track (cached in <output>/mane/). Writes <output>/clinvar_annotated.tsv (one row per record in a window)
and <output>/clinvar_counts.tsv (per element, >= 1-star P/LP and B/LB before and after the protein-change filter).
Usage: python tools/proxy_benign/clinvar_noncoding.py --root results/proxy_benign_20261106 \
           --output results/proxy_benign_20261106/v3_noncoding_clinvar
"""
import argparse
import json
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / 'promoter_variant_census'))
import census_lib as cl  # noqa: E402
import truth_sets as ts  # noqa: E402

cl.WINDOWS_TSV = HERE / 'elements.tsv'
PROTEIN = ('missense', 'nonsense', 'stop gained', 'stop_gained', 'stop lost', 'stop_lost', 'frameshift', 'inframe',
           'initiator', 'start lost', 'start_lost')
PAD = 2000


def stars(review):
    r = str(review)
    if r in ('practice guideline', 'reviewed by expert panel', 'criteria provided, multiple submitters, no conflicts'):
        return 2
    return 1 if r.startswith('criteria provided') else 0


def mane(out, w):
    p = out / 'mane' / f'{w.window}.json'
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        url = (f'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=mane;chrom={w.chrom};'
               f'start={w.start - 1 - PAD};end={w.end + PAD}')
        with urllib.request.urlopen(url, timeout=60) as r:
            p.write_text(r.read().decode())
        time.sleep(0.5)
    return [t for t in json.loads(p.read_text()).get('mane', []) if t['maneStat'] == 'MANE Select']


def region(t, pos):
    """Region of 1-based pos on MANE transcript t ('' if outside the transcript and its 2 kb upstream)."""
    s0, e0, plus = t['chromStart'], t['chromEnd'], t['strand'] == '+'
    q = pos - 1
    if s0 <= q < e0:
        bs = [int(v) for v in t['blockSizes'].rstrip(',').split(',')]
        st = [int(v) for v in t['chromStarts'].rstrip(',').split(',')]
        if not any(s0 + o <= q < s0 + o + b for b, o in zip(bs, st)):
            return 'intron'
        if t['thickStart'] <= q < t['thickEnd']:
            return 'CDS'
        return "5'UTR" if (q < t['thickStart']) == plus else "3'UTR"
    tss = s0 + 1 if plus else e0
    if (plus and tss - PAD <= pos < tss) or (not plus and tss < pos <= tss + PAD):
        return 'upstream (<= 2 kb)'
    return ''


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, out = Path(a.root), Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    cv = ts.clinvar_table(root)
    cv = cv[cv.in_window.astype(bool)].copy()
    cons = {}
    for w in cl.windows().window:
        r = json.load(open(root / 'clinvar' / f'{w}.json'))['result']
        for u in r['uids']:
            cons[r[u]['accession']] = ';'.join(sorted(r[u].get('molecular_consequence_list', [])))
    cv['consequence'] = cv.accession.map(cons).fillna('')
    cv['protein_changing'] = cv.consequence.str.lower().apply(lambda c: any(k in c for k in PROTEIN))
    cv['synonymous'] = cv.consequence.str.contains('synonymous') & ~cv.protein_changing
    cv['stars'] = cv.review.map(stars)
    W = cl.windows().set_index('window')
    lab, oth = [], []
    for x in cv.itertuples():
        tx = mane(out, W.loc[x.window].to_frame().T.assign(window=x.window).iloc[0])
        own = [t for t in tx if t['geneName2'] == W.loc[x.window, 'gene']]
        r0 = region(own[0], int(x.pos)) if own and pd.notna(x.pos) else ''
        lab.append(r0 or 'other (> 2 kb from the MANE gene)')
        oth.append(';'.join(f"{t['geneName2']}:{region(t, int(x.pos))}" for t in tx
                            if t['geneName2'] != W.loc[x.window, 'gene'] and pd.notna(x.pos) and region(t, int(x.pos))))
    cv['mane_region'], cv['other_mane_genes'] = lab, oth
    cv['benchmark'] = cv.cls.isin(['plp', 'blb']) & (cv.stars >= 1) & ~cv.protein_changing
    cv.to_csv(out / 'clinvar_annotated.tsv', sep='\t', index=False)
    rows = []
    for w, c in cv.groupby('window', sort=False):
        r = dict(element=w)
        for cls in ('plp', 'blb'):
            x = c[(c.cls == cls) & (c.stars >= 1)]
            r[f'{cls}_1star'] = len(x)
            r[f'{cls}_protein_changing_dropped'] = int(x.protein_changing.sum())
            r[f'{cls}_noncoding'] = int((~x.protein_changing).sum())
            r[f'{cls}_noncoding_snv'] = int((~x.protein_changing & (x.kind == 'SNV')).sum())
            r[f'{cls}_synonymous_kept'] = int(x.synonymous.sum())
        rows.append(r)
    counts = pd.DataFrame(rows)
    counts.to_csv(out / 'clinvar_counts.tsv', sep='\t', index=False)
    print(counts.to_string(index=False))
    b = cv[cv.benchmark]
    print('\nbenchmark records by MANE region:')
    print(pd.crosstab(b.window, b.mane_region).to_string())
    print('\ndropped (protein-changing, >= 1 star P/LP or B/LB):')
    d = cv[cv.cls.isin(['plp', 'blb']) & (cv.stars >= 1) & cv.protein_changing]
    print(d[['window', 'cls', 'stars', 'title', 'consequence']].to_string(index=False))


if __name__ == '__main__':
    main()

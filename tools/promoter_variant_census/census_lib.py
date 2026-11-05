"""Shared helpers for the promoter variant census: windows, reference sequence, gnomAD parsing, window membership.

Coordinates are 1-based inclusive hg38 throughout. Only pandas and the standard library, so the AlphaGenome venv
can import it too.
"""
import json
import time
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
RULES = json.load(open(HERE / 'rules.json'))
GROUP_MIN_AN = 2000
GROUPS_EXCLUDED = {'remaining'}
# AlphaGenome scoring context: one 1 Mb interval per context; SLC39A4 W1/W2 share one (union midpoint)
CONTEXT = {'LDLR': 'LDLR', 'HBB': 'HBB', 'GATA2_IS': 'GATA2_IS', 'GATA2_IG': 'GATA2_IG',
           'SLC39A4_W1': 'SLC39A4', 'SLC39A4_W2': 'SLC39A4', 'GATA2_E95_W1': 'GATA2_E95', 'GATA2_E95_W2': 'GATA2_E95'}


def windows():
    return pd.read_csv(HERE / 'windows.tsv', sep='\t')


def contexts():
    """context -> (chrom, gene, lo, hi) spanning its windows."""
    w = windows().assign(context=lambda d: d.window.map(CONTEXT))
    return {c: (g.chrom.iloc[0], g.gene.iloc[0], int(g.start.min()), int(g.end.max())) for c, g in w.groupby('context')}


def reference(out, w):
    """Window sequence on the + strand (UCSC hg38), cached under out/reference/."""
    p = Path(out) / 'reference' / f'{w.window}.txt'
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        url = f'https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={w.chrom};start={w.start - 1};end={w.end}'
        with urllib.request.urlopen(url, timeout=60) as r:
            p.write_text(json.loads(r.read())['dna'].upper())
        time.sleep(1)
    s = p.read_text().strip()
    assert len(s) == w.end - w.start + 1, (w.window, len(s))
    return s


def all_snvs(out):
    """Every possible SNV in every window: window, chrom, pos, ref, alt, key."""
    rows = []
    for w in windows().itertuples():
        for i, ref in enumerate(reference(out, w)):
            for alt in 'ACGT':
                if alt != ref:
                    rows.append(dict(window=w.window, gene=w.gene, chrom=w.chrom, pos=w.start + i, ref=ref, alt=alt))
    return pd.DataFrame(rows).assign(key=lambda d: d.chrom + ':' + d.pos.astype(str) + ':' + d.ref + ':' + d.alt)


def in_window(pos, ref, alt, lo, hi):
    """SNV: pos in window. Indel (VCF left-anchored): any deleted base pos+1..pos+len(ref)-1, or for an insertion the
    base after the insertion point (pos+1), lies in [lo, hi]. Complex ref/alt of equal length > 1: any ref base."""
    if len(ref) == 1 and len(alt) == 1:
        return lo <= pos <= hi
    if len(ref) == len(alt):
        return pos <= hi and pos + len(ref) - 1 >= lo
    if len(ref) > 1:
        return pos + 1 <= hi and pos + len(ref) - 1 >= lo
    return lo <= pos + 1 <= hi


def gnomad_table(out):
    """One row per (window, gnomAD record in the window). Genome and exome fields kept separately."""
    rows = []
    for w in windows().itertuples():
        for v in json.load(open(Path(out) / 'gnomad' / f'{w.window}.json'))['data']['region']['variants']:
            if not in_window(v['pos'], v['ref'], v['alt'], w.start, w.end):
                continue
            g, e = v['genome'], v['exome']
            row = dict(window=w.window, gene=w.gene, chrom=w.chrom, pos=v['pos'], ref=v['ref'], alt=v['alt'],
                       gnomad_id=v['variant_id'], rsid=';'.join(v['rsids']),
                       kind='SNV' if len(v['ref']) == 1 and len(v['alt']) == 1 else 'indel',
                       g_ac=g['ac'] if g else 0, g_an=g['an'] if g else 0, g_af=(g['af'] or 0) if g else 0,
                       g_filters=';'.join(g['filters']) if g else '', g_hom=g['homozygote_count'] if g else 0,
                       e_ac=e['ac'] if e else 0, e_an=e['an'] if e else 0, e_filters=';'.join(e['filters']) if e else '')
            grp = [(p['id'], p['ac'], p['an']) for p in (g['populations'] if g else [])
                   if p['id'] not in GROUPS_EXCLUDED and '_' not in p['id'] and p['id'] not in ('XX', 'XY')]
            big = [(i, ac / an) for i, ac, an in grp if an >= GROUP_MIN_AN]
            row['g_grpmax_af'] = max((af for _, af in big), default=0.0)
            row['g_grpmax_group'] = max(big, key=lambda x: x[1])[0] if big and row['g_grpmax_af'] > 0 else ''
            row['g_groups_af_gt_1pct'] = ';'.join(i for i, af in big if af > 0.01)
            rows.append(row)
    t = pd.DataFrame(rows)
    return t.assign(key=t.chrom + ':' + t.pos.astype(str) + ':' + t.ref + ':' + t.alt,
                    gtex_id=t.chrom + '_' + t.pos.astype(str) + '_' + t.ref + '_' + t.alt + '_b38')

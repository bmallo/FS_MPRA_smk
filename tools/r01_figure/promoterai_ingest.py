"""Pull PromoterAI precomputed scores for the 633 designed LDLR promoter SNVs (chr19:11,089,285-11,089,495, hg38).

Reads a TSV/CSV (optionally .gz/.bgz) in chunks, keeps rows in the element window, and matches them to designed SNVs
by position, ref and alt (ref must equal our reference base). Column names are detected case-insensitively
(chrom/chr/#chrom, pos/position, ref, alt, score/promoterai/promoterai_score, optional strand, gene, transcript); pass
--cols chrom,pos,ref,alt,score to override. If several rows match one SNV (e.g. several TSSs/genes), rows for gene
LDLR on the + strand are preferred; all matches are kept in matches.tsv. A file with a tabix index is read by region.
Writes scores.tsv (snv, score[, gene, strand, transcript]) for merge_evidence.py --promoterai, matches.tsv,
summary.json.
Usage: python tools/r01_figure/promoterai_ingest.py --scores <file> --root results/analysis_v2_20261028 --output <dir>
"""
import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
CHROM, LO, HI, OFFSET = 'chr19', 11089285, 11089495, 11089264
ALIASES = {'chrom': ['chrom', 'chr', '#chrom', 'chromosome', '#chr'], 'pos': ['pos', 'position', 'start'],
           'ref': ['ref', 'reference'], 'alt': ['alt', 'alternate'],
           'score': ['score', 'promoterai', 'promoterai_score', 'promoterai score'], 'tss_pos': ['tss_pos', 'tss'],
           'strand': ['strand'], 'gene': ['gene', 'gene_name', 'symbol'], 'transcript': ['transcript', 'transcript_id', 'tx']}


def detect(cols, override=None):
    if override:
        keys = ['chrom', 'pos', 'ref', 'alt', 'score']
        return dict(zip(keys, override.split(',')))
    low = {c.lower().strip(): c for c in cols}
    out = {}
    for k, al in ALIASES.items():
        for a in al:
            if a in low:
                out[k] = low[a]
                break
    missing = [k for k in ('chrom', 'pos', 'ref', 'alt', 'score') if k not in out]
    if missing:
        raise SystemExit(f'cannot find columns {missing} among {list(cols)}; use --cols')
    return out


def read_window(path, override):
    p = str(path)
    if Path(p + '.tbi').exists() or Path(p + '.csi').exists():
        hdr = subprocess.run(f'tabix -H {p}', shell=True, capture_output=True, text=True).stdout.strip().splitlines()
        txt = subprocess.run(['tabix', p, f'{CHROM}:{LO}-{HI}'], capture_output=True, text=True).stdout
        if not txt.strip():
            txt = subprocess.run(['tabix', p, f'19:{LO}-{HI}'], capture_output=True, text=True).stdout
        names = hdr[-1].lstrip('#').split('\t') if hdr else None
        from io import StringIO
        df = pd.read_csv(StringIO(txt), sep='\t', header=None, names=names)
        return df, detect(df.columns, override)
    sep = ',' if '.csv' in p else '\t'
    if p.endswith('.gz') and sep == '\t':
        # plain gzip (no index): prefilter the window with awk (keeps the header), far faster than pandas chunks
        cmd = (f"zcat '{p}' | awk -F'\\t' 'NR==1 || (($1==\"{CHROM}\" || $1==\"19\") && $2>={LO} && $2<={HI})'")
        txt = subprocess.run(['bash', '-c', cmd], capture_output=True, text=True, check=True).stdout
        from io import StringIO
        df = pd.read_csv(StringIO(txt), sep='\t')
        return df, detect(df.columns, override)
    it = pd.read_csv(p, sep=sep, chunksize=2_000_000, low_memory=False)
    keep, cmap = [], None
    for ch in it:
        cmap = cmap or detect(ch.columns, override)
        c = ch[cmap['chrom']].astype(str).str.replace('^chr', '', regex=True)
        pos = pd.to_numeric(ch[cmap['pos']], errors='coerce')
        keep.append(ch[(c == '19') & pos.between(LO, HI)])
    return pd.concat(keep), cmap


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scores', required=True)
    ap.add_argument('--root', required=True)
    ap.add_argument('--cols')
    ap.add_argument('--prefer-transcript', default='ENST00000558518',
                    help='transcript (prefix match) to prefer among LDLR rows; default MANE Select (NM_000527.5)')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    s = importlib.util.spec_from_file_location('clinical', HERE / 'clinical.py')
    cl = importlib.util.module_from_spec(s)
    sys.modules['clinical'] = cl
    s.loader.exec_module(cl)
    ref, vis = cl.reference(Path(a.root))
    snvs = set(vis[vis.cls == 'SNV'].source_variant)

    df, cm = read_window(a.scores, a.cols)
    df = df.rename(columns={v: k for k, v in cm.items()})
    df['oligo'] = df.pos.astype(int) - OFFSET
    df['snv'] = df.oligo.astype(str) + df.ref.astype(str).str.upper() + '>' + df.alt.astype(str).str.upper()
    df['ref_ok'] = [ref[o - 21] == r.upper() if 21 <= o <= 231 and len(r) == 1 else False
                    for o, r in zip(df.oligo, df.ref.astype(str))]
    m = df[df.snv.isin(snvs)].copy()
    m.to_csv(out / 'matches.tsv', sep='\t', index=False)
    pref = m.copy()
    pref['rank'] = 0
    if 'gene' in pref:
        pref['rank'] += (pref.gene.astype(str).str.upper() != 'LDLR').astype(int) * 2
    if 'strand' in pref:
        pref['rank'] += (~pref.strand.astype(str).isin(['+', '1', '+1'])).astype(int) * 4
    if 'transcript' in pref and a.prefer_transcript:
        pref['rank'] += (~pref.transcript.astype(str).str.startswith(a.prefer_transcript)).astype(int)
    best = pref.sort_values(['snv', 'rank']).drop_duplicates('snv')
    cols = ['snv', 'score'] + [c for c in ('gene', 'strand', 'transcript', 'tss_pos') if c in best]
    best[cols].to_csv(out / 'scores.tsv', sep='\t', index=False)
    json.dump(dict(file=str(a.scores), columns=cm, rows_in_window=len(df), ref_mismatch=int((~df.ref_ok).sum()),
                   designed_matched=int(best.snv.nunique()), designed_total=len(snvs),
                   snvs_with_multiple_rows=int((m.snv.value_counts() > 1).sum()),
                   transcripts=m.get('transcript', pd.Series(dtype=str)).astype(str).value_counts().to_dict(),
                   chosen_transcripts=best.get('transcript', pd.Series(dtype=str)).astype(str).value_counts().to_dict()), open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

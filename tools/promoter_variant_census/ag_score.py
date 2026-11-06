"""Score the census variants with the AlphaGenome API (settings: rules.json "alphagenome").

Variants: every possible SNV in every window plus gnomAD indels with genome AC >= 1 in a window. One 1 Mb interval per
context (census_lib.CONTEXT; SLC39A4 W1/W2 share one), recommended RNA_SEQ / CAGE / DNASE scorers. Writes, resumably:
  alphagenome/tracks.tsv  per variant x track, for the kept biosamples (KEEP; RNA_SEQ rows for the target gene only)
  alphagenome/maxabs.tsv  per variant x output type, max |quantile_score| and |raw_score| over ALL tracks (RNA_SEQ:
                          target gene only) with the track holding the max
Runs in the alphagenome venv; the key is read from ~/.config/alphagenome/key and never printed. Requests sequential.
Usage: python tools/promoter_variant_census/ag_score.py --output <dir> [--windows-tsv <file>]
"""
import argparse
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import census_lib as cl  # noqa: E402

KEEP = {
    'EFO:0001187': 'HepG2', 'UBERON:0002107': 'liver', 'UBERON:0001114': 'liver', 'UBERON:0001115': 'liver',
    'EFO:0001182': 'HEK293', 'EFO:0002067': 'K562', 'UBERON:0013756': 'venous blood', 'UBERON:0000178': 'blood',
    'UBERON:0002371': 'bone marrow', 'UBERON:0012168': 'umbilical cord blood', 'CL:0002092': 'bone marrow cell',
    'CL:0001059': 'CMP CD34+', 'CL:0000837': 'hematopoietic multipotent progenitor',
    'UBERON:0002108': 'small intestine', 'UBERON:0002114': 'duodenum', 'EFO:0001099': 'Caco-2',
    'CL:0002563': 'intestinal epithelial cell',
    # immune promoters and TERT (tools/immune_promoter_scan/rules.json "alphagenome")
    'CL:0000084': 'T-cell', 'CL:0000624': 'CD4 T cell', 'CL:0000625': 'CD8 T cell', 'CL:0000895': 'naive CD4 T cell',
    'CL:0000236': 'B cell', 'CL:0000623': 'NK cell', 'CL:0000576': 'monocyte', 'CL:0001054': 'CD14 monocyte',
    'CL:2000001': 'PBMC', 'CLO:0007045': 'Jurkat E6-1', 'EFO:0002796': 'Jurkat', 'EFO:0002784': 'GM12878',
    'UBERON:0002106': 'spleen', 'UBERON:0002370': 'thymus', 'UBERON:0000029': 'lymph node',
    'EFO:0003042': 'H1', 'EFO:0003045': 'H9', 'EFO:0002791': 'HeLa-S3', 'EFO:0001086': 'A549', 'EFO:0001203': 'MCF-7',
    'EFO:0003072': 'SK-N-SH'}
TYPES = ('RNA_SEQ', 'CAGE', 'DNASE')
COLS = ['key', 'context', 'output_type', 'gene_name', 'track_name', 'track_strand', 'ontology_curie', 'biosample_name',
        'raw_score', 'quantile_score']


def variants(out):
    s = cl.all_snvs(out)[['window', 'gene', 'chrom', 'pos', 'ref', 'alt', 'key']]
    g = cl.gnomad_table(out)
    ind = g[(g.kind == 'indel') & (g.g_ac > 0)][s.columns]
    v = pd.concat([s, ind]).assign(context=lambda d: d.window.map(cl.context))
    return v.drop_duplicates(['context', 'key'])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output', required=True)
    ap.add_argument('--windows-tsv', help='window table (default windows.tsv next to this script)')
    ap.add_argument('--limit', type=int, default=0, help='score only the first N pending variants (testing)')
    a = ap.parse_args()
    if a.windows_tsv:
        cl.WINDOWS_TSV = Path(a.windows_tsv)
    out = Path(a.output)
    d = out / 'alphagenome'
    d.mkdir(parents=True, exist_ok=True)
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers as vs
    m = dna_client.create(Path('~/.config/alphagenome/key').expanduser().read_text().strip())
    L = dna_client.SEQUENCE_LENGTH_1MB
    ctx = cl.contexts()
    scorers = [vs.RECOMMENDED_VARIANT_SCORERS[k] for k in TYPES]
    v = variants(out)
    tp, mp = d / 'tracks.tsv', d / 'maxabs.tsv'
    done = set(zip(*pd.read_csv(mp, sep='\t', usecols=['context', 'key']).values.T)) if mp.exists() else set()
    todo = [r for r in v.itertuples() if (r.context, r.key) not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(len(v), 'variant-contexts;', len(todo), 'to score', flush=True)
    t0 = time.time()
    for i, r in enumerate(todo):
        chrom, gene, lo, hi = ctx[r.context]
        mid = (lo + hi) // 2
        iv = genome.Interval(chrom, mid - L // 2, mid + L // 2)
        for attempt in range(5):
            try:
                s = m.score_variant(interval=iv, variant=genome.Variant(chrom, int(r.pos), r.ref, r.alt),
                                    variant_scorers=scorers)
                break
            except Exception as e:                       # transient API errors: back off and retry
                print(f'{r.key}: attempt {attempt + 1} failed: {type(e).__name__}', flush=True)
                time.sleep(15 * (attempt + 1))
        else:
            raise RuntimeError(f'{r.key} failed 5 times')
        df = vs.tidy_scores([s])
        df = df[~(df.output_type == 'RNA_SEQ') | (df.gene_name == gene)]
        df = df.assign(key=r.key, context=r.context, track_strand=df.get('track_strand', ''))
        df[df.ontology_curie.isin(KEEP)][COLS].to_csv(tp, sep='\t', index=False, mode='a', header=not tp.exists(),
                                                     float_format='%.6g')
        mx = []
        for ot, g in df.groupby('output_type'):
            j = g.quantile_score.abs().idxmax()
            mx.append(dict(context=r.context, key=r.key, output_type=ot, n_tracks=len(g),
                           max_abs_q=abs(g.quantile_score[j]), q_at_max=g.quantile_score[j],
                           track_at_max=g.track_name[j], biosample_at_max=g.biosample_name[j],
                           max_abs_raw=g.raw_score.abs().max()))
        pd.DataFrame(mx).to_csv(mp, sep='\t', index=False, mode='a', header=not mp.exists(), float_format='%.6g')
        if i % 100 == 0:
            print(f'{i + 1}/{len(todo)} {time.time() - t0:.0f} s', flush=True)
    print('done', time.time() - t0, 's', flush=True)


if __name__ == '__main__':
    main()

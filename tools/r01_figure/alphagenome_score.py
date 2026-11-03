"""Score the 633 designed LDLR promoter SNVs with the AlphaGenome API (settings: alphagenome_rules.json).

Runs in the alphagenome venv (/mmfs1/gscratch/stergachislab/bmallo/large_home/envs/alphagenome). The API key is read
from ~/.config/alphagenome/key and never printed. Requests are sequential. Per-variant rows for the kept biosamples
are appended to tracks.tsv (resumable: variants already present are skipped); summarise() then writes scores.tsv
(one row per SNV, columns used by merge_evidence.py --alphagenome).
Usage: python tools/r01_figure/alphagenome_score.py --snvs <evidence.tsv> --output <dir> [--summarise-only]
"""
import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

OLIGO_OFFSET, MID = 11089264, 11089390
KEEP = {'EFO:0001187': 'hepg2', 'UBERON:0002107': 'liver', 'UBERON:0001114': 'liver', 'UBERON:0001115': 'liver',
        'EFO:0001182': 'hek293', 'EFO:0001184': 'hek293'}
TF_SETS = {'SP': ('SP1', 'SP2', 'SP3', 'SP4'), 'KLF': ('KLF',), 'SREBF': ('SREBF1', 'SREBF2'),
           'CREB_ATF': ('CREB1', 'CREM', 'ATF1', 'ATF2', 'ATF3', 'ATF4', 'ATF7')}
COLS = ['snv', 'output_type', 'gene_name', 'track_name', 'ontology_curie', 'biosample_name', 'transcription_factor',
        'raw_score', 'quantile_score']


def score(snvs, out):
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers as vs
    m = dna_client.create(Path('~/.config/alphagenome/key').expanduser().read_text().strip())
    L = dna_client.SEQUENCE_LENGTH_1MB
    iv = genome.Interval('chr19', MID - L // 2, MID + L // 2)
    scorers = [vs.RECOMMENDED_VARIANT_SCORERS[k] for k in ('RNA_SEQ', 'CAGE', 'DNASE', 'CHIP_TF')]
    path = out / 'tracks.tsv'
    done = set(pd.read_csv(path, sep='\t', usecols=['snv']).snv) if path.exists() else set()
    t0 = time.time()
    for i, v in enumerate(snvs):
        if v in done:
            continue
        k, ref, alt = int(v[:-3]), v[-3], v[-1]
        for attempt in range(4):
            try:
                s = m.score_variant(interval=iv, variant=genome.Variant('chr19', OLIGO_OFFSET + k, ref, alt),
                                    variant_scorers=scorers)
                break
            except Exception as e:                       # transient API errors: back off and retry
                print(f'{v}: attempt {attempt + 1} failed: {type(e).__name__}', flush=True)
                time.sleep(10 * (attempt + 1))
        else:
            raise RuntimeError(f'{v} failed 4 times')
        df = vs.tidy_scores([s])
        df = df[df.ontology_curie.isin(KEEP)]
        df = df[~df.output_type.isin(['RNA_SEQ']) | (df.gene_name == 'LDLR')]
        df = df.assign(snv=v)[COLS]
        df.to_csv(path, sep='\t', index=False, mode='a', header=not path.exists(), float_format='%.6g')
        if i % 50 == 0:
            print(f'{i + 1}/{len(snvs)} {time.time() - t0:.0f} s', flush=True)


def summarise(out):
    t = pd.read_csv(out / 'tracks.tsv', sep='\t')
    t['cell'] = t.ontology_curie.map(KEEP)
    rows = {}
    for (v, ot, cell), g in t.groupby(['snv', 'output_type', 'cell']):
        r = rows.setdefault(v, {})
        if ot == 'CHIP_TF':
            tf = g.transcription_factor.fillna('').str.upper()
            for name, pref in TF_SETS.items():
                h = g[tf.str.startswith(pref)]
                if len(h):
                    r[f'chiptf_{name}_{cell}'] = h.raw_score.mean()
            r[f'chiptf_maxabs_{cell}'] = g.raw_score.loc[g.raw_score.abs().idxmax()]
            continue
        r[f'{ot.lower()}_{cell}'] = g.raw_score.mean()
        r[f'{ot.lower()}_{cell}_q'] = g.quantile_score.mean()
    s = pd.DataFrame.from_dict(rows, orient='index')
    s.index.name = 'snv'
    s.to_csv(out / 'scores.tsv', sep='\t', float_format='%.6g')
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--snvs', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--summarise-only', action='store_true')
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    if not a.summarise_only:
        snvs = pd.read_csv(a.snvs, sep='\t').snv.tolist()
        score(snvs, out)
    s = summarise(out)
    print(f'scores for {len(s)} SNVs, {s.shape[1]} columns', flush=True)


if __name__ == '__main__':
    main()

"""Replace a site extraction's footprint coverage bits with plasmid-LLR calls.

Writes a new extraction directory: scores/known/qualities, meta.npz, variants.json and
molecules.tsv.gz are symlinked from the source extraction; coverage.npy is rebuilt from LLR
calls (global mode, one min_llr) in the same bit layout site_effects reads:
  1 TF-sized call (10-60 bp), 2 MSP, 4 nucleosome-sized call (85-250 bp), 8 call > 250 bp.
MSP = element positions not inside any call longer than 60 bp, so TF-sized footprints sit inside
MSPs as they do for FiberHMM. Call lengths are classified before clipping to the element.
reverse=True maps an element on the reference minus strand into oligo order (HEK293T).

Usage: python -m fsmpra.llr_coverage --extract <dir> --calls <a.parquet> [<b.parquet> ...] \
           --element 1234 1445 --min-llr 5 --output <new dir> [--reverse]
"""
import argparse
import csv
import os
import time
from pathlib import Path

import numpy as np

from .io import atomic_output, open_text, read_json, source_info, write_json

LINKED = ('scores.npy', 'known.npy', 'qualities.npy', 'meta.npz', 'variants.json', 'molecules.tsv.gz')


def call_bits(length):
    return (1 if 10 <= length <= 60 else 0) | (4 if 85 <= length <= 250 else 0) | (8 if length > 250 else 0)


def coverage_from_calls(rows, starts, ends, n_rows, element, reverse=False):
    """rows: extraction row per call; starts/ends: 0-based half-open reference intervals."""
    a, b = element
    width = b - a
    cov = np.zeros((n_rows, width), dtype=np.uint8)
    blocked = np.zeros((n_rows, width), dtype=bool)
    for r, s, e in zip(rows, starts, ends):
        lo, hi = max(s, a) - a, min(e, b) - a
        if lo >= hi:
            continue
        length = e - s
        cov[r, lo:hi] |= call_bits(length)
        if length > 60:
            blocked[r, lo:hi] = True
    cov[~blocked] |= 2
    return cov[:, ::-1].copy() if reverse else cov


def run(extract, calls, element, min_llr, output, reverse=False):
    import pandas as pd
    started = time.monotonic()
    extract = Path(extract).resolve()
    audit = read_json(extract / 'audit.json')
    with open_text(extract / 'molecules.tsv.gz') as h:
        index = {r['read_id']: int(r['row']) for r in csv.DictReader(h, delimiter='\t')}
    frames = []
    for p in calls:
        df = pd.read_parquet(p, columns=['min_llr', 'read_name', 'start', 'end'])
        frames.append(df[np.isclose(df.min_llr, min_llr)])
    df = pd.concat(frames, ignore_index=True)
    called_reads = set(df.read_name.unique())
    df = df[df.read_name.isin(index.keys())]
    missing = [r for r in index if r not in called_reads]
    rows = df.read_name.map(index).to_numpy()
    cov = coverage_from_calls(rows, df.start.to_numpy() - 1, df.end.to_numpy(), len(index), element, reverse)
    with atomic_output(output) as out:
        for name in LINKED:
            os.symlink(extract / name, out / name)
        np.save(out / 'coverage.npy', cov)
        write_json(out / 'audit.json', dict(
            audit, sources=audit['sources'] + [source_info(Path(p)) for p in calls] + [source_info(Path(__file__), True)],
            source_extraction=str(extract),
            llr_coverage=dict(min_llr=min_llr, element_reference_interval0=list(element), reversed=reverse,
                              calls=int(len(df)), molecules=len(index), molecules_without_calls_file_entry=len(missing),
                              msp_fraction_mean=float(((cov & 2) > 0).mean()),
                              nuc_85_250_fraction_mean=float(((cov & 4) > 0).mean()),
                              tf_fraction_mean=float(((cov & 1) > 0).mean()),
                              definition='MSP = not inside an LLR call > 60 bp; bits 1 TF 10-60, 4 nuc 85-250, 8 > 250; '
                                         'lengths before clipping to the element.')))
    print(f'{len(index):,} molecules, {len(df):,} calls, {len(missing):,} reads with no calls in any parquet '
          f'(called fully accessible), {time.monotonic() - started:.0f}s', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--extract', required=True)
    ap.add_argument('--calls', nargs='+', required=True)
    ap.add_argument('--element', nargs=2, type=int, required=True)
    ap.add_argument('--min-llr', type=float, default=5.0)
    ap.add_argument('--output', required=True)
    ap.add_argument('--reverse', action='store_true')
    a = ap.parse_args()
    run(a.extract, a.calls, tuple(a.element), a.min_llr, a.output, a.reverse)


if __name__ == '__main__':
    main()

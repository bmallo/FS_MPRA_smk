"""Step A8 (part 1): read-name CRC32 and exact barcode for every A6 molecule (settings: rules.json call_confidence).

Walks the same reads with the same filters as a6_extract.py (a6x.molecule, same subsample) and, for each kept
molecule, records CRC32(read name) (emission fold = CRC32 % 2) and the reference-strand bases aligned to the barcode
interval (empty if any base there is unaligned or an insertion sits inside it). With --flank-regex (HEK293T, whose
rebuilt reference holds a 14-N placeholder for a 20-bp barcode), the barcode is instead the 17-23 bp found between
the 12-bp reference flanks in the read sequence, the rule of fsmpra/site_extract_hek.py. Read names are not stored.
The kept count and per-molecule distal fraction must equal the A6 extraction's, so rows line up with molecules.npz.

Usage: python tools/analysis_v2/a8_ids.py --a6 <a6 extract dir> --output <new dir>
"""
import argparse
import importlib.util
import json
import re
import time
import zlib
from pathlib import Path

import numpy as np
import pysam

_spec = importlib.util.spec_from_file_location('a6x', Path(__file__).with_name('a6_extract.py'))
a6x = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a6x)
a5x = a6x.a5x


def barcode(rec, start, end):
    """Query bases aligned to reference [start, end), or '' if any position is unaligned or an insertion is inside."""
    seq = rec.query_sequence
    got = []
    inside = False
    for q, r in rec.get_aligned_pairs():
        if r is not None and start <= r < end:
            inside = True
            if q is None:
                return ''
            got.append(seq[q])
        elif inside and r is None and q is not None and len(got) < end - start:
            return ''                                    # insertion between barcode bases
        elif r is not None and r >= end:
            break
    return ''.join(got).upper() if len(got) == end - start else ''


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--a6', required=True, help='A6 extraction dir (its extract.json gives the BAM and settings)')
    ap.add_argument('--flank-regex', action='store_true', help='barcode between reference flanks (HEK293T rule)')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    cfg = json.load(open(Path(a.a6) / 'extract.json'))
    element, bc, wt_region = tuple(cfg['element']), tuple(cfg['barcode']), tuple(cfg['wt_region'])
    ref = a5x.ref_seq(cfg['reference'])
    L = len(ref)
    ref_arr = np.frombuffer(ref.encode(), np.uint8)
    at_idx = np.flatnonzero((ref_arr == ord('A')) | (ref_arr == ord('T')))
    dmask = a5x.distal_mask(L, element, bc)
    flank = re.compile(ref[bc[0] - 12:bc[0]] + '([ACGT]{17,23})' + ref[bc[1]:bc[1] + 12]) if a.flank_regex else None
    crc, bcs, distal, t0 = [], [], [], time.time()
    with pysam.AlignmentFile(cfg['bam']) as f:
        for rec in f.fetch(cfg['contig']):
            if rec.is_secondary or rec.is_supplementary or not a5x.keep_read(rec.query_name, cfg['fraction']):
                continue
            st, x = a6x.molecule(rec, ref_arr, at_idx, L, dmask, wt_region)
            if x is None:
                continue
            crc.append(zlib.crc32(rec.query_name.encode()))
            if flank is not None:
                m = flank.search(rec.query_sequence.upper())
                bcs.append(m.group(1) if m else '')
            else:
                bcs.append(barcode(rec, *bc))
            distal.append(x['distal_frac'])
            if len(crc) % 20000 == 0:
                print(f'{len(crc)} molecules, {time.time() - t0:.0f} s', flush=True)
    z = np.load(Path(a.a6) / 'molecules.npz')
    if len(crc) != z['distal_frac'].size or not np.allclose(distal, z['distal_frac']):
        raise ValueError(f'sidecar does not line up with A6: {len(crc)} vs {z["distal_frac"].size} molecules')
    bcs = np.array(bcs)
    np.savez_compressed(out / 'ids.npz', name_crc=np.array(crc, np.uint32), barcode=bcs)
    s = dict(a6=a.a6, molecules=len(crc), barcode_readable=int((bcs != '').sum()),
             distinct_barcodes=int(np.unique(bcs[bcs != '']).size), barcode_interval=bc, flank_regex=a.flank_regex,
             seconds=round(time.time() - t0))
    json.dump(s, open(out / 'ids.json', 'w'), indent=1)
    print(json.dumps(s))


if __name__ == '__main__':
    main()

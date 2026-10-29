"""HEK293T (older LDLR library) site extraction in the replicate-3 array contract.

Writes the same outputs as site_extract.py (scores/known/qualities/coverage .npy, n x 211 in the
HepG2 oligo frame, plus meta.npz, variants.json, molecules.tsv.gz, audit.json) so site_effects
machinery can read either cell. The 211-bp element sits reverse-complemented at HEK ref0
[3236, 3447): arrays are reversed into oligo order and alleles complemented onto HepG2 IDs.

Selection: molecules marked eligible by tools/hek_compare/extract_hek.py whose 318-bp mutable
promoter (ref0 3184-3502) is exact WT or carries exactly one substitution inside the shared
element, with no indel and every base >= Q20. Coverage bits are left zero here and filled from
plasmid-LLR calls afterwards (fsmpra.llr_coverage).

Clone QC: the 20-bp barcode is read between its flanks; per clone the script records molecule
count, genotype purity and backbone positions where most of the clone's reads differ from the
reference (candidate clone-specific mutations outside the promoter).
"""
import argparse
import csv
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from .io import atomic_output, open_text, read_json, source_info, write_json, write_table
from .raw_candidate_extract import project_read

ARRAYS = [('scores', np.uint8), ('known', np.bool_), ('qualities', np.uint8), ('coverage', np.uint8)]
MUT = (3184, 3502)            # older library mutable promoter, ref0 half-open
ELEM = (3236, 3447)           # shared 211-bp element, reverse strand; oligo = 3446 - ref0 + 21
BARCODE = (1343, 1357)        # N run in the rebuilt reference (the barcode itself is 20 bp)
COMP = str.maketrans('ACGT', 'TGCA')
META = ['passes', 'read_quality', 'mapq', 'whole_m6a', 'distal_m6a', 'distal_nucleosome', 'nucleosome_segments',
        'element_msp', 'element_tf_small']


def reference(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>')).upper()


def span(read, a, b):
    """Query [qa, qb) aligned to reference [a, b) with no indel inside, else None."""
    rp, qp, qa = read.reference_start, 0, None
    for op, size in read.cigartuples:
        if op in (0, 7, 8):
            if qa is None and rp <= a < rp + size:
                qa = qp + a - rp
            if qa is not None and rp < b <= rp + size:
                qb = qp + b - rp
                return (qa, qb) if qb - qa == b - a else None
            rp += size
            qp += size
        elif op in (1, 4):
            if qa is not None:
                return None
            qp += size
        elif op in (2, 3):
            if rp + size > a and rp < b:
                return None
            rp += size
    return None


def ref_aligned(read, length):
    """Reference-indexed read bases (0 = not aligned); matches only."""
    out = np.zeros(length, dtype=np.uint8)
    seq = np.frombuffer(read.query_sequence.encode(), dtype=np.uint8)
    rp, qp = read.reference_start, 0
    for op, size in read.cigartuples:
        if op in (0, 7, 8):
            out[rp:rp + size] = seq[qp:qp + size]
            rp += size
            qp += size
        elif op in (1, 4):
            qp += size
        elif op in (2, 3):
            rp += size
    return out


def oligo_id(ref0, ref_base, alt):
    return f'{ELEM[1] - 1 - ref0 + 21}{ref_base.translate(COMP)}>{alt.translate(COMP)}'


def run(bam_path, molecules, ref_path, variants_path, output, limit=None):
    import pysam
    started = time.monotonic()
    ref = reference(ref_path)
    L = len(ref)
    rs = ref[MUT[0]:MUT[1]]
    variants = read_json(variants_path)
    if [v['class'] for v in variants].count('WT') != 1 or len(variants) != 634:
        raise ValueError('Expected WT plus 633 SNVs')
    wt_seq = next(v['sequence'] for v in variants if v['class'] == 'WT')
    if ref[ELEM[0]:ELEM[1]].translate(COMP)[::-1] != wt_seq:
        raise ValueError('HEK element does not reverse-complement onto the HepG2 WT element')
    byid = {('WT' if v['class'] == 'WT' else f"{v['oligo_pos1']}{v['ref']}>{v['alt']}"): i for i, v in enumerate(variants)}
    wt_i = byid['WT']
    flank = re.compile(ref[BARCODE[0] - 12:BARCODE[0]] + '([ACGT]{17,23})' + ref[BARCODE[1]:BARCODE[1] + 12])
    mol = {}
    with open_text(molecules) as h:
        for r in csv.DictReader(h, delimiter='\t'):
            if r['status'] == 'eligible':
                mol[r['read_id']] = r
    rows, meta = [], {k: [] for k in ['variant_index', 'barcode', 'supported', 'eligible', 'eligible_with_barcode', 'circular',
                                      'minus', 'evaluated'] + META}
    barcodes = {}
    genotype_counts = Counter()
    clone_genotypes = defaultdict(Counter)
    clone_mismatch = defaultdict(lambda: np.zeros(L, dtype=np.int32))
    clone_depth = defaultdict(lambda: np.zeros(L, dtype=np.int32))
    refarr = np.frombuffer(ref.encode(), dtype=np.uint8)
    chunks = {name: [] for name, _ in ARRAYS}
    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            r = mol.get(read.query_name)
            if r is None:
                continue
            sp = span(read, *MUT)
            if sp is None:
                genotype_counts['promoter_indel'] += 1
                continue
            seq = read.query_sequence
            if min(read.query_qualities[sp[0]:sp[1]]) < 20:
                genotype_counts['low_quality'] += 1
                continue
            seg = seq[sp[0]:sp[1]]
            diffs = [i for i in range(len(seg)) if seg[i] != rs[i]]
            m = flank.search(seq)
            bc = m.group(1) if m else ''
            if len(diffs) == 0:
                key = 'WT'
            elif len(diffs) == 1 and ELEM[0] <= MUT[0] + diffs[0] < ELEM[1]:
                p = MUT[0] + diffs[0]
                key = oligo_id(p, ref[p], seg[diffs[0]])
            else:
                key = 'flank_snv' if len(diffs) == 1 else 'multi'
            genotype_counts[key if key in ('flank_snv', 'multi') else 'WT' if key == 'WT' else 'shared_snv'] += 1
            if bc:
                clone_genotypes[bc][key] += 1
                al = ref_aligned(read, L)
                clone_depth[bc] += al > 0
                clone_mismatch[bc] += (al > 0) & (al != refarr)
            if key not in byid:
                continue
            expected = ref[ELEM[0]:ELEM[1]] if key == 'WT' else None
            if expected is None:
                p = MUT[0] + diffs[0]
                expected = ref[ELEM[0]:p] + seg[diffs[0]] + ref[p + 1:ELEM[1]]
            values = project_read(read, expected, start=ELEM[0], end=ELEM[1])
            for (name, _), value in zip(ARRAYS, values):
                chunks[name].append(np.asarray(value)[::-1].copy())
            k = len(rows)
            rows.append(dict(row=k, read_id=read.query_name, observed_variant_id=variants[byid[key]]['variant_id'],
                             source_variant=key, barcode=bc, strand='-' if read.is_reverse else '+'))
            meta['variant_index'].append(byid[key])
            meta['barcode'].append(barcodes.setdefault(bc, len(barcodes)) if bc else -1)
            for flag in ('supported', 'eligible_with_barcode'):
                meta[flag].append(bool(bc))
            meta['eligible'].append(True)
            meta['circular'].append(False)
            meta['minus'].append(read.is_reverse)
            meta['evaluated'].append(True)
            meta['passes'].append(float(r['passes']))
            meta['read_quality'].append(float(r['read_quality']))
            meta['mapq'].append(float(r['mapq']))
            meta['whole_m6a'].append(float(r['whole_m6a_fraction']) if r['whole_m6a_fraction'] else np.nan)
            meta['distal_m6a'].append(float(r['distal_m6a_fraction']) if r['distal_m6a_fraction'] else np.nan)
            for name in ('distal_nucleosome', 'nucleosome_segments', 'element_msp', 'element_tf_small'):
                meta[name].append(np.nan)
            if limit and len(rows) >= limit:
                break
    n = len(rows)
    dtypes = dict(variant_index=np.int16, barcode=np.int32, supported=bool, eligible=bool, eligible_with_barcode=bool,
                  circular=bool, minus=bool, evaluated=bool)
    meta = {k: np.asarray(v, dtype=dtypes.get(k, np.float32)) for k, v in meta.items()}
    cluster = meta['barcode'].astype(np.int64)
    lone = cluster < 0
    cluster[lone] = len(barcodes) + np.arange(lone.sum())
    meta['cluster'] = cluster
    backbone = np.ones(L, bool)
    backbone[MUT[0]:MUT[1]] = False
    backbone[BARCODE[0] - 3:BARCODE[1] + 3] = False
    clones = []
    for bc, g in clone_genotypes.items():
        total = sum(g.values())
        if total < 2:
            continue
        top, top_n = g.most_common(1)[0]
        depth, mism = clone_depth[bc], clone_mismatch[bc]
        cons = np.flatnonzero(backbone & (depth >= 5) & (mism > 0.5 * np.maximum(depth, 1)))
        clones.append(dict(barcode=bc, molecules=total, top_genotype=top, top_fraction=top_n / total,
                           genotypes=len(g), backbone_consensus_mismatches=len(cons),
                           backbone_mismatch_ref0=','.join(map(str, cons[:20]))))
    clones.sort(key=lambda c: -c['molecules'])
    with atomic_output(output) as out:
        for name, dt in ARRAYS:
            np.save(out / f'{name}.npy', np.asarray(chunks[name], dtype=dt).reshape(n, 211))
        np.savez(out / 'meta.npz', **meta)
        write_json(out / 'variants.json', variants)
        write_table(out / 'molecules.tsv.gz', rows, list(rows[0]))
        write_table(out / 'clones.tsv', clones, list(clones[0]))
        per_variant = Counter(r['source_variant'] for r in rows)
        write_json(out / 'audit.json', dict(
            schema_version='fsmpra.site_extract.v1', cell='HEK293T', library='LDLR older 318-bp promoter',
            complete=limit is None, molecule_limit=limit, molecules=n, barcodes=len(barcodes),
            sources=[source_info(Path(p), Path(p).suffix not in ('.bam', '.gz')) for p in (bam_path, molecules, ref_path, variants_path)]
                    + [source_info(Path(__file__), True)],
            genotype_counts=dict(genotype_counts), variants_observed=len([k for k in per_variant if k != 'WT']),
            element_reference_interval0=list(ELEM), orientation='arrays reversed into HepG2 oligo order; alleles complemented',
            seconds=time.monotonic() - started,
            arrays='scores: retained ML byte; known: A/T opportunity; qualities: base quality; coverage: zero here, '
                   'filled from plasmid-LLR calls (bits 1 TF 10-60, 2 MSP, 4 nuc 85-250, 8 nuc >250).',
            limitations='Existing plasmid-aligned BAM (fibertools 0.5.4 m6A, kinetics stripped); reference rebuilt by consensus. '
                        'Most HEK293T SNVs are a single plasmid clone.'))
    print(f'Complete: {n:,} molecules, {len(barcodes):,} barcodes, {time.monotonic() - started:.0f}s', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--molecules', required=True, help='molecules.tsv.gz from tools/hek_compare/extract_hek.py')
    ap.add_argument('--reference', required=True)
    ap.add_argument('--variants', required=True, help='HepG2 site extraction variants.json (WT + 633 SNVs)')
    ap.add_argument('--output', required=True)
    ap.add_argument('--limit', type=int)
    a = ap.parse_args()
    run(a.bam, a.molecules, a.reference, a.variants, a.output, a.limit)


if __name__ == '__main__':
    main()

"""Step A9 (part 1): per-molecule m6A calls, genotype and clone for every analysis-ready WT and SNV molecule.

Read set (rules.json statistics.molecules):
  --ledger (HepG2): fsmpra molecule ledger, exact design, sequence/length/methylation QC passed, class WT or SNV.
  --hek-molecules (HEK293T): fsmpra site_extract_hek molecules.tsv.gz (WT promoter or one element substitution).
Clones (rules.json clones): key = barcode + genotype; a read is dropped when its barcode has >= 3 full-length
qualified reads with a majority genotype >= 0.9 that the read contradicts; no readable barcode = its own clone.
HepG2 WT keeps all molecules of 25% of WT clones (crc32(key) % 4 == 0); every SNV molecule is kept.

m6A calls and distal labeling come from a6_extract.molecule (A6 read filters, no WT check) on the same reads,
joined by read name. The BAM is walked by --procs workers, each taking crc32(read name) % procs.
Writes molecules.npz (A6 format) and meta.npz (genotype label, clone id, name CRC, read quality) in table order,
clones.tsv and extract.json. Read names are not stored.

Usage: python tools/analysis_v2/a9_extract.py --a6 <A6 extract dir of the run> \
           (--ledger <ledger dir> --catalog <catalog.json> | --hek-molecules <molecules.tsv.gz>) \
           --procs 4 --output <new dir>
"""
import argparse
import importlib.util
import json
import multiprocessing as mp
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

_spec = importlib.util.spec_from_file_location('a6x', Path(__file__).with_name('a6_extract.py'))
a6x = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a6x)
a5x = a6x.a5x

WT_KEEP_MOD = 4
DROP_MIN_READS, DROP_MAJORITY = 3, 0.9
LEDGER_COLS = ['read_id', 'observed_variant_id', 'genotype_class', 'genotype_status', 'barcode', 'barcode_status',
               'sequence_qc', 'length_qc', 'methylation_qc', 'read_quality']


def crc(s):
    return zlib.crc32(s.encode())


def clone_keys(read_id, genotype, barcode, qualified):
    """Clone key per read and the reads to drop. qualified marks the reads that vote on a barcode's majority."""
    t = pd.DataFrame(dict(g=genotype, b=barcode, q=qualified))
    v = t[t.q & (t.b != '')].groupby(['b', 'g']).size().rename('n').reset_index()
    v = v.sort_values(['b', 'n'], ascending=[True, False])
    tot = v.groupby('b').n.sum()
    top = v.drop_duplicates('b').set_index('b')
    maj = pd.DataFrame(dict(reads=tot, majority=top.g, share=top.n / tot))
    j = t.join(maj, on='b')
    drop = ((j.reads >= DROP_MIN_READS) & (j.share >= DROP_MAJORITY) & (j.g != j.majority)).to_numpy()
    key = np.where(t.b != '', t.b + '|' + t.g, '#' + pd.Series(read_id).map(lambda r: str(crc(r))).to_numpy())
    return key, drop


def ledger_table(ledger_dir, catalog_path):
    d = pd.read_csv(Path(ledger_dir) / 'ledger/molecule_ledger.tsv.gz', sep='\t', usecols=LEDGER_COLS,
                    keep_default_na=False, dtype={'barcode': str, 'observed_variant_id': str})
    ok = (d.genotype_status == 'exact_design') & (d.sequence_qc == 1)
    d['barcode'] = np.where(d.barcode_status == 'ok', d.barcode, '')
    qualified = (ok & (d.length_qc == 1) & (d.barcode != '')).to_numpy()
    key, drop = clone_keys(d.read_id.to_numpy(), d.observed_variant_id.to_numpy(), d.barcode.to_numpy(), qualified)
    d['key'], d['drop'] = key, drop
    ready = ok & (d.length_qc == 1) & (d.methylation_qc == 1) & d.genotype_class.isin(['WT', 'SNV'])
    d = d[ready].copy()
    cat = pd.DataFrame(json.load(open(catalog_path))['library']['variants']).set_index('variant_id')
    d['label'] = d.observed_variant_id.map(cat.source_variant)
    d['cls'] = d.genotype_class
    d['rq'] = pd.to_numeric(d.read_quality, errors='coerce')
    return d[['read_id', 'label', 'cls', 'key', 'drop', 'rq']].reset_index(drop=True)


def hek_table(path):
    d = pd.read_csv(path, sep='\t', keep_default_na=False, dtype={'barcode': str})
    key, drop = clone_keys(d.read_id.to_numpy(), d.observed_variant_id.to_numpy(), d.barcode.to_numpy(),
                           (d.barcode != '').to_numpy())
    d['key'], d['drop'] = key, drop
    d['label'] = d.source_variant
    d['cls'] = np.where(d.source_variant == 'WT', 'WT', 'SNV')
    d['rq'] = np.nan
    return d[['read_id', 'label', 'cls', 'key', 'drop', 'rq']].reset_index(drop=True)


_G = {}


def _work(w):
    cfg, want, procs = _G['cfg'], _G['want'], _G['procs']
    ref = a5x.ref_seq(cfg['reference'])
    L = len(ref)
    ref_arr = np.frombuffer(ref.encode(), np.uint8)
    at_idx = np.flatnonzero((ref_arr == ord('A')) | (ref_arr == ord('T')))
    dmask = a5x.distal_mask(L, tuple(cfg['element']), tuple(cfg['barcode']))
    rows, vecs, vals, status, t0 = [], [], [], {}, time.time()
    with pysam.AlignmentFile(cfg['bam']) as f:
        for rec in f.fetch(cfg['contig']):
            if rec.is_secondary or rec.is_supplementary:
                continue
            name = rec.query_name
            if crc(name) % procs != w:
                continue
            row = want.get(name)
            if row is None:
                continue
            st, x = a6x.molecule(rec, ref_arr, at_idx, L, dmask, None)
            status[st] = status.get(st, 0) + 1
            if x is None:
                continue
            rows.append(row)
            vecs.append(np.packbits(x['vec'] >= 0))
            vecs.append(np.packbits(x['vec'] == 1))
            vals.append((x['distal_frac'], x['passes'], x['strand']))
            if w == 0 and len(rows) % 50000 == 0:
                print(f'worker 0: {len(rows)} molecules, {time.time() - t0:.0f} s', flush=True)
    nb = (at_idx.size + 7) // 8
    v = np.array(vecs, np.uint8).reshape(-1, 2, nb)
    return np.array(rows, np.int64), v, np.array(vals, float).reshape(-1, 3), status, at_idx


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--a6', required=True, help="A6 extract dir of the run (its extract.json gives BAM and layout)")
    ap.add_argument('--ledger')
    ap.add_argument('--catalog')
    ap.add_argument('--hek-molecules')
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    t0 = time.time()
    cfg = json.load(open(Path(a.a6) / 'extract.json'))
    if a.ledger:
        t = ledger_table(a.ledger, a.catalog)
        wt_sub = True
    else:
        t = hek_table(a.hek_molecules)
        wt_sub = False
    n_ready = t.cls.value_counts().to_dict()
    n_drop = t[t['drop']].cls.value_counts().to_dict()
    t = t[~t['drop']]
    if wt_sub:
        keep_wt = t.key.map(lambda k: crc(k) % WT_KEEP_MOD == 0)
        t = t[(t.cls != 'WT') | keep_wt]
    t = t.reset_index(drop=True)
    print(f'table: {len(t)} reads ({t.cls.value_counts().to_dict()}), {time.time() - t0:.0f} s', flush=True)
    _G.update(cfg=cfg, want=dict(zip(t.read_id, range(len(t)))), procs=a.procs)
    with mp.get_context('fork').Pool(a.procs) as pool:
        parts = pool.map(_work, range(a.procs))
    rows = np.concatenate([p[0] for p in parts])
    order = np.argsort(rows)
    rows = rows[order]
    V = np.concatenate([p[1] for p in parts])[order]
    vals = np.concatenate([p[2] for p in parts])[order]
    status = {}
    for p in parts:
        for k, v in p[3].items():
            status[k] = status.get(k, 0) + v
    at_idx = parts[0][4]
    s = t.iloc[rows].reset_index(drop=True)
    uniq, clone = np.unique(s.key.to_numpy(), return_inverse=True)
    np.savez_compressed(out / 'molecules.npz', pos=at_idx, observed=V[:, 0], methylated=V[:, 1],
                        distal_frac=vals[:, 0], passes=vals[:, 1].astype(np.int32), strand=vals[:, 2].astype(np.int8))
    np.savez_compressed(out / 'meta.npz', label=s.label.to_numpy().astype(str), cls=s.cls.to_numpy().astype(str),
                        clone=clone.astype(np.int64), name_crc=np.array([crc(r) for r in s.read_id], np.uint32),
                        read_quality=s.rq.to_numpy(float))
    c = pd.DataFrame(dict(clone=clone, label=s.label, readable=~s.key.str.startswith('#')))
    c.groupby(['clone', 'label', 'readable']).size().rename('molecules').reset_index() \
        .to_csv(out / 'clones.tsv', sep='\t', index=False)
    info = dict(a6=a.a6, ledger=a.ledger, catalog=a.catalog, hek_molecules=a.hek_molecules, bam=cfg['bam'],
                ready_reads=n_ready, dropped_contradicting=n_drop, wt_clone_fraction=1 / WT_KEEP_MOD if wt_sub else 1.0,
                table_reads=int(len(t)), status=status, kept=int(len(s)), kept_by_class=s.cls.value_counts().to_dict(),
                clones=int(uniq.size), snvs=int(s[s.cls == 'SNV'].label.nunique()), seconds=round(time.time() - t0))
    json.dump(info, open(out / 'extract.json', 'w'), indent=1)
    print(json.dumps(info))


if __name__ == '__main__':
    main()

"""R01 figure: is there any TF motif under footprint F04 (oligo 60-86, protection ~66-81)?

Scans the WT LDLR element (oligo 21-231) and SNV alleles with every motif in JASPAR 2026 CORE vertebrates
(non-redundant) and HOCOMOCO v14 H14CORE, both strands. Scoring follows FIMO's defaults: frequencies with a 0.1
pseudocount spread by the background, log2-odds against the element's base composition, exact P-value of each
score from the motif's score distribution under that background (scores binned at 0.01 bits), q-value = BH over
every motif x position x strand test on the WT element (FIMO's q). Reported hits use FIMO's default P < 1e-4.

Positive controls: the literature Sp1 repeats (R1 89-104, R3 140-155) and SRE-1 (124-139) should be found.

Writes to --output:
  wt_hits.tsv       every WT match with P < 1e-3 (motif, db, TF, family, start/end oligo, strand, score, P, q)
  windows.tsv       per window (F04 core 66-81, F04 span 60-86, R1, SRE-1, R3, every other 16-bp tile): number of
                    motifs with P < 1e-4 overlapping by >= half their length, best P, best motif
  alleles.tsv       for every SNV at oligo 55-91: per motif, best P over matches overlapping the SNV, WT vs variant;
                    rows kept when either is < 1e-4
  summary.json      versions, md5, counts
Usage: python tools/r01_figure/motif_scan.py --root results/analysis_v2_20261028 --db <motif_db dir> --output <dir>
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('clinical', Path(__file__).with_name('clinical.py'))
cl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cl)

BASES = 'ACGT'
PSEUDO = 0.1                # FIMO --motif-pseudo default
BIN = 0.01                  # bits per score bin for the exact P-value
P_HIT, P_KEEP = 1e-4, 1e-3
JASPAR = 'JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt'
HOCO = 'H14CORE_jaspar_format.txt'
HOCO_ANN = 'H14CORE_annotation.jsonl'
WINDOWS = {'F04 core 66-81': (66, 81), 'F04 span 60-86': (60, 86), 'R1 Sp1 89-104': (89, 104),
           'R2 SRE-1 124-139': (124, 139), 'R3 Sp1 140-155': (140, 155)}


def read_motifs(path, db):
    """JASPAR-style file: '>ID[ \\t]NAME' then 4 rows A C G T (with or without 'A [ ... ]')."""
    lines = [x.strip() for x in open(path) if x.strip()]
    out = []
    for i, x in enumerate(lines):
        if not x.startswith('>'):
            continue
        h = x[1:].split(None, 1)
        rows = []
        for r in lines[i + 1:i + 5]:
            r = r.replace('[', ' ').replace(']', ' ').split()
            if r[0] in BASES:
                r = r[1:]
            rows.append([float(v) for v in r])
        out.append(dict(id=h[0], name=h[1] if len(h) > 1 else h[0].split('.')[0], db=db, counts=np.array(rows).T))
    return out


def log_odds(counts, bg):
    f = counts / counts.sum(1, keepdims=True)
    f = (f + PSEUDO * bg) / (1 + PSEUDO)
    return np.log2(f / bg)


def pvalue_table(lo, bg):
    """Exact survival function of the integer-binned score under an iid background: (min score bin, P(S >= s))."""
    s = np.round(lo / BIN).astype(int)
    mn = s.min(1).sum()
    dist = np.ones(1)                             # P(score = mn + i) after each position
    for row in s:
        rmin = row.min()
        new = np.zeros(len(dist) + row.max() - rmin)
        for b in range(4):
            k = row[b] - rmin
            new[k:k + len(dist)] += dist * bg[b]
        dist = new
    sf = np.cumsum(dist[::-1])[::-1]
    return mn, sf, s


def scan(seq, lo_int, mn, sf):
    """Best-strand integer scores and P per start position. Returns arrays (fwd P, rev P) of length len(seq)-L+1."""
    idx = np.array([BASES.index(c) for c in seq])
    L = lo_int.shape[0]
    n = len(seq) - L + 1
    if n <= 0:
        return np.ones(0), np.ones(0)
    win = np.lib.stride_tricks.sliding_window_view(idx, L)
    fwd = lo_int[np.arange(L), win].sum(1)
    rc = lo_int[::-1, ::-1]                       # reverse complement of the motif
    rev = rc[np.arange(L), win].sum(1)
    pf = sf[np.clip(fwd - mn, 0, len(sf) - 1)]
    pr = sf[np.clip(rev - mn, 0, len(sf) - 1)]
    return pf, pr


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--db', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, db, out = Path(a.root), Path(a.db), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    ref, vis = cl.reference(root)
    bg = np.array([ref.count(b) for b in BASES], float)
    bg = (bg + bg[::-1]) / (2 * len(ref))            # strand-symmetric element composition

    fam = {}
    for line in open(db / HOCO_ANN):
        r = json.loads(line)
        fam[r['name']] = (r['masterlist_info']['tf'], r['masterlist_info'].get('tfclass_family', ''))
    motifs = read_motifs(db / JASPAR, 'JASPAR2026') + read_motifs(db / HOCO, 'HOCOMOCOv14')
    for m in motifs:
        lo = log_odds(m['counts'], bg)
        m['mn'], m['sf'], m['int'] = pvalue_table(lo, bg)
        if m['db'] == 'HOCOMOCOv14':
            m['tf'], m['family'] = fam.get(m['id'], (m['name'], ''))
        else:
            m['tf'], m['family'] = m['name'], ''

    # WT scan
    rows = []
    n_tests = 0
    for m in motifs:
        pf, pr = scan(ref, m['int'], m['mn'], m['sf'])
        n_tests += 2 * len(pf)
        L = m['int'].shape[0]
        for strand, pv in (('+', pf), ('-', pr)):
            for i in np.flatnonzero(pv < P_KEEP):
                rows.append(dict(db=m['db'], motif=m['id'], tf=m['tf'], family=m['family'], start=cl.LO + i,
                                 end=cl.LO + i + L - 1, strand=strand, p=pv[i],
                                 site=ref[i:i + L] if strand == '+' else ref[i:i + L][::-1].translate(
                                     str.maketrans('ACGT', 'TGCA'))))
    wt = pd.DataFrame(rows).sort_values('p')
    # BH over all tests (untested remainder have P >= 1e-3, so ranks of kept rows are exact)
    wt['q'] = np.minimum.accumulate((wt.p * n_tests / np.arange(1, len(wt) + 1))[::-1])[::-1].clip(upper=1)
    wt.to_csv(out / 'wt_hits.tsv', sep='\t', index=False, float_format='%.4g')

    # windows
    hits = wt[wt.p < P_HIT]
    wins = dict(WINDOWS)
    for s in range(cl.LO, cl.HI - 15 + 1, 8):
        wins.setdefault(f'tile {s}-{s + 15}', (s, s + 15))
    wrows = []
    for name, (s, e) in wins.items():
        ov = (np.minimum(hits.end, e) - np.maximum(hits.start, s) + 1)
        h = hits[ov >= (hits.end - hits.start + 1) / 2]
        best = h.iloc[0] if len(h) else None
        # chance expectation: (motif, start, strand) placements overlapping the window by >= half the motif, x 1e-4
        n_place = 0
        for m in motifs:
            L = m['int'].shape[0]
            st = np.arange(cl.LO, cl.HI - L + 2)
            n_place += 2 * int(((np.minimum(st + L - 1, e) - np.maximum(st, s) + 1) >= L / 2).sum())
        wrows.append(dict(window=name, start=s, end=e, sequence=ref[s - cl.LO:e - cl.LO + 1],
                          expected_matches_p1e4=n_place * P_HIT, matches_p1e4=len(h), motifs_p1e4=h.motif.nunique(),
                          jaspar=int((h.drop_duplicates('motif').db == 'JASPAR2026').sum()),
                          hocomoco=int((h.drop_duplicates('motif').db == 'HOCOMOCOv14').sum()),
                          best_p=best.p if best is not None else np.nan,
                          best_q=best.q if best is not None else np.nan,
                          best_motif=f'{best.tf} ({best.motif})' if best is not None else '',
                          best_site=f'{best.start}-{best.end}{best.strand} {best.site}' if best is not None else '',
                          tfs=';'.join(h.drop_duplicates('tf').tf.head(12))))
    pd.DataFrame(wrows).to_csv(out / 'windows.tsv', sep='\t', index=False, float_format='%.4g')

    # alleles at 55-91
    snvs = vis[(vis.cls == 'SNV') & vis.oligo_pos1.between(55, 91)].source_variant.tolist()
    arows = []
    for v in snvs:
        k, alt = int(v[:-3]), v[-1]
        seq = cl.edit(ref, k, k, alt)
        for m in motifs:
            L = m['int'].shape[0]
            best = []
            for sq in (ref, seq):
                pf, pr = scan(sq, m['int'], m['mn'], m['sf'])
                st = np.arange(len(pf)) + cl.LO
                cov = (st <= k) & (st + L - 1 >= k)
                best.append(min(pf[cov].min(), pr[cov].min()) if cov.any() else 1.0)
            if min(best) < P_HIT:
                arows.append(dict(snv=v, oligo=k, db=m['db'], motif=m['id'], tf=m['tf'], family=m['family'],
                                  p_wt=best[0], p_variant=best[1],
                                  log10_ratio=np.log10(best[1]) - np.log10(best[0])))
    pd.DataFrame(arows).sort_values(['oligo', 'snv', 'p_variant']).to_csv(
        out / 'alleles.tsv', sep='\t', index=False, float_format='%.4g')

    json.dump(dict(
        databases={JASPAR: md5(db / JASPAR), HOCO: md5(db / HOCO), HOCO_ANN: md5(db / HOCO_ANN)},
        motifs=dict(JASPAR2026=sum(m['db'] == 'JASPAR2026' for m in motifs),
                    HOCOMOCOv14=sum(m['db'] == 'HOCOMOCOv14' for m in motifs)),
        background=dict(zip(BASES, bg.round(4).tolist())), wt_tests=n_tests,
        wt_hits_p1e4=int((wt.p < P_HIT).sum()), wt_hits_q05=int((wt.q < 0.05).sum()),
        settings=dict(pseudocount=PSEUDO, bin_bits=BIN, p_hit=P_HIT)), open(out / 'summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

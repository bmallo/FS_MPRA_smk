"""R01 figure, step 3 (clinical variants): map reported LDLR promoter variants onto the library and collect coverage,
activity and the stored analysis v2 results. No new statistics.

Sources (downloaded 2026-11-02 into <output>/../clinical_raw, treated as data):
  ClinVar   esearch 19[chr] AND 11089285:11089495[chrpos38] + esummary (GRCh38 SPDI, classification, traits)
  LOVD      shared LDLR variants (databases.lovd.nl REST API, all LDLR entries; hg38 or hg19 position)
  paper     plasmid Fiber-seq paper, Fig. 6b / Extended Data Fig. 12 (PMC12154926, preprint version):
            c.-223A>G, c.-214T>C, c.-214T>A, c.-156C>T, c.-193C>T, c.-79C>G (sixth variant in Fig. 6b unnamed)

Every variant is turned into the element sequence it produces (oligo 21-231) and matched to designed SNVs and
1-bp indels by sequence, so homopolymer/repeat indels match whatever position the design used.
Coordinates: oligo k = chr19:11,089,264 + k (hg38); NM_000527.4 c.N (N < 0) = oligo 285 + N (checked against
ClinVar SNVs in summary.json).

Writes clinical_variants.tsv (one row per distinct variant sequence) and summary.json to --output.
Usage: python tools/r01_figure/clinical.py --root results/analysis_v2_20261028 --raw <clinical_raw> --output <dir>
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

OLIGO_OFFSET = 11089264
LO, HI = 21, 231
C_OFFSET = 285                  # oligo = 285 + c for c < 0 (NM_000527.4)
PAPER = [('c.-223A>G', 'Fig. 6b'), ('c.-214T>C', 'Fig. 6b'), ('c.-214T>A', 'Ext. Data Fig. 12'),
         ('c.-156C>T', 'Fig. 6b'), ('c.-193C>T', 'Fig. 6b'), ('c.-79C>G', 'Fig. 6b')]


def reference(root):
    v = pd.read_csv(root / 'a4_informative/visibility/variants_visibility.tsv', sep='\t')
    s = v[v.cls == 'SNV'].drop_duplicates('oligo_pos1').set_index('oligo_pos1')['ref']
    assert list(s.index) == list(range(LO, HI + 1))
    return ''.join(s.values), v


def edit(ref, start, end, ins):
    """Replace oligo start..end (1-based inclusive; end = start - 1 for a pure insertion after start-1) by ins."""
    i, j = start - LO, end - LO + 1
    return ref[:i] + ins + ref[j:]


def design_sequences(ref, vis):
    out = {}
    for name, cls in vis[['source_variant', 'cls']].itertuples(index=False):
        m = re.fullmatch(r'(\d+)([ACGT])>([ACGT])', name)
        if m:
            k = int(m[1]); assert ref[k - LO] == m[2]
            out[edit(ref, k, k, m[3])] = (name, 'SNV')
            continue
        m = re.fullmatch(r'(\d+)del([ACGT])', name)
        if m:
            k = int(m[1]); assert ref[k - LO] == m[2]
            out[edit(ref, k, k, '')] = (name, '1-bp deletion')
            continue
        m = re.fullmatch(r'(\d+)_(\d+)ins([ACGT])', name)
        if m:
            out[edit(ref, int(m[2]), int(m[1]), m[3])] = (name, '1-bp insertion')
    return out


def from_spdi(spdi, ref):
    """NC_000019.10:pos0:del:ins -> (oligo start, end, ins) after trimming shared bases; None outside the element."""
    _, p0, d, i = spdi.split(':')
    p = int(p0) + 1
    while d and i and d[0] == i[0]:
        d, i, p = d[1:], i[1:], p + 1
    while d and i and d[-1] == i[-1]:
        d, i = d[:-1], i[:-1]
    s = p - OLIGO_OFFSET
    e = s + len(d) - 1
    if s < LO or e > HI or (not d and s > HI + 1):
        return None
    if d and ref[s - LO:e - LO + 1] != d:
        raise ValueError(f'reference mismatch {spdi}')
    return s, e, i


def from_hgvs_g(g, ref):
    """chr19:g.<hg38 HGVS> simple forms -> (start, end, ins) in oligo; None if unparsed or outside."""
    g = g.split('g.', 1)[1]
    m = re.fullmatch(r'(\d+)([ACGT])>([ACGT])', g)
    if m:
        k = int(m[1]) - OLIGO_OFFSET
        if not LO <= k <= HI or ref[k - LO] != m[2]:
            return None
        return k, k, m[3]
    m = re.fullmatch(r'(\d+)(?:_(\d+))?(del|dup)', g)
    if m:
        a = int(m[1]) - OLIGO_OFFSET
        b = int(m[2] or m[1]) - OLIGO_OFFSET
        if a < LO or b > HI:
            return None
        return (a, b, '') if m[3] == 'del' else (b + 1, b, ref[a - LO:b - LO + 1])
    m = re.fullmatch(r'(\d+)_(\d+)ins([ACGT]+)', g)
    if m:
        a = int(m[1]) - OLIGO_OFFSET
        if not LO <= a + 1 <= HI + 1:
            return None
        return a + 1, a, m[3]
    m = re.fullmatch(r'(\d+)(?:_(\d+))?delins([ACGT]+)', g)
    if m:
        a = int(m[1]) - OLIGO_OFFSET
        b = int(m[2] or m[1]) - OLIGO_OFFSET
        if a < LO or b > HI:
            return None
        return a, b, m[3]
    return None


def kind(s, e, ins):
    n_del, n_ins = e - s + 1, len(ins)
    if n_del == 1 and n_ins == 1:
        return 'SNV'
    if n_del == 1 and n_ins == 0:
        return '1-bp deletion'
    if n_del == 0 and n_ins == 1:
        return '1-bp insertion'
    return f'other ({n_del} del, {n_ins} ins)'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', required=True)
    ap.add_argument('--raw', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    root, raw, out = Path(a.root), Path(a.raw), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    ref, vis = reference(root)
    designs = design_sequences(ref, vis)
    summ = dict(element_length=len(ref), designs=len(designs))

    recs = []   # (source, label, start, end, ins, extra dict)
    # paper
    for c, where in PAPER:
        m = re.fullmatch(r'c\.(-\d+)([ACGT])>([ACGT])', c)
        k = C_OFFSET + int(m[1])
        assert ref[k - LO] == m[2], (c, k, ref[k - LO])
        recs.append(('paper', c, k, k, m[3], dict(paper_figure=where)))
    # ClinVar
    cv = json.load(open(raw / 'clinvar/esummary.json'))['result']
    n_cv, cv_out, c_check = 0, [], []
    for uid in cv['uids']:
        r = cv[uid]
        n_cv += 1
        vs = r['variation_set'][0]
        spdi = vs.get('canonical_spdi', '')
        x = from_spdi(spdi, ref) if spdi else None
        if x is None:
            cv_out.append(r['title'])
            continue
        gc = r.get('germline_classification', {})
        traits = ';'.join(sorted({t['trait_name'] for t in gc.get('trait_set', [])}))
        recs.append(('ClinVar', r['title'], *x, dict(
            clinvar=r['accession'], clinvar_class=gc.get('description', ''), clinvar_review=gc.get('review_status', ''),
            clinvar_traits=traits)))
        m = re.search(r'c\.(-\d+)([ACGT])>([ACGT])', r['title'])
        if m and kind(*x) == 'SNV':
            c_check.append(x[0] - int(m[1]))
    summ['clinvar'] = dict(records=n_cv, outside_or_unparsed=cv_out,
                           oligo_minus_c=sorted(set(c_check)), snvs_checked=len(c_check))
    # LOVD
    lv = json.load(open(raw / 'lovd/ldlr_all.json'))
    pairs = [(int(v['position_genomic']['hg19'].split(':')[1].split('_')[0]),
              int(v['position_genomic']['hg38'].split(':')[1].split('_')[0]))
             for v in lv if 'hg38' in v['position_genomic'] and v['position_genomic']['hg19'] != 'chr19:0']
    off = sorted({h19 - h38 for h19, h38 in pairs if 11_190_000 < h19 < 11_210_000})
    n_lv, lv_skip = 0, 0
    for v in lv:
        vg = v['variant_genomic']
        if 'hg38' in vg:
            g = vg['hg38']
        elif vg.get('hg19', '').startswith('chr19:g.') and len(off) == 1:
            g19 = vg['hg19'].split('g.', 1)[1]
            g = 'chr19:g.' + re.sub(r'\d+', lambda m: str(int(m[0]) - off[0]), g19, count=2)
        else:
            continue
        try:
            x = from_hgvs_g(g, ref)
        except (IndexError, ValueError):
            x = None
        if x is None:
            continue
        n_lv += 1
        recs.append(('LOVD', ';'.join(v['Variant/DNA']), *x, dict(
            lovd=v['id'], lovd_effect=';'.join(v.get('effect_concluded') or v.get('effect_reported') or []),
            lovd_times_reported=v.get('Times_reported', ''))))
    summ['lovd'] = dict(records=len(lv), hg19_minus_hg38=off, in_element=n_lv)

    # collapse by resulting sequence
    act = pd.concat([_mpra(root, f, sep, el) for el, f, sep in
                     (('', 'GRCh38_LDLR.tsv', '\t'), ('2', 'GRCh38_LDLR.2.csv', ','))], axis=1)
    hep = _clones(root / 'a9_statistics/extract/HepG2_run3/clones.tsv', 'hepg2')
    hek = _clones(root / 'a9_statistics/extract/HEK293T_run1/clones.tsv', 'hek')
    s7 = pd.read_csv(root / 'a10_effects/HepG2_run3/snvs.tsv', sep='\t').query('endpoint == "E1"').set_index('snv')
    s5 = pd.read_csv(root / 'a10_effects/HepG2_run3_pm5/snvs.tsv', sep='\t').query('endpoint == "E1"').set_index('snv')
    fn = pd.read_csv(root / 'a11_motifs_function/function_snvs.tsv', sep='\t').set_index('snv')
    fe = pd.read_parquet(root / 'fp_features/effects.parquet').query('mask == "A7" and calls == "confident"')
    fe5 = pd.read_parquet(root / 'fp_features/effects.parquet').query('mask == "pm5" and calls == "confident"')
    fe5 = fe5.set_index(['snv', 'footprint'])
    vmol = vis.set_index('source_variant')['molecules']

    groups = {}
    for src, label, s, e, ins, extra in recs:
        seq = edit(ref, s, e, ins)
        g = groups.setdefault(seq, dict(start=s, end=e, ins=ins, sources=set(), labels=set(), extra={}))
        g['sources'].add(src)
        g['labels'].add(label)
        for k, val in extra.items():
            g['extra'].setdefault(k, set()).add(str(val))
    rows = []
    for seq, g in groups.items():
        name, dkind = designs.get(seq, ('', ''))
        s = g['start']
        row = dict(oligo=s, c_dot=s - C_OFFSET, hg38=OLIGO_OFFSET + s, type=kind(g['start'], g['end'], g['ins']),
                   designed=name or 'not designed', in_paper='paper' in g['sources'],
                   sources=';'.join(sorted(g['sources'])), names=' | '.join(sorted(g['labels'])))
        row.update({k: ';'.join(sorted(v)) for k, v in g['extra'].items()})
        if name:
            row.update(act.loc[name].to_dict() if name in act.index else {})
            if dkind == 'SNV':
                row.update(hep.loc[name].to_dict() if name in hep.index else dict(hepg2_molecules=0, hepg2_clones=0))
                row.update(hek.loc[name].to_dict() if name in hek.index else dict(hek_molecules=0, hek_clones=0))
                if name in s7.index:
                    row.update(a10_simes_p=s7.loc[name, 'simes_p'], a10_selected=bool(s7.loc[name, 'selected']),
                               a10_hits=int(s7.loc[name, 'hits']), a10_group=s7.loc[name, 'group'])
                if name in s5.index:
                    row.update(a10pm5_selected=bool(s5.loc[name, 'selected']), a10pm5_hits=int(s5.loc[name, 'hits']))
                if name in fn.index:
                    row.update(C1=fn.loc[name, 'C1'], C2=fn.loc[name, 'C2'])
                f = fe[(fe.snv == name)]
                sig = f[f.sig]
                row['fp_sig_A7_conf'] = ';'.join(
                    f'{x.footprint}:{100 * x.D:+.2f}' + ('*' if fe5.loc[(name, x.footprint), 'sig'] else '')
                    for x in sig.itertuples())
                if 'F04' in set(f.footprint):
                    x = f[f.footprint == 'F04'].iloc[0]
                    row.update(F04_D_pts=100 * x.D, F04_p=x.p)
            else:
                row['hepg2_molecules'] = int(vmol.get(name, 0))
                row['hepg2_note'] = 'indel: molecules from A4 design table; never analysed'
        rows.append(row)
    t = pd.DataFrame(rows).sort_values(['oligo', 'type', 'names'])
    t.to_csv(out / 'clinical_variants.tsv', sep='\t', index=False, float_format='%.4g')
    cls = t.get('clinvar_class', pd.Series(dtype=str)).fillna('')
    summ['variants'] = dict(
        distinct=len(t), by_type=t.type.value_counts().to_dict(),
        designed=int((t.designed != 'not designed').sum()),
        designed_by_type=t[t.designed != 'not designed'].type.value_counts().to_dict(),
        clinvar_plp=int((cls.str.contains('athogenic') & ~cls.str.contains('onflicting')).sum()),
        clinvar_classes=cls[cls != ''].value_counts().to_dict())
    json.dump(summ, open(out / 'summary.json', 'w'), indent=1, default=str)


def _mpra(root, f, sep, el):
    t = pd.read_csv(Path('/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/mpra_published') / f, sep=sep)
    t['oligo'] = t['Position'] - OLIGO_OFFSET
    t = t[(t.oligo >= LO) & (t.oligo <= HI)]
    t['name'] = np.where(t.Alt == '-', t.oligo.astype(str) + 'del' + t.Ref, t.oligo.astype(str) + t.Ref + '>' + t.Alt)
    return t.set_index('name')[['Value', 'P-Value', 'Tags']].rename(
        columns={'Value': f'act{el}', 'P-Value': f'act{el}_p', 'Tags': f'act{el}_tags'})


def _clones(path, tag):
    c = pd.read_csv(path, sep='\t')
    return c[c.label != 'WT'].groupby('label').agg(**{f'{tag}_molecules': ('molecules', 'sum'),
                                                      f'{tag}_clones': ('clone', 'size')})


if __name__ == '__main__':
    main()

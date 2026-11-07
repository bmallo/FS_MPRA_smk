"""Data for the standalone Fig 3b (fig3_rules.json panel_b_isolated_draft15): the whole Kircher 2019 LDLR element.

Steps (each writes into --output):
  snvs      every SNV of chr19:11,089,231-11,089,548 from the published MPRA table (954), reference bases checked
            against hg38 (results/proxy_benign_20261106/reference/LDLR.txt, 11,088,549-11,089,548); snvs.tsv and
            snvs_extra.tsv (the 321 outside the plasmid library, input to alphagenome_score.py)
  controls  ClinVar germline P/LP and B/LB SNVs in the region (proxy_benign clinvar/LDLR.json); controls.tsv;
            every ClinVar SNV record in the region with its class: clinvar_all.tsv (panel a)
  gtex      GTEx V10 variant-table lookup per control (present, maf01) and significant eQTLs; gtex.tsv
  primate   Ensembl REST EPO primate alignment of the region; per control, non-human primates whose aligned base
            equals the human alt allele; primate.tsv (+ raw alignment JSON)
  assemble  per SNV and readout table (scores.tsv) joining MPRA, AlphaGenome, PromoterAI, Fiber-seq (HepG2, HEK293T),
            gnomAD, GTEx, primate
Network steps are sequential with pauses (klone policy). SNV names follow the library: f'{oligo}{ref}>{alt}',
oligo = hg38 - 11,089,264 (negative upstream of the library); HGVS c. = oligo - 285.
Usage: python tools/r01_figure/fig3b_data.py <step> --output results/r01_figure_20261102/fig3b_data
"""
import argparse
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

LO, HI = 11089231, 11089548
OLIGO_OFFSET, C0 = 11089264, 285
REF_LO = 11088549
ROOT = Path('/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA')
MPRA = ROOT / 'mpra_published/GRCh38_LDLR.tsv'
PAI = ROOT / 'promoterai/promoterAI_tss500.tsv.gz'
PB = Path('results/proxy_benign_20261106')
W = Path('results/r01_figure_20261102')
GTEX = 'https://gtexportal.org/api/v2'
ENSEMBL = 'https://rest.ensembl.org'
TOP = 0.8
MIN_MOLECULES = 100


def name(pos, ref, alt):
    return f'{pos - OLIGO_OFFSET}{ref}>{alt}'


def get(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'fs-mpra-fig3b'})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception as e:                                   # transient: back off
            print(f'attempt {i + 1} failed: {type(e).__name__}', flush=True)
            time.sleep(5 * (i + 1))
    raise RuntimeError(url)


def snvs(out):
    m = pd.read_csv(MPRA, sep='\t')
    m = m[(m.Alt != '-') & m.Position.between(LO, HI)]
    ref = (PB / 'reference/LDLR.txt').read_text().strip()
    assert len(ref) == 1000
    bad = [(p, r) for p, r in zip(m.Position, m.Ref) if ref[p - REF_LO] != r]
    assert not bad, bad[:5]
    t = pd.DataFrame(dict(pos=m.Position, ref=m.Ref, alt=m.Alt, mpra=m.Value, mpra_p=m['P-Value'],
                          mpra_tags=m.Tags))
    t['snv'] = [name(*x) for x in zip(t.pos, t.ref, t.alt)]
    t['oligo'] = t.pos - OLIGO_OFFSET
    t['c_dot'] = t.oligo - C0
    t['in_library'] = t.oligo.between(21, 231)
    assert len(t) == 954 and t.snv.is_unique and t.in_library.sum() == 633
    t.to_csv(out / 'snvs.tsv', sep='\t', index=False)
    t[~t.in_library][['snv']].to_csv(out / 'snvs_extra.tsv', sep='\t', index=False)
    print(f'{len(t)} SNVs, {t.in_library.sum()} in the library, {(~t.in_library).sum()} extra; hg38 ref bases all match')


def cv_class(desc):                       # same mapping as the census truth_sets.py / panel b drafts
    d = desc.lower()
    if 'conflicting' in d:
        return 'conflicting'
    if 'pathogenic' in d and 'benign' not in d and 'uncertain' not in d:
        return 'plp'
    if 'benign' in d and 'pathogenic' not in d and 'uncertain' not in d:
        return 'blb'
    return 'vus' if 'uncertain' in d else 'other'


STARS = {'practice guideline': 4, 'reviewed by expert panel': 3,
         'criteria provided, multiple submitters, no conflicts': 2, 'criteria provided, single submitter': 1,
         'criteria provided, conflicting classifications': 1}


def controls(out):
    res = json.load(open(PB / 'clinvar/LDLR.json'))['result']
    rows = []
    for uid in res['uids']:
        r = res[uid]
        spdi = r['variation_set'][0].get('canonical_spdi', '')
        gc = r.get('germline_classification', {})
        if not spdi:
            continue
        _, p0, d, i = spdi.split(':')
        if len(d) != 1 or len(i) != 1:
            continue
        p = int(p0) + 1
        if not LO <= p <= HI:
            continue
        rows.append(dict(snv=name(p, d, i), pos=p, ref=d, alt=i, accession=r['accession'], title=r['title'],
                         clinvar_desc=gc.get('description', ''), review=gc.get('review_status', ''),
                         stars=STARS.get(gc.get('review_status', ''), 0), cls=cv_class(gc.get('description', ''))))
    t = pd.DataFrame(rows)
    print('ClinVar SNV records in the region by class:', t.cls.value_counts().to_dict())
    assert t.snv.is_unique
    t[['snv', 'pos', 'cls', 'clinvar_desc', 'stars', 'accession']].to_csv(out / 'clinvar_all.tsv', sep='\t', index=False)
    c = t[t.cls.isin(['plp', 'blb'])].sort_values(['cls', 'pos'], ascending=[False, True])
    c.to_csv(out / 'controls.tsv', sep='\t', index=False)
    print(c[['snv', 'cls', 'clinvar_desc', 'stars', 'title']].to_string(index=False))


def gtex(out):
    c = pd.read_csv(out / 'controls.tsv', sep='\t')
    rows = []
    for v in c.itertuples():
        vid = f'chr19_{v.pos}_{v.ref}_{v.alt}_b38'
        r = get(f'{GTEX}/dataset/variant?' + urllib.parse.urlencode({'variantId': vid, 'datasetId': 'gtex_v10'}))
        r = r['data']
        row = dict(snv=v.snv, gtex_id=vid, gtex_present=bool(r), gtex_maf01=bool(r and r[0].get('maf01')),
                   gtex_sig_genes='')
        time.sleep(0.5)
        if row['gtex_present']:
            e = get(f'{GTEX}/association/singleTissueEqtl?' + urllib.parse.urlencode(
                {'variantId': vid, 'datasetId': 'gtex_v10', 'itemsPerPage': 250}))['data']
            row['gtex_sig_genes'] = ';'.join(sorted({x['geneSymbol'] for x in e}))
            time.sleep(0.5)
        rows.append(row)
    t = pd.DataFrame(rows)
    t['gtex_tested'] = t.gtex_maf01
    t.to_csv(out / 'gtex.tsv', sep='\t', index=False)
    print(t.to_string(index=False))


def primate(out):
    url = (f'{ENSEMBL}/alignment/region/homo_sapiens/19:{LO}-{HI}:1?' +
           urllib.parse.urlencode({'method': 'EPO', 'species_set_group': 'primates', 'content-type': 'application/json'}))
    cache = out / 'ensembl_epo_primates.json'
    if cache.exists():
        raw = json.load(open(cache))
    else:
        raw = get(url)
        json.dump(raw, open(cache, 'w'))
    blocks = []
    for b in raw:
        seqs = {s['species']: s for s in b['alignments']}
        print('genomes:', sorted(k for k in seqs if '[' not in k), flush=True)
        h = next(s for k, s in seqs.items() if k == 'homo_sapiens')
        blocks.append((h, seqs))
        print('block human', h['start'], h['end'], h['strand'], 'species', len(seqs), flush=True)
    c = pd.read_csv(out / 'controls.tsv', sep='\t')
    rows = []
    for v in c.itertuples():
        hit, seen, alts = None, [], []
        for h, seqs in blocks:
            if not h['start'] <= v.pos <= h['end']:
                continue
            assert h['strand'] == 1
            col, p = None, h['start'] - 1
            for j, ch in enumerate(h['seq']):
                if ch != '-':
                    p += 1
                    if p == v.pos:
                        col = j
                        break
            assert h['seq'][col].upper() == v.ref, (v.snv, h['seq'][col])
            hit = True
            for sp, s in seqs.items():
                if sp == 'homo_sapiens' or '[' in sp:            # '[n]' = inferred ancestor, not a genome
                    continue
                b = s['seq'][col].upper()
                if b in 'ACGT':
                    seen.append(f'{sp}:{b}')
                    if b == v.alt:
                        alts.append(sp)
        rows.append(dict(snv=v.snv, aligned=bool(hit), primates_with_base=len(seen), primate_alt=bool(alts),
                         primate_alt_species=';'.join(alts), primate_bases=';'.join(seen)))
    t = pd.DataFrame(rows)
    t.to_csv(out / 'primate.tsv', sep='\t', index=False)
    print(t.drop(columns='primate_bases').to_string(index=False))


def gnomad_status():
    """Per SNV in gnomAD v4.1 (genomes or exomes, AC >= 1, PASS), from the proxy_benign 1-kb gnomAD query."""
    rows = []
    for v in json.load(open(PB / 'gnomad/LDLR.json'))['data']['region']['variants']:
        if len(v['ref']) != 1 or len(v['alt']) != 1 or not LO <= v['pos'] <= HI:
            continue
        g, e = v['genome'], v['exome']
        gp = bool(g and g['ac'] > 0 and not g['filters'])
        ep = bool(e and e['ac'] > 0 and not e['filters'])
        rows.append(dict(snv=name(v['pos'], v['ref'], v['alt']), gnomad_genome=gp, gnomad_exome=ep,
                         gnomad=gp or ep, g_ac=g['ac'] if g else 0, g_an=g['an'] if g else 0,
                         e_ac=e['ac'] if e else 0, e_an=e['an'] if e else 0,
                         g_hom=g['homozygote_count'] if g else 0))
    return pd.DataFrame(rows)


def assemble(out):
    """scores.tsv: one row per SNV of the region with every readout (variant - WT) and gnomAD; controls_full.tsv: the
    18 controls with GTEx and primate columns. Library SNVs are checked against the draft-14 inputs (evidence_v3)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location('fig3', Path(__file__).with_name('fig3.py'))
    f3 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(f3)
    t = pd.read_csv(out / 'snvs.tsv', sep='\t')
    ag = pd.concat([pd.read_csv(W / 'alphagenome/scores.tsv', sep='\t'),
                    pd.read_csv(out / 'alphagenome/scores.tsv', sep='\t')])
    assert ag.snv.is_unique and set(t.snv) <= set(ag.snv), 'AlphaGenome scores missing'
    t = t.merge(ag[['snv', 'rna_seq_hepg2', 'dnase_hepg2']].rename(
        columns={'rna_seq_hepg2': 'ag_rna', 'dnase_hepg2': 'ag_dnase'}), on='snv', how='left', validate='1:1')
    pr = out / 'promoterai_region.tsv'                         # region rows of the genome-wide table (awk filter)
    if not pr.exists():
        import subprocess
        subprocess.run(f"zcat {PAI} | awk -F'\\t' 'NR==1 || ($1==\"chr19\" && $2>={LO} && $2<={HI})' > {pr}",
                       shell=True, check=True)
    pai = pd.read_csv(pr, sep='\t')
    pai = pai[(pai.chrom == 'chr19') & pai.pos.between(LO, HI) & (pai.transcript_id == 'ENST00000558518.6')]
    pai['snv'] = [name(*x) for x in zip(pai.pos, pai.ref, pai.alt)]
    t = t.merge(pai[['snv', 'promoterAI']].rename(columns={'promoterAI': 'pai'}), on='snv', how='left', validate='1:1')
    S = f3.readout_scores(W)                                   # draft-14 inputs, library SNVs only
    t = t.merge(S[['snv', 'Fiber-seq|signed']].rename(columns={'Fiber-seq|signed': 'fs_hepg2'}), on='snv', how='left',
                validate='1:1')
    H = f3.hek_scores(W / 'hek_fiberseq_20261104')
    t = t.merge(H[['snv', 'hek_signed', 'hek_n', 'hek_clones']].rename(columns={'hek_signed': 'fs_hek293t'}), on='snv',
                how='left', validate='1:1')
    g = gnomad_status()
    t = t.merge(g, on='snv', how='left', validate='1:1')
    t['gnomad'] = t.gnomad.fillna(False).astype(bool)
    # checks against draft 14 (evidence_v3) for the 633 library SNVs
    e = pd.read_csv(W / 'evidence_v3/evidence.tsv', sep='\t').set_index('snv')
    lib = t[t.in_library].set_index('snv')
    chk = {}
    for col, ecol in (('mpra', 'act'), ('ag_rna', 'ag_rna_seq_hepg2'), ('ag_dnase', 'ag_dnase_hepg2'), ('pai', 'pai_score')):
        d = (lib[col] - e.loc[lib.index, ecol]).abs()
        chk[col] = float(d.max())
        assert d.max() < 1e-3, (col, d.max())
    chk['gnomad_vs_evidence'] = int((lib.gnomad != (e.loc[lib.index, 'gnomad_ac'].fillna(0) > 0)).sum())
    for c in ('mpra', 'ag_rna', 'ag_dnase', 'pai'):
        assert t[c].notna().all(), c
    t.to_csv(out / 'scores.tsv', sep='\t', index=False, float_format='%.6g')
    c = pd.read_csv(out / 'controls.tsv', sep='\t').merge(pd.read_csv(out / 'gtex.tsv', sep='\t'), on='snv',
                                                          validate='1:1')
    c = c.merge(pd.read_csv(out / 'primate.tsv', sep='\t'), on='snv', validate='1:1')
    c = c.merge(t.drop(columns=['pos', 'ref', 'alt']), on='snv', how='left', validate='1:1')
    c.to_csv(out / 'controls_full.tsv', sep='\t', index=False, float_format='%.6g')
    print('max |difference| vs evidence_v3 on library SNVs:', chk)
    print('readout coverage (non-missing SNVs):', {k: int(t[k].notna().sum()) for k in
                                                   ('mpra', 'ag_rna', 'ag_dnase', 'pai', 'fs_hepg2', 'fs_hek293t')})
    print('gnomAD SNVs in region:', int(t.gnomad.sum()), '(library', int(t[t.in_library].gnomad.sum()), ')')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('step', choices=['snvs', 'controls', 'gtex', 'primate', 'assemble'])
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    dict(snvs=snvs, controls=controls, gtex=gtex, primate=primate, assemble=assemble)[a.step](out)


if __name__ == '__main__':
    main()

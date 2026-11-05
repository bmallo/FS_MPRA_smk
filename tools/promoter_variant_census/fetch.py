"""Retrieve the raw inputs of the promoter variant census (settings: rules.json, windows: windows.tsv).

Subcommands, each writing under --output (one network request at a time, klone policy):
  gnomad      gnomAD v4.1 GraphQL region query per window (padded 50 bp) -> gnomad/<window>.json
  gtex        GTEx Portal API v2: significant pairs per target gene (V10, V8), V10/V8 variant-table lookup and
              per-variant eQTLs (any gene) for the gnomAD variants, V8 fine-mapping per gene -> gtex/*.tsv
  promoterai  rows of the precomputed PromoterAI file inside the windows -> promoterai/rows.tsv
  clinvar     NCBI E-utilities esearch '<chr>[chr] AND <start>:<end>[chrpos38]' + esummary per window
              -> clinvar/<window>.json (truth_set_rules.json)
Usage: python tools/promoter_variant_census/fetch.py {gnomad,gtex,promoterai,clinvar} --output <dir> [--scores <file>]
"""
import argparse
import json
import subprocess
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
PAD = 50
LOOKUP_COLS = ['variant_id', 'gtex_v10_present', 'gtex_v10_maf01', 'gtex_v8_present', 'gtex_v8_maf01', 'v10_sig_genes',
               'v10_sig_pairs']                  # gtex/variant_lookup.tsv column order (appends must match)
GNOMAD = 'https://gnomad.broadinstitute.org/api'
GTEX = 'https://gtexportal.org/api/v2'
QUERY = """{ region(chrom: "%s", start: %d, stop: %d, reference_genome: GRCh38) {
  variants(dataset: gnomad_r4) { variant_id pos ref alt rsids
    genome { ac an af filters homozygote_count populations { id ac an } }
    exome { ac an af filters } } } }"""


def windows():
    return pd.read_csv(HERE / 'windows.tsv', sep='\t')


def get(url, data=None, tries=4):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json'})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())
        except Exception as e:                       # transient API errors: back off and retry
            print(f'{url[:90]}: attempt {attempt + 1} failed: {type(e).__name__} {e}', flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f'failed {tries} times: {url}')


def gnomad(out, only=None):
    d = out / 'gnomad'
    d.mkdir(parents=True, exist_ok=True)
    for w in windows().itertuples():
        if only and w.window not in only:
            continue
        q = QUERY % (w.chrom.removeprefix('chr'), w.start - PAD, w.end + PAD)
        r = get(GNOMAD, json.dumps({'query': q}).encode())
        if r.get('errors'):
            raise RuntimeError(r['errors'])
        json.dump(r, open(d / f'{w.window}.json', 'w'))
        print(w.window, len(r['data']['region']['variants']), 'variant records (padded region)', flush=True)
        time.sleep(1)


def gtex_paged(path, **params):
    rows, page = [], 0
    while True:
        r = get(f'{GTEX}/{path}?' + urllib.parse.urlencode({**params, 'page': page, 'itemsPerPage': 250}))
        rows += r['data']
        page += 1
        if page >= r['paging_info']['numberOfPages']:
            return rows
        time.sleep(0.5)


def gtex(out):
    d = out / 'gtex'
    d.mkdir(parents=True, exist_ok=True)
    # target-gene GENCODE ids per release (V10 = GENCODE v39, V8 = v26); gene-level tables are fetched once
    have = set()
    if (d / 'significant_target_gene.tsv').exists() and (d / 'finemapping_v8_target_gene.tsv').exists():
        have = set(pd.read_csv(d / 'significant_target_gene.tsv', sep='\t').target_gene)
    if not set(windows().gene) <= have:
        genes = {}
        for g in sorted(set(windows().gene)):
            for ds, gv in (('gtex_v10', 'v39'), ('gtex_v8', 'v26')):
                r = get(f'{GTEX}/reference/gene?' + urllib.parse.urlencode({'geneId': g, 'gencodeVersion': gv,
                                                                           'genomeBuild': 'GRCh38/hg38'}))
                genes[(g, ds)] = [x['gencodeId'] for x in r['data'] if x['geneSymbol'] == g][0]
                time.sleep(0.5)
        json.dump({f'{g}|{ds}': v for (g, ds), v in genes.items()}, open(d / 'gencode_ids.json', 'w'), indent=1)
        sig, fm = [], []
        for (g, ds), gid in genes.items():
            sig += [dict(x, target_gene=g) for x in gtex_paged('association/singleTissueEqtl', gencodeId=gid, datasetId=ds)]
            if ds == 'gtex_v8':
                fm += [dict(x, target_gene=g) for x in gtex_paged('association/fineMapping', gencodeId=gid, datasetId=ds)]
            print(g, ds, gid, 'significant rows so far', len(sig), flush=True)
        pd.DataFrame(sig).to_csv(d / 'significant_target_gene.tsv', sep='\t', index=False)
        pd.DataFrame(fm).to_csv(d / 'finemapping_v8_target_gene.tsv', sep='\t', index=False)
    # variant table lookups for every gnomAD record with genome AC >= 1 in the padded regions
    ids = set()
    for f in sorted((out / 'gnomad').glob('*.json')):
        for v in json.load(open(f))['data']['region']['variants']:
            if v['genome'] and v['genome']['ac'] > 0:
                c, p, ref, alt = v['variant_id'].split('-')
                ids.add(f'chr{c}_{p}_{ref}_{alt}_b38')
    path = d / 'variant_lookup.tsv'
    done = set(pd.read_csv(path, sep='\t').variant_id) if path.exists() else set()
    todo = sorted(ids - done)
    print(len(ids), 'variants;', len(todo), 'to look up', flush=True)
    t0 = time.time()
    for i, vid in enumerate(todo):
        row = dict(variant_id=vid, v10_sig_genes='', v10_sig_pairs=0)
        for ds in ('gtex_v10', 'gtex_v8'):
            r = get(f'{GTEX}/dataset/variant?' + urllib.parse.urlencode({'variantId': vid, 'datasetId': ds}))['data']
            row[f'{ds}_present'] = bool(r)
            row[f'{ds}_maf01'] = bool(r and r[0].get('maf01'))
            time.sleep(0.3)
        if row['gtex_v10_present']:                     # eQTLs for any gene; only variants in the V10 table can have one
            e = get(f'{GTEX}/association/singleTissueEqtl?' + urllib.parse.urlencode(
                {'variantId': vid, 'datasetId': 'gtex_v10', 'itemsPerPage': 250}))['data']
            row['v10_sig_genes'] = ';'.join(sorted({x['geneSymbol'] for x in e}))
            row['v10_sig_pairs'] = len(e)
            time.sleep(0.3)
        pd.DataFrame([row])[LOOKUP_COLS].to_csv(path, sep='\t', index=False, mode='a', header=not path.exists())
        if i % 50 == 0:
            print(f'{i + 1}/{len(todo)} {time.time() - t0:.0f} s', flush=True)


def promoterai(out, scores):
    d = out / 'promoterai'
    d.mkdir(parents=True, exist_ok=True)
    w = windows()
    cond = ' || '.join(f'($1=="{r.chrom}" && $2>={r.start} && $2<={r.end})' for r in w.itertuples())
    cmd = f"zcat {scores} | awk -F'\\t' 'NR==1 || {cond}' > {d / 'rows.tsv'}"
    subprocess.run(cmd, shell=True, check=True)
    t = pd.read_csv(d / 'rows.tsv', sep='\t')
    print(len(t), 'rows;', t.groupby('gene').size().to_dict(), flush=True)


def clinvar(out, only=None):
    d = out / 'clinvar'
    d.mkdir(parents=True, exist_ok=True)
    eu = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
    for w in windows().itertuples():
        if only and w.window not in only:
            continue
        term = f'{w.chrom.removeprefix("chr")}[chr] AND {w.start}:{w.end}[chrpos38]'
        ids = get(f'{eu}/esearch.fcgi?' + urllib.parse.urlencode({'db': 'clinvar', 'term': term, 'retmax': 5000,
                                                                   'retmode': 'json'}))['esearchresult']['idlist']
        res = {'uids': []}
        for i in range(0, len(ids), 150):
            time.sleep(0.5)
            r = get(f'{eu}/esummary.fcgi?' + urllib.parse.urlencode({'db': 'clinvar', 'id': ','.join(ids[i:i + 150]),
                                                                      'retmode': 'json'}))['result']
            res['uids'] += r.pop('uids')
            res.update(r)
        json.dump(dict(term=term, result=res), open(d / f'{w.window}.json', 'w'))
        print(w.window, term, len(ids), 'records', flush=True)
        time.sleep(0.5)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('step', choices=['gnomad', 'gtex', 'promoterai', 'clinvar'])
    ap.add_argument('--output', required=True)
    ap.add_argument('--windows', nargs='*', help='gnomad/clinvar: only these windows (default all)')
    ap.add_argument('--scores', default='/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/promoterai/'
                                        'promoterAI_tss500.tsv.gz')
    a = ap.parse_args()
    out = Path(a.output)
    if a.step == 'gnomad':
        gnomad(out, a.windows)
    elif a.step == 'gtex':
        gtex(out)
    elif a.step == 'clinvar':
        clinvar(out, a.windows)
    else:
        promoterai(out, a.scores)


if __name__ == '__main__':
    main()

"""Step A1: inventory of every sequencing run (provenance, chemistry, software, read counts).

Reads the run registry (config/analysis_v2/runs.json), then for each run:
  - @RG: movie, chemistry (CM), barcode, sample; movie date from the movie name
  - @PG: CCS, Jasmine (and whether it called 6mA), fibertools (subcommands), pbmm2 versions,
    and the m6A caller actually used (jasmine --6ma / on-instrument vs ft predict-m6a)
  - read counts with samtools: raw HiFi reads (primary records of the raw BAM, or of the aligned BAM
    when it keeps unmapped reads), primary plasmid alignments and supplementary plasmid records
Writes inventory.tsv, inventory.json and differences.tsv (which attributes differ for each run pair).

Usage: python tools/analysis_v2/a1_inventory.py --runs config/analysis_v2/runs.json \
           --output results/analysis_v2_20261028/a1_inventory --threads 3
"""
import argparse
import datetime as dt
import itertools
import json
import re
import subprocess
from pathlib import Path

import pandas as pd


def sh(*cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def header(bam):
    rg, pg = {}, []
    for line in sh('samtools', 'view', '-H', bam).splitlines():
        f = dict(x.split(':', 1) for x in line.split('\t')[1:] if ':' in x)
        if line.startswith('@RG') and not rg:
            rg = f
        elif line.startswith('@PG'):
            pg.append(f)
    return rg, pg


def software(pg):
    out = {}
    for p in pg:
        pid, vn, cl = p.get('ID', ''), p.get('VN', '').split(' ')[0], p.get('CL', '')
        name = pid.split('.')[0]
        if name in ('ccs', 'jasmine', 'pbmm2', 'lima'):
            out[name] = vn
        if name == 'jasmine':
            out['jasmine_6ma_flag'] = '--6ma' in cl
            out['jasmine_desc'] = p.get('DS', '')
        if name == 'ft':
            out['fibertools'] = vn
            sub = re.search(r'\bft\s+([a-z0-9-]+)', cl)
            out.setdefault('ft_steps', []).append(sub.group(1) if sub else '?')
    return out


def m6a_caller(sw, chem):
    if 'predict-m6a' in sw.get('ft_steps', []):
        return f"fibertools predict-m6a {sw.get('fibertools')}"
    # SPRQ chemistry (P2-C3 here) calls 6mA on instrument without an explicit flag
    if sw.get('jasmine_6ma_flag') or 'C3' in chem:
        return f"Jasmine {sw.get('jasmine')} (on instrument)"
    return 'unknown'


def count(bam, flags, region=None, threads=1):
    cmd = ['samtools', 'view', '-c', '-@', str(threads), *flags, bam] + ([region] if region else [])
    return int(sh(*cmd).strip())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--runs', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--threads', type=int, default=1)
    a = ap.parse_args()
    reg = json.load(open(a.runs))
    root = Path(reg['project_root'])
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for r in reg['runs']:
        aligned = str(root / r['bam_aligned'])
        rg, pg = header(aligned)
        sw = software(pg)
        movie = rg.get('PU', '')
        m = re.match(r'm\d+_(\d{6})_', movie)
        date = dt.datetime.strptime(m.group(1), '%y%m%d').date().isoformat() if m else None
        chem = rg.get('CM', '')
        if r['bam_raw']:
            raw = count(str(root / r['bam_raw']), ['-F', '0x900'], threads=a.threads)
        elif r['aligned_includes_unmapped']:
            raw = count(aligned, ['-F', '0x900'], threads=a.threads)
        else:
            raw = None
        prim = count(aligned, ['-F', '0x904'], r['plasmid_contig'], a.threads)
        supp = count(aligned, ['-f', '0x800'], r['plasmid_contig'], a.threads)
        rows.append(dict(run=r['id'], cell=r['cell'], movie=movie, sequenced=date, chemistry=chem,
                         reaction_min=r['reaction_min'], hia5_ul=r['hia5_ul'], cells_per_reaction=r['cells_per_reaction'],
                         construct=r['construct'], ccs=sw.get('ccs'), jasmine=sw.get('jasmine'), fibertools=sw.get('fibertools'),
                         ft_steps='+'.join(sw.get('ft_steps', [])), pbmm2=sw.get('pbmm2'), m6a_caller=m6a_caller(sw, chem),
                         raw_hifi_reads=raw, plasmid_primary=prim, plasmid_supplementary=supp,
                         supplementary_per_primary=round(supp / prim, 4) if prim else None,
                         reaction_source=r['reaction_source'], bam_aligned=r['bam_aligned']))
        print(rows[-1]['run'], raw, prim, supp, flush=True)
    inv = pd.DataFrame(rows)
    inv.to_csv(out / 'inventory.tsv', sep='\t', index=False)
    attrs = ['cell', 'construct', 'reaction_min', 'hia5_ul', 'chemistry', 'ccs', 'm6a_caller', 'fibertools', 'sequenced']
    diff = []
    for x, y in itertools.combinations(rows, 2):
        d = [k for k in attrs if x[k] != y[k]]
        diff.append(dict(run_a=x['run'], run_b=y['run'], n_differences=len(d), differs_in=', '.join(d)))
    pd.DataFrame(diff).to_csv(out / 'differences.tsv', sep='\t', index=False)
    json.dump(dict(runs=rows, differences=diff, activity=reg['activity']), open(out / 'inventory.json', 'w'), indent=1, default=str)
    print(inv.drop(columns=['reaction_source', 'bam_aligned']).to_string(index=False))
    print(pd.DataFrame(diff).to_string(index=False))


if __name__ == '__main__':
    main()

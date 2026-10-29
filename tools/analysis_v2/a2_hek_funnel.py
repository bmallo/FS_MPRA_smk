"""Step A2 (HEK293T): read funnel from the existing HEK293T selection and extraction audits.

HEK293T was processed by a different pipeline (tools/hek_compare + fsmpra/site_extract_hek.py, construct
with a 318-bp mutable region), so this funnel is reported on its own, not merged with HepG2.
Order: primary plasmid alignments -> passes < 10 -> not full-length -> no m6A tags -> MAPQ -> eligible,
then the eligible molecules by genotype class.

Usage: python tools/analysis_v2/a2_hek_funnel.py --selection results/hek293t_llr_20261028/config/selection.json \
           --audit results/variant_effects_20261028/hek293t/extract/audit.json \
           --output results/analysis_v2_20261028/a2_accounting/HEK293T_run1
"""
import argparse
import json
from pathlib import Path

import pandas as pd

LABELS = {'passes': 'low sequencing quality (passes < 10)', 'not_full_length': 'not full-length', 'no_m6a_tags': 'no m6A data',
          'mapq': 'low MAPQ', 'promoter_indel': 'indel in the mutable region', 'low_quality': 'element base quality too low',
          'multi': 'two or more SNVs', 'flank_snv': 'single SNV outside the 211-bp element', 'shared_snv': 'eligible single SNV (in HepG2 design)',
          'WT': 'eligible WT'}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--selection', required=True)
    ap.add_argument('--audit', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    f = json.load(open(a.selection))['funnel']
    g = json.load(open(a.audit))['genotype_counts']
    rows = [dict(step=LABELS[k], molecules=f['status'][k]) for k in ('not_full_length', 'passes', 'mapq', 'no_m6a_tags')]
    rows += [dict(step=LABELS[k], molecules=g[k]) for k in ('promoter_indel', 'low_quality', 'multi', 'flank_snv', 'shared_snv', 'WT')]
    t = pd.DataFrame(rows)
    assert t.molecules.sum() == f['primary_records'], (t.molecules.sum(), f['primary_records'])
    t['fraction_of_plasmid'] = t.molecules / f['primary_records']
    t.insert(0, 'run', 'HEK293T_run1')
    t.to_csv(out / 'funnel.tsv', sep='\t', index=False)
    json.dump(dict(plasmid_primary=f['primary_records'], eligible_status=f['eligible'], eligible_wt=g['WT'], eligible_shared_snv=g['shared_snv'],
                   note='HiFi read total not available (aligned BAM keeps plasmid-mapped reads only); different pipeline from HepG2.'),
              open(out / 'summary.json', 'w'), indent=1)
    print(t.to_string(index=False))


if __name__ == '__main__':
    main()

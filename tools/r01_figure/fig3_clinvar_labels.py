"""Figure 3 panel a with every ClinVar-classified SNV labelled by its ClinVar record, for mapping to ClinVar inputs
(reference panel, not a figure draft).

Panel a data (fig_A_v9/panel_a.tsv: HepG2 MPRA log2 effect and ClinVar class per designed SNV) joined to the ClinVar
esummary records (clinical_raw/clinvar/esummary.json) by canonical SPDI (NC_000019.10, 0-based position; single-base
substitutions only). Labels give the HGVS c. (NM_000527.5) and the VCV accession; a TSV gives the full ClinVar title,
VCV accession, classification and review status. ClinVar titles come as NM_000527.4/.5 c., NC_000019.10 g. or NR_
names; any classified SNV without a record, an NM title whose c. differs from the panel c., or a g. title whose
position differs from hg38 is reported.
Usage: python tools/r01_figure/fig3_clinvar_labels.py --work results/r01_figure_20261102 --output <dir>
"""
import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HPRC = Path('/gscratch/stergachislab/bmallo/large_home/projects/HPRC_driver_variants')
for f in ('Arimo-Regular.ttf', 'Arimo-Bold.ttf'):
    font_manager.fontManager.addfont(str(HPRC / 'scripts/fonts' / f))
plt.rcParams.update({'font.family': 'Arimo', 'svg.fonttype': 'none', 'font.size': 6.5})
sys.path.insert(0, str(Path(__file__).parent))
from fig3 import CLS, C0, INK, INK2  # noqa: E402

HG38_OFF = 11089264                    # hg38 = oligo + HG38_OFF (evidence_v3)
ORDER = ('vus', 'conf', 'blb', 'plp')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--work', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    w, out = Path(a.work), Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    pa = pd.read_csv(w / 'fig_A_v9/panel_a.tsv', sep='\t')
    pa['hg38'] = pa.oligo + HG38_OFF
    pa['alt'] = pa.snv.str[-1]
    cv = json.load(open(w / 'clinical_raw/clinvar/esummary.json'))['result']
    recs = []
    for uid in cv['uids']:
        r = cv[uid]
        spdi = r['variation_set'][0].get('canonical_spdi', '')
        m = re.fullmatch(r'NC_000019\.10:(\d+):([ACGT]):([ACGT])', spdi)
        if not m:
            continue
        gc = r.get('germline_classification', {})
        recs.append(dict(hg38=int(m[1]) + 1, ref=m[2], alt=m[3], clinvar_title=r['title'], vcv=r['accession'],
                         clinvar_classification=gc.get('description', ''), review_status=gc.get('review_status', '')))
    R = pd.DataFrame(recs)
    t = pa[pa.cls != 'none'].merge(R, on=['hg38', 'alt'], how='left')
    t['c_panel'] = 'c.' + t.c_dot.astype(str) + t.snv.str[-3:]
    nm = t.clinvar_title.str.startswith('NM_000527', na=False)
    gt = t.clinvar_title.str.startswith('NC_000019.10:g.', na=False)
    c_title = t.clinvar_title.str.extract(r'(c\.[^ ]+?)(?: \(|$)')[0]
    g_pos = pd.to_numeric(t.clinvar_title.str.extract(r'g\.(\d+)')[0], errors='coerce')
    miss = t[t.vcv.isna()]
    diff = t[(nm & (c_title != t.c_panel)) | (gt & (g_pos != t.hg38))]
    print(f'classified SNVs {len(t)}; with a ClinVar record {t.vcv.notna().sum()}; without {len(miss)}; titles: '
          f'NM {nm.sum()}, NC g. {gt.sum()}, other {(~nm & ~gt).sum()}; position disagreements {len(diff)}')
    if len(diff):
        print(diff[['snv', 'c_panel', 'clinvar_title']].to_string(index=False))
    assert t.groupby('snv').size().max() == 1, 'more than one ClinVar record per SNV'
    t = t.sort_values(['oligo', 'snv'])
    cols = ['snv', 'oligo', 'hg38', 'c_panel', 'clinvar_title', 'vcv', 'clinvar_classification', 'review_status',
            'cls', 'act']
    t[cols].rename(columns={'c_panel': 'hgvs_c_NM_000527.5', 'act': 'mpra_log2_hepg2', 'cls': 'figure_class'}).to_csv(
        out / 'clinvar_variants.tsv', sep='\t', index=False)

    # plot: all SNVs as in panel a; labels for the classified ones spread evenly above, with leader lines
    fig = plt.figure(figsize=(17, 6.2))
    ax = fig.add_axes([0.05, 0.12, 0.93, 0.48])
    for k in ('none',) + ORDER:
        m = pa.cls == k
        lab, kw = CLS[k]
        ax.scatter(pa.c_dot[m], pa.act[m], label=f'{lab} ({m.sum()})', **kw)
    ax.axhline(0, color=INK2, lw=0.5, zorder=1)
    xlo, xhi = pa.c_dot.min() - 2, pa.c_dot.max() + 2
    ax.set_xlim(xlo, xhi)
    ylo, yhi = pa.act.min() - 0.2, pa.act.max() + 0.2
    ax.set_ylim(ylo, yhi)
    ax.set_xticks(range(-260, -50, 10))
    ax.set_xlabel('Position relative to the start codon (HGVS c., NM_000527.5)')
    ax.set_ylabel('Expression change (log2, HepG2 MPRA)')
    ax.spines[['top', 'right']].set_visible(False)
    n = len(t)
    lx = np.linspace(xlo + 1, xhi - 1, n)
    for (_, r), x_lab in zip(t.iterrows(), lx):
        col = CLS[r.cls][1]['color']
        lab = f'{r.c_panel}  {r.vcv}' if isinstance(r.vcv, str) else r.c_panel + ' (no record)'
        ax.plot([r.c_dot, r.c_dot, x_lab], [r.act, yhi + 0.05, yhi + 0.3], color=col, lw=0.35, alpha=0.7,
                clip_on=False, zorder=1)
        ax.text(x_lab, yhi + 0.33, lab, rotation=90, ha='center', va='bottom', fontsize=5.2,
                color=INK if r.cls == 'vus' else col, fontweight='normal' if r.cls == 'vus' else 'bold', clip_on=False)
    hs = [Line2D([], [], ls='none', marker=CLS[k][1].get('marker', 'o'), ms={'none': 2.5, 'vus': 3}.get(k, 4.5),
                 color=CLS[k][1]['color'], mec='white', mew=0.3) for k in ('none',) + ORDER]
    ax.legend(hs, [f'{CLS[k][0]} ({(pa.cls == k).sum()})' for k in ('none',) + ORDER], loc='upper center',
              bbox_to_anchor=(0.5, -0.13), ncol=5, frameon=False, fontsize=7)
    fig.text(0.05, 0.985, f'LDLR promoter (chr19:{pa.hg38.min():,}-{pa.hg38.max():,}, hg38): {n} ClinVar-classified '
             'SNVs labelled with HGVS c. (NM_000527.5) and ClinVar VCV accession (full ClinVar titles in clinvar_variants.tsv). '
             'Bold = P/LP, B/LB or conflicting; plain = VUS.', fontsize=7, va='top')
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(out / f'panel_a_clinvar_labels.{ext}', dpi=250 if ext == 'png' else None)


if __name__ == '__main__':
    main()

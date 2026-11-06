"""Figure 2b source (draft 10+): read the HEK293T luciferase points and box geometry from the vector Fig. S15A
(Supplemental.pdf page 27), since the raw firefly/renilla values are not on klone (fig2_rules.json revision_draft10).

The S15A HEK293T half (PDF x 56-196, y 68-284 pt, y down) holds 8 WT and 8 c.-169A>G points as filled circles, one box
(rect), median line and two whiskers per group. y is calibrated on the 0.0 and 2.0 tick marks. Checks: 8 points per
group; WT median 1.000 (S15A is normalised to WT); the median of the points equals the drawn median; points outside
the whiskers are reported (drawn as outliers).
Needs PyMuPDF (installed outside the conda envs; run with PYTHONPATH pointing at it).
Writes points.tsv (group, value, x_offset = jitter in box widths) and box.tsv (group, q1, median, q3, whisker_lo,
whisker_hi).
Usage: PYTHONPATH=<pymupdf dir> python tools/r01_figure/fig2_luc_digitize.py --pdf <Supplemental.pdf> --output <dir>
"""
import argparse
import csv
from pathlib import Path

import pymupdf

CROP = (56, 68, 196, 284)                     # S15A HEK293T half, PDF pt
TICKS = {0.0: 252.829, 2.0: 89.824}           # y tick marks (PDF pt) -> value; checked against the tick labels below
COLS = {'WT': (0.182, 0.190, 0.573), 'Var': (0.930, 0.111, 0.141)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pdf', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    pg = pymupdf.open(a.pdf)[26]
    x0, y0, x1, y1 = CROP
    dr = [d for d in pg.get_drawings() if x0 <= d['rect'].x0 and d['rect'].x1 <= x1 and y0 <= d['rect'].y0
          and d['rect'].y1 <= y1]
    tick_y = sorted(round(d['rect'].y0, 3) for d in dr if d['rect'].height == 0 and 3.7 < d['rect'].width < 3.9)
    assert tick_y[0] == TICKS[2.0] and tick_y[-1] == TICKS[0.0] and len(tick_y) == 5, tick_y
    labels = sorted((s['bbox'][1], s['text']) for b in pg.get_text('dict')['blocks'] for ln in b.get('lines', [])
                    for s in ln['spans'] if x0 <= s['bbox'][0] <= 100 and y0 <= s['bbox'][1] <= y1
                    and s['text'] in ('0.0', '0.5', '1.0', '1.5', '2.0'))
    assert [t for _, t in labels] == ['2.0', '1.5', '1.0', '0.5', '0.0'], labels
    per = (TICKS[0.0] - TICKS[2.0]) / 2.0                          # pt per unit
    val = lambda y: (TICKS[0.0] - y) / per                          # noqa: E731
    group = lambda c: next((g for g, rgb in COLS.items() if c and all(abs(u - v) < 0.01 for u, v in zip(c, rgb))),
                           None)                                    # noqa: E731
    box, pts = {}, []
    for g in COLS:
        rect = [d for d in dr if group(d.get('fill')) == g and ''.join(i[0] for i in d['items']) == 're']
        assert len(rect) == 1
        r = rect[0]['rect']
        med = [d for d in dr if d['rect'].height == 0 and abs(d['rect'].x0 - r.x0) < 0.01 and d.get('width') == 1.5]
        caps = sorted(d['rect'].y0 for d in dr if group(d.get('color')) == g and d['rect'].height == 0)
        assert len(med) == 1 and len(caps) == 2
        cx, bw = (r.x0 + r.x1) / 2, r.x1 - r.x0
        box[g] = dict(group=g, q1=val(r.y1), median=val(med[0]['rect'].y0), q3=val(r.y0), whisker_lo=val(caps[1]),
                      whisker_hi=val(caps[0]))
        circ = [d['rect'] for d in dr if group(d.get('fill')) == g and ''.join(i[0] for i in d['items']) == 'cccccccc']
        assert len(circ) == 8, (g, len(circ))
        pts += [dict(group=g, value=val((c.y0 + c.y1) / 2), x_offset=((c.x0 + c.x1) / 2 - cx) / bw) for c in circ]
    for g in COLS:
        v = sorted(p['value'] for p in pts if p['group'] == g)
        pmed = (v[3] + v[4]) / 2
        out_w = [x for x in v if not box[g]['whisker_lo'] - 0.01 <= x <= box[g]['whisker_hi'] + 0.01]
        print(f"{g}: n={len(v)} values={[round(x, 3) for x in v]} point median={pmed:.3f} drawn median="
              f"{box[g]['median']:.3f} mean={sum(v) / len(v):.3f} outside whiskers={[round(x, 3) for x in out_w]}")
        assert abs(pmed - box[g]['median']) < 0.01, 'point median differs from the drawn median'
    assert abs(box['WT']['median'] - 1.0) < 0.002, box['WT']['median']
    for name, rows in (('points.tsv', pts), ('box.tsv', list(box.values()))):
        with open(out / name, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
            w.writeheader()
            w.writerows({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()} for r in rows)


if __name__ == '__main__':
    main()

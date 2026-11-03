"""Figure 2b source: export Supplemental.pdf page 27 (Figure S15, vector: luciferase A + footprint B) as SVG with live
text, and list text positions so the HEK293T luciferase half can be clipped out by fig2.py.

Needs PyMuPDF (installed outside the conda envs; run with PYTHONPATH pointing at it, in its own process).
Writes S15_page.svg and S15_text.tsv (text spans with bounding boxes in PDF points, y down).
Usage: PYTHONPATH=<pymupdf dir> python tools/r01_figure/fig2_s15_svg.py --pdf <Supplemental.pdf> --page 27 --output <dir>
"""
import argparse
import csv
import re
from pathlib import Path

import pymupdf

FONT_SIZE = re.compile(r'font-size="([0-9.]+)"')


def split_runs(svg):
    """PyMuPDF writes one tspan per text run with a per-glyph x list, merging separate labels (e.g. the tick labels
    '2.01.51.0'). Renderers that honour only the first x (rsvg) then misplace them. Split each tspan into word-level
    tspans wherever the next glyph does not follow on (goes back, or jumps > 1.2 x font size), each with one x; the
    glyph advances inside a word come from the font (Arial / metric-compatible Arimo). Subset font names
    (ABCDEF+ArialMT) become Arial."""
    def fix_text(m):
        el = m.group(0)
        size = float(FONT_SIZE.search(el).group(1))

        def fix_tspan(t):
            y, xs, txt = t.group(1), [float(v) for v in t.group(2).split()], t.group(3)
            if len(xs) != len(txt):
                return t.group(0)
            parts, start = [], 0
            for i in range(1, len(txt) + 1):
                if i == len(txt) or xs[i] < xs[i - 1] or xs[i] - xs[i - 1] > 1.2 * size:
                    seg = txt[start:i]
                    if seg.strip():
                        lead = len(seg) - len(seg.lstrip())
                        parts.append(f'<tspan y="{y}" x="{xs[start + lead]}">{seg.strip()}</tspan>')
                    start = i
            return ''.join(parts)
        return re.sub(r'<tspan y="([^"]+)" x="([^"]+)">([^<]*)</tspan>', fix_tspan, el)
    svg = re.sub(r'<text[^>]*>.*?</text>', fix_text, svg, flags=re.S)
    return re.sub(r'font-family="[A-Z]{6}\+Arial[^"]*"', 'font-family="Arial, Helvetica, sans-serif"', svg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pdf', required=True)
    ap.add_argument('--page', type=int, required=True, help='1-based page number')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    page = pymupdf.open(a.pdf)[a.page - 1]
    raw = page.get_svg_image(text_as_path=False)
    (out / 'S15_page.svg').write_text(raw)
    (out / 'S15_page_runs.svg').write_text(split_runs(raw))
    with open(out / 'S15_text.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['text', 'x0', 'y0', 'x1', 'y1', 'size', 'font'])
        for b in page.get_text('dict')['blocks']:
            for ln in b.get('lines', []):
                for sp in ln['spans']:
                    w.writerow([sp['text'], *[round(v, 2) for v in sp['bbox']], round(sp['size'], 2), sp['font']])
    print('page size (pt):', page.rect)


if __name__ == '__main__':
    main()

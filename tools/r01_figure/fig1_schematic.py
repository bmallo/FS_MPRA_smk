"""Figure 1a: the plasmid Fiber-seq workflow row from the HPRC main figure (vector), adapted for the grant.

Source: projects/HPRC_driver_variants/figure/main/panels_ab_paste_in_place.svg (panel a, 'Experimental Outline' row:
plasmid library, transfection, Hia5-marked plasmid, linearize, PacBio). Changes:
  - the row heading and the arrow into the HPRC read-analysis step are removed; labels: '5,204 alleles' -> 'reporter
    library' (second line of 'Plasmid library' replaced), 'HEK293T' -> 'cells';
  - the variant-allele plasmids (lib-plasmid / lib-insert / lib-variant groups) are replaced by reporter plasmids drawn
    here: grey backbone, a coloured regulatory element (colour = a different element) and the reporter gene as a dark
    arrow, at the same seven positions.
Writes schematic.svg: the full source page with these edits, viewBox set to CROP (Illustrator-editable vector), and
schematic.png (600 dpi preview). fig1.py embeds schematic.svg's content clipped to CROP.
Usage: python tools/r01_figure/fig1_schematic.py --output <dir>
"""
import argparse
import math
import re
import subprocess
from pathlib import Path

SRC = Path('/gscratch/stergachislab/bmallo/large_home/projects/HPRC_driver_variants/figure/main/panels_ab_paste_in_place.svg')
RSVG = '/mmfs1/home/bmallo/miniconda3/envs/fiberseq-qc/bin/rsvg-convert'
OX, OY = 46.63 + 2.0, 46.63               # panel-a origin on the page (translates in the source)
EX, EY = 10, 82                           # experiment box in panel coordinates (s46_panel_a_mockup.py)
CROP = (OX + EX + 1, OY + EY + 1.5, 185.0, 54)  # x, y, w, h (page pt): library .. PacBio, inside the old frame
REPLACE = {'>Plasmid library<': '>Reporter<', '>5,204 alleles<': '>plasmid library<', '>HEK293T<': '>cells<',
           '>Experimental Outline<': '><'}
# library layout as in s46 (panel coordinates): hexagonal cluster centred (EX+25, EY+25.5), spacing 9.6, radius 4.1
LC, SP, R = (EX + 25, EY + 25.5), 9.6, 4.1
CENTRES = [(LC[0] - SP / 2, LC[1] - 0.87 * SP), (LC[0] + SP / 2, LC[1] - 0.87 * SP), (LC[0] - SP, LC[1]), LC,
           (LC[0] + SP, LC[1]), (LC[0] - SP / 2, LC[1] + 0.87 * SP), (LC[0] + SP / 2, LC[1] + 0.87 * SP)]
ELEMENT_COLOURS = ['#E77C34', '#338B8C', '#755B9D', '#D0312D', '#2A78D6', '#1BAF7A', '#EDA100']
BACKBONE, REPORTER = '#8B9094', '#3A3A3A'


def arc(cx, cy, r, a0, a1):
    """SVG path for a circular arc from angle a0 to a1 (degrees, clockwise in screen space)."""
    x0, y0 = cx + r * math.cos(math.radians(a0)), cy + r * math.sin(math.radians(a0))
    x1, y1 = cx + r * math.cos(math.radians(a1)), cy + r * math.sin(math.radians(a1))
    large = 1 if (a1 - a0) % 360 > 180 else 0
    return f'M {x0:.3f} {y0:.3f} A {r} {r} 0 {large} 1 {x1:.3f} {y1:.3f}'


def reporter_plasmid(i, cx, cy, a0):
    """Backbone circle, regulatory element (70 deg, coloured) and reporter gene arrow (95 deg, dark) clockwise after it."""
    x, y = cx + OX, cy + OY                      # page coordinates
    col = ELEMENT_COLOURS[i % len(ELEMENT_COLOURS)]
    e0, e1 = a0, a0 + 70
    g0, g1 = e1 + 8, e1 + 8 + 95
    th = math.radians(g1)
    tip = (x + R * math.cos(th), y + R * math.sin(th))
    tx, ty = -math.sin(th), math.cos(th)         # clockwise tangent
    nx, ny = math.cos(th), math.sin(th)
    head = [(tip[0] + 1.6 * tx, tip[1] + 1.6 * ty), (tip[0] + 1.1 * nx, tip[1] + 1.1 * ny),
            (tip[0] - 1.1 * nx, tip[1] - 1.1 * ny)]
    pts = ' '.join(f'{px:.3f},{py:.3f}' for px, py in head)
    return (f'<g id="reporter-plasmid-{i}">'
            f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{R}" fill="#ffffff" stroke="{BACKBONE}" stroke-width="0.7"/>'
            f'<path d="{arc(x, y, R, e0, e1)}" fill="none" stroke="{col}" stroke-width="1.9" stroke-linecap="butt"/>'
            f'<path d="{arc(x, y, R, g0, g1)}" fill="none" stroke="{REPORTER}" stroke-width="1.1" stroke-linecap="butt"/>'
            f'<polygon points="{pts}" fill="{REPORTER}"/></g>')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    s = SRC.read_text()
    for k, v in REPLACE.items():
        assert s.count(k) >= 1, k
        s = s.replace(k, v)
    s, n = re.subn(r'<g id="step-arrow-3">.*?</g>', '', s, count=1, flags=re.S)   # arrow into the HPRC analysis step
    assert n == 1
    removed = 0
    for kind in ('lib-plasmid', 'lib-insert', 'lib-variant', 'lib-bg'):
        s, k = re.subn(rf'<g id="{kind}-[0-9-]+">.*?</g>', '', s, flags=re.S)
        removed += k
    assert removed >= 21, removed
    starts = [200, 320, 30, 250, 140, 290, 80]   # element start angles, varied like the source cluster
    lib = '<g id="reporter-library">' + ''.join(reporter_plasmid(i, cx, cy, a0)
                                                for i, ((cx, cy), a0) in enumerate(zip(CENTRES, starts))) + '</g>'
    s = s.replace('</svg>', lib + '\n</svg>')
    x, y, w, h = CROP
    s = re.sub(r'<svg ([^>]*?)width="[^"]*" height="[^"]*" viewBox="[^"]*"',
               rf'<svg \1width="{w}pt" height="{h}pt" viewBox="{x:.2f} {y:.2f} {w} {h}"', s, count=1)
    (out / 'schematic.svg').write_text(s)
    subprocess.run([RSVG, '-d', '600', '-p', '600', '-f', 'png', '-o', str(out / 'schematic.png'),
                    str(out / 'schematic.svg')], check=True)
    print(f'removed {removed} library groups')


if __name__ == '__main__':
    main()

"""Render the promoter explorer offline with headless Firefox and save screenshots (for checking visuals).

Builds a self-contained test copy of the page in --work: Plotly from the Python plotly package instead of the CDN,
and the data (base.json, chunks/) served from inline scripts instead of fetch. Then runs firefox --headless
--screenshot for each requested variant/hash. Nothing is published; the real page is not modified.

Usage: python tools/analysis_v2/render_page.py --page tools/analysis_v2/explorer.html \
           --data results/analysis_v2_20261028/explorer_v6 --work <scratch dir> --width 1400 --height 5200 \
           [--variant 146C>T] [--hash hek293t] [--dark]
"""
import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path

import plotly

PLOTLY_JS = Path(plotly.__file__).parent / 'package_data' / 'plotly.min.js'
CDN = 'https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.27.0/plotly.min.js'


def build(page, data, work, variant, dark):
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy(PLOTLY_JS, work / 'plotly.min.js')
    files = {'base.json': json.load(open(data / 'base.json'))}
    for p in sorted((data / 'chunks').glob('*.json')):
        files[f'chunks/{p.name}'] = json.load(open(p))
    with open(work / 'data.js', 'w') as f:
        f.write('window.__FILES = ' + json.dumps(files, separators=(',', ':')) + ';')
    shim = """<script src="data.js"></script>
<script>
  window.fetch = url => Promise.resolve(window.__FILES[url]
    ? { ok: true, json: () => Promise.resolve(window.__FILES[url]) } : { ok: false, status: 404 });
  try { localStorage.clear(); %s } catch (e) {}
  %s
</script>
""" % (f"localStorage.setItem('ldlr-explorer-3', JSON.stringify({{sel: {json.dumps(variant)}}}));" if variant else '',
       "document.documentElement.setAttribute('data-theme', 'dark');" if dark else '')
    html = open(page).read().replace(CDN, 'plotly.min.js')
    head, _, rest = html.partition('<style>')
    (work / 'test.html').write_text(head + shim + '<style>' + rest)
    return work / 'test.html'


def shoot(test, out, width, height, hash_):
    prof = test.parent / 'ff-profile'
    prof.mkdir(exist_ok=True)
    url = test.resolve().as_uri() + (f'#{hash_}' if hash_ else '')
    env = dict(os.environ, MOZ_HEADLESS='1', HOME=str(test.parent))
    if out.exists():
        out.unlink()
    # headless Firefox often crashes on exit after writing the screenshot, so the exit code is not checked
    subprocess.run(['/usr/bin/firefox', '--headless', '--no-remote', '--profile', str(prof), '--screenshot', str(out),
                    f'--window-size={width},{height}', url], env=env, timeout=240,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    if not out.exists():
        raise RuntimeError('firefox wrote no screenshot')
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--page', required=True)
    ap.add_argument('--data', required=True)
    ap.add_argument('--work', required=True)
    ap.add_argument('--width', type=int, default=1400)
    ap.add_argument('--height', type=int, default=5200)
    ap.add_argument('--variant')
    ap.add_argument('--hash')
    ap.add_argument('--dark', action='store_true')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    test = build(Path(a.page), Path(a.data), Path(a.work), a.variant, a.dark)
    print(shoot(test, Path(a.out), a.width, a.height, a.hash))


if __name__ == '__main__':
    main()

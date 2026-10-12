#!/usr/bin/env python3
"""Build the single-file page from template.html, solver.js and data/cloud_points.csv.

    python3 build.py                       -> index.html, a complete page for GitHub Pages
    python3 build.py --artifact OUT URL    -> page body only, with the home link pointing to URL
"""
import csv
import json
import re
import sys
from pathlib import Path

here = Path(__file__).resolve().parent


def cloud_points():
    pts = {}
    with open(here / 'data' / 'cloud_points.csv', encoding='utf-8') as f:
        rows = [r for r in f if not r.startswith('#')]
    for r in csv.DictReader(rows):
        pts.setdefault(r['solvent'], []).append([float(r['nonsolvent_wt']), float(r['solvent_wt']), float(r['polymer_wt'])])
    return pts


def build(home, full):
    src = (here / 'solver.js').read_text(encoding='utf-8')
    solver = re.sub(r"\nif \(typeof module !== 'undefined'\) \{.*?\n\}\n", "\n", src, flags=re.S)
    assert solver != src, 'no export block found in solver.js'
    page = (here / 'template.html').read_text(encoding='utf-8')
    page = page.replace('/*SOLVER*/', solver.strip())
    page = page.replace('/*CLOUDPOINTS*/', json.dumps(cloud_points()))
    if not home.startswith('http'):
        page = page.replace('href="/*HOME*/" target="_blank" rel="noopener"', f'href="{home}"')
    page = page.replace('/*HOME*/', home)
    left = re.findall(r'/\*[A-Z]+\*/', page)
    assert not left, f'placeholders left: {left}'
    for ch in ('\u2013', '\u2014'):
        assert ch not in page, f'dash {ch!r} found'
    if full:
        i = page.index('</style>') + len('</style>')
        page = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
                '<meta name="description" content="Flory-Huggins ternary phase diagrams for polymer, solvent and non-solvent systems, with binodal, spinodal and tie-lines.">\n'
                + page[:i].strip() + '\n</head>\n<body>\n' + page[i:].strip() + '\n</body>\n</html>\n')
    return page


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--artifact':
        out, home, full = Path(sys.argv[2]), sys.argv[3], False
    else:
        out, home, full = here / 'index.html', '/', True
    page = build(home, full)
    out.write_text(page, encoding='utf-8')
    print(f'wrote {out} ({len(page)} bytes)')

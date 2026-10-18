#!/usr/bin/env python3
"""Build the single-file page from template.html, solver.js and data/cloud_points.csv.

    python3 build.py                       -> index.html, a complete page for GitHub Pages
    python3 build.py --artifact OUT URL    -> page body only, with the home link pointing to URL
"""
import csv
import html
import json
import re
import sys
from pathlib import Path

here = Path(__file__).resolve().parent

# ---- search engines and link previews (full page only) ----
META = {
    'repo': 'ternary-phase-diagrams',
    'name': 'Ternary phase diagrams',
    'title': 'Flory-Huggins ternary phase diagram calculator | Hashim Ishaque',
    'headline': 'Ternary phase diagrams for polymer, solvent and non-solvent systems',
    'description': ('Ternary phase diagrams for membrane casting: Flory-Huggins binodal, spinodal, tie-lines, '
                    'dope cloud point and coagulation bath mixing line. CSV download.'),
    'image_alt': 'Ternary phase diagram of Torlon, NMP and 2-propanol with the binodal, the spinodal and tie-lines',
    'keywords': ('ternary phase diagram, Flory-Huggins, binodal, spinodal, tie-line, critical point, cloud point, '
                 'phase inversion, membrane casting, polymer solution, casting solution, coagulation bath, Hereon, PAN, PES, PVDF, PEI, CSV'),
    'languages': ['JavaScript', 'Python'],
    'about': [
        ('Flory-Huggins solution theory', 'https://en.wikipedia.org/wiki/Flory%E2%80%93Huggins_solution_theory'),
        ('Ternary plot', 'https://en.wikipedia.org/wiki/Ternary_plot'),
        ('Binodal', 'https://en.wikipedia.org/wiki/Binodal'),
        ('Spinodal decomposition', 'https://en.wikipedia.org/wiki/Spinodal_decomposition'),
        ('Synthetic membrane', 'https://en.wikipedia.org/wiki/Synthetic_membrane'),
    ],
}
SITE = 'https://hashimishaque.com'
PERSON = {
    '@type': 'Person', '@id': SITE + '/#person', 'name': 'Hashim Ishaque', 'url': SITE + '/',
    'sameAs': ['https://www.linkedin.com/in/hashimishaque/', 'https://github.com/hashimishaque-platform'],
}
WEBSITE = {'@type': 'WebSite', '@id': SITE + '/#website', 'url': SITE + '/', 'name': 'Hashim Ishaque'}
MIT = 'https://opensource.org/licenses/MIT'


def head(meta, extra=()):
    """Title, description, canonical address, icons, link preview and structured data (schema.org)."""
    url = f"{SITE}/{meta['repo']}/"
    image = url + 'social.png'
    page = {
        '@type': 'WebPage', '@id': url, 'url': url, 'name': meta['headline'], 'description': meta['description'],
        'inLanguage': 'en', 'isPartOf': {'@id': WEBSITE['@id']}, 'author': {'@id': PERSON['@id']},
        'primaryImageOfPage': {'@type': 'ImageObject', 'url': image, 'width': 1200, 'height': 630},
        'about': [{'@type': 'Thing', 'name': n, 'sameAs': u} for n, u in meta['about']],
        'keywords': meta['keywords'], 'breadcrumb': {'@id': url + '#breadcrumb'},
    }
    code = {
        '@type': 'SoftwareSourceCode', '@id': url + '#code', 'name': meta['name'], 'description': meta['description'],
        'url': url, 'codeRepository': f"https://github.com/hashimishaque-platform/{meta['repo']}",
        'programmingLanguage': meta['languages'], 'runtimePlatform': 'Web browser', 'license': MIT,
        'author': {'@id': PERSON['@id']},
    }
    crumbs = {'@type': 'BreadcrumbList', '@id': url + '#breadcrumb', 'itemListElement': [
        {'@type': 'ListItem', 'position': 1, 'name': 'Hashim Ishaque', 'item': SITE + '/'},
        {'@type': 'ListItem', 'position': 2, 'name': meta['name'], 'item': url}]}
    graph = [page, code, crumbs, WEBSITE, PERSON, *extra]
    ld = json.dumps({'@context': 'https://schema.org', '@graph': graph}, ensure_ascii=False, indent=1).replace('</', '<\\/')
    e = lambda s: html.escape(s, quote=True)
    return '\n'.join([
        f"<title>{e(meta['title'])}</title>",
        f'<meta name="description" content="{e(meta["description"])}">',
        '<meta name="author" content="Hashim Ishaque">',
        f'<link rel="canonical" href="{url}">',
        '<link rel="icon" href="/favicon.ico" sizes="32x32">',
        '<link rel="icon" href="/favicon.svg" type="image/svg+xml">',
        '<link rel="apple-touch-icon" href="/apple-touch-icon.png">',
        '<meta property="og:type" content="website">',
        '<meta property="og:site_name" content="Hashim Ishaque">',
        f'<meta property="og:title" content="{e(meta["headline"])}">',
        f'<meta property="og:description" content="{e(meta["description"])}">',
        f'<meta property="og:url" content="{url}">',
        f'<meta property="og:image" content="{image}">',
        '<meta property="og:image:width" content="1200">',
        '<meta property="og:image:height" content="630">',
        f'<meta property="og:image:alt" content="{e(meta["image_alt"])}">',
        '<meta name="twitter:card" content="summary_large_image">',
        f'<script type="application/ld+json">\n{ld}\n</script>',
    ])


THESIS = {
    '@type': 'Thesis', '@id': 'https://publications.hereon.de/id/38756/', 'url': 'https://publications.hereon.de/id/38756/',
    'name': 'Development of novel polymeric/inorganic composite membranes for the application in gas phase membrane reactors',
    'author': {'@id': PERSON['@id']}, 'datePublished': '2020', 'inSupportOf': "Master's degree",
}
DATASET = {
    '@type': 'Dataset', '@id': f"{SITE}/{META['repo']}/#cloud-points",
    'name': 'Cloud points of Torlon 4000HV in NMP and in DMF with 2-propanol',
    'description': ('Nine cloud points of the polyamide-imide Torlon 4000HV in N-methyl-2-pyrrolidone (NMP) and in '
                    'N,N-dimethylformamide (DMF), with 2-propanol as the non-solvent, given as weight fractions of '
                    'polymer, solvent and non-solvent. Read off Figure 36 of Ishaque (2020); reading error about 0.5 wt %.'),
    'url': f"{SITE}/{META['repo']}/", 'creator': {'@id': PERSON['@id']},
    'isBasedOn': {'@id': THESIS['@id']}, 'citation': {'@id': THESIS['@id']},
    'keywords': ['cloud point', 'Torlon 4000HV', 'polyamide-imide', 'NMP', 'DMF', '2-propanol', 'phase inversion', 'membrane'],
    'variableMeasured': ['polymer weight fraction', 'solvent weight fraction', 'non-solvent weight fraction'],
    'distribution': {'@type': 'DataDownload', 'encodingFormat': 'text/csv',
                     'contentUrl': f"{SITE}/{META['repo']}/data/cloud_points.csv"},
}
EXTRA = (THESIS, DATASET)


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
        page = re.sub(r'<title>.*?</title>\n', '', page, count=1)
        i = page.index('</style>') + len('</style>')
        page = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
                + head(META, EXTRA) + '\n'
                + page[:i].strip() + '\n</head>\n<body>\n' + page[i:].strip() + '\n</body>\n</html>\n')
        for ch in ('\u2013', '\u2014'):
            assert ch not in page, f'dash {ch!r} found'
    return page


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--artifact':
        out, home, full = Path(sys.argv[2]), sys.argv[3], False
    else:
        out, home, full = here / 'index.html', '/', True
    page = build(home, full)
    out.write_text(page, encoding='utf-8')
    print(f'wrote {out} ({len(page)} bytes)')

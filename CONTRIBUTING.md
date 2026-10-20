# Contributing

Corrections, new models and new checks are welcome. Please open an issue first for anything
larger than a small fix.

## Ground rules

- **Every formula and every literature value has a source:** full reference with DOI, and the
  equation, table or page number, in the code comment and in the page's References.
- **Every new part has a check** in `checks.js` against an exact solution or a published value,
  and the script must exit with an error when a check fails.
- **Data.** The only measured data are the nine published cloud points in `data/cloud_points.csv`, cited as Ishaque (2020). Other measured data can be added only if it is published and cited.
- **Keep the style of the page:** plain HTML, Times New Roman, tables with grey borders,
  "Figure N." captions, short sentences, no frameworks and no web fonts, and no horizontal
  scrolling at 400 px width.

## How to work on it

```
node checks.js        # the checks
python3 build.py      # builds index.html from template.html and the code
```

| File | Contents |
|---|---|
| `solver.js` | Flory-Huggins solver: binodal, spinodal, tie-lines, cloud points |
| `template.html` | Text, figures and form of the page |
| `data/cloud_points.csv` | The nine published cloud points |

Do not edit `index.html` by hand. Commit it after `python3 build.py`; the workflow in
`.github/workflows/checks.yml` runs the checks, rebuilds the page and fails if the committed
page is not reproduced.

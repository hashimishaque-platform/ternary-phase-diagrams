# Roadmap

Each step is developed on its own branch and merged into `main` once `node checks.js` passes.
A step is ticked here with the date it was merged.

## Membrane fabrication advisor
Guidance for choosing a dope and a coagulation bath. Every rule is cited, and results are labelled as trends where the model cannot give numbers.
- [x] Casting box: dope and bath compositions, the mixing line between them on the diagram, the cloud point of the dope and the most solvent the bath can hold before the mixing line misses the two-phase region. (2026-10-18)
- [x] Instantaneous or delayed demixing: a qualitative indicator from the dope cloud-point margin and the bath solvent content, with the structure it tends to give (macrovoids, sponge-like sublayer, dense or porous skin), each rule cited (Reuvers and Smolders 1987; Wijmans, Baaij and Smolders 1983; Smolders et al. 1992). (2026-10-18)
- [x] Recipe suggestions: for a target structure, suggest ranges for the polymer content, non-solvent in the dope and solvent in the bath from the diagram of the chosen system, with the reasoning shown. (2026-10-18)
- [x] Porosity estimates only where a published correlation exists for that polymer, solvent and bath, with its source and range of validity. The page says so where none exists. (2026-10-18)
  No predictive correlation was found for the systems on the page, so the page says so and gives the gravimetric formula (Ali et al. 2018) to measure porosity instead.
- [x] Pore-forming additives (PVP, PEG): the published trends and their sources, without a four-component model. (2026-10-18)
- [x] A cited table of the polymers used for membranes at Hereon, with the casting solution and bath from the source and whether the calculator covers each one. (2026-10-18)
- [ ] PEI (Ultem 1000) in the calculator. Hansen parameters alone give chi23 = 0.69 with NMP, which wrongly says NMP does not dissolve Ultem, so this waits for published interaction parameters or cloud points.
- [ ] 1,4-dioxane and GBL as solvents (used at Hereon for CA, PES and PEI), with cited Hansen parameters and chi23 values for each polymer, or a note where none are published.
- [ ] Replace the rough-estimate parameters of PAN, PES, PVDF, CA and the polyimide with published values where they exist (for example Romay, Diban and Urtiaga 2021 for PES and PVDF with NMP, DMAc and water).
- [ ] More polymers (for example PBI, PVC) only where interaction parameters or cloud points are published.
  Checked 2026-10-18: for PVC only one tabulated cloud point (Padilha and Borges 2019) and for PBI only Hansen parameters were found in open sources. Not enough to set three chi values.

## Planner for students
Turn the page into a starting point for choosing a casting recipe. Every value is cited; nothing is invented.
- [ ] Recipe library: a table of published casting recipes (polymer, dope, bath, form, reported structure, porosity or flux, source), about 20 to start for PES, PSf, PVDF and PAN, with filters and a link that shows a recipe on the diagram. Results are filled only from the source.
- [ ] Confidence label on every system (fitted to data, published parameters, rough estimate), shown next to the diagram and on every recipe.
- [ ] Validation page: the model against published cloud points without refitting, one figure and one RMS per system. Systems that miss keep the rough-estimate label.
- [ ] Weigh-out calculator: masses of polymer, additive, solvent and non-solvent for a chosen dope mass, and of solvent and non-solvent for the bath.
- [ ] Lab sheet: one printable page with the recipe, the weighing plan, the expected tendency, the references and blank fields for bench conditions and results.
- [ ] Shareable link: the address holds the system and the settings.
- [ ] One clear step order on the page: system, diagram, target structure, recipe, weigh-out, lab sheet.
- [ ] Uncertainty band on the binodal over the plausible range of chi.
- [ ] Bath temperature and time in air as cited trends, before temperature-dependent chi.
- [ ] Composition-dependent g12 for water systems (see the Composition-dependent section below).
- [ ] Green solvents (Cyrene, PolarClean, GVL) only where published parameters exist.
- [ ] CITATION.cff and a versioned release, so the tool can be cited.

## Export
- [x] Download the binodal, the spinodal and the tie-lines of the current diagram as CSV files, in wt % and in vol %. (2026-10-18)
- [ ] Download the table of measured and calculated cloud points as a CSV file.

## Composition-dependent interaction parameter
- [ ] Add the Altena and Smolders form g12(u2) = α + β / (1 − γ u2) to the solver, with the chemical potentials derived from it, and a check that they match numerical derivatives of the free energy.
- [ ] Use it in the tie-line solver and the spinodal, with a check that β = 0 reproduces the constant-χ results exactly.
- [ ] Form fields for α, β and γ, off by default, and the text in Method and Limitations.

## Temperature
- [ ] Temperature-dependent χ = A + B / T for each pair, with defaults that reproduce the present values at 25 °C.
- [ ] A figure with the cloud-point curve at two temperatures.

## Titration path
- [ ] Draw the cloud-point titration path for a chosen polymer concentration and mark where it meets the binodal.

## Python cross-check
- [ ] `python/ternary.py`: an independent implementation of the binodal with SciPy root finding.
- [ ] `python/compare.py`: compare tie-lines and critical points with the JavaScript solver to 10⁻⁶.

## Robustness
- [ ] Checks for extreme parameters (χ13 up to 3.5, n3 = 50 and 2000, χ12 near 2): the binodal closes, or the page says why not.
- [ ] Move the test composition with the arrow keys.

## Documentation
- [ ] `docs/derivation.md`: the chemical potentials derived step by step from the free energy of mixing.
- [ ] A worked example in the README: using the diagram to choose a casting solution.

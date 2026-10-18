# Roadmap

Each step is developed on its own branch and merged into `main` once `node checks.js` passes.
A step is ticked here with the date it was merged.

## Membrane fabrication advisor
Guidance for choosing a dope and a coagulation bath. Every rule is cited, and results are labelled as trends where the model cannot give numbers.
- [x] Casting box: dope and bath compositions, the mixing line between them on the diagram, the cloud point of the dope and the most solvent the bath can hold before the mixing line misses the two-phase region. (2026-10-18)
- [ ] Instantaneous or delayed demixing: a qualitative indicator from the dope cloud-point margin and the bath solvent content, with the structure it tends to give (macrovoids, sponge-like sublayer, dense or porous skin), each rule cited (Reuvers and Smolders 1987; Wijmans, Baaij and Smolders 1983; Mulder 1996).
- [ ] Recipe suggestions: for a target structure, suggest ranges for the polymer content, non-solvent in the dope and solvent in the bath from the diagram of the chosen system, with the reasoning shown.
- [ ] Porosity estimates only where a published correlation exists for that polymer, solvent and bath, with its source and range of validity. The page says so where none exists.
- [ ] Pore-forming additives (PVP, PEG): the published trends and their sources, without a four-component model.
- [x] A cited table of the polymers used for membranes at Hereon, with the casting solution and bath from the source and whether the calculator covers each one. (2026-10-18)
- [ ] PEI (Ultem 1000) in the calculator. Hansen parameters alone give chi23 = 0.69 with NMP, which wrongly says NMP does not dissolve Ultem, so this waits for published interaction parameters or cloud points.
- [ ] 1,4-dioxane and GBL as solvents (used at Hereon for CA, PES and PEI), with cited Hansen parameters and chi23 values for each polymer, or a note where none are published.
- [ ] Replace the rough-estimate parameters of PAN, PES, PVDF, CA and the polyimide with published values where they exist (for example Romay, Diban and Urtiaga 2021 for PES and PVDF with NMP, DMAc and water).
- [ ] More polymers (for example PBI, PVC) only where interaction parameters or cloud points are published.

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

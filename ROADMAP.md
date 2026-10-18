# Roadmap

Each step is developed on its own branch and merged into `main` once `node checks.js` passes.
A step is ticked here with the date it was merged.

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

# Relative Hodge spectra of multiscale conifold smoothings

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 12 November 2026. This version: 12 November 2026. 41 pages.

[Read the paper](relative-hodge-conifold-smoothings.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/relative-hodge-conifold-smoothings/) · [Zenodo:24175848](https://zenodo.org/records/24175848) · [DOI](https://doi.org/10.5281/zenodo.24175848)

We determine the complete small two-form spectrum of the Ricci-flat metrics on an explicit compact polarized Calabi–Yau conifold smoothing with eight nodes. The estimates hold on its full six-dimensional period base, with no restriction on the relative sizes of the nonzero local smoothing parameters. Exactly two positive eigenvalues tend to zero. A fixed exterior cohomological map identifies their spectral space with the two-dimensional quotient of divisor classes by those surviving the smoothing, and compares its entire energy form relatively with the capacity-weighted restriction form. The harmonic mass and complementary spectral gaps are controlled as well. The proof constructs a uniform Ricci-flat correction through independently scaled Stenzel regions, retaining their quadratic harmonic modes by exact complex automorphisms, and combines weighted capacity bounds with compactness in degrees one and two. The resulting energy and mass forms give uniform asymptotics for both eigenlines in these cohomological coordinates and yield separated logarithmic and iterated-logarithmic spectral scales on the same algebraic family.

## Build and reproduce

Run `latexmk -pdf main.tex` in this directory. The descriptive PDF is the released manuscript; rebuilding writes `main.pdf`.

Figure programs are supplied in `figures/code/` or `scripts/` where applicable. Run a supplied program with `uv run --no-project --with-requirements requirements.txt python <program-path>`, or run `python3 reproduce.py` at the repository root. The integral reconstruction procedure in Paper 31 is proved to terminate; a numerical execution of its full integral matrix is not included. These programs verify finite identities and draw the accompanying figures; the analytic assertions rely on the manuscript's proofs.

[The programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

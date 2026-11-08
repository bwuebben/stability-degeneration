# Light vectors and scalar interactions near heterotic K3 gauge enhancement

By [Bernd Johannes Wuebben](https://bwuebben.github.io/),
[ORCID](https://orcid.org/0009-0002-4588-3284).

[Current version (PDF)](heterotic-k3-enhancement.pdf) ·
[Abstract and citation](https://bwuebben.github.io/papers/heterotic-k3-enhancement/) ·
[Zenodo:24063938](https://zenodo.org/records/24063938) · [DOI](https://doi.org/10.5281/zenodo.24063938)

An SO(32) heterotic compactification on a smooth K3 has three light vector
species near a repeated bundle factor. Their leading squared masses have a
relative geometric estimate independent of their hierarchy. The first tree
vector pencil includes the complete leading metric and two-form towers.
The same family has a nonzero leading neutral-scalar interaction with a
quadratic relative geometric error at fixed angular ratio. The vector bound
concerns the defined first-order pencil; the scalar statement takes the
leading large-volume coefficient before the geometric limit.

## Build and reproduce

With a LaTeX installation, run `latexmk -pdf main.tex` in this directory.
The official JHEP style, bibliography style and generated bibliography are
included. The archived PDF is the checked-in file named above; a local build
produces `main.pdf`.

The four `figures/code/compute*.py` programs use the Python standard library
and reproduce finite matrix, projection, cohomology and curvature calculations.
They do not replace the analytic proofs or the perturbative physical derivation.
The figure program `figures/code/fig01_scales.py` requires NumPy and Matplotlib
and plots exact finite-model formulas, not numerical HE solutions.

[Programme and companion papers](../../README.md) ·
[Physics guide](https://bwuebben.github.io/research/physics/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved.
The JHEP support files retain their own notices. No blanket code license
applies to the manuscript. arXiv submission is pending.

# Light vectors and scalar interactions near heterotic K3 gauge enhancement

By [Bernd Johannes Wuebben](https://bwuebben.github.io/),
[ORCID](https://orcid.org/0009-0002-4588-3284).

[Current version (PDF)](heterotic-k3-enhancement.pdf) ·
[Abstract and citation](https://bwuebben.github.io/papers/heterotic-k3-enhancement/) ·
[Zenodo:24072269](https://zenodo.org/records/24072269) · [DOI](https://doi.org/10.5281/zenodo.24072269)

Original version: 8 November 2026. This version: 8 November 2026.

Two actual massless bundle scalars in a smooth SO(32) heterotic K3 family
have a complete leading Yang–Mills tree amplitude that approaches an explicit
light-vector exchange formula, uniformly as the Higgs scales separate.
The comparison includes the full internal vector and scalar spectrum, below
the heavier Higgs masses and away from light poles. Harmonic orthogonality
controls the vertices and an exact scalar contact–exchange cancellation
removes inverse-mass losses. The same estimates give a selected curvature
coefficient with uniform inverse-square growth.

The defined first-order coupled vector pencil has a ratio-uniform relative
correction bound; a fixed-ray coefficient satisfies a Pauli selection rule.
At the next heterotic tree order, quadratic potential mixing with K3 metric
deformations and universal graviton exchange are controlled separately.
The complete coefficient at that order, finite-string-scale remainders and
quantum effects remain open. The physical curvature statement extracts its
leading large-volume coefficient before applying the geometric estimate.

## Build and reproduce

With a LaTeX installation, run `latexmk -pdf main.tex` in this directory.
The official JHEP support files and generated bibliography are included.
`JHEP-doi.bst` is a local derivative that adds DOI hyperlinks for preprints;
the original `JHEP.bst` is retained. The archived PDF is the checked-in file named above; a local build
produces `main.pdf`.

The six `figures/code/compute*.py` programs reproduce finite matrix, projection,
cohomology, curvature, spectral-amplitude and graviton calculations. The first
four use the Python standard library; the last two require SymPy.
They do not replace the analytic proofs or the perturbative physical derivation.
The two `figures/code/fig*.py` programs require NumPy and Matplotlib and plot
exact finite-model formulas. Pinned dependencies are in `requirements.txt`.
For example, use `uv run --with-requirements requirements.txt python
figures/code/compute06_gravity.py`. These are not numerical HE solutions.

[Programme and companion papers](../../README.md) ·
[Physics guide](https://bwuebben.github.io/research/physics/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved.
The JHEP support files retain their own notices. No blanket code license
applies to the manuscript. arXiv submission is pending.

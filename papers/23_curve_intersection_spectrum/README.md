# Small eigenvalues on degenerating curves and the intersection form

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 10 November 2026. This version: 10 November 2026. 11 pages.

[Read the paper](curve-intersection-spectrum.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/curve-intersection-spectrum/) · [Zenodo:24120959](https://zenodo.org/records/24120959) · [DOI](https://doi.org/10.5281/zenodo.24120959)

Let a holomorphic family of compact Riemann surfaces over a disc, with smooth total space, degenerate to a reduced curve with N irreducible components and arbitrary singularities, and let the fibres carry the metrics induced by a Riemannian metric on the total space that is Hermitian near the singular points of the central fibre. We prove that the N−1 eigenvalues of the Laplace–Beltrami operator on functions that tend to zero, multiplied by log(1/|s|), where s is the parameter, converge to 2π times the eigenvalues of the negative intersection matrix of the components relative to the diagonal matrix of their areas, and that the remaining positive eigenvalues stay bounded away from zero. Equivalently, the small eigenvalues are asymptotic to those of a finite weighted graph on the components whose edge weights are intersection numbers, in analogy with Burger's theorem for hyperbolic surfaces. This answers questions of Dai and Yoshikawa on the existence and geometric meaning of these limits, and evaluates the constant in their asymptotic formula for the product of the small eigenvalues as a weighted count of spanning trees. For Kähler metrics the limits can also be derived from recent work of Cao on the Archimedean height pairing combined with the eigenvalue bounds of Dai and Yoshikawa; the proof given here is direct. It uses the logarithms of local equations of the components as explicit harmonic functions on the smooth fibres, whose Dirichlet energies are computed by winding numbers.

## Build

With an AMS LaTeX installation, run `latexmk -pdf main.tex` in this directory.
The descriptive PDF is the released manuscript; rebuilding writes `main.pdf`.
The exact released checksum is recorded in `../../manuscripts.json`.

[The programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

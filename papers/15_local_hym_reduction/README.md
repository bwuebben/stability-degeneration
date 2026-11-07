# Quantitative local models for Hermitian–Einstein metrics

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).
Paper 15. Preprint, 7 November 2026. 22 pages.

[Read the paper](local-hym-reduction.pdf) · [Full abstract and citation](https://bwuebben.github.io/papers/local-hym-reduction/) · [Zenodo:24029090](https://zenodo.org/records/24029090) · [DOI](https://doi.org/10.5281/zenodo.24029090)

Near any fixed locally free polystable bundle on a compact Kähler manifold,
the paper constructs a common quantitative local model for Hermitian–Einstein
metrics, the complete small endomorphism spectrum and harmonic metrics on
smooth pieces of the reduced moduli space. It retains the full convergent
Kuranishi obstruction. The comparisons are relative and uniform through
stabilizer changes and arbitrary separation of positive finite eigenvalues.
Canonical harmonic coordinates give quadratic error at fixed polarization;
moving harmonic backgrounds preserve that order for varying polarizations.

A stable rank-two bundle on a smooth bicubic Calabi–Yau threefold has analytic
deformation ring C{a,b}/(a²b,ab²). The theorem applies to the full deformation
problem for three copies, with eighteen harmonic coordinates and eight small
trace-free modes. A diagonal subfamily illustrates arbitrarily separated
positive scales and collision-kernel dimensions two, four and eight. The full
matrix-valued obstruction is retained; the scalar cubic equations are not
replaced by noncommutative substitution.

The results concern local stationary bundle moduli. HYM evolution, singular
graded replacement laws and global metric completion remain separate questions.

## Build and reproduce

With an AMS LaTeX installation, run `latexmk -pdf paper.tex` in this directory.
The source directly includes the supplied preamble and sections. The released
PDF is `local-hym-reduction.pdf`; a local build produces `paper.pdf`.

Run `python3 figures/code/compute06_bicubic.py`. The script uses only the Python
standard library and prints exact checks of the cubic coefficients, monomial
ideal decomposition, degree polynomials, spectral multiplicities and collision
kernels. The expected output is in `reproducibility/finite_identities.json`.
These finite checks corroborate the printed identities; they do not replace
the geometric or elliptic proofs.

[Programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved.

# Mass increase and acyclic cancellation for phase flows on quiver complexes

By [Bernd Johannes Wuebben](https://bwuebben.github.io/),
[ORCID](https://orcid.org/0009-0002-4588-3284).
Paper 11 in the Stability and Geometric Degeneration series.
Preprint. Original date: 2 November 2026. This version: 3 November 2026. 27 pages.

[Read the paper](paper.pdf) ·
[Abstract and citation](https://bwuebben.github.io/papers/mass-increase-acyclic-cancellation/) ·
[Zenodo:23917733](https://zenodo.org/records/23917733) ·
[DOI](https://doi.org/10.5281/zenodo.23917733).

For two-term complexes over a two-vertex parallel-arrow quiver, trace-norm mass
can increase along the positivity-defined phase flow. Nevertheless, an invertible
source differential gives exponential acyclic cancellation for arbitrary positive
initial metrics and upper-half-plane charges. For charges −1+i and 1+i, initially
balanced cohomology and a sufficiently large least acyclic singular value give
quantitative cancellation in arbitrary surviving dimension and acyclic rank,
with arbitrary bounded mixing. Constants are independent of the number of
arrows, the surviving stabilizer and the acyclic condition number. The surviving
moment-map defects tend to zero and the mass tends to twice the surviving
dimension. The exponential comparison is with instantaneous cohomology; no
uniform final relaxation rate or general quotient metric is asserted. The
rank-two Pauli family is a corollary exhibiting a stabilizer change.

## Build and reproduce

The main [LaTeX source](paper.tex) includes `sections/*.tex` and `references.tex`.
From this directory, run `latexmk -pdf paper.tex` with a standard AMS LaTeX
installation. The checked-in PDF is the archived manuscript; its checksum is in
the repository's `manuscripts.json`.

Install the pinned versions in `requirements.txt` in a virtual environment.
Run `figures/code/compute05_interval.py` and `figures/code/compute05_independent.py`
to reproduce two independent interval enclosures of the limiting mass derivative.
The second program uses only standard-library integer arithmetic. Run
`figures/code/fig07_mass_profile.py` to regenerate the exact static mass profile
and its limiting curve. These finite calculations support the stated examples;
the dynamical and uniform estimates rely on the proofs in the manuscript.

[Programme and companion papers](../../README.md) ·
[Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved.
No blanket code license applies to the manuscript.

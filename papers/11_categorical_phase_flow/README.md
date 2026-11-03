# Mass increase and acyclic cancellation for phase flows on quiver complexes

By [Bernd Johannes Wuebben](https://bwuebben.github.io/),
[ORCID](https://orcid.org/0009-0002-4588-3284).
Paper 11 in the Stability and Geometric Degeneration series.
Preprint. Original date: 2 November 2026. This version: 3 November 2026. 29 pages.

[Read the paper](paper.pdf) ·
[Abstract and citation](https://bwuebben.github.io/papers/mass-increase-acyclic-cancellation/) ·
[Zenodo:23921919](https://zenodo.org/records/23921919) ·
[DOI](https://doi.org/10.5281/zenodo.23921919).

Trace-norm mass can increase along the positivity-defined phase flow on the
specified two-term quiver complexes. An invertible source differential nevertheless
gives exponential acyclic cancellation from arbitrary positive metrics, for
strict upper-half-plane charges and without a large initial-scale assumption.

For a surjective source differential with surviving cohomology at both vertices,
a sufficiently large least initial acyclic singular value gives exponential
comparison with instantaneous cohomology. Theorem 1.3 allows arbitrary surviving
tuples, unequal dimensions, arbitrary finite arrow count and fixed strict
upper-half-plane charges. For fixed dimensions, acyclic rank, charges and total
initial norm bounds, constants are independent of arrow count, stabilizer and
acyclic condition number. The initial threshold is also independent of acyclic
rank. Theorem 1.4 separately proves relaxation for equal dimensions and charges
−1+i, 1+i when the cohomology Kempf–Ness function has finite infimum, without
initial balance or attainment. Its energy estimate retains the initial excess.
The balanced and Pauli cases are Corollaries 1.5–1.6.

Static coefficient and mass estimates hold for arbitrary fixed invertible
external differentials. An explicit scalar-scale family converges in C2 while
every positive Hölder seminorm of its Hessian diverges. No uniform final
cohomology relaxation rate, normalized mixed-metric limit, general arrow-space
C2 estimate or unrestricted quotient metric is asserted.

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
and its analytically differentiated Hessian. The independent high-precision
check is `validation/check_revision.py`. These finite calculations support the stated examples;
the dynamical and uniform estimates rely on the proofs in the manuscript.

[Programme and companion papers](../../README.md) ·
[Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved.
No blanket code license applies to the manuscript.

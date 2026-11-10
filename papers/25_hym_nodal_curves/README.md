# Effective Hermitian–Yang–Mills flow on nodal degenerations of curves

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 10 November 2026. This version: 10 November 2026. 16 pages.

[Read the paper](effective-hym-nodal-curves.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/effective-hym-nodal-curves/) · [Zenodo:24120966](https://zenodo.org/records/24120966) · [DOI](https://doi.org/10.5281/zenodo.24120966)

We prove a finite-dimensional limit for Hermitian–Yang–Mills metric flow on the trivial bundle over a compact curve approaching a reduced nodal fibre. The total space is smooth and Kähler, the central components are smooth, and each node joins distinct components. On the time scale log(1/|t|), where t is the degeneration parameter, the limiting component metrics follow the gradient flow of the sum of their squared symmetric-space distances, with component areas as masses and one coupling for each node. For uniformly bounded determinant-one initial metrics with prescribed local L² limits away from the nodes, we prove uniform global L² tracking on bounded rescaled time intervals and pointwise tracking away from the nodes after a fixed physical-time transient. The construction combines component Green functions, exact stationary neck metrics and the Donaldson trace inequality. At the constant flat metric, we also identify the rescaled linearized heat evolution on the small spectral subspace and prove a uniform gap on its orthogonal complement. A rank-two example on a smoothing of a triangle of projective lines has component metrics that admit no simultaneous diagonalization by a fixed holomorphic frame.

## Build

With an AMS LaTeX installation, run `latexmk -pdf main.tex` in this directory.
The descriptive PDF is the released manuscript; rebuilding writes `main.pdf`.
The exact released checksum is recorded in `../../manuscripts.json`.

## Reproduce the finite example

Run `uv run --with sympy python figures/code/compute_initial_model.py`.
This checks the displayed triangle matrices, commutator, graph spectrum and
mode count. It does not certify the analytic flow estimates.

[The programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

# All-time approximation of Hermitian–Yang–Mills flow on degenerating curves

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 11 November 2026. This version: 11 November 2026. 29 pages.

[Read the paper](all-time-hym-curves.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/all-time-hym-curves/) · [Zenodo:24149212](https://zenodo.org/records/24149212) · [DOI](https://doi.org/10.5281/zenodo.24149212)

We prove an all-time approximation of Hermitian–Yang–Mills flow by a finite matrix gradient flow on a degenerating family of compact genus-three curves. Four independent node parameters and two component-holonomy parameters determine flat unitary bundles passing from stable to strictly polystable. With fixed infinitesimal holonomy generators and node identifications, the spatially integrated squared error is bounded, uniformly in time, by the initial integrated trace discrepancy plus an error vanishing as all six parameters tend independently to zero. Actual initial metrics are uniformly bounded and finite initial data lie in a small ball; no lower bound on the weakest decay rate is required. An open spatial initial class enters after preparation negligible on every finite relaxation scale. The proof writes the first-order component-frame error as a global divergence using its zero angular mean, and combines HYM energy with time-integrated bounds from the finite gradient flow. We also prove relative asymptotics for every positive small eigenvalue of the endomorphism connection Laplacian, with exact correspondence of kernels and varying node identifications. One stable family has relaxation scales of orders h⁻¹, h⁻² and h⁻⁴ as h→0, all covered by the nonlinear theorem.

## Build and reproduce

Run `latexmk -pdf main.tex` in this directory. The descriptive PDF is the released manuscript; rebuilding writes `main.pdf`.

The exact computations and figure sources are in `figures/code/`. Run each Python program with `uv run --no-project --with-requirements requirements.txt python figures/code/<program>.py`, or run `python3 reproduce.py` at the repository root. These programs verify finite identities and draw the accompanying figures; the analytic assertions rely on the manuscript's proofs.

[The programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

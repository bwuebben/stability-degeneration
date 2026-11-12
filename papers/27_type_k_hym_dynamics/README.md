# Hermitian–Yang–Mills flow on degenerating Calabi–Yau threefolds of type K

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 11 November 2026. This version: 11 November 2026. 28 pages.

[Read the paper](type-k-hym-dynamics.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/type-k-hym-dynamics/) · [Zenodo:24149210](https://zenodo.org/records/24149210) · [DOI](https://doi.org/10.5281/zenodo.24149210)

We construct a rank-four holomorphic bundle family on a compact degeneration, with smoothing parameter t, of Enriques Calabi–Yau threefolds. The threefolds are free quotients of products of a K3 surface and an elliptic curve degenerating to a cycle of four rational curves. On the unitary gluing locus, the bundles have nonflat Hermitian–Einstein connections for prescribed product-induced Kähler metrics, which need not be Ricci-flat. We identify all five small eigenvalues of the trace-free endomorphism Laplacian, prove relative asymptotics uniform up to strict polystability, and establish a uniform complementary gap. The effective quadratic energy lies on an interval with one-dimensional endpoint constraints whose relative angle records holonomy around the elliptic central-fibre cycle. We derive a finite-dimensional gradient flow of noncommuting positive Hermitian matrices and reconstruct comparison metrics approximating actual Hermitian–Yang–Mills trajectories in global L² affine distance. With determinant fixed to that of the Hermitian–Einstein reference metric, an open class of initial metrics whose relative logarithms are globally bounded and whose component-exterior restrictions admit uniformly bounded C⁶,ᵝ extensions to the normalization components, 0<β<1, enters this description after physical time s=o(L), where L=log(1/|t|). Tracking is uniform on bounded s/L intervals and, where the least effective eigenvalue is bounded below by a positive constant, for all subsequent times. Within the specified family, the five eigenvalues, elliptic component areas and known normalization L asymptotically recover the projective conjugacy class of the cycle transport. Infinite physical-time limits need not commute with degeneration. The rescaled limits commute under this gap bound when the selected parallel limits of the component flows and the gluing parameters converge, but can fail near reducibility.

## Build and reproduce

Run `latexmk -pdf main.tex` in this directory. The descriptive PDF is the released manuscript; rebuilding writes `main.pdf`.

The exact computations and figure sources are in `figures/code/`. Run each Python program with `uv run --no-project --with-requirements requirements.txt python figures/code/<program>.py`, or run `python3 reproduce.py` at the repository root. These programs verify finite identities and draw the accompanying figures; the analytic assertions rely on the manuscript's proofs.

[The programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

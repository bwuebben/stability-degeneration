# Integral vanishing cycles from periods on conifold smoothings

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 12 November 2026. This version: 12 November 2026. 36 pages.

[Read the paper](conifold-integral-reconstruction.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/conifold-integral-reconstruction/) · [Zenodo:24175968](https://zenodo.org/records/24175968) · [DOI](https://doi.org/10.5281/zenodo.24175968)

We give a terminating construction of the marked integral vanishing-cycle homomorphism for a six-parameter family of Calabi–Yau conifold smoothings in ℙ³ × ℙ³. Its eight local spheres span a rank-six lattice in free middle homology of rank 134. Starting from algebraically selected equations, we construct a known neighbourhood in vanishing-period coordinates, a geometric symplectic marking, and absolute error bounds for transport and integration of represented cycles. Picard–Lefschetz period jumps recover the integral columns, their saturation and monodromy without a height prior or a restriction on nonzero period ratios. We verify the hypotheses of established analytic results on a possibly ineffective smaller neighbourhood. There the capacity-weighted squared norm of restriction of divisor classes from the small resolution to its exceptional curves governs the two small positive Ricci-flat two-form eigenvalues; the image of this restriction is the kernel of the recovered vanishing-cycle homomorphism. For nonflat bundles on the same fibres, a finite-dimensional invariant manifold of Hermitian–Yang–Mills metric flow is tangent at each constructed Hermitian–Einstein metric to the full low eigenspace on trace-free Hermitian endomorphisms.

## Build and reproduce

Run `latexmk -pdf main.tex` in this directory. The descriptive PDF is the released manuscript; rebuilding writes `main.pdf`.

Figure programs are supplied in `figures/code/` or `scripts/` where applicable. Run a supplied program with `uv run --no-project --with-requirements requirements.txt python <program-path>`, or run `python3 reproduce.py` at the repository root. The integral reconstruction procedure in Paper 31 is proved to terminate; a numerical execution of its full integral matrix is not included. These programs verify finite identities and draw the accompanying figures; the analytic assertions rely on the manuscript's proofs.

[The programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

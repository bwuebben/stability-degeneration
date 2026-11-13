# Integral vanishing cycles from periods on conifold smoothings

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 12 November 2026. This version: 12 November 2026. 37 pages.

[Read the paper](conifold-integral-reconstruction.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/conifold-integral-reconstruction/) · [Zenodo:24182113](https://zenodo.org/records/24182113) · [DOI](https://doi.org/10.5281/zenodo.24182113)

We give a terminating construction of the marked integral vanishing-cycle homomorphism for a six-parameter family of Calabi–Yau conifold smoothings in ℙ³ × ℙ³. Its eight local spheres span a primitive rank-six lattice in free middle homology of rank 134. Integral topology determines its symplectic normal form and torsion-free quotient. The period construction recovers its coordinates in an independently constructed geometric marking. Starting from algebraically selected equations, we construct a known neighbourhood in vanishing-period coordinates, a geometric symplectic marking, and absolute error bounds for transport and integration of represented cycles. Picard–Lefschetz period jumps recover the integral columns and monodromy without an a priori bound on the integral coordinates or a restriction on nonzero period ratios. The recovered matrix is constant in the transported marking, so recovery at one reference fibre suffices throughout the smooth family. On a possibly ineffective smaller neighbourhood, we verify applicability of established Ricci-flat Hodge, Hermitian–Yang–Mills flow and harmonic moduli results, and identify the kernel of the recovered map with the divisor restriction image governing the weak Hodge energy.

## Build and reproduce

Run `latexmk -pdf main.tex` in this directory. The descriptive PDF is the released manuscript; rebuilding writes `main.pdf`.

Figure programs are supplied in `figures/code/` or `scripts/` where applicable. Run a supplied program with `uv run --no-project --with-requirements requirements.txt python <program-path>`, or run `python3 reproduce.py` at the repository root. The integral reconstruction procedure in Paper 31 is proved to terminate; a numerical execution of its full integral matrix is not included. These programs verify finite identities and draw the accompanying figures; the analytic assertions rely on the manuscript's proofs.

[The programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

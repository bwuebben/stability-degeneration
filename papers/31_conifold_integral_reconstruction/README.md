# Integral vanishing cycles from periods on conifold smoothings

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 12 November 2026. This version: 13 November 2026. 37 pages.

[Read the paper](conifold-integral-reconstruction.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/conifold-integral-reconstruction/) · [Zenodo:24185494](https://zenodo.org/records/24185494) · [DOI](https://doi.org/10.5281/zenodo.24185494)

We give a terminating construction of the marked integral vanishing-cycle map for a six-parameter family of Calabi–Yau conifold smoothings in ℙ³ × ℙ³. The eight vanishing spheres span a primitive rank-six lattice whose symplectic normal form is determined by topology. We recover their coordinates in an independent geometric marking by constructing an effective neighbourhood in vanishing-period coordinates and absolute error bounds for cycle transport and integration. Picard–Lefschetz period jumps determine the integral columns and monodromy without an a priori bound on the integer entries or restrictions on nonzero period ratios. Recovery at one reference fibre suffices throughout the smooth locus. As applications, we verify the hypotheses of established analytic results on the same fibres in a possibly ineffective smaller neighbourhood. Over ℝ, the kernel of the recovered map is the image of divisor restriction to the exceptional curves; its capacity-weighted quadratic form and limiting harmonic mass give relative asymptotics for the two small positive Ricci-flat two-form eigenvalues. For the constructed nonflat polystable rank-four bundles, a local invariant manifold of Hermitian–Yang–Mills metric flow is tangent at the Hermitian–Einstein metric to the complete low spectral subspace of the Laplacian on trace-free Hermitian endomorphisms, including its kernel.

## Build and reproduce

Run `latexmk -pdf main.tex` in this directory. The descriptive PDF is the released manuscript; rebuilding writes `main.pdf`.

Figure programs are supplied in `figures/code/` or `scripts/` where applicable. Run a supplied program with `uv run --no-project --with-requirements requirements.txt python <program-path>`, or run `python3 reproduce.py` at the repository root. The integral reconstruction procedure in Paper 31 is proved to terminate; a numerical execution of its full integral matrix is not included. These programs verify finite identities and draw the accompanying figures; the analytic assertions rely on the manuscript's proofs.

[The programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

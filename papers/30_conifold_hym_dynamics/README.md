# Hermitian–Yang–Mills dynamics and harmonic moduli near a conifold

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 12 November 2026. This version: 12 November 2026. 51 pages.

[Read the paper](conifold-hym-dynamics.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/conifold-hym-dynamics/) · [Zenodo:24175924](https://zenodo.org/records/24175924) · [DOI](https://doi.org/10.5281/zenodo.24175924)

We construct a nonflat rank-four bundle family on compact Calabi–Yau threefolds smoothing eight ordinary double points and compare its Hermitian–Einstein metrics, metric flow and local moduli geometry uniformly through degeneration. A two-arrow matrix form determines all seven small trace-free endomorphism eigenvalues with relative error, including changing kernels and arbitrarily separated positive rates. On each polystable member, an open spatial C¹ basin admits a finite invariant graph and an initial phase for which the actual Hermitian–Yang–Mills flow tracks an exact graph trajectory with exponentially decaying C¹ error for all time. At the Hermitian–Einstein equilibrium, the Hessian of twice Donaldson's functional restricted to the graph is relatively equivalent to the same matrix form in the original metric coordinates. We also identify the actual local bundle classes and compare their harmonic metric on every stabilizer stratum. For paths confined to this locus, the unique metric tangent cone at the split bundle is (ℂ², 4|du|²)/((ℤ/2)² ⋊ S₂). All comparisons allow independent node and bundle scales.

## Build and reproduce

Run `latexmk -pdf main.tex` in this directory. The descriptive PDF is the released manuscript; rebuilding writes `main.pdf`.

Figure programs are supplied in `figures/code/` or `scripts/` where applicable. Run a supplied program with `uv run --no-project --with-requirements requirements.txt python <program-path>`, or run `python3 reproduce.py` at the repository root. The integral reconstruction procedure in Paper 31 is proved to terminate; a numerical execution of its full integral matrix is not included. These programs verify finite identities and draw the accompanying figures; the analytic assertions rely on the manuscript's proofs.

[The programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

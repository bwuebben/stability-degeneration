# The L² metric at the Uhlenbeck boundary

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 13 November 2026. This version: 13 November 2026. 22 pages.

[Read the paper](l2-metric-uhlenbeck-boundary.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/l2-metric-uhlenbeck-boundary/) · [Zenodo:24211932](https://zenodo.org/records/24211932) · [DOI](https://doi.org/10.5281/zenodo.24211932)

We construct local Kähler potentials for the L² metric on moduli spaces of Hermitian–Yang–Mills connections over polarized projective surfaces that extend continuously to ideal connections. The potential combines a zeta determinant on a curve avoiding the bubble points with a weighted curvature integral and an explicit term for the defect cycle. On a projective K3 surface with a generic polarization and its Ricci-flat metric, let the Mukai vector be primitive, of rank at least two and positive square, with nonempty stable bundle locus. The potential identifies the L² metric with the canonical singular Ricci-flat metric on the normalized Uhlenbeck compactification. Its completion is the Uhlenbeck compactification, and in complex dimension 2N its volume is (16π⁴ Vol X)ᴺ/N!. In rank two, the strata with stable residual bundle and distinct simple bubble points are totally geodesic hyperkähler orbifold strata, with the base-surface metric scaled by 4π². At stable-residual boundary points, we identify the affine tangent cones and their Reeb actions and compute their volume densities from instanton partition functions.

The full identification of the tangent metric with the flat ADHM quotient remains conditional outside the separately proved rank-two distinct-point orbifold case.

## Build and reproduce

Run `latexmk -pdf main.tex` in this directory. The descriptive PDF is the released manuscript; rebuilding writes `main.pdf`.

The two programs in `figures/code/` reproduce the normalization constants, volume and curvature coefficients, examples and instanton densities printed in the manuscript. Run each with `uv run --no-project --with-requirements requirements.txt python <program-path>`. They use exact SymPy calculations, with cancellation before specialization in the localization formula. Their expected output is in `figures/output/`. These finite computations complement the analytic and geometric proofs.

[Programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/#global-geometry).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

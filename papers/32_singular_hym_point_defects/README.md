# Point defects and Hermitian–Yang–Mills flow on surfaces

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 12 November 2026. This version: 12 November 2026. 33 pages.

[Read the paper](point-defects-hym-flow.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/point-defects-hym-flow/) · [Zenodo:24175884](https://zenodo.org/records/24175884) · [DOI](https://doi.org/10.5281/zenodo.24175884)

We determine the leading metric and concentration scales of Hermitian–Yang–Mills flow for locally free extensions 0 → L → V → 𝓘_W M → 0 on a compact Kähler surface, where L, M are distinct line bundles of equal degree and W is a nonempty reduced finite scheme. Fix Hermitian–Einstein metrics on the line factors and their product as determinant. For every smooth initial metric with this determinant, the integrated squared norm Z(T) of the fixed inclusion satisfies Z(T) ∼ v₀²/(4π²σT), where v₀ is the volume and σ = ∑ₚ∈W (det Qₚ)⁻¹ᐟ²; Qₚ is the quadratic form of the inclusion norm at p in orthonormal coordinates. On the radius {Z(T)/(v₀√(det Qₚ))}¹ᐟ², the actual metrics converge smoothly to an explicit charge-one anti-self-dual metric. The inclusion half-height radius and the half-defect-energy radius have this same leading value. The proof corrects a compact family of metrics orthogonally to its small connection eigenspace and compares the resulting trajectory with every actual flow. It also gives a global inclusion profile and quantitative convergence outside shrinking balls. For arbitrary finite defects we prove Z(T) = O(T⁻¹); if a defect is nonreduced, then TZ(T) → 0.

## Build and reproduce

Run `latexmk -pdf main.tex` in this directory. The descriptive PDF is the released manuscript; rebuilding writes `main.pdf`.

Figure programs are supplied in `figures/code/` or `scripts/` where applicable. Run a supplied program with `uv run --no-project --with-requirements requirements.txt python <program-path>`, or run `python3 reproduce.py` at the repository root. The integral reconstruction procedure in Paper 31 is proved to terminate; a numerical execution of its full integral matrix is not included. These programs verify finite identities and draw the accompanying figures; the analytic assertions rely on the manuscript's proofs.

[The programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

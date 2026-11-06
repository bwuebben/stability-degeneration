# Metric wall crossing in rank-three K3 bundle moduli

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).
Paper 13. Preprint, 6 November 2026. 21 pages.

[Read the paper](k3-metric-wall-crossing.pdf) · [Full abstract and citation](https://bwuebben.github.io/papers/k3-metric-wall-crossing/) · [Zenodo:24021822](https://zenodo.org/records/24021822) · [DOI](https://doi.org/10.5281/zenodo.24021822)

The paper compares actual harmonic L² metrics on a complex twelve-dimensional
rank-three K3 bundle-moduli neighborhood with the hyperkähler quotient metrics
on T*Gr(2,5) and T*Gr(3,5). The relative error is quadratic in deformation size,
uniform through stabilizer changes and unrelated rates of approach to the wall
and its singular strata. Consequences include exceptional-fibre collapse,
complete pointed Gromov–Hausdorff limits after rescaling, quantitative matching
to the wall metric and intrinsic distance asymptotics to the contracted loci.
The theorem is local and concerns locally free bundles.

Projection-relative estimates and finite balancing from the repeated-factor
papers, together with complementary-gauge minimization from the wall-moduli
paper, supply the methodological starting point. The necessary arguments are
developed here on the exact obstructed K3 deformation locus. The result does
not assert global metric completion or describe bubbling.

## Build and reproduce

With an AMS LaTeX installation, run `latexmk -pdf paper.tex` in this directory.
The source directly includes the supplied section files and figure. The released
PDF is named `k3-metric-wall-crossing.pdf`; a local build produces `paper.pdf`.

With [uv](https://docs.astral.sh/uv/), run:

    uv run --no-project reproducibility/k3_model.py
    uv run --no-project reproducibility/plot_distance_profile.py

Dependencies are pinned in the scripts. The [calculation guide](reproducibility/README.md)
explains the polynomial calculations and the distance-integral illustration.
These computations support the explicit examples; the geometric and analytic
results depend on the manuscript's proofs.

[Programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved.

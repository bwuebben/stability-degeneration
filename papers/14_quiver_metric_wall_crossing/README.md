# Quantitative metric wall crossing for K3 bundle moduli

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).
Paper 14. Preprint, 6 November 2026. 30 pages.

[Read the paper](quiver-metric-wall-crossing.pdf) · [Full abstract and citation](https://bwuebben.github.io/papers/quiver-metric-wall-crossing/) · [Zenodo:24024998](https://zenodo.org/records/24024998) · [DOI](https://doi.org/10.5281/zenodo.24024998)

The paper proves a uniform local comparison between actual harmonic L² metrics
and finite quiver quotient metrics near polystable bundles on a projective K3
surface. The locally free factors satisfy the stated slope-stability, common
reduced Hilbert polynomial and central Ricci-flat background hypotheses. The
comparison remains relative through stabilizer changes and arbitrary separation
of positive degeneration scales. The spatial error is linear, improving to
quadratic under an additional central sign symmetry.

A separate intrinsic theorem proves pointed Gromov–Hausdorff continuity of the
complete stratified quotient metrics at arbitrary converging real parameters
when the complex zero moment fiber admits a module-simple representation.
Complete actual rescaled limits and finite successive limits retain convergence
of the rescaled parameters and balanced normal base representatives at each stage.

The geometrically realized rank-five example has three spherical line factors,
multiplicities (2,2,1), a doubled-triangle Ext-quiver and two independent
polarization parameters. Its local moduli have complex dimension 48. Five
realized walls divide the parameter plane into ten chambers; selected transverse
contractions have models T*P⁹, T*P¹¹ and T*P¹⁷. The paper determines contraction
scales, distances to the corresponding strata and matching to the surrounding
wall metric. All conclusions concern local stationary bundle moduli; global
metric completion, bubbling and HYM evolution are separate questions.

## Build and reproduce

With an AMS LaTeX installation and TikZ, run `latexmk -pdf paper.tex` in this
directory. The source includes the supplied preamble, sections and wall figure.
The released PDF is named `quiver-metric-wall-crossing.pdf`; a local build
produces `paper.pdf`.

Run `python3 figures/code/check_quiver_geometry.py`. Only the Python standard
library is required. The script regenerates `checks/quiver_geometry.json` and
checks the lattice embedding, Chern and Mukai data, geometric parameter signs,
five wall normals, table dimensions and eleven simple balanced representations.
The printed formulas and matrices supply the inputs. These exact finite checks
support the example; they do not replace its geometric and analytic proofs.

[Programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved.

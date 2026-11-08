# Integral attachments and wall correspondences for K3 bundle moduli

By [Bernd Johannes Wuebben](https://bwuebben.github.io/),
[ORCID](https://orcid.org/0009-0002-4588-3284).
Preprint. Original date: 8 November 2026. This version: 8 November 2026. 24 pages.

[Current version (PDF)](integral-attachments.pdf) ·
[Abstract and citation](https://bwuebben.github.io/papers/integral-attachments/) ·
[Zenodo:24054004](https://zenodo.org/records/24054004) · [DOI](https://doi.org/10.5281/zenodo.24054004)

For a compact six-chamber family of rank-three K3 bundle moduli, the integral
complement ring and all boundary restrictions are calculated in an explicit
Hilbert-square marking. Four adjacent crossings are elementary Mukai flops with
their corrected correspondences; two are isomorphisms over divisorial targets.
Opposite divisor signs determine the reflection action. The calculation includes
all normalized topological universal-index Chern classes, the mod-eight boundary
restriction, signed connecting maps, and the determinant-line correction under
wall transport. The resulting classes have exact periods and volumes in the
established harmonic metric degeneration of this family.

The construction uses the actual marked K3 surface and compact moduli geometry.
Local singularity types alone do not supply the global multiplication. The
cohomological path relation makes no assertion about Chow or derived compositions,
genuine smoothing monodromy, or reconstruction of a manifold or metric tensor.

## Build and reproduce

The [LaTeX source](paper.tex) inputs the supplied sections and TikZ figures.
With an AMS LaTeX installation and TikZ, run the following in this directory:

    latexmk -pdf paper.tex

The archived PDF is integral-attachments.pdf; a local build creates paper.pdf.
Released hashes are recorded in the repository's manuscripts.json.
Run the four exact arithmetic checks with:

    python3 reproducibility/run_all.py

They use the Python standard library. [Calculation guide](reproducibility/README.md).
Rebuild the standalone figure PDF/PNG exports with:

    python3 figures/code/build_figures.py

This requires LaTeX and Poppler.

[Programme and companion papers](../../README.md) ·
[Stability and Geometric Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. No blanket code
license applies to the manuscript.

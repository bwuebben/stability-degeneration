# Harmonic paths and honeycombs

By [Bernd Johannes Wuebben](https://bwuebben.github.io/),
[ORCID](https://orcid.org/0009-0002-4588-3284).
Paper 12 in the Stability and Geometric Degeneration series.
Preprint. Original date: 4 November 2026. 45 pages.

[Read the paper](paper.pdf) ·
[Abstract and citation](https://bwuebben.github.io/papers/harmonic-paths-and-honeycombs/) ·
[Zenodo:23947189](https://zenodo.org/records/23947189) ·
[DOI](https://doi.org/10.5281/zenodo.23947189).

For the identity representation of the A₂ quiver over the specified real-exponent
Novikov field, every triple of full diagonalizable boundary norms in every finite
rank admits a harmonic norm path. The charges are symmetric and the plane carries
the fixed radial area form. The path has finite apartment presentations, exact
normalized ends, a finite balanced honeycomb support, actual semistable edge
labels and compatible local residue filtrations at the vertices.

Uniqueness holds in the finite regular principal-velocity class with intrinsic
first-order matching at joins. This is an additional condition in the comparison
with HKKP Conjecture 7.10. General boundary continuity and a global
Fukaya–Maurer–Cartan correspondence remain open. A worked collapsing-face family
has its own full-path Lipschitz estimate.

## Build and reproduce

The [LaTeX source](paper.tex) includes `sections/*.tex` and the TikZ diagram
`figures/diamond-collapse.tex`. With a standard AMS LaTeX installation and TikZ,
run `latexmk -pdf -no-shell-escape paper.tex` in this directory. The checked-in
PDF is the archived manuscript; its checksum is recorded in `../../manuscripts.json`.

Run `uv run reproducibility/check_examples.py` to reproduce the exact valuation,
minor and finite-graph calculations. SymPy 1.14.0 is pinned in the script.
The [reproduction guide](reproducibility/README.md) describes the calculations
and their scope. They verify the worked examples, not the analytic existence
or uniqueness theorems.

[Programme and companion papers](../../README.md) ·
[Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved.
No blanket code license applies to the manuscript.

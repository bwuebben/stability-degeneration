# Integral reconstruction and quantitative metric degeneration of K3 bundle moduli

By [Bernd Johannes Wuebben](https://bwuebben.github.io/),
[ORCID](https://orcid.org/0009-0002-4588-3284).
Preprint. Original manuscript date: 7 November 2026.
This version: 7 November 2026.

[Read the paper](integral-reconstruction.pdf) ·
[Abstract and citation](https://bwuebben.github.io/papers/integral-reconstruction/) ·
[Zenodo:24047889](https://zenodo.org/records/24047889) · [DOI](https://doi.org/10.5281/zenodo.24047889)

A compact moduli eightfold of rank-seven bundles on a K3 surface has moving
projective-line fibres that specialize to Schubert lines in a deepest
Grassmannian Gr(2,4). Both middle Schubert classes extend integrally, every
homology inclusion from the fibre is primitive, and its self-intersection six
gives an explicitly calculated degree-eight attachment modulo six. These same
cycles have exact harmonic periods and volumes, ambient diameter asymptotics,
and complete cotangent-Grassmannian rescaled limits.

For compact K3 bundle-moduli contractions without bubbling, complete convergent
quiver charts with the simple-zero property give a variational distortion bound
for the actual contraction graph. Whole-fibre coverage, principal-stratum path
approximation and both ambient distance comparisons identify the central harmonic
completion and prove compact convergence through arbitrary adjacent chambers.
Integral cochains retain the actual geometric overlaps and monodromy; the
middle-degree attachment calculation retains a marked complement ring.
The conclusions do not reconstruct an unknown complement from local fibre types.

## Build and reproduce

The [LaTeX source](paper.tex) directly includes its sections, bibliography and two
TikZ figures. With an AMS LaTeX installation and TikZ, run from this directory:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error paper.tex
python3 figures/code/check_rank_partitions.py
python3 figures/code/check_integral_data.py
```

The Python programs use only the standard library. They reproduce the finite
rank and charge tables, lattice witnesses, Schubert and index identities,
boundary data and Göttsche-product Betti numbers in the supplied JSON files.
The manuscript proves existence, compactness, all-integer formulas and metric
asymptotics; finite checks are supplementary. The checked-in
`integral-reconstruction.pdf` is the released manuscript; rebuilding produces
`paper.pdf`. Released checksums are in the repository's `manuscripts.json`.

[Programme and companion papers](../../README.md) ·
[Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved.

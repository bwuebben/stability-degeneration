# Exceptional cycles and harmonic metrics in K3 bundle moduli

By [Bernd Johannes Wuebben](https://bwuebben.github.io/),
[ORCID](https://orcid.org/0009-0002-4588-3284).
Preprint. Original manuscript date: 7 November 2026.

[Read the paper](exceptional-cycles.pdf) ·
[Abstract and citation](https://bwuebben.github.io/papers/exceptional-cycles/) ·
[Zenodo:24046092](https://zenodo.org/records/24046092) · [DOI](https://doi.org/10.5281/zenodo.24046092)

Compact K3 bundle-moduli fourfolds contain exceptional cycles whose integral
attachments and harmonic metric degeneration can be calculated together. In the
main two-parameter family, a moving projective line specializes with multiplicity
one to a line in each of two Lagrangian planes. The primitive plane classes have
intersection matrix [[3,1],[1,3]], and their degree-four attachment is an explicit
congruence modulo eight. An isolated-plane example gives the corresponding
congruence modulo three.

The same cycles have exact harmonic areas and volumes, ambient diameter estimates
uniform at unrelated positive scales, and complete Calabi rescaled limits. All six
chambers have a common compact harmonic completion. The proof includes exhaustive
factor and lost-charge classifications, complete exceptional fibres, a controlled
convergent common-slope K3 normal form, and an open-stratum Brauer calculation.
The marked complement remains input to the integral attachment; the constructed
compact moduli spaces contain bundles only and have no bubbling boundary.

## Build and reproduce

The [LaTeX source](paper.tex) directly includes its sections, bibliography and two
TikZ figures. With an AMS LaTeX installation and TikZ, run from this directory:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error paper.tex
python3 figures/code/check_isolated.py
python3 figures/code/check_interacting.py
```

The Python programs use only the standard library and write their exact JSON
results beside the scripts. They check lattice and Mukai identities, integral
witnesses, index normalizations, intersection matrices, boundary cokernels,
period coefficients and chamber supports. They support the printed proofs;
they do not prove existence, gauge compactness or metric asymptotics.
The checked-in `exceptional-cycles.pdf` is the released manuscript; rebuilding
produces `paper.pdf`. Released checksums are in the repository's `manuscripts.json`.

[Programme and companion papers](../../README.md) ·
[Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved.

# Uniform asymptotic phase for Hermitian–Yang–Mills flow near a polystable bundle

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).
Paper 16. Preprint, 7 November 2026.

[Read the paper](uniform-hym-phase.pdf) · [Full abstract and citation](https://bwuebben.github.io/papers/uniform-hym-phase/) · [Zenodo:24038504](https://zenodo.org/records/24038504) · [DOI](https://doi.org/10.5281/zenodo.24038504)

Confined Hermitian–Yang–Mills trajectories near a fixed locally free polystable
bundle admit an exact finite-dimensional comparison flow with continuously
differentiable initial phase and exponentially vanishing relative metric error.
The full convergent Kuranishi obstruction, arbitrary multiplicities and
noncommuting automorphisms are retained. The constants are controlled by the
complementary elliptic gap, independently of weak positive finite eigenvalues,
their ratios and stabilizer changes. Reconstructing complex gauges may be unbounded.

The theorem assumes confinement of the actual harmonic coefficients and small
initial normal displacement. The graph and phase are chosen; the exact effective
field is distinguished from its leading moment-map approximation. The result does
not establish universal confinement or classify the corrected flow's algebraic
asymptotics. A stable bundle on a smooth bicubic supplies an obstructed
eighteen-coordinate repeated-factor application, an explicit noncommuting
integrable family and its finite commutator spectrum. Those eigenvalues are
distinguished from actual HYM relaxation rates.

## Build and reproduce

Run `sh build.sh` or `latexmk -pdf paper.tex` in this directory with an AMS LaTeX
installation. The source includes the supplied preamble, sections and bibliography.
A build produces `paper.pdf`; the released PDF is `uniform-hym-phase.pdf`.

Run `python3 figures/code/compute06_spectrum.py`. It uses the Python standard
library to check exact coefficient Gram matrices, nilpotent matrix products and
complex and real Hermitian kernel ranks. Its expected output is
`checks/finite_spectrum.json`. These finite calculations corroborate the printed
identities; the geometric and parabolic results rely on their proofs.

[Programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved.

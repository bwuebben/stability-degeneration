# L² metrics on Hermitian–Yang–Mills moduli spaces near a stability wall

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).
Preprint. Original date: 31 October 2026. This version: 31 October 2026. 34 pages.

[Read the paper](paper.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/wall-moduli-metrics/) · [Zenodo:23850992](https://zenodo.org/records/23850992) · [DOI](https://doi.org/10.5281/zenodo.23850992).

Uniform relative harmonic L² metrics on full obstructed local quotients near a four-line stability wall, with quantitative changes of chart and intrinsic distance estimates.

## Build and reproduce

The assembled [LaTeX source](paper.tex) includes the bibliography. With an AMS LaTeX installation and TikZ, run `latexmk -pdf paper.tex` in this directory. The checked-in PDF is the archived manuscript; its checksum is recorded in the repository's `manuscripts.json`.

Install the versions in `requirements.txt` in a virtual environment. Run the 5 programs named `figures/code/compute*.py` from this directory to reproduce the finite identities and examples. The `fig*.py` programs regenerate the associated figures. These calculations do not replace the geometric, elliptic or parabolic proofs in the manuscript.

[The programme and companion papers](../../README.md) · [Geometry & moduli reading guide](https://bwuebben.github.io/research/geometry/#stability-degeneration).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. No blanket code license applies to the manuscript. arXiv submission is pending.

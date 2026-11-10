# Integral topology and small eigenvalues in Calabi–Yau conifold transitions

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 9 November 2026. This version: 9 November 2026. 41 pages.

[Read the paper](cy-conifold-hodge-spectrum.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/cy-conifold-hodge-spectrum/) · [Zenodo:24096887](https://zenodo.org/records/24096887) · [DOI](https://doi.org/10.5281/zenodo.24096887)

For a fixed polarized Calabi–Yau conifold smoothing, approached along a fixed ray and smoothing every node to first order, the lost divisor directions determine the entire small positive two-form spectrum of the actual Ricci-flat metrics. The leading inverse-logarithmic coefficients are generalized eigenvalues of the exceptional-period form relative to the full quotient harmonic norm. Global intersection theory computes that norm. The theorem assumes a simply connected projective small resolution; its constants and leading node-scale ratios belong to the fixed family.

The applications include the sixteen-node split quintic, the eight-node toric transition from the stability programme, and the original-side four-node family established in Calabi–Yau smoothing Papers V–VI. The split quintic also has explicit integral rings, attachment and specialization maps, vanishing lattice and monodromy image. Classical surgery and metric existence are inputs; the leading spectral theorem for the Ricci-flat metrics is the new analytic contribution. Independently separated node scales, arbitrary gauge bundles and HYM flow on these smoothing families remain further questions.

The leading spectral law is proved on exact conic models and transferred using a pointwise metric comparison. A separate argument establishes that comparison for the prescribed Ricci-flat family, including finitely many nodes. The differentiated estimate is used separately for the sharper rank-one upper remainder. The full surviving harmonic quotient and the finite-cap norm compactness argument are explicit. Supporting integral topology and tangent-sheaf calculations appear in Appendices A–B.

## Build and reproduce

With an AMS LaTeX installation and TikZ, run `latexmk -pdf main.tex` in this directory. The checked-in descriptive PDF is the released manuscript; rebuilding writes `main.pdf`. Its release checksum is recorded in `../../manuscripts.json`.

The first two exact computations use the Python standard library. The third requires SageMath:

```sh
python3 figures/code/compute01_split_quintic.py
python3 figures/code/compute02_eight_node.py
sage -python figures/code/compute03_x9.py --output figures/data/x9_exact.json
```

The JSON outputs record the intersection, period and spectral calculations. These computations supplement the proofs and do not certify the analytic estimates. To regenerate the diagram, run `pdflatex -output-directory=figures/output figures/code/fig01_transition.tex`.

[The programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

# Hermitian–Yang–Mills flow and gluing on degenerating curves

By [Bernd Johannes Wuebben](https://bwuebben.github.io/), [ORCID](https://orcid.org/0009-0002-4588-3284).

Preprint. Original date: 11 November 2026. This version: 11 November 2026. 39 pages.

[Read the paper](global-gluing-hym-curves.pdf) · [Abstract and citation](https://bwuebben.github.io/papers/global-gluing-hym-curves/) · [Zenodo:24149197](https://zenodo.org/records/24149197) · [DOI](https://doi.org/10.5281/zenodo.24149197)

For a six-parameter degeneration of plane quartics, we approximate Hermitian–Yang–Mills metric flow on flat unitary bundles with trivial component connections and noncommuting gluing by a finite gradient flow. Its energy is a sum of squared distances between positive Hermitian matrices after unitary transport, weighted by neck capacities. The integrated squared-distance comparison holds on time intervals proportional to the shortest logarithmic neck length, uniformly in all node-length ratios and as the gluing representation becomes reducible. We prove relative asymptotics for every small linearized eigenvalue and reach the comparison from an open class of bounded spatial initial metrics with uniform component regularity after a negligible initial layer. The comparison holds for all time when the smallest effective eigenvalue times the shortest logarithmic neck length stays bounded away from zero. A second compact family permits independently vanishing component holonomy and neck size. On the quartic family with equal node parameters we determine the leading Hodge metric on the unitary Schottky locus and asymptotically recover projective gluing from low eigenvalues and component averages of orthonormal eigensections in fixed unitary frames. On a finite class of first-order rational-component smoothings, these observations and the base-loop orientation determine unknown adjacency and the integral symplectic conjugacy class of monodromy.

## Build and reproduce

Run `latexmk -pdf main.tex` in this directory. The descriptive PDF is the released manuscript; rebuilding writes `main.pdf`.

The exact computations and figure sources are in `figures/code/`. Run each Python program with `uv run --no-project --with-requirements requirements.txt python figures/code/<program>.py`, or run `python3 reproduce.py` at the repository root. These programs verify finite identities and draw the accompanying figures; the analytic assertions rely on the manuscript's proofs.

[The programme and companion papers](../../README.md) · [Stability & Degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/).

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved. Public access does not grant a blanket reuse license.

# Nonsimplicial local stability chambers of toric tangent bundles

By Bernd Johannes Wuebben, 27 October 2026.

[Read the paper](paper.pdf). The self-contained LaTeX source is [paper.tex](paper.tex). Compile with `pdflatex paper.tex` twice, or `latexmk -pdf paper.tex`. An AMS LaTeX installation with TikZ is required.

The Python programs in `figures/code/` reproduce the finite calculations printed in the paper. They use only the Python standard library and write their records to `figures/data/`. Run `python3 -B figures/code/<program>.py` from this directory. The TeX figure sources are supplied alongside the computations; the manuscript includes its figures directly.

These computations check finite identities and examples. The geometric, elliptic and parabolic theorems rely on the proofs in the manuscript.

Archived version: [doi:10.5281/zenodo.23717722](https://doi.org/10.5281/zenodo.23717722).
[Abstract and citation](https://bwuebben.github.io/papers/nonsimplicial-toric-stability-chambers/) on the author's website.

# Stability walls and geometric degeneration

Three papers by [Bernd Johannes Wuebben](https://bwuebben.github.io/), connecting slope-stability chambers, small Dolbeault eigenvalues, Hermitian–Einstein metrics and Hermitian Yang–Mills flow.

| Paper | Manuscript | Contribution |
| --- | --- | --- |
| III | [Uniform Hermitian–Einstein spectra for bipartite two-step extensions](papers/iii/paper.pdf) | Uniform relative asymptotics for every small eigenvalue, including cyclic supports and arbitrarily separated positive scales. |
| IV | [Nonsimplicial local stability chambers of toric tangent bundles](papers/iv/paper.pdf) | Nonsimplicial chambers, positive-spanning slope access without local freeness, a Bott-filtration obstruction and a Calabi–Yau conifold transition. Original version 27 October; this version 29 October 2026. |
| V | [Metric response and Hermitian Yang–Mills flow for two-step extensions](papers/v/paper.pdf) | Intrinsic metric response, inverse estimates and actual-flow comparisons, with an integrated theorem on the same tangent family. |

III and IV have independent proofs. V combines their inputs with its own analytic arguments. The internal labels III–V refer to the research programme; each paper has its own title and bibliography. Earlier work on [equilibrium spectra near stability walls](https://bwuebben.github.io/papers/small-dolbeault-stability-walls/), [HYM flow spectra](https://bwuebben.github.io/papers/small-dolbeault-eigenvalues/), and [toric local stability chambers](https://bwuebben.github.io/papers/toric-local-stability-chambers/) gives entry points.

The spectral and analytic results concern specified two-step extensions with pairwise nonisomorphic stable factors. The tangent examples are smooth projective weak-Fano varieties, not Fano varieties. Cycles can retain harmonic extension data, and flows can retain initial-metric information. Repeated-factor actual metric calibration, higher-step extensions, singular graded objects and moduli-metric comparisons require further results.

## Reproduce the finite calculations

All included Python programs use only the Python standard library. From the repository root:

```sh
python3 reproduce.py
```

The programs check graph spectra and spanning-tree identities, toric intersections and chamber counts, matrix and metric-response identities, confinement estimates and scalar profiles. They do not numerically establish the elliptic or parabolic theorems. Their output records and native figure sources are included with each paper.

## Build the manuscripts

Each `paper.tex` is self-contained, including figures and bibliography. With a standard TeX Live or MiKTeX installation providing AMS packages and TikZ, run `latexmk -pdf paper.tex` inside the selected paper directory. The checked-in PDFs are the released versions; their checksums are recorded in `manuscripts.json`.

The [Geometry & moduli research guide](https://bwuebben.github.io/research/geometry/) explains the programme and its relation to the other papers.

Copyright © 2026 Bernd Johannes Wuebben. Existing author rights are reserved; this repository does not apply a blanket code license to the manuscripts.

## Archived manuscripts

- III: [doi:10.5281/zenodo.23717660](https://doi.org/10.5281/zenodo.23717660).
- IV: [doi:10.5281/zenodo.23774770](https://doi.org/10.5281/zenodo.23774770).
- V: [doi:10.5281/zenodo.23717742](https://doi.org/10.5281/zenodo.23717742).

The Zenodo records preserve these released PDFs. arXiv submissions remain pending.

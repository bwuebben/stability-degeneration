# Stability walls and geometric degeneration

Four papers by [Bernd Johannes Wuebben](https://bwuebben.github.io/), connecting slope-stability chambers, small Dolbeault eigenvalues, Hermitian–Einstein metrics and Hermitian Yang–Mills flow.

| Paper | Manuscript | Contribution |
| --- | --- | --- |
| III | [Uniform Hermitian–Einstein spectra for bipartite two-step extensions](papers/iii/paper.pdf) | Uniform relative asymptotics for every small eigenvalue, including cyclic supports and arbitrarily separated positive scales. |
| IV | [Nonsimplicial local stability chambers of toric tangent bundles](papers/iv/paper.pdf) | Nonsimplicial chambers, positive-spanning slope access without local freeness, a Bott-filtration obstruction and a Calabi–Yau conifold transition. Original version 27 October; this version 29 October 2026. |
| V | [Metric response and Hermitian Yang–Mills flow for two-step extensions](papers/v/paper.pdf) | Intrinsic metric response, inverse estimates and actual-flow comparisons, with an integrated theorem on the same tangent family. |
| Repeated factors | [Uniform Hermitian–Einstein spectra with repeated stable factors](papers/repeated_factor_spectra/main.pdf) | Actual metric comparison and every small matrix mode, with uniform relative accuracy at arbitrarily separated positive scales. |

III and IV have independent proofs. V combines their inputs with its own analytic arguments. The internal labels III–V refer to the research programme; each paper has its own title and bibliography. Earlier work on [equilibrium spectra near stability walls](https://bwuebben.github.io/papers/small-dolbeault-stability-walls/), [HYM flow spectra](https://bwuebben.github.io/papers/small-dolbeault-eigenvalues/), and [toric local stability chambers](https://bwuebben.github.io/papers/toric-local-stability-chambers/) gives entry points.

The analytic theorems concern specified locally free two-step extensions. The repeated-factor paper extends stationary HE metric and spectral comparison to arbitrary fixed multiplicities within either layer, retaining distinct stable types across the layers. The response and actual-flow results retain their distinct-factor hypotheses. The tangent examples are smooth projective weak-Fano varieties, not Fano varieties. Cycles can retain harmonic extension data, and flows can retain initial-metric information. Repeated-factor response and actual flow, higher-step extensions, singular graded objects and higher-dimensional bundle moduli-metric comparisons remain open. An abstract moment-map correction and spectral theory with new geometric applications is a further direction; it is not a theorem of these papers.

## Reproduce the finite calculations

The III–V calculations use the Python standard library. The repeated-factor calculations use the versions pinned in that paper’s `requirements.txt`; the runner requires `uv` to supply them. From the repository root:

```sh
python3 reproduce.py
```

The programs check graph spectra and spanning-tree identities, toric intersections and chamber counts, matrix and metric-response identities, confinement estimates and scalar profiles. They do not numerically establish the elliptic or parabolic theorems. Their output records and native figure sources are included with each paper.

## Build the manuscripts

For III–V, each `paper.tex` is self-contained, including figures and bibliography. With a standard TeX Live or MiKTeX installation providing AMS packages and TikZ, run `latexmk -pdf paper.tex` inside the selected paper directory. For the repeated-factor paper, run `sh scripts/build.sh` in `papers/repeated_factor_spectra`; `main.tex` includes its section sources, bibliography and generated figure. The checked-in PDFs are the released versions; their checksums are recorded in `manuscripts.json`.

The [Geometry & moduli research guide](https://bwuebben.github.io/research/geometry/) explains the programme and its relation to the other papers.

Copyright © 2026 Bernd Johannes Wuebben. Existing author rights are reserved; this repository does not apply a blanket code license to the manuscripts.

## Archived manuscripts

- Uniform Hermitian–Einstein spectra for bipartite two-step extensions: [doi:10.5281/zenodo.23717660](https://doi.org/10.5281/zenodo.23717660).
- Nonsimplicial local stability chambers of toric tangent bundles: [doi:10.5281/zenodo.23774770](https://doi.org/10.5281/zenodo.23774770).
- Metric response and Hermitian Yang–Mills flow for two-step extensions: [doi:10.5281/zenodo.23717742](https://doi.org/10.5281/zenodo.23717742).
- Uniform Hermitian–Einstein spectra with repeated stable factors: [doi:10.5281/zenodo.23780241](https://doi.org/10.5281/zenodo.23780241).

The Zenodo records preserve these released PDFs. arXiv submissions remain pending.

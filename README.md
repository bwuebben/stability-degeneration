# Stability walls and geometric degeneration

Ten papers by [Bernd Johannes Wuebben](https://bwuebben.github.io/) connect stability chambers, small Dolbeault eigenvalues, Hermitian–Einstein metrics, local moduli geometry and Hermitian–Yang–Mills flow.

| No. | Manuscript | Status | Resources |
| ---: | --- | --- | --- |
| 01 | [Uniform Hermitian–Einstein spectra for bipartite two-step extensions](papers/01_two_step_spectra/paper.pdf) | Preprint | [Abstract and citation](https://bwuebben.github.io/papers/uniform-hermitian-einstein-spectra/) · [Source](papers/01_two_step_spectra/) · [doi:10.5281/zenodo.23717660](https://doi.org/10.5281/zenodo.23717660) |
| 02 | [Nonsimplicial local stability chambers of toric tangent bundles](papers/02_nonsimplicial_chambers/paper.pdf) | Preprint | [Abstract and citation](https://bwuebben.github.io/papers/nonsimplicial-toric-stability-chambers/) · [Source](papers/02_nonsimplicial_chambers/) · [doi:10.5281/zenodo.23774770](https://doi.org/10.5281/zenodo.23774770) |
| 03 | [Metric response and Hermitian Yang–Mills flow for two-step extensions](papers/03_metric_response_flow/paper.pdf) | Preprint | [Abstract and citation](https://bwuebben.github.io/papers/metric-response-hym-flow/) · [Source](papers/03_metric_response_flow/) · [doi:10.5281/zenodo.23717742](https://doi.org/10.5281/zenodo.23717742) |
| 04 | [Uniform Hermitian–Einstein spectra with repeated stable factors](papers/04_repeated_factor_spectra/main.pdf) | Preprint | [Abstract and citation](https://bwuebben.github.io/papers/repeated-factor-he-spectra/) · [Source](papers/04_repeated_factor_spectra/) · [doi:10.5281/zenodo.23780241](https://doi.org/10.5281/zenodo.23780241) |
| 05 | [Hermitian–Einstein metrics and small spectra for three-step extensions](papers/05_three_step_metrics/paper.pdf) | Preprint | [Abstract and citation](https://bwuebben.github.io/papers/three-step-he-metrics/) · [Source](papers/05_three_step_metrics/) · [doi:10.5281/zenodo.23850990](https://doi.org/10.5281/zenodo.23850990) |
| 06 | [L² metrics on Hermitian–Yang–Mills moduli spaces near a stability wall](papers/06_wall_moduli_metrics/paper.pdf) | Preprint | [Abstract and citation](https://bwuebben.github.io/papers/wall-moduli-metrics/) · [Source](papers/06_wall_moduli_metrics/) · [doi:10.5281/zenodo.23850992](https://doi.org/10.5281/zenodo.23850992) |
| 07 | [Hermitian–Einstein metrics and moduli of three-step extensions with repeated factors](papers/07_repeated_factor_moduli/paper.pdf) | Preprint | [Abstract and citation](https://bwuebben.github.io/papers/repeated-factor-moduli/) · [Source](papers/07_repeated_factor_moduli/) · [doi:10.5281/zenodo.23850994](https://doi.org/10.5281/zenodo.23850994) |
| 08 | [Initial conditions and asymptotic phase for Hermitian–Yang–Mills flow](papers/08_flow_phase/paper.pdf) | Preprint | [Abstract and citation](https://bwuebben.github.io/papers/hym-initial-condition-phase/) · [Source](papers/08_flow_phase/) · [doi:10.5281/zenodo.23850996](https://doi.org/10.5281/zenodo.23850996) |
| 09 | [Hermitian–Yang–Mills flow at and beyond a stability wall](papers/09_wall_unstable_flow/paper.pdf) | Preprint | [Abstract and citation](https://bwuebben.github.io/papers/hym-wall-unstable-flow/) · [Source](papers/09_wall_unstable_flow/) · [doi:10.5281/zenodo.23850998](https://doi.org/10.5281/zenodo.23850998) |
| 10 | [Singular graded limits and metric growth under Hermitian–Yang–Mills flow](papers/10_singular_hym_asymptotics/paper.pdf) | Preprint | [Abstract and citation](https://bwuebben.github.io/papers/singular-hym-metric-growth/) · [Source](papers/10_singular_hym_asymptotics/) · [doi:10.5281/zenodo.23893608](https://doi.org/10.5281/zenodo.23893608) |

## Reading the programme

Directories 01–04 contain the two-step foundations and their repeated-factor extension. Directories 05–07 develop three-step metrics and local moduli geometry; 08–09 form the dynamical pair. Directory 10 proves a general singular-graded obstruction on compact Kähler surfaces. This numbering follows the development and mathematical dependencies of the programme and supplies a reading order; it is not part of the manuscript titles.

The bipartite spectral theorem and the toric chamber constructions give independent analytic and geometric starting points. The two-step metric-response paper joins them on the same tangent bundles. The repeated-factor spectral paper replaces scalar graph directions by matrix modes within its two-step class.

The three-step metric paper treats composable extensions with integrability relations and compares actual metrics and every small spectral scale with an independently balanced graph. The wall-moduli paper retains the full analytic obstruction and controls the harmonic L² metric on smooth pieces of the local quotient. The repeated-factor moduli paper develops the corresponding matrix estimates in a specified rank-five deformation problem. These are the stationary and moduli-geometric branches.

For dynamics, the initial-condition paper proves differentiable phase selection relative to a chosen invariant graph, with uniform exponential tracking near a stable family. The wall-and-unstable paper proves confinement in moving coefficients even when diagonal metrics diverge, and identifies a residual carrying initial conditions. The singular-graded paper proves a general rank-two surface obstruction: the squared norm on an equal-slope line decays strictly faster than a bounded-remainder extension of the categorical weight model would permit. The theorem includes nonreduced singularities and repeated line factors, with toric and K3 examples. Its stationary curvature theorem concerns the toric family.

These theorems retain their stated bundle, polarization and local-domain hypotheses. A general theory for arbitrary filtrations and stable factors remains open, as do repeated-factor dynamics and a replacement asymptotic law for singular graded sheaves. The surface obstruction does not contradict the weight-filtration theorem with locally free graded object.

Earlier work on [equilibrium spectra](https://bwuebben.github.io/papers/small-dolbeault-stability-walls/), [flow spectra](https://bwuebben.github.io/papers/small-dolbeault-eigenvalues/), [toric local chambers](https://bwuebben.github.io/papers/toric-local-stability-chambers/) and [relative metrics near reducible Yang–Mills connections](https://bwuebben.github.io/papers/relative-metrics-reducible-yang-mills/) provides further entry points. The [Stability walls and geometric degeneration guide](https://bwuebben.github.io/research/geometry/quantitative-moduli/) explains the reading routes.

## Build and reproduce

Each paper has a source and reproduction guide in its directory. The `paper.tex` files are assembled, self-contained manuscripts apart from their supplied figures. Compile with `latexmk -pdf paper.tex`. The earlier repeated-factor spectral paper uses `main.tex` and `sh scripts/build.sh`.

Run `python3 reproduce.py` from the repository root to regenerate the finite computations and figures. Programs with third-party dependencies use `uv` and the paper's pinned `requirements.txt`; the others use the Python standard library. The finite calculations establish the printed identities and examples, while the geometric, elliptic and parabolic assertions rely on the proofs. Released PDF hashes and archive identities are in `manuscripts.json`.

The Zenodo records preserve the released manuscripts. arXiv submissions remain pending.

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved; this repository applies no blanket code license to the manuscripts.

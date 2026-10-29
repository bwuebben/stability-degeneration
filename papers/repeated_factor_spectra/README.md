# Uniform Hermitian–Einstein spectra with repeated stable factors

By [Bernd Johannes Wuebben](https://bwuebben.github.io/).
Initial version: 29 October 2026. Preprint, 18 pages.

[Read the paper](main.pdf) ·
[Abstract and citation](https://bwuebben.github.io/papers/repeated-factor-he-spectra/) ·
[Programme and reading routes](https://bwuebben.github.io/research/geometry/#stability-degeneration).

For fixed simple locally free two-step extensions with repeated stable factors
within the layers, the paper compares the actual Hermitian–Einstein metric with
an independently balanced matrix model. It determines every small complex
endomorphism eigenvalue with uniform relative error and a fixed complementary
gap, allowing arbitrarily separated positive scales. A fixed rank-five bundle
on a product of four projective lines has six small modes and realizes the weak
scale explicitly. The result concerns stationary metrics; repeated-factor flow,
higher-step extensions and singular graded objects require additional theorems.

## Build

`main.tex` includes the seven section sources, `references.tex` and the figure.
With a standard AMS LaTeX installation and `latexmk`, run:

```sh
sh scripts/build.sh
```

## Reproduce

The dependency versions are pinned in `requirements.txt`. With `uv`:

```sh
uv run --no-project --with-requirements requirements.txt python figures/code/verify_star_model.py
uv run --no-project --with-requirements requirements.txt python figures/code/verify_geometry_independent.py
uv run --no-project --with-requirements requirements.txt python figures/code/plot_weak_modes.py --recompute
```

The two verifiers independently check the finite algebra and geometric model.
The plot uses the printed matrices; omit `--recompute` to redraw the supplied
data. Outputs are in `figures/data/` and `figures/output/`. These computations
illustrate and check the finite model; the metric and elliptic assertions are
proved in the manuscript. No HE partial differential equation is solved
numerically.

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved.

Archived version: [doi:10.5281/zenodo.23780241](https://doi.org/10.5281/zenodo.23780241).

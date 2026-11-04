# Exact example calculations

From the Paper12 directory run:

```bash
uv run reproducibility/check_examples.py
```

The inline dependency is SymPy 1.14.0. The script checks rational-function
orders at b=0, so with b=T^s its affine valuation comparisons hold for
every 0<s<1, not a sample of values. It verifies 74 nonzero pairwise minor
candidates and 14 triple determinant candidates, all six oriented pairwise
spectra, D111=6-s, and the six vertices/six finite edges/one bounded face
of the normalized diamond. Exact rational coordinates produce area 1/2
in (x/p,sqrt(3)t/p), hence p^2/(2sqrt(3)) in the original plane.

The support figure is the hand-readable coordinate construction in
../figures/diamond-collapse.tex and builds with the manuscript. These finite
checks do not prove the analytic path reconstruction or its continuity.
The mathematical derivations are printed in Section8.

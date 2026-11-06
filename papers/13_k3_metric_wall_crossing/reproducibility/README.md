# Finite calculations

The polynomial calculation verifies the Segre product, five transverse intersections,
unit affine Jacobian ideals modulo 101, and intersection and Mukai arithmetic in
Section 2. The figure evaluates the distance integral in Section 3 with 256- and
512-point Gauss–Legendre quadrature. These finite calculations do not replace the
geometric or analytic proofs.

From this directory, run:

    uv run --no-project k3_model.py
    uv run --no-project plot_distance_profile.py

The first command writes k3-model.json. The second writes distance-profile.json
and regenerates the figure in ../figures. Timing and PDF creation metadata may vary.

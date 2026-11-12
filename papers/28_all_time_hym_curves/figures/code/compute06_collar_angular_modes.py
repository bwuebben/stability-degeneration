"""Exact local checks for R33; not a proof of its global nonlinear theorem.

Run from the project root with uv run --with sympy python <this path>.
The noncommuting matrix check reconstructs the two polar Chern derivatives.
The scalar example checks both Neumann boundary fluxes and exact rho scaling.
"""

import json
from pathlib import Path

import sympy as sp


t, theta = sp.symbols("t theta", real=True)
p = sp.symbols("p", real=True, nonzero=True)
c = sp.Function("chi")(t)
T = sp.Matrix([[0, 1], [1, 0]])
H = sp.diag(sp.exp(p * t), sp.exp(-p * t))
C = H.inv() * H.diff(t)
# The derivative at alpha=0 of exp(-i alpha T F) H exp(i alpha T F),
# multiplied by H^{-1}; chi models the radial cutoff in the first variation.
M = sp.I * c * (T - H.inv() * T * H)


def comm(a, b):
    return a * b - b * a


matrix_checks = []
for n in (-1, 1):
    W = sp.exp(t + sp.I * n * theta) * M
    # A_z = C/(2z); first form delta A_z, then take 4 partial_bar_z.
    z = sp.exp(t + sp.I * theta)
    dz_W = (W.diff(t) - sp.I * W.diff(theta)) / (2 * z)
    delta_A = dz_W + comm(C / (2 * z), W)
    four_dbar = 2 * sp.exp(-t + sp.I * theta) * (
        delta_A.diff(t) + sp.I * delta_A.diff(theta)
    )
    w = sp.exp(t) * M
    expected = sp.exp(-2 * t + sp.I * n * theta) * (
        w.diff(t, 2) - n * n * w + comm(C.diff(t), w)
        + comm(C, w.diff(t)) - n * comm(C, w)
    )
    assert all(sp.simplify(v) == 0 for v in four_dbar - expected)
    coefficient = sp.simplify(sp.exp(-sp.I * n * theta) * four_dbar)
    assert all(sp.simplify(v.diff(theta)) == 0 for v in coefficient)
    assert comm(C, M) != sp.zeros(2)
    matrix_checks.append({"mode": n, "polar_identity": True,
                          "coefficient_independent_of_angle": True,
                          "noncommuting_background": True})

x = sp.symbols("x", positive=True)
y = x - 1
chi = 6 * y**5 - 15 * y**4 + 10 * y**3
# u(r,theta)=rho*f(r/rho)*cos(theta).  The source is
# -Delta(chi(r/rho)*r*cos(theta)); harmonic additions set both normal fluxes.
f = sp.expand(chi * x - sp.Rational(4, 3) * x - sp.Rational(4, 3) / x)
radial_laplacian = lambda g: sp.diff(g, x, 2) + sp.diff(g, x) / x - g / x**2
assert sp.simplify(radial_laplacian(f) - radial_laplacian(chi * x)) == 0
assert sp.diff(f, x).subs(x, 1) == 0
assert sp.diff(f, x).subs(x, 2) == 0
assert sp.integrate(sp.cos(theta), (theta, 0, 2 * sp.pi)) == 0
energy = sp.simplify(sp.integrate(x * sp.diff(f, x)**2 + f**2 / x, (x, 1, 2)))
source_norm = sp.simplify(sp.integrate(x * radial_laplacian(f)**2, (x, 1, 2)))
assert energy > 0 and source_norm > 0

result = {
    "scope": "Exact local matrix/Fourier and scalar Neumann checks only; no global or asymptotic PDE certification.",
    "matrix_checks": matrix_checks,
    "neumann_example": {
        "cutoff_in_scaled_radius": str(sp.expand(chi)),
        "radial_solution_f": str(f),
        "source_equals_minus_laplacian_of_cutoff_first_mode": True,
        "inner_normal_derivative_zero": True,
        "outer_normal_derivative_zero": True,
        "source_integral_zero_by_first_angular_mode": True,
        "dirichlet_energy_divided_by_pi_rho_squared": str(energy),
        "source_L2_squared_divided_by_pi": str(source_norm),
    },
    "result": "PASS",
}
out = Path(__file__).with_name("collar_angular_modes.json")
out.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))

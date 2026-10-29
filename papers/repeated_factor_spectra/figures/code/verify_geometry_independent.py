"""Independent exact audit of geometric_star.md (no imports from author code).

Portable edition: run from this manuscript directory.
uv run --no-project --with sympy==1.14.0 python \
  figures/code/verify_geometry_independent.py

Uses the unsquared arrow map in a complex Pauli basis, square-free intersection
arithmetic, and an explicit two-source boundary solution. Writes figures/data/geometry_independent_results.json. No PDE claims.
"""

from pathlib import Path
import hashlib
import json
import sys

import sympy as s


HERE = Path(__file__).resolve().parent
ROUND = HERE.parent
I = s.eye(2)
R = s.Rational
lam = s.symbols("lambda")


def clean(M):
    return M.applyfunc(s.simplify)


def squarefree_product(p, q):
    result = {}
    for a, x in p.items():
        for b, y in q.items():
            if a & b == 0:
                result[a | b] = result.get(a | b, 0) + x * y
    return {key: s.expand(value) for key, value in result.items()}


def p1_cohomology(k):
    return max(k + 1, 0), max(-k - 1, 0)


def product_h1(degrees):
    hs = [p1_cohomology(k) for k in degrees]
    return sum(hs[j][1] * s.prod(hs[k][0] for k in range(4) if k != j)
               for j in range(4))


h = s.symbols("h1:5", positive=True)
omega = {1 << j: h[j] for j in range(4)}
powers = [{0: s.Integer(1)}]
for _ in range(4):
    powers.append(squarefree_product(powers[-1], omega))
volume = s.factor(powers[4][15] / s.factorial(4))
assert volume == s.prod(h)
degrees = [(0, 0, 0, 0)]
degrees.extend(tuple(2 * (k == j) - 2 * (k == 3) for k in range(4))
               for j in range(3))
curvatures = []
for d in degrees:
    c1 = {1 << j: d[j] for j in range(4) if d[j]}
    degree = squarefree_product(c1, powers[3]).get(15, 0) / s.factorial(3)
    curvatures.append(s.simplify(2 * s.pi * degree / volume))
assert curvatures[0] == 0
for j in range(3):
    assert s.simplify(curvatures[j + 1] - 4 * s.pi * (1 / h[j] - 1 / h[3])) == 0
ext_dimensions = [product_h1(tuple(-x for x in d)) for d in degrees[1:]]
assert ext_dimensions == [3, 3, 3]
a = s.symbols("a1:4", real=True)
A = sum(a)
chart = {h[3]: 1, **{h[j]: 1 / (1 + (A / 2 + a[j]) / (4 * s.pi))
                          for j in range(3)}}
kappas = [s.simplify(k.subs(chart)) for k in curvatures]
kE = s.simplify(sum(kappas) / 5)
assert kE == A / 2
demands = s.Matrix([s.simplify(k - kE) for k in kappas[1:]])
assert demands == s.Matrix(a)
assert demands.jacobian(a) == s.eye(3)

# The imaginary Pauli matrix is included: this is a basis of all M_2(C).
pauli = [I / s.sqrt(2), s.diag(1, -1) / s.sqrt(2),
         s.Matrix([[0, 1], [1, 0]]) / s.sqrt(2),
         s.Matrix([[0, -s.I], [s.I, 0]]) / s.sqrt(2)]
basis = [(p, (0, 0, 0)) for p in pauli]
basis.extend((s.zeros(2), tuple(int(j == k) for j in range(3))) for k in range(3))
Z = [s.Matrix([[s.sqrt(R(3, 2)), 0, 0], [0, 1 / s.sqrt(2), 0]]),
     s.Matrix([[1, 0, 0], [R(1, 2), s.sqrt(3) / 2, 0]]),
     s.Matrix([[1 / s.sqrt(2), 0, 0], [-1 / s.sqrt(2), 1, 0]])]
B = [clean(z * z.adjoint()) for z in Z]
assert sum(B, s.zeros(2)) == 3 * I
assert [z.rank() for z in Z] == [2, 2, 2]
assert [b.det() for b in B] == [R(3, 4), R(3, 4), R(1, 2)]
assert B[0] * B[1] - B[1] * B[0] == s.Matrix([[0, R(1, 2)], [-R(1, 2), 0]])
columns = []
for target, sources in basis:
    column = []
    for j, z in enumerate(Z):
        column.extend((target - sources[j] * I) * z)
    columns.append(s.Matrix(column))
D = s.Matrix.hstack(*columns)
M = clean(D.adjoint() * D)
identity = s.Matrix([s.sqrt(2), 0, 0, 0, 1, 1, 1])
assert M * identity == s.zeros(7, 1)
assert D.rank() == 6
assert M == M.adjoint()
coupling = -M[:4, 4:]
assert M[:4, :4] == 3 * s.eye(4)
assert M[4:, 4:] == 2 * s.eye(3)
assert set((coupling.adjoint() * coupling).eigenvals()) == {6, R(3, 2), R(1, 2)}
expected = (lam - 5) * (lam - 3) * (lam**2 - 5 * lam + R(9, 2)) * (lam**2 - 5 * lam + R(11, 2))
full_polynomial = s.factor(M.charpoly(lam).as_expr())
assert s.expand(full_polynomial - lam * expected) == 0
scalar = M.extract([0, 4, 5, 6], [0, 4, 5, 6])
scalar_polynomial = s.factor(scalar.charpoly(lam).as_expr())
assert scalar_polynomial == lam * (lam - 5) * (lam - 2)**2

# Exact endpoint with a=(0,1,1). This is an auxiliary balanced limit,
# not a claim that the original nonsplit boundary bundle has an HE metric.
B2, B3 = B[1:]
S = B2 / s.sqrt(B2.det()) + B3 / s.sqrt(B3.det())
dS = s.simplify(S.det())
assert s.simplify(dS - (2 + 5 * s.sqrt(6) / 3)) == 0
H = clean(s.sqrt(dS) * S.inv())
assert s.simplify(H.det() - 1) == 0
traces = [s.simplify(s.trace(H * b)) for b in (B2, B3)]
for trace, b in zip(traces, (B2, B3)):
    assert s.simplify(trace**2 - dS * b.det()) == 0
assert clean(B2 / traces[0] + B3 / traces[1] - H.inv()) == s.zeros(2)
rho2 = s.simplify(1 - 4 / dS)
assert s.simplify(rho2 - (5 * s.sqrt(6) - 6) / (5 * s.sqrt(6) + 6)) == 0
r = s.symbols("rho", positive=True)
# A unitary change diagonalizes C2, C3=I-C2; it preserves the full spectrum.
boundary_covariances = [s.zeros(2), s.diag(1 + r, 1 - r) / 2,
                        s.diag(1 - r, 1 + r) / 2]
def covariance_gram(covariances):
    return s.Matrix(7, 7, lambda u, v: s.simplify(sum(
        s.trace((basis[u][0] - basis[u][1][j] * I).adjoint()
                * (basis[v][0] - basis[v][1][j] * I) * covariances[j])
        for j in range(3))))
boundary = covariance_gram(boundary_covariances)
boundary_polynomial = s.factor(boundary.charpoly(lam).as_expr())
assert s.expand(boundary_polynomial - lam**2 * (lam - 2) * (lam - 1)**2
                * ((lam - 1)**2 - r**2)) == 0
# The new kernel, orthogonal to the global identity, has squared mass 20.
weak = s.Matrix([s.sqrt(2), 0, 0, 0, -4, 1, 1])
assert (identity.adjoint() * weak)[0] == 0
assert boundary * weak == s.zeros(7, 1)
weak_mass = (weak.adjoint() * weak)[0]
assert weak_mass == 20
u, v, w = s.symbols("u v w", real=True)
unit_trace = s.Matrix([[u, v + s.I * w], [v - s.I * w, 1 - u]])
arbitrary_weak_derivative = covariance_gram([unit_trace, s.zeros(2), s.zeros(2)])
weak_energy_derivative = s.simplify((weak.adjoint() * arbitrary_weak_derivative * weak)[0])
assert weak_energy_derivative == 25
weak_coefficient = weak_energy_derivative / weak_mass
assert weak_coefficient == R(5, 4)

pins = {Path(__file__).name: {'sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}
result = {
    "scope": "Independent exact finite/geometric identities; no verification of elliptic or HE transfer arguments.",
    "versions": {"python": sys.version.split()[0], "sympy": s.__version__},
    "method": "Square-free intersection ring; raw 18-by-7 complex arrow map in Pauli basis; exact two-source endpoint; universal unit-trace weak perturbation.",
    "volume": str(volume), "curvatures": [str(k) for k in curvatures],
    "ext1_dimensions": ext_dimensions, "polarization_demands": list(map(str, demands)),
    "polarization_jacobian": str(demands.jacobian(a)),
    "covariances": list(map(str, B)), "operator_pauli_basis": str(M),
    "centralizer_dimension": 7 - D.rank(), "complex_trace_free_dimension": 6,
    "coupling_gram_eigenvalues": ["6", "3/2", "1/2"],
    "trace_free_characteristic_polynomial": str(s.factor(expected)),
    "scalar_compression_characteristic_polynomial_with_zero": str(scalar_polynomial),
    "boundary_H": str(H), "boundary_det_S": str(dS),
    "boundary_rho_squared": str(rho2), "boundary_gap": str(s.N(1 - s.sqrt(rho2), 30)),
    "boundary_trace_free_polynomial": str(s.factor(boundary_polynomial / lam)),
    "boundary_positive_eigenvalues": ["1-rho", "1", "1", "1+rho", "2"],
    "weak_kernel_mass": str(weak_mass), "weak_energy_derivative": str(weak_energy_derivative),
    "weak_eigenvalue_coefficient": str(weak_coefficient), "source_pins": pins,
    "all_exact_assertions_passed": True,
}
output = HERE.parent / "data" / "geometry_independent_results.json"
output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"output": str(output), "all_exact_assertions_passed": True,
                  "trace_free_characteristic_polynomial": result["trace_free_characteristic_polynomial"],
                  "boundary_gap": result["boundary_gap"], "weak_coefficient": str(weak_coefficient)}, indent=2))

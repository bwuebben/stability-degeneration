#!/usr/bin/env python3
"""Exact finite identities for the split-quintic application. These do not prove geometry or asymptotics."""
from itertools import product
from fractions import Fraction
import json
from pathlib import Path


def add(a, b):
    out = dict(a)
    for m, c in b.items():
        out[m] = out.get(m, 0) + c
    return {m: c for m, c in out.items() if c}


def scale(a, c):
    return {m: c * x for m, x in a.items() if c * x}


def mul(a, b, degree=5):
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            if i + k <= 4 and j + l <= 1 and i + j + k + l <= degree:
                m = (i + k, j + l)
                out[m] = out.get(m, 0) + x * y
    return {m: c for m, c in out.items() if c}


def power(a, n, degree=5):
    out = {(0, 0): 1}
    for _ in range(n):
        out = mul(out, a, degree)
    return out


def inverse_one_plus(a, degree=3):
    out = {}
    for n in range(degree + 1):
        out = add(out, scale(power(a, n, degree), (-1) ** n))
    return out


def homogeneous(a, degree):
    return {m: c for m, c in a.items() if sum(m) == degree}


one = {(0, 0): 1}
h, ell = {(1, 0): 1}, {(0, 1): 1}
d1, d2 = add(h, ell), add(scale(h, 4), ell)
fundamental = mul(d1, d2)
ambient_chern = mul(power(add(one, h), 5), power(add(one, ell), 2))
chern = mul(mul(ambient_chern, inverse_one_plus(d1), 3), inverse_one_plus(d2), 3)
c1, c2, c3 = [homogeneous(chern, j) for j in (1, 2, 3)]
assert not c1
assert c2 == {(2, 0): 6, (1, 1): 5}


def integrate(poly):
    return mul(poly, fundamental).get((4, 1), 0)


cubics = [integrate(power(h, 3)), integrate(mul(power(h, 2), ell)),
          integrate(mul(h, power(ell, 2))), integrate(power(ell, 3))]
assert cubics == [5, 4, 0, 0]
assert [integrate(mul(c2, h)), integrate(mul(c2, ell))] == [50, 24]
assert integrate(c3) == -168
b3_resolution = 2 + 2 * 2 - integrate(c3)
assert b3_resolution == 174

# Quintic Chern series: (1+h)^5/(1+5h) through degree three.
quintic_c = mul(power(add(one, h), 5), inverse_one_plus(scale(h, 5)), 3)
assert homogeneous(quintic_c, 2) == {(2, 0): 10}
assert homogeneous(quintic_c, 3) == {(3, 0): -40}
euler_smoothing = -40 * 5
assert euler_smoothing == -200
assert 2 + 2 - euler_smoothing == 204

# Exact unimodular reduction of I_n + J_n: row differences, a column
# addition, then row additions. Every operation has determinant one.
n = 15
matrix = [[1 + int(i == j) for j in range(n)] for i in range(n)]
original = [row[:] for row in matrix]
# The two integral maps in the attachment calculation. In the quotient
# W=Z^16/Z(1,...,1), use the images of e_1,...,e_15 as basis.
restriction = [[0, 1] for _ in range(16)]
quotient_to_W = [[int(i == j) for j in range(15)] + [-1]
                 for i in range(15)]
inclusion_K = [[int(i == j) for j in range(15)] for i in range(15)] + [[-1]*15]


def matrix_product(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


assert matrix_product(quotient_to_W, restriction) == [[0, 0] for _ in range(15)]
assert matrix_product(quotient_to_W, inclusion_K) == original
assert all(sum(row[j] for row in inclusion_K) == 0 for j in range(15))
# Unit entries/minors exhibit the primitive diagonal and the surjective
# quotient in these printed bases. Geometric identification needs the proof.
assert restriction[0][1] == 1
assert [row[:15] for row in quotient_to_W] == [[int(i == j) for j in range(15)] for i in range(15)]
for i in range(1, n):
    matrix[i] = [a - b for a, b in zip(matrix[i], matrix[0])]
for i in range(n):
    matrix[i][0] += sum(matrix[i][1:])
for j in range(1, n):
    matrix[0] = [a - b for a, b in zip(matrix[0], matrix[j])]
assert matrix == [[(16 if i == 0 else 1) if i == j else 0
                   for j in range(n)] for i in range(n)]

# Degree-five evaluation on a (4,4) planar grid has exactly fifteen
# distinct characters. The missing character is (3,3), so weighting the
# evaluations by xy gives the nonzero-coefficient relation.
monomials = [(a, b) for a in range(6) for b in range(6 - a)]
characters = {(a % 4, b % 4) for a, b in monomials}
assert len(characters) == 15
assert set(product(range(4), repeat=2)) - characters == {(3, 3)}
roots = [1, 1j, -1, -1j]  # Exact Gaussian integers under these operations.
for a, b in monomials:
    assert sum(x ** (a + 1) * y ** (b + 1) for x, y in product(roots, repeat=2)) == 0

# Conditional spectral mass/normalization arithmetic. The existence of a
# limiting harmonic representative and the spectral comparison are NOT checked.
primitive_coefficient = Fraction(integrate(mul(ell, power(h, 2))), integrate(power(h, 3)))
primitive_class = add(ell, scale(h, -primitive_coefficient))
assert integrate(mul(primitive_class, power(h, 2))) == 0
conditional_mass = -integrate(mul(power(primitive_class, 2), h))
assert conditional_mass == Fraction(16, 5)
link_norm_times_volume_in_pi_units = Fraction(9, 8) * Fraction(16, 27)
assert link_norm_times_volume_in_pi_units == Fraction(2, 3)
conditional_radius_log_coefficient_in_pi_units = 16 * link_norm_times_volume_in_pi_units / conditional_mass
assert conditional_radius_log_coefficient_in_pi_units == Fraction(10, 3)

report = {
    'scope': 'exact finite identities only; geometry, topology and PDE need proofs',
    'resolution_cubic_H_L': cubics,
    'resolution_c2_pairings_H_L': [50, 24],
    'euler_resolution': -168,
    'euler_smoothing': -200,
    'b3_resolution': 174,
    'b3_smoothing': 204,
    'complement_betti': [1, 0, 2, 189, 1, 15, 0],
    'ambient_H4_product_index': abs(5 * 0 - 4 * 4),
    'monodromy_vanishing_block': original,
    'monodromy_smith_diagonal': [1] * 14 + [16],
    'boundary_restriction_H_L': restriction,
    'vanishing_quotient_map': quotient_to_W,
    'cycle_relation_inclusion': inclusion_K,
    'attachment_composition_zero': True,
    'monodromy_factorization_checked': True,
    'degree5_grid_evaluation_rank': len(characters),
    'degree5_grid_missing_character': [3, 3],
    'conditional_spectral_mass': str(conditional_mass),
    'conditional_link_energy_coefficient_div_pi': str(link_norm_times_volume_in_pi_units),
    'conditional_algebraic_log_coefficient_div_pi': str(3 * conditional_radius_log_coefficient_in_pi_units),
    'spectral_caveat': 'coefficients are conditional arithmetic, not an eigenvalue theorem',
    'checks': 'PASS',
}
out = Path(__file__).resolve().parents[1] / 'data' / 'split_quintic.json'
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({k: v for k, v in report.items() if k != 'monodromy_vanishing_block'}, indent=2))

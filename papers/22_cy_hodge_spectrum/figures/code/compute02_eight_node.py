#!/usr/bin/env python3
"""Exact intersection and quotient-mass checks; not a geometric or PDE proof."""
from fractions import Fraction as F
from pathlib import Path
import json


def add(a, b):
    out = dict(a)
    for exp, value in b.items():
        out[exp] = out.get(exp, 0) + value
    return {e: v for e, v in out.items() if v}


def mul(a, b):
    out = {}
    for ea, va in a.items():
        for eb, vb in b.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            out[e] = out.get(e, 0) + va * vb
    return out


def integral_monomial(exp):
    # f^2=g^2=0, xi_i^2=(f-g)xi_i, integral(f*g*xi1*xi2)=1.
    if exp[0] > 1 or exp[1] > 1:
        return 0
    for i in (2, 3):
        if exp[i] > 1:
            a, b = list(exp), list(exp)
            a[i] -= 1
            b[i] -= 1
            a[0] += 1
            b[1] += 1
            return integral_monomial(tuple(a)) - integral_monomial(tuple(b))
    return int(exp == (1, 1, 1, 1))


def integral(poly):
    return sum(v * integral_monomial(e) for e, v in poly.items())


def transpose(a):
    return list(map(list, zip(*a)))


def product(a, b):
    return [[sum(F(x) * F(y) for x, y in zip(row, col))
             for col in transpose(b)] for row in a]


def minus(a, b):
    return [[x-y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def inverse2(a):
    det = a[0][0]*a[1][1]-a[0][1]*a[1][0]
    return [[a[1][1]/det, -a[0][1]/det],
            [-a[1][0]/det, a[0][0]/det]]


basis = [{tuple(int(j == i) for j in range(4)): 1} for i in range(4)]
H = {}
for coefficient, poly in zip((0, 4, 2, 2), basis):
    H = add(H, {e: coefficient*v for e, v in poly.items()})
H2 = mul(H, H)
B = [[F(integral(mul(mul(di, dj), H2))) for dj in basis] for di in basis]
assert B == [[0, 8, 4, 4], [8, 0, 12, 12], [4, 12, -8, 0], [4, 12, 0, -8]]
assert integral(mul(H2, H2)) == 320

K = [[0, 0], [1, 1], [1, 0], [0, 1]]
D = [[1, 0], [0, 1], [0, 0], [0, 0]]
R = [[1, 0, 0, 0]]*4 + [[0, 1, -1, -1]]*4
BKK = product(product(transpose(K), B), K)
BKD = product(product(transpose(K), B), D)
P = minus(D, product(K, product(inverse2(BKK), BKD)))
assert product(R, K) == [[0, 0]]*8
assert product(R, P) == product(R, D)
assert product(product(transpose(K), B), P) == [[0, 0], [0, 0]]
mass = [[-x for x in row] for row in product(product(transpose(P), B), P)]
periods = product(R, P)
energy = product(transpose(periods), periods)
assert mass == [[F(36, 5), F(-4, 5)], [F(-4, 5), F(36, 5)]]
assert energy == [[4, 0], [0, 4]]

modes = []
for vec, mu in [([[1], [-1]], F(1, 2)), ([[1], [1]], F(5, 8))]:
    norm = product(product(transpose(vec), mass), vec)[0][0]
    period_norm = product(product(transpose(vec), energy), vec)[0][0]
    assert product(energy, vec) == [[mu*x[0]] for x in product(mass, vec)]
    assert period_norm/norm == mu
    modes.append({"quotient_vector": vec, "squared_mass": str(norm),
                  "squared_period_norm": str(period_norm), "mu": str(mu),
                  "coefficient_in_units_of_pi": str(2*mu)})

report = {
    "status": "PASS: exact finite arithmetic only",
    "limit": "No smoothing-existence, metric, PDE, integral-saturation or asymptotic proof.",
    "source": "Wuebben, Nonsimplicial local stability chambers of toric tangent bundles, Appendix A, Proposition A.1",
    "source_sha256_at_verification": "e5078553511b2a06d65397130e1701e74e84793082b0aa29e4fa432069046727",
    "polarization": [0, 4, 2, 2],
    "volume_intersection": 320,
    "intersection_matrix": [[str(x) for x in row] for row in B],
    "quotient_mass": [[str(x) for x in row] for row in mass],
    "period_energy": [[str(x) for x in row] for row in energy],
    "modes": modes,
}
out = Path(__file__).resolve().parents[1]/"data/eight_node.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(report, indent=2)+"\n")
print(json.dumps(report, indent=2))

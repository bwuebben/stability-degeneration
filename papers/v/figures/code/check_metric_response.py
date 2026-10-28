#!/usr/bin/env python3
"""Exact rational block algebra only; no elliptic or actual HE verification.

The printed model uses K_2,2, ranks (1,2,3,5), and orbit-compatible weights
(q^2,q^(2+2h),q^(2+4h),q^(2+6h)).  Its low/high symmetric form is
[[G,B^T],[B,I]], with G=(1+2q^2)Q and B^T B=q^2 Q.  All identities use
the nontrivial mass matrix, not an accidental Euclidean-coordinate inverse.
"""

from fractions import Fraction as F
import json
from pathlib import Path
import time


def zero(n, m):
    return [[F(0) for _ in range(m)] for _ in range(n)]


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def trans(a):
    return [list(row) for row in zip(*a)]


def mul(a, b):
    bt = trans(b)
    return [[sum((x * y for x, y in zip(row, col)), F(0))
             for col in bt] for row in a]


def add(a, b):
    return [[x + y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def scale(a, s):
    return [[s * x for x in row] for row in a]


def inverse(a):
    n = len(a)
    w = [list(row) + e for row, e in zip(a, eye(n))]
    for j in range(n):
        pivot = next(i for i in range(j, n) if w[i][j])
        w[j], w[pivot] = w[pivot], w[j]
        v = w[j][j]
        w[j] = [x / v for x in w[j]]
        for i in range(n):
            if i != j:
                v = w[i][j]
                w[i] = [x - v * y for x, y in zip(w[i], w[j])]
    return [row[n:] for row in w]


def positive_definite(a):
    assert a == trans(a)
    w = [row[:] for row in a]
    for j in range(len(a)):
        pivot = w[j][j]
        if pivot <= 0:
            return False
        for i in range(j + 1, len(a)):
            for k in range(j + 1, len(a)):
                w[i][k] -= w[i][j] * w[j][k] / pivot
    return True


def diagonal(values):
    return [[x if i == j else F(0) for j in range(len(values))]
            for i, x in enumerate(values)]


def block(a, b, c, d):
    return [ra + rb for ra, rb in zip(a, b)] + [rc + rd for rc, rd in zip(c, d)]


def verify_model(q, h):
    ranks = [F(1), F(2), F(3), F(5)]
    incidence = [[F(0) for _ in range(4)] for _ in range(4)]
    for e, (target, source) in enumerate([(0, 2), (0, 3), (1, 2), (1, 3)]):
        incidence[target][e] = F(1)
        incidence[source][e] = F(-1)
    exponents = [2, 2 + 2 * h, 2 + 4 * h, 2 + 6 * h]
    weights = [q ** e for e in exponents]
    square_roots = [q ** (e // 2) for e in exponents]
    sigma = sum(weights, F(0))
    assert weights[0] * weights[3] == weights[1] * weights[2]
    z = [[F(i == j) / ranks[i] for j in range(3)] for i in range(3)]
    z += [[-F(1) / ranks[3] for _ in range(3)]]
    assert mul([ranks], z) == zero(1, 3)
    mass = mul(mul(trans(z), diagonal(ranks)), z)
    graph = mul(mul(incidence, diagonal(weights)), trans(incidence))
    form = mul(mul(trans(z), graph), z)
    coupling = scale(mul(mul(diagonal(square_roots), trans(incidence)), z), q)
    bt_b = mul(trans(coupling), coupling)
    assert bt_b == scale(form, q * q)
    g = scale(form, 1 + 2 * q * q)
    schur = add(g, scale(bt_b, -1))
    assert schur == scale(form, 1 + q * q)
    assert positive_definite(mass) and positive_definite(form)
    assert positive_definite(add(scale(form, sigma), scale(bt_b, -1)))
    assert positive_definite(add(scale(mass, 2 * sigma * sigma), scale(bt_b, -1)))
    hilbert_mass = block(mass, zero(3, 4), zero(4, 3), eye(4))
    bilinear = block(g, trans(coupling), coupling, eye(4))
    operator = mul(inverse(hilbert_mass), bilinear)
    assert positive_definite(bilinear)
    j = eye(3) + scale(coupling, -1)
    jstar = mul(mul(inverse(mass), trans(j)), hilbert_mass)
    schur_operator_inverse = mul(inverse(schur), mass)
    high_inverse = block(zero(3, 3), zero(3, 4), zero(4, 3), eye(4))
    decomposed = add(mul(mul(j, schur_operator_inverse), jstar), high_inverse)
    assert decomposed == inverse(operator)
    assert mul(operator, decomposed) == eye(7)
    for shift in [F(1, 1000), F(1), F(1000)]:
        cz_inverse = F(1) / (1 + shift)
        jz = eye(3) + scale(coupling, -cz_inverse)
        jzstar = mul(mul(inverse(mass), trans(jz)), hilbert_mass)
        sz = add(add(g, scale(mass, shift)), scale(bt_b, -cz_inverse))
        high_z = scale(high_inverse, cz_inverse)
        rhs = add(mul(mul(jz, mul(inverse(sz), mass)), jzstar), high_z)
        assert rhs == inverse(add(operator, scale(eye(7), shift)))
        assert positive_definite(add(sz, scale(add(form, scale(mass, shift)), -1)))
        upper = add(scale(form, 1 + 2 * sigma), scale(mass, shift))
        assert positive_definite(add(upper, scale(sz, -1)))
    # An ordinary operator perturbation of size sigma^2 can be much larger
    # than the energy of the cut isolating vertex 1 (rank 2).
    x = [-F(2, 11), F(9, 11), -F(2, 11), -F(2, 11)]
    xmass = sum((r * v * v for r, v in zip(ranks, x)), F(0))
    cut_energy = weights[2] + weights[3]
    false_relative_ratio = sigma * sigma * xmass / cut_energy
    assert false_relative_ratio > 1
    return {
        "q": str(q), "h": h, "weight_exponents": exponents,
        "weakest_to_strongest_weight": str(weights[3] / weights[0]),
        "cycle_orbit_identity": True,
        "mass_correct_inverse_identity": True,
        "positive_resolvent_identities": 3,
        "relative_schur_and_almost_isometry_bounds": True,
        "ordinary_sigma_squared_error_over_weak_cut_energy": str(false_relative_ratio),
    }


def two_factor_controls():
    count = 0
    for r1 in range(1, 5):
        for r2 in range(1, 5):
            for b in [F(1, 10), F(1, 10**20)]:
                dx1 = F(r2, r1 + r2) / b
                dx2 = -F(r1, r1 + r2) / b
                output_squared = r1 * dx1 * dx1 + r2 * dx2 * dx2
                input_squared = F(1, r1) + F(1, r2)
                eigenvalue = b * input_squared
                assert output_squared / input_squared == 1 / (eigenvalue * eigenvalue)
                count += 1
    return count


def main():
    started = time.monotonic()
    controls = two_factor_controls()
    records = [verify_model(q, h) for q in [F(1, 10), F(1, 100)]
               for h in [1, 3, 7, 15]]
    elapsed = time.monotonic() - started
    assert elapsed < 60, "Finite computation exceeded its 60-second cap"
    output = {
        "scope": "Exact rational finite block algebra; not geometric verification",
        "ranks": [1, 2, 3, 5], "graph": "K_2,2", "models": records,
        "two_factor_exact_normalization_cases": controls,
        "all_checks_passed": True,
        "limits": ["Does not test elliptic estimates or actual HE metrics",
                   "Does not prove harmonic gauge or parameter regularity"],
    }
    destination = (Path(__file__).resolve().parents[1] / "data/metric_response_checks.json")
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"all_checks_passed": True, "models": len(records),
                      "two_factor_cases": controls, "elapsed_seconds": round(elapsed, 3)}))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact finite data printed in Sections 1, 2 and 6; no PDE simulation."""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
import json
import math


def laplacian(n, edges):
    a = [[F(0) for _ in range(n)] for _ in range(n)]
    for u, v, weight in edges:
        a[u][u] += weight
        a[v][v] += weight
        a[u][v] -= weight
        a[v][u] -= weight
    return a


def action(a, x):
    return [sum((v * y for v, y in zip(row, x)), F(0)) for row in a]


def determinant(a):
    n = len(a)
    result = F(0)
    for p in permutations(range(n)):
        sign = (-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        result += sign * math.prod(a[i][p[i]] for i in range(n))
    return result


def poly_mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def determinant_shift(a, masses):
    """Coefficients of det(a + z diag(masses)), computed independently."""
    n = len(a)
    out = [F(0)] * (n + 1)
    for p in permutations(range(n)):
        term = [F((-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n)))]
        for i in range(n):
            term = poly_mul(term, [a[i][p[i]], F(masses[i]) if p[i] == i else F(0)])
        for i, value in enumerate(term):
            out[i] += value
    return out


def fraction_list(x):
    return [str(v) for v in x]


def main():
    records = {"scope": "exact finite algebra; does not verify HE equations",
               "path_controls": [], "four_cycles": []}
    for q in (F(1, 10), F(1, 100), F(1, 10**8)):
        a = laplacian(3, [(0, 1, q), (1, 2, q**3)])
        potential = [-1 / (3*q), -1 / (3*q), 2 / (3*q)]
        forcing = [F(0), -q*q, q*q]
        assert action(a, potential) == forcing
        assert max(potential) - min(potential) == 1 / q
        assert q*q <= (q + q**3)**2
        records["path_controls"].append({
            "q": str(q), "forcing": fraction_list(forcing),
            "potential": fraction_list(potential), "oscillation": str(1/q),
            "weak_cut_capacity": str(q**3)})

    vectors = ([1, -1, -1, 1], [1, -1, 1, -1], [1, 1, -1, -1])
    for rho in (F(1), F(1, 7)):
        for sqrt_theta in (F(1, 2), F(1), F(3)):
            alpha = sqrt_theta / (1 + sqrt_theta)
            weights = [rho*alpha, rho*(1-alpha), rho*(1-alpha), rho*alpha]
            matrix = laplacian(4, [(0, 2, weights[0]), (0, 3, weights[1]),
                                    (1, 2, weights[2]), (1, 3, weights[3])])
            eigenvalues = [2*rho*alpha, 2*rho*(1-alpha), 2*rho]
            for vector, eigenvalue in zip(vectors, eigenvalues):
                assert action(matrix, vector) == [eigenvalue*v for v in vector]
            assert action(matrix, [1, 1, 1, 1]) == [0, 0, 0, 0]
            cofactor = determinant([row[1:] for row in matrix[1:]])
            assert math.prod(eigenvalues) == 4*cofactor
            assert weights[0]*weights[3] / (weights[1]*weights[2]) == sqrt_theta**2
            records["four_cycles"].append({
                "rho": str(rho), "theta": str(sqrt_theta**2),
                "weights": fraction_list(weights),
                "eigenvalues_in_display_order": fraction_list(eigenvalues),
                "cofactor": str(cofactor)})

    unit_cycle = laplacian(4, [(0, 2, F(1, 2)), (0, 3, F(1, 2)),
                               (1, 2, F(1, 2)), (1, 3, F(1, 2))])
    masses = [1, 2, 3, 5]
    coefficients = determinant_shift(unit_cycle, masses)
    product = coefficients[1] / math.prod(masses)
    assert coefficients[0] == 0
    assert coefficients[-1] == math.prod(masses)
    assert product == F(11, 60)
    assert product != F(1, 60)  # Reject a missing total-rank factor.
    records["nonunit_mass"] = {
        "masses": masses, "det_L_plus_zM": fraction_list(coefficients),
        "positive_spectral_product": str(product),
        "missing_total_rank_factor_rejected": True}
    target = Path(__file__).resolve().parents[1] / "data" / "graph_values.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(records, indent=2) + "\n")
    print("PASS: exact path control, six four-cycle spectra, nonunit-rank product 11/60.")


if __name__ == "__main__":
    main()

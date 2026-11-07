#!/usr/bin/env python3
"""Exact finite identities only; this does not check any elliptic/PDE proof."""
from fractions import Fraction
import json


def degree(a, b, N):
    # c1 = a*w1+b*w2; mixed triple intersections on a bicubic are 3.
    return 3 * (a * (2 * N + 1) + b * (N * N + 2 * N))


def monomial_contains(gens, mon):
    return any(all(a <= b for a, b in zip(g, mon)) for g in gens)


def pair_eigenvalues(z):
    # Unit harmonic direction; mass 2*Vol, with Vol=1 for this finite check.
    return sorted((z[i] - z[j]) ** 2 / 2
                  for i in range(3) for j in range(3) if i != j)


def main():
    # Exact restriction of x^3+V*y^3+U*z^3 to x=a*y+b*z:
    # Exponent order: a,b,y,z. Multiply rather than assume the expansion.
    linear = {(1, 0, 1, 0): 1, (0, 1, 0, 1): 1}
    power = {(0, 0, 0, 0): 1}
    for _ in range(3):
        product = {}
        for ex, cx in power.items():
            for ey, cy in linear.items():
                mon = tuple(a + b for a, b in zip(ex, ey))
                product[mon] = product.get(mon, 0) + cx * cy
        power = product
    assert power == {(3, 0, 3, 0): 1, (2, 1, 2, 1): 3,
                     (1, 2, 1, 2): 3, (0, 3, 0, 3): 1}
    # Eliminating y^3,z^3 coefficients leaves these actual monomials.
    residual_generators = [(2, 1), (1, 2)]
    # Verify the monomial ideal identity via minimal generators of intersections.
    intersect = [(max(a[0], b[0]), max(a[1], b[1]))
                 for a in [(1, 1)] for b in [(3, 0), (2, 1), (1, 2), (0, 3)]]
    minimal = [g for g in intersect
               if not any(h != g and all(x <= y for x, y in zip(h, g))
                          for h in intersect)]
    assert sorted(minimal) == sorted(residual_generators)
    assert not any(monomial_contains(residual_generators, (i, 2 - i))
                   for i in range(3))
    assert monomial_contains(residual_generators, (2, 1))

    # Compare the coefficient polynomials, not only values at a sample N.
    # degree(a,b,N)=3*b*N^2+6*(a+b)*N+3*a.
    coeff = lambda a, b: (3 * b, 6 * (a + b), 3 * a)
    assert coeff(-2, 2) == (6, 0, -6)
    assert coeff(2, 0) == (0, 12, 6)  # E(2,-1): 6*(2*N+1).
    assert coeff(4, -2) == (-6, 12, 12)
    assert [degree(*c, 2) for c in [(-2, 2), (2, 0), (4, -2)]] == [18, 30, 12]

    # Arbitrarily separated rational scales: no floating-point eigenvalue loss.
    s, eps = Fraction(1, 10), Fraction(1, 10**30)
    vals = pair_eigenvalues([0, s, s + eps])
    assert vals == sorted([eps**2 / 2] * 2 + [s**2 / 2] * 2
                          + [(s + eps)**2 / 2] * 2)
    # Include two diagonal trace-free kernel directions.
    assert 2 + sum(v == 0 for v in pair_eigenvalues([0, s, s])) == 4
    assert 2 + sum(v == 0 for v in pair_eigenvalues([0, 0, 0])) == 8

    print(json.dumps({
        "finite_identities": "PASS",
        "hilbert_ideal": ["a^2*b", "a*b^2"],
        "degrees_at_N_2": [18, 30, 12],
        "parallel_tracefree_dimension": 8,
        "kernel_dimensions": [2, 4, 8],
        "limitations": "Exact polynomial and matrix-unit identities only; no PDE or geometric smoothness proof."
    }, indent=2))


if __name__ == "__main__":
    main()

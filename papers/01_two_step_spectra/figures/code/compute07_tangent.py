#!/usr/bin/env python3
"""Exact fourfold intersections and leading graph coefficients in Section 7."""
from fractions import Fraction as F
from pathlib import Path
import json


class Quadratic:
    """The field Q(sqrt(21)), with exact rational coordinates."""
    def __init__(self, a=0, b=0):
        if isinstance(a, Quadratic):
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = F(a), F(b)

    def __add__(self, other):
        other = Quadratic(other)
        return Quadratic(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Quadratic(-self.a, -self.b)

    def __sub__(self, other):
        return self + -Quadratic(other)

    def __mul__(self, other):
        other = Quadratic(other)
        return Quadratic(self.a*other.a + 21*self.b*other.b,
                         self.a*other.b + self.b*other.a)

    __rmul__ = __mul__

    def __eq__(self, other):
        other = Quadratic(other)
        return self.a == other.a and self.b == other.b

    def record(self):
        return {"rational": str(self.a), "sqrt21_coefficient": str(self.b)}


def ring_add(a, b):
    out = dict(a)
    for key, value in b.items():
        out[key] = out.get(key, Quadratic()) + value
    return out


def ring_mul(a, b):
    out = {}
    for p, ap in a.items():
        for q, bq in b.items():
            degrees = [x+y for x, y in zip(p, q)]
            value = ap*bq
            for i in (0, 1):
                while degrees[i] >= 2:
                    degrees[i] -= 2
                    degrees[2] += 1
                    degrees[3] += 1
                    value *= F(-1, 2)
            if degrees[2] > 1 or degrees[3] > 1:
                continue
            key = tuple(degrees)
            out[key] = out.get(key, Quadratic()) + value
    return out


def ring_power(a, n):
    result = {(0, 0, 0, 0): Quadratic(1)}
    for _ in range(n):
        result = ring_mul(result, a)
    return result


def integral(a):
    return a.get((1, 1, 1, 1), Quadratic())


ORDER = 9


def series(*coefficients):
    return [F(x) for x in coefficients] + [F(0)]*(ORDER-len(coefficients))


def add(a, b):
    return [x+y for x, y in zip(a, b)]


def scale(a, value):
    return [value*x for x in a]


def multiply(a, b):
    return [sum((a[j]*b[i-j] for j in range(i+1)), F(0)) for i in range(ORDER)]


def invert(a):
    assert a[0]
    out = series(1/a[0])
    for n in range(1, ORDER):
        out[n] = -sum((a[j]*out[n-j] for j in range(1, n+1)), F(0))/a[0]
    return out


def main():
    eta1, eta2, h1, h2 = [
        {tuple(int(i == j) for i in range(4)): Quadratic(1)}
        for j in range(4)]
    r = Quadratic(F(1, 2), F(1, 6))
    wall = ring_add(ring_add(eta1, eta2),
                    {key: r*value for key, value in ring_add(h1, h2).items()})
    volume4 = integral(ring_power(wall, 4))
    assert volume4 == Quadratic(16, 4)
    degrees = [integral(ring_mul(
        {key: 2*value for key, value in divisor.items()}, ring_power(wall, 3)))
        for divisor in (eta1, eta2, h1, h2)]
    assert all(x == Quadratic(6, 2) for x in degrees)
    boundary = ring_add(ring_add(eta1, eta2), ring_add(h1, h2))
    assert integral(ring_power(boundary, 4)) == Quadratic(20)
    # 3r^2 - 3r - 1 = 0, the equal-slope equation; sqrt(21)>3 gives r>1.
    assert 3*r*r - 3*r - 1 == 0

    # Factor out eps powers before formal division:
    # sigma=eps(1+eps^2), Sx=eps^2(1+eps^4),
    # Sy=eps^2(1-2eps+4eps^2-2eps^3+eps^4).
    one_plus = series(1, 0, 1)
    sx = series(1, 0, 0, 0, 1)
    sy = series(1, -2, 4, -2, 1)
    e2_bracket = add(
        add(scale(multiply(one_plus, one_plus), 4), scale(add(sx, sy), -1)),
        scale(multiply(multiply(sx, sy), invert(multiply(one_plus, one_plus))), -2))
    e2_series = scale(e2_bracket, F(1, 2))
    assert e2_series[0] == 0
    assert e2_series[1] == 3  # e2=eps^2(3eps+...).
    e3_factor = scale(multiply(series(1, -1, 1), invert(one_plus)), 4)
    assert e3_factor[0] == 4
    records = {
        "scope": "exact intersections and finite series, not a numerical HE computation",
        "wall_r": r.record(), "H0_power4": volume4.record(),
        "factor_raw_degrees": [d.record() for d in degrees],
        "boundary_H1_power4": 20,
        "e2_over_eps2_series": [str(x) for x in e2_series],
        "e3_over_eps6_series": [str(x) for x in e3_factor],
        "symmetric_leading_terms": {"e1": "2*eps", "e2": "3*eps^3", "e3": "4*eps^6"},
        "predicted_ordered_graph_eigenvalues": ["(4/3)*eps^3", "(3/2)*eps^2", "2*eps"],
        "actual_elliptic_transfer": "requires the proved analytic theorem"}
    target = Path(__file__).resolve().parents[1] / "data" / "tangent_values.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(records, indent=2) + "\n")
    print("PASS: Q(sqrt21) wall/intersections and finite three-scale coefficients.")


if __name__ == "__main__":
    main()

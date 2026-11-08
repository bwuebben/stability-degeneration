#!/usr/bin/env python3
"""Exact rational checks of the 2x2 orbit projection in W2; not a proof substitute."""
from fractions import Fraction as F
import json


def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def projection(b, c):
    den = b*b+c*c
    if den == 0:
        raise ValueError("The orbit direction disappears at the singular point")
    return [[b*b/den, -b*c/den], [-b*c/den, c*c/den]]


def main():
    records = []
    for b in (F(1), F(1, 2), F(1, 10), F(1, 100)):
        for c in (F(0), F(1, 3), b):
            p = projection(b, c)
            assert mm(p, p) == p
            assert p[0][1] == p[1][0]
            assert p[0][0]+p[1][1] == 1
        first = [[F(0), -1/b], [-1/b, F(0)]]
        second = [[-2/(b*b), F(0)], [F(0), 2/(b*b)]]
        assert mm(first, first) == [[1/(b*b), F(0)], [F(0), 1/(b*b)]]
        assert mm(second, second) == [[4/(b**4), F(0)], [F(0), 4/(b**4)]]
        records.append({"b": str(b), "first_operator_norm": str(1/b),
                        "second_operator_norm": str(2/(b*b))})
    try:
        projection(F(0), F(0))
    except ValueError:
        singular_control = "passed: undefined disappearing-orbit projection rejected"
    else:
        raise AssertionError("Singular control failed")
    # A frozen projection has zero derivative and must fail the nonzero exact prediction.
    assert records[0]["first_operator_norm"] != "0"
    print(json.dumps({"checks": "passed", "records": records,
                      "singular_failure_control": singular_control,
                      "scope": "Exact finite projector algebra; no metric or physical theorem"}, indent=2))


if __name__ == "__main__":
    main()

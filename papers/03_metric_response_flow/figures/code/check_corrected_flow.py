#!/usr/bin/env python3
"""Exact finite matrix-series checks for the corrected radial flow.

Scope: algebra of logarithmic metric velocity and conjugation only.
This is not an elliptic or parabolic verification.
"""

from fractions import Fraction as Q
from pathlib import Path
import json
import math


ORDER = 12
ZERO = ((Q(0), Q(0)), (Q(0), Q(0)))
IDENTITY = ((Q(1), Q(0)), (Q(0), Q(1)))


def ma(a, b):
    return tuple(tuple(a[i][j] + b[i][j] for j in range(2)) for i in range(2))


def ms(a, s):
    return tuple(tuple(s * a[i][j] for j in range(2)) for i in range(2))


def mm(a, b):
    return tuple(
        tuple(sum((a[i][k] * b[k][j] for k in range(2)), Q(0)) for j in range(2))
        for i in range(2)
    )


def series(matrix=ZERO):
    return [matrix] + [ZERO for _ in range(ORDER)]


def add(a, b):
    return [ma(x, y) for x, y in zip(a, b)]


def scale(a, c):
    return [ms(x, c) for x in a]


def multiply(a, b):
    result = series()
    for n in range(ORDER + 1):
        total = ZERO
        for i in range(n + 1):
            total = ma(total, mm(a[i], b[n - i]))
        result[n] = total
    return result


def exponential(u, coefficient):
    result = series(IDENTITY)
    power = series(IDENTITY)
    for degree in range(1, ORDER + 1):
        power = multiply(power, u)
        result = add(result, scale(power, coefficient**degree / math.factorial(degree)))
    return result


def t_derivative(u):
    """Return t d/dt of a power series."""
    return [ms(value, Q(n)) for n, value in enumerate(u)]


def is_zero(u):
    return all(value == ZERO for value in u)


def first_nonzero(u):
    return next((i for i, value in enumerate(u) if value != ZERO), None)


def conjugate(left, value, right):
    return multiply(multiply(left, value), right)


def run():
    cases = 0
    noncommuting_cases = 0
    wrong_speed_detected = True
    missing_dexp_detected = True
    for a, b in [(1, 1), (1, 3), (2, 5), (7, 2)]:
        q = ((Q(b, a + b), Q(0)), (Q(0), -Q(a, a + b)))
        for offdiag in [Q(0), Q(1, 3), -Q(2, 5)]:
            u = series()
            u[2] = ((Q(1, 4), offdiag), (offdiag, -Q(1, 4)))
            u[3] = ((Q(0), Q(2, 7)), (Q(2, 7), Q(0)))
            u[4] = ((Q(3, 11), -Q(1, 5)), (-Q(1, 5), -Q(3, 11)))
            em = exponential(u, Q(-1))
            ep = exponential(u, Q(1))
            ehm = exponential(u, Q(-1, 2))
            ehp = exponential(u, Q(1, 2))
            qs = series(q)
            tq = conjugate(ehm, qs, ehp)
            t2q = conjugate(em, qs, ep)
            bracket = add(add(t2q, qs), scale(tq, Q(-2)))
            square = add(
                add(conjugate(ehm, tq, ehp), scale(tq, Q(-2))), qs
            )
            assert bracket == square
            assert first_nonzero(bracket) is None or first_nonzero(bracket) >= 4

            dexp = multiply(em, t_derivative(ep))
            residual_over_k = add(scale(bracket, Q(-1)), scale(dexp, Q(-1)))
            assert first_nonzero(residual_over_k) == 2
            assert residual_over_k[2] == ms(u[2], Q(-2))

            # Incorrect A=-(k/2)q leaves the nonzero constant term kq.
            wrong_speed = add(
                add(scale(t2q, Q(-1, 2)), scale(qs, Q(-1, 2))),
                add(scale(tq, Q(2)), scale(dexp, Q(-1, 2))),
            )
            wrong_speed_detected &= wrong_speed[0] == q and wrong_speed[0] != ZERO

            # Omitting e^-u (d/dtau)e^u loses the genuine quadratic term.
            without_dexp = scale(bracket, Q(-1))
            missing_dexp_detected &= without_dexp[2] == ZERO
            missing_dexp_detected &= residual_over_k[2] != ZERO

            if offdiag:
                assert mm(u[2], q) != mm(q, u[2])
                noncommuting_cases += 1
            cases += 1

    assert wrong_speed_detected
    assert missing_dexp_detected
    data = {
        "scope": "Exact finite matrix-series algebra only; no elliptic or PDE verification",
        "truncation_order": ORDER,
        "rational_matrix_cases": cases,
        "noncommuting_quadratic_correction_cases": noncommuting_cases,
        "conjugation_square_identity": True,
        "residual_over_k_starts_at_t_squared": True,
        "wrong_radial_time_factor_detected": wrong_speed_detected,
        "omitted_exponential_derivative_detected": missing_dexp_detected,
    }
    output = (Path(__file__).resolve().parents[1] / "data/corrected_flow_checks.json")
    output.write_text(json.dumps(data, indent=2) + "\n")
    print(json.dumps(data, indent=2))


if __name__ == "__main__":
    run()

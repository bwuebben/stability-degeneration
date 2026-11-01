#!/usr/bin/env python3
"""Exact matrix-series checks of the sliced metric velocity, not a PDE test."""
from fractions import Fraction as F
from math import factorial
import json
from pathlib import Path

ORDER = 10


def zero(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def eye(n):
    a = zero(n)
    for i in range(n):
        a[i][i] = F(1)
    return a


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(a, q):
    return [[q * x for x in row] for row in a]


def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def comm(a, b):
    return add(mul(a, b), scale(mul(b, a), -1))


def tr(a):
    return sum(a[i][i] for i in range(len(a)))


def is_zero(a):
    return all(x == 0 for row in a for x in row)


def series_const(a):
    return [a] + [zero(len(a)) for _ in range(ORDER)]


def series_add(a, b):
    return [add(x, y) for x, y in zip(a, b)]


def series_scale(a, q):
    return [scale(x, q) for x in a]


def series_mul(a, b):
    out = []
    for n in range(ORDER + 1):
        c = zero(len(a[0]))
        for j in range(n + 1):
            c = add(c, mul(a[j], b[n - j]))
        out.append(c)
    return out


def exp_series(u, q):
    power = eye(len(u))
    out = []
    for j in range(ORDER + 1):
        out.append(scale(power, q ** j / factorial(j)))
        power = mul(power, u)
    return out


def dexp_series(u, v):
    powers = [eye(len(u))]
    for _ in range(ORDER):
        powers.append(mul(powers[-1], u))
    out = []
    for n in range(ORDER + 1):
        a = zero(len(u))
        for j in range(n + 1):
            a = add(a, mul(mul(powers[j], v), powers[n - j]))
        out.append(scale(a, F(1, factorial(n + 1))))
    return out


def ad_power(u, a, n):
    for _ in range(n):
        a = comm(u, a)
    return a


C_COEF = [F(0) for _ in range(ORDER + 1)]
S_COEF = C_COEF.copy()
for n in range(0, ORDER + 1, 2):
    C_COEF[n] = F(1, 2 ** n * factorial(n))
    S_COEF[n] = F(1, 2 ** n * factorial(n + 1))
B_COEF = []
for n in range(ORDER + 1):
    B_COEF.append(C_COEF[n] - sum(S_COEF[j] * B_COEF[n-j]
                                      for j in range(1, n+1)))


def op_series(u, values, coef):
    out = []
    for n in range(ORDER + 1):
        a = zero(len(u))
        for j in range(n + 1):
            if coef[j]:
                a = add(a, scale(ad_power(u, values[n-j], j), coef[j]))
        out.append(a)
    return out


def solve_s(u, rhs):
    out = []
    for n in range(ORDER + 1):
        a = rhs[n]
        for j in range(1, n+1):
            if S_COEF[j]:
                a = add(a, scale(ad_power(u, out[n-j], j), -S_COEF[j]))
        out.append(a)
    return out


def blocks(ranks):
    start = 0
    out = []
    for r in ranks:
        out.append(list(range(start, start+r)))
        start += r
    return out


def h_projection(a, ranks):
    n = len(a)
    out = zero(n)
    for indices in blocks(ranks):
        mean = sum(a[i][i] for i in indices) / len(indices) - tr(a) / n
        for i in indices:
            out[i][i] = mean
    return out


def h_coords(a, ranks):
    values = [sum(a[i][i] for i in indices) / len(indices)
              for indices in blocks(ranks)]
    return [x-values[-1] for x in values[:-1]]


def h_basis(ranks):
    n = sum(ranks)
    out = []
    for indices in blocks(ranks)[:-1]:
        a = scale(eye(n), -F(len(indices), n))
        for i in indices:
            a[i][i] += 1
        out.append(a)
    return out


def from_coords(values, ranks):
    out = zero(sum(ranks))
    for q, a in zip(values, h_basis(ranks)):
        out = add(out, scale(a, q))
    return out


def sliced_solve(u, curvature, ranks, wrong_factor=False):
    invsk = solve_s(u, curvature)
    basis = h_basis(ranks)
    dim = len(basis)
    rseries = [[[F(0) for _ in range(dim)] for _ in range(dim)]
               for _ in range(ORDER+1)]
    for j, a in enumerate(basis):
        ba = op_series(u, series_const(a), B_COEF)
        for n in range(ORDER+1):
            for i, v in enumerate(h_coords(h_projection(ba[n], ranks), ranks)):
                rseries[n][i][j] = v
    assert rseries[0] == eye(dim)
    coords = []
    for n in range(ORDER+1):
        rhs = [-v for v in h_coords(h_projection(invsk[n], ranks), ranks)]
        for j in range(1, n+1):
            for i in range(dim):
                rhs[i] -= sum(rseries[j][i][k] * coords[n-j][k]
                              for k in range(dim))
        coords.append(rhs)
    a = [from_coords(v, ranks) for v in coords]
    if wrong_factor:
        a = series_scale(a, 2)
    udot = series_add(series_scale(invsk, -2),
                     series_scale(op_series(u, a, B_COEF), -2))
    return a, udot


def build_case(ranks, seed):
    n = sum(ranks)
    raw = [[F(((i+1)*(j+1)+seed*(i+j+1)) % 7 - 3, seed+2)
            for j in range(n)] for i in range(n)]
    u = add(raw, scale(h_projection(raw, ranks), -1))
    # Remove the scalar trace separately; P only removes trace-free block scalars.
    u = add(u, scale(eye(n), -tr(u)/n))
    v = [[F(((i+2)*(j+2)+seed*(i+j+2)) % 11 - 5, seed+3)
          for j in range(n)] for i in range(n)]
    v = add(v, scale(h_projection(v, ranks), -1))
    v = add(v, scale(eye(n), -tr(v)/n))
    a = from_coords([F(seed+j+1, seed+2) for j in range(len(ranks)-1)], ranks)
    k = series_const(scale(a, -1))
    k[1] = v
    k[2] = add(u, scale(a, F(1, 3)))
    return u, v, a, k


def run():
    total = 0
    noncommuting = 0
    wrong_factor_detected = False
    omitted_centering_detected = False
    for ranks in [(1, 1), (1, 1, 1), (2, 1), (1, 2, 1)]:
        for seed in range(1, 4):
            u, v, a0, curvature = build_case(ranks, seed)
            assert is_zero(h_projection(u, ranks)) and tr(u) == 0
            ep = exp_series(u, F(1, 2))
            em = exp_series(u, -F(1, 2))
            direct = series_add(
                series_add(series_mul(series_mul(em, series_const(a0)), ep),
                           series_mul(series_mul(ep, series_const(a0)), em)),
                series_mul(series_mul(em, dexp_series(u, v)), em))
            formula = series_add(series_scale(op_series(u, series_const(a0), C_COEF), 2),
                                 op_series(u, series_const(v), S_COEF))
            assert direct == formula
            a, udot = sliced_solve(u, curvature, ranks)
            lhs = series_add(series_scale(op_series(u, a, C_COEF), 2),
                             op_series(u, udot, S_COEF))
            assert lhs == series_scale(curvature, -2)
            assert all(is_zero(h_projection(vn, ranks)) and tr(vn)==0 for vn in udot)
            assert solve_s(u, op_series(u, series_const(v), S_COEF)) == series_const(v)
            wrong_a, wrong_u = sliced_solve(u, curvature, ranks, wrong_factor=True)
            wrong_factor_detected |= any(not is_zero(h_projection(vn, ranks)) for vn in wrong_u)
            omitted_centering_detected |= direct != series_add(
                series_scale(series_const(a0), 2), series_const(v))
            noncommuting += not is_zero(comm(u, a0))
            total += 1
    assert wrong_factor_detected and omitted_centering_detected
    result = {
        'scope': 'Exact finite matrix-series algebra; no elliptic, parabolic or invariant-manifold verification',
        'truncation_order': ORDER,
        'rational_block_cases': total,
        'noncommuting_cases': noncommuting,
        'centered_velocity_identity': True,
        'sliced_solution_and_normal_projection': True,
        'principal_dexp_coefficient_cancellation': True,
        'wrong_factor_in_slice_detected': wrong_factor_detected,
        'omitted_centering_detected': omitted_centering_detected,
    }
    print(json.dumps(result, indent=2))
    (Path(__file__).resolve().parents[1] / "data/sliced_equations_checks.json").write_text(
        json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    run()

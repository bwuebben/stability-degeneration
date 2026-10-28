#!/usr/bin/env python3
"""Exact finite checks of the tangent-family spectral and response formulas.

Uses only the standard library and writes tangent_bridge_checks.json.
The geometric and analytic assertions are proved in the paper."""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path


def det(a):
    a = [list(row) for row in a]
    n = len(a)
    if n == 0:
        return F(1)
    out = F(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            out = -out
        v = a[j][j]
        out *= v
        for i in range(j + 1, n):
            ratio = a[i][j] / v
            for k in range(j + 1, n):
                a[i][k] -= ratio * a[j][k]
    return out


def lap(k):
    d = len(k[0])
    a = [[F(0) for _ in range(d + 2)] for _ in range(d + 2)]
    for i in range(2):
        for j in range(d):
            v = k[i][j]
            a[i][i] += v
            a[j + 2][j + 2] += v
            a[i][j + 2] -= v
            a[j + 2][i] -= v
    return a


def submat(a, inds):
    return [[a[i][j] for j in inds] for i in inds]


def char_coeff(a):
    """det(lambda I - A), ascending, via independent principal minors."""
    n = len(a)
    out = [F(0)] * (n + 1)
    for size in range(n + 1):
        out[n - size] = (-1) ** size * sum(
            (det(submat(a, inds)) for inds in combinations(range(n), size)), F(0)
        )
    return out


def poly_mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def product(a):
    out = F(1)
    for x in a:
        out *= x
    return out


def check_graph(d, x, y):
    x, y = list(map(F, x)), list(map(F, y))
    sig = sum(x)
    assert sig == sum(y) and min(x + y) > 0
    k = [[xx * yy / sig for yy in y] for xx in x]
    assert [sum(row) for row in k] == x
    assert [sum(k[i][j] for i in range(2)) for j in range(d)] == y
    # Independently chosen row/column coefficients, including uneven scales.
    aa = [F(2, 3), F(7, 5)]
    qq = [F(j + 2, j + 5) for j in range(d)]
    wr = [x[i] / aa[i] for i in range(2)]
    wq = [sig * qq[j] / y[j] for j in range(d)]
    for i in range(2):
        for j in range(d):
            assert aa[i] * qq[j] * wr[i] / wq[j] == k[i][j]
    for j, ell in combinations(range(d), 2):
        assert k[0][j] * k[1][ell] == k[0][ell] * k[1][j]
    a = lap(k)
    assert all(sum(row) == 0 for row in a)
    cofactor = det(submat(a, [0] + list(range(2, d + 2))))
    assert cofactor == product(x + y) / sig
    # Nonzero eigenvalue product is n times a cofactor, not the cofactor.
    cp = char_coeff(a)
    n = d + 2
    spec_prod = (-1) ** (n - 1) * cp[1]
    assert spec_prod == n * product(x + y) / sig
    assert spec_prod != cofactor
    delta = min(x + y)
    # On mean-zero vectors u = Bz, B_i = e_i - e_(n-1).
    # All principal minors prove this particular rational restricted form PSD.
    restricted = [[
        a[i][j] - a[i][n - 1] - a[n - 1][j] + a[n - 1][n - 1]
        - delta * (int(i == j) + 1)
        for j in range(n - 1)] for i in range(n - 1)]
    count = 0
    for size in range(1, n):
        for inds in combinations(range(n - 1), size):
            assert det(submat(restricted, inds)) >= 0
            count += 1
    # Trial e_v - 1/n verifies the upper Rayleigh quotient exactly.
    v = (x + y).index(delta)
    trial = [F(int(i == v)) - F(1, n) for i in range(n)]
    energy = sum(trial[i] * a[i][j] * trial[j] for i in range(n) for j in range(n))
    mass = sum(t * t for t in trial)
    assert energy / mass == F(n, n - 1) * delta
    return {
        "d": d, "x": x, "y": y, "sigma": sig, "minimum_margin": delta,
        "cofactor": cofactor, "positive_eigenvalue_product": spec_prod,
        "restricted_lower_bound_principal_minors": count,
    }


class Qr:
    """Exact Q[r]/(r^2-r-1/3), the fourfold ample wall field."""
    def __init__(self, a=0, b=0):
        self.a, self.b = F(a), F(b)

    @staticmethod
    def coerce(x):
        return x if isinstance(x, Qr) else Qr(x)

    def __add__(self, x):
        x = self.coerce(x)
        return Qr(self.a + x.a, self.b + x.b)

    __radd__ = __add__

    def __neg__(self):
        return Qr(-self.a, -self.b)

    def __sub__(self, x):
        return self + -self.coerce(x)

    def __rsub__(self, x):
        return self.coerce(x) - self

    def __mul__(self, x):
        x = self.coerce(x)
        return Qr(self.a * x.a + self.b * x.b / 3,
                  self.a * x.b + self.b * x.a + self.b * x.b)

    __rmul__ = __mul__

    def inverse(self):
        norm = self.a * (self.a + self.b) - self.b * self.b / 3
        return Qr((self.a + self.b) / norm, -self.b / norm)

    def __truediv__(self, x):
        return self * self.coerce(x).inverse()

    def __bool__(self):
        return bool(self.a or self.b)

    def __eq__(self, x):
        x = self.coerce(x)
        return self.a == x.a and self.b == x.b


def chow_mul(a, b):
    # eta_1^2 = eta_2^2 = -h_1 h_2/2; h_j^2 = 0.
    out = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            mon = [ma[i] + mb[i] for i in range(4)]
            coef = ca * cb
            for i in range(2):
                while mon[i] >= 2:
                    mon[i] -= 2
                    mon[2] += 1
                    mon[3] += 1
                    coef *= F(-1, 2)
            if mon[2] > 1 or mon[3] > 1:
                continue
            key = tuple(mon)
            out[key] = out.get(key, Qr()) + coef
    return out


def chow_pow(a, p):
    out = {(0, 0, 0, 0): Qr(1)}
    for _ in range(p):
        out = chow_mul(out, a)
    return out


def integral(a):
    return a.get((1, 1, 1, 1), Qr())


def fixed_volume_fourfold():
    # Computes derivatives directly in Chow ring, not by differentiating F_k.
    e = []
    for j in range(4):
        key = tuple(int(i == j) for i in range(4))
        e.append({key: Qr(1)})
    r = Qr(0, 1)
    h = {next(iter(e[i])): [Qr(1), Qr(1), r, r][i] for i in range(4)}
    h2, h3 = chow_pow(h, 2), chow_pow(h, 3)
    # C_p = 2e_p; K_p = 2e_p - 1/2 sum e_q, n=4.
    deg = [2 * integral(chow_mul(ei, h3)) for ei in e]
    assert all(v == deg[0] for v in deg)
    dz = []
    for p in range(4):
        kp = {next(iter(e[q])): Qr(F(2 if p == q else 0) - F(1, 2)) for q in range(4)}
        dz.append([3 * integral(chow_mul(chow_mul(kp, ei), h2)) for ei in e])
    dv = [integral(chow_mul(ei, h3)) / 6 for ei in e]
    hv = [Qr(1), Qr(1), r, r]
    assert all(sum(dz[p][i] * hv[i] for i in range(4)) == 0 for p in range(4))
    volume = integral(chow_pow(h, 4)) / 24
    assert sum(dv[i] * hv[i] for i in range(4)) == 4 * volume
    augmented_det = det(dz[:3] + [dv])
    assert augmented_det
    # Missing volume row leaves the radial direction in the kernel.
    assert det(dz) == 0
    return {"ring": "r^2=r+1/3; r=(3+sqrt(21))/6", "common_C_degree": deg[0],
            "volume": volume, "dz": dz, "dV": dv,
            "det_three_centered_rows_and_volume": augmented_det,
            "radial_kernel_verified": True}


def encode(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, Qr):
        return {"1": str(x.a), "r": str(x.b)}
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def main():
    cases, spectra = [], []
    for d in [2, 4, 6]:
        for tiny in [F(1, 7), F(1, 10**12), F(1, 10**40)]:
            x = [tiny, 1 - tiny]
            y = [F(1, d)] * d
            cases.append(check_graph(d, x, y))
            # Intersecting faces with arbitrarily separated positive orders.
            y = [tiny**2] + [(1 - tiny**2) / (d - 1)] * (d - 1)
            cases.append(check_graph(d, x, y))
        sig = F(7, 5)
        x, y = [sig / 2] * 2, [sig / d] * d
        a = lap([[xx * yy / sig for yy in y] for xx in x])
        eigs = [F(0), sig / 2] + [sig / d] * (d - 1) + [sig * (F(1, 2) + F(1, d))]
        expected = [F(1)]
        for ev in eigs:
            expected = poly_mul(expected, [-ev, F(1)])
        assert char_coeff(a) == expected
        spectra.append({"d": d, "sigma": sig, "eigenvalues": eigs,
                        "characteristic_polynomial_ascending": expected})
    # Frozen rank-one coefficients do not imply actual rank-one nearby array.
    h = F(1, 10**20)
    c_actual = [[F(2), F(3)], [F(2), 3 * (1 + h)]]
    actual_cycle = c_actual[0][0] * c_actual[1][1] / (c_actual[0][1] * c_actual[1][0])
    assert actual_cycle == 1 + h and actual_cycle != 1
    # A diagonal-orbit trajectory has rank-one frozen weights even without balance.
    g = [F(2), F(5), F(3), F(7)]
    frozen = [[F(j + 2) * (g[i] / g[j + 2])**2 for j in range(2)] for i in range(2)]
    xx = [sum(row) for row in frozen]
    yy = [sum(frozen[i][j] for i in range(2)) for j in range(2)]
    total = sum(xx)
    assert all(frozen[i][j] == xx[i] * yy[j] / total for i in range(2) for j in range(2))
    assert xx != [F(1), F(1)]  # Instantaneous margins are not prescribed demands.
    record = {
        "status": "PASS finite exact identities only; no PDE or symmetry proof by computation",
        "graph_cases": cases, "equal_margin_spectra": spectra,
        "fixed_volume_fourfold": fixed_volume_fourfold(),
        "controls": {"missing_mass_n_factor_rejected": True,
                     "full_centered_derivative_singular_without_volume": True,
                     "nearby_coefficient_rank_one_not_assumed": actual_cycle,
                     "instantaneous_margins_not_demands": True,
                     "frozen_diagonal_orbit_rank_one": True},
    }
    path = (Path(__file__).resolve().parents[1] / "data/tangent_bridge_checks.json")
    path.write_text(json.dumps(encode(record), indent=2, sort_keys=True) + "\n")
    print(f"PASS: {len(cases)} exact separated-margin cases; {len(spectra)} full characteristic polynomials; exact fourfold fixed-volume rank; controls.")
    print(path)


if __name__ == "__main__":
    main()

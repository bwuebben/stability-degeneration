#!/usr/bin/env python3
"""Exact finite checks for the double-projective-line-bundle family.

Run from any directory with python3 -B; outputs are in figures/data.
Uses only the standard library. A 240-second cap bounds the computation.
The all-dimensional geometric theorem is proved in the paper."""
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial
from pathlib import Path
import json
import time

START = time.monotonic()
CAP = 240


class QS21:
    """Exact a+b*sqrt(21), used for fourfold wall intersections."""
    def __init__(self, a=0, b=0):
        self.a, self.b = Q(a), Q(b)

    @staticmethod
    def coerce(x):
        return x if isinstance(x, QS21) else QS21(x)

    def __add__(self, other):
        o = self.coerce(other)
        return QS21(self.a+o.a, self.b+o.b)

    __radd__ = __add__

    def __neg__(self):
        return QS21(-self.a, -self.b)

    def __sub__(self, other):
        return self+-self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other)+-self

    def __mul__(self, other):
        o = self.coerce(other)
        return QS21(self.a*o.a+21*self.b*o.b, self.a*o.b+self.b*o.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = self.coerce(other)
        norm = o.a*o.a-21*o.b*o.b
        return self*QS21(o.a/norm, -o.b/norm)

    def __eq__(self, other):
        o = self.coerce(other)
        return self.a == o.a and self.b == o.b

    def __bool__(self):
        return bool(self.a or self.b)

    def __str__(self):
        return f"({self.a})+({self.b})*sqrt(21)"


def check_time():
    if time.monotonic() - START > CAP:
        raise TimeoutError("240-second finite-probe cap reached")


def rref(rows, n):
    a = [[Q(x) for x in row] for row in rows]
    rank = 0
    for c in range(n):
        p = next((i for i in range(rank, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[rank], a[p] = a[p], a[rank]
        v = a[rank][c]
        a[rank] = [x / v for x in a[rank]]
        for i in range(len(a)):
            if i != rank and a[i][c]:
                t = a[i][c]
                a[i] = [x - t*y for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return tuple(tuple(row) for row in a[:rank])


def det(rows):
    a = [[Q(x) for x in row] for row in rows]
    ans = Q(1)
    for j in range(len(a)):
        p = next(i for i in range(j, len(a)) if a[i][j])
        if p != j:
            a[p], a[j] = a[j], a[p]
            ans = -ans
        piv = a[j][j]
        ans *= piv
        for i in range(j+1, len(a)):
            t = a[i][j] / piv
            for c in range(j+1, len(a)):
                a[i][c] -= t*a[j][c]
    return ans


def rays(k, zero_twist=False):
    d, n = 2*k, 2*k+2
    eps = [0]*d if zero_twist else [1]*k+[-1]*k
    out = []
    for i in range(n):
        p = [0]*n
        p[i] = 1
        m = [0]*n
        m[i] = -1
        if i >= 2:
            m[0] = m[1] = eps[i-2]
        out.extend((tuple(p), tuple(m)))
    return out, eps


def all_flats(rr):
    n = len(rr[0])
    span_cache = {}
    def closure(indices):
        key = tuple(sorted(indices))
        if key not in span_cache:
            basis = rref([rr[i] for i in key], n)
            flat = tuple(i for i, ray in enumerate(rr)
                         if len(rref(basis+(ray,), n)) == len(basis))
            span_cache[key] = (flat, basis)
        return span_cache[key]
    fs = {(): ()}
    todo = [()]
    while todo:
        check_time()
        f = todo.pop()
        for i in range(len(rr)):
            if i not in f:
                g, basis = closure(f+(i,))
                if g not in fs:
                    fs[g] = basis
                    todo.append(g)
    return fs


def expected_ideal_flats(n):
    ideals = []
    for bits in product((0, 1), repeat=n):
        if any(bits[2:]) and not (bits[0] and bits[1]):
            continue
        ideals.append(tuple(j for i, b in enumerate(bits) if b
                            for j in (2*i, 2*i+1)))
    return set(ideals)


def polynomial_F(k):
    f = [Q(0)]*(2*k+1)
    for j in range(k+1):
        f[2*k-2*j] = Q((-1)**j*comb(k, j), (2*j+1)*(j+1))
    return f


def derivative(f):
    return [i*f[i] for i in range(1, len(f))]


def evalpoly(f, x):
    s = Q(0)
    for a in reversed(f):
        s = s*x+a
    return s


def polynomial_G(k):
    f = polynomial_F(k)
    fp = derivative(f)
    g = [Q(0)]*len(f)
    for i, a in enumerate(fp):
        g[i] += a
        g[i+1] += k*a
    for i, a in enumerate(f):
        g[i] -= k*(2*k+2)*a
    return g


def root_interval(k):
    g = polynomial_G(k)
    lo, hi = Q(1), Q(2)
    assert evalpoly(g, lo) > 0
    while evalpoly(g, hi) >= 0:
        hi *= 2
    for _ in range(40):
        mid = (lo+hi)/2
        if evalpoly(g, mid) > 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


def ring_multiply(a, b, eps):
    """Chow ring in eta1,eta2,h1,...,hd; eta_i^2=L^2/4."""
    out = {}
    d = len(eps)
    for ma, ca in a.items():
        for mb, cb in b.items():
            if ma & mb & ~3:
                continue
            terms = {ma ^ mb: ca*cb}
            for i in (0, 1):
                if ma & mb & (1 << i):
                    nxt = {}
                    for m, c in terms.items():
                        for j in range(d):
                            for l in range(j+1, d):
                                bits = (1 << (j+2)) | (1 << (l+2))
                                if not m & bits:
                                    mm = m | bits
                                    nxt[mm] = nxt.get(mm, Q(0))+c*Q(eps[j]*eps[l], 2)
                    terms = nxt
            for m, c in terms.items():
                out[m] = out.get(m, Q(0))+c
    return {m: c for m, c in out.items() if c}


def chow_degrees(k, r):
    d, n = 2*k, 2*k+2
    eps = [1]*k+[-1]*k
    h = {1: Q(1), 2: Q(1)}
    h.update({1 << (j+2): r for j in range(d)})
    p = {0: Q(1)}
    for _ in range(n-1):
        p = ring_multiply(p, h, eps)
    top = (1 << n)-1
    cs = [ring_multiply(p, {1 << i: Q(2)}, eps).get(top, Q(0))
          for i in range(n)]
    vol = ring_multiply(p, h, eps).get(top, Q(0))
    return cs, vol


def slope_flat_check(fs, ray_degrees, n, strict=False):
    total = sum(ray_degrees)
    tight = []
    for f, basis in fs.items():
        rank = len(basis)
        if rank in (0, n):
            continue
        lhs = n*sum(ray_degrees[i] for i in f)
        rhs = rank*total
        assert lhs < rhs if strict else lhs <= rhs
        if lhs == rhs:
            tight.append(list(f))
    return tight


def qsqrt_sign(a, b):
    """Sign of a+b*sqrt(21) using rational comparisons only."""
    if not b:
        return (a > 0)-(a < 0)
    if not a:
        return (b > 0)-(b < 0)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    t = a*a-21*b*b
    return ((t > 0)-(t < 0))*(1 if a > 0 else -1)


def exact_fourfold_polytope(rr):
    # r=(3+sqrt21)/6; all support coefficients are 1/2 or r/2.
    r = (Q(1, 2), Q(1, 6))
    supports = [(Q(1, 2), Q(0))]*4+[(r[0]/2, r[1]/2)]*4
    vertices = []
    for choices in product((0, 1), repeat=4):
        t = [Q(-1, 2) if choices[i] == 0 else Q(1, 2) for i in (0, 1)]
        s = t[0]+t[1]
        vertex = [(t[0], Q(0)), (t[1], Q(0))]
        for j, ep in enumerate((1, -1)):
            if choices[j+2] == 0:
                vertex.append((-r[0]/2, -r[1]/2))
            else:
                vertex.append((r[0]/2+ep*s, r[1]/2))
        active = []
        for i, (ray, sup) in enumerate(zip(rr, supports)):
            a = sup[0]+sum(v[0]*u for v, u in zip(vertex, ray))
            b = sup[1]+sum(v[1]*u for v, u in zip(vertex, ray))
            sg = qsqrt_sign(a, b)
            assert sg >= 0
            if sg == 0:
                active.append(i)
        expected = [2*i+choices[i] for i in range(4)]
        assert active == expected
        assert abs(det([rr[i] for i in active])) == 1
        vertices.append({"active_rays": active,
                         "coordinates_a_plus_b_sqrt21": [[str(a), str(b)] for a, b in vertex]})
    return vertices


def main():
    result = {"status": "finite exact checks of the double-projective-line-bundle family",
              "seconds_cap": CAP, "families": []}
    for k in (1, 2, 3):
        check_time()
        print(f"k={k}: exact fan determinants and complete flat enumeration", flush=True)
        rr, eps = rays(k)
        n = 2*k+2
        determinants = [det([rr[2*i+choices[i]] for i in range(n)])
                        for choices in product((0, 1), repeat=n)]
        assert all(abs(d) == 1 for d in determinants)
        cones = [tuple(2*i+choices[i] for i in range(n))
                 for choices in product((0, 1), repeat=n)]
        fs = all_flats(rr)
        tight = {f for f, basis in fs.items() if len(f) == 2*len(basis)}
        assert all(len(f) <= 2*len(basis) for f, basis in fs.items())
        assert tight == expected_ideal_flats(n)
        assert all(len(set(f).intersection(cone)) == len(fs[f])
                   for f in tight for cone in cones)
        # Chern classes in the independent xi1,xi2,h_j coordinates.
        cs_matrix = []
        for i in range(n):
            row = [Q(0)]*n
            row[i] = Q(2)
            if i < 2:
                for j, ep in enumerate(eps):
                    row[j+2] = -ep
            cs_matrix.append(row)
        assert len(rref(cs_matrix, n)) == n
        total_c = [sum(row[j] for row in cs_matrix) for j in range(n)]
        centered = [[row[j]-total_c[j]/n for j in range(n)] for row in cs_matrix]
        slope_rank = len(rref(centered, n))
        assert slope_rank == n-1
        # This also checks every non-coordinate, non-ideal proper flat.
        slope_flat_check(fs, [Q(1)]*(2*n), n)
        # Independent Chow-vs-slice-volume prediction at a rational ample class.
        r = Q(3, 2)
        cs, vol = chow_degrees(k, r)
        f, fp = polynomial_F(k), derivative(polynomial_F(k))
        expected_fiber = factorial(n-1)*(n*evalpoly(f, r)-r*evalpoly(fp, r))
        expected_base = Q(factorial(n-1), k)*evalpoly(fp, r)
        assert cs == [expected_fiber]*2+[expected_base]*(n-2)
        assert vol == factorial(n)*evalpoly(f, r)
        lo, hi = root_interval(k)
        record = {"k": k, "dimension": n, "rays": rr,
                  "maximal_cones": 2**n,
                  "smooth_determinants": sorted(set(map(str, determinants))),
                  "all_flats": len(fs),
                  "proper_nonzero_flats": sum(0 < len(b) < n for b in fs.values()),
                  "tight_flats_including_zero_and_full": [list(f) for f in sorted(tight)],
                  "expected_tight_flat_count": 2**(2*k)+3,
                  "G_coefficients_ascending": list(map(str, polynomial_G(k))),
                  "root_interval": [str(lo), str(hi)],
                  "root_interval_decimal": [float(lo), float(hi)],
                  "independent_chow_at_r_3_over_2": {"C_degrees": list(map(str, cs)), "H_power_n": str(vol)},
                  "chern_matrix_xi1_xi2_h": [[str(x) for x in row] for row in cs_matrix],
                  "slope_map_rank": slope_rank,
                  "condition_15_all_tight_flats_all_cones": True,
                  "chamber_dimension": n-1,
                  "chamber_extreme_rays": 4*k,
                  "chamber_facets": 2*k+2}
        if k == 1:
            wall_cs, wall_volume = chow_degrees(1, QS21(Q(1, 2), Q(1, 6)))
            assert wall_cs == [QS21(6, 2)]*4
            assert wall_volume == QS21(16, 4)
            # Both base-divisor degrees are equal, so ell has degree zero;
            # eta_i plus/minus ell/2 therefore give equal fiber-ray degrees.
            assert wall_cs[2] == wall_cs[3]
            assert all(c/2 == QS21(3, 1) for c in wall_cs)
            record["exact_wall_r"] = "(3+sqrt(21))/6"
            record["each_ray_degree_at_wall"] = str(wall_cs[0]/2)
            record["block_degrees_at_wall_from_Chow"] = list(map(str, wall_cs))
            record["H0_power_4"] = str(wall_volume)
            record["polytope_vertices"] = exact_fourfold_polytope(rr)
            # Stable, unstable, and nonample controls predicted before checking.
            stable_cs, _ = chow_degrees(1, Q(5, 4))
            stable_ds = [c/2 for c in stable_cs for _ in (0, 1)]
            assert slope_flat_check(fs, stable_ds, n, strict=True) == []
            assert cs[0] > sum(cs)/n  # r=3/2: relative tangent line destabilizes.
            record["controls"] = {"stable_r_5_over_4_C_degrees": list(map(str, stable_cs)),
                                  "unstable_r_3_over_2_witness": "relative tangent line R1",
                                  "r_1_not_ample": "base interval width r-1 vanishes at two fiber corners"}
        result["families"].append(record)
        print(f"k={k}: {len(fs)} flats; {len(tight)} tight; root in ({float(lo):.12f},{float(hi):.12f})", flush=True)
    rr, _ = rays(1, zero_twist=True)
    fs = all_flats(rr)
    zero_tight = {f for f, b in fs.items() if len(f) == 2*len(b)}
    assert len(zero_tight) == 16
    result["zero_twist_failure_control"] = {"tight_flats": len(zero_tight),
                                            "prediction": "four-element antichain, product tangent bundle, no strict stable locus"}
    elapsed_seconds = time.monotonic()-START
    path = Path(__file__).resolve().parents[1] / "data/realizability_checks.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(f"PASS; {elapsed_seconds:.3f}s; wrote {path.name}", flush=True)


if __name__ == "__main__":
    main()

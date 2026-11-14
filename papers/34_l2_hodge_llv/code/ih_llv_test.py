#!/usr/bin/env python3
"""Exact character calculations for Section 10 of Paper 34.

Uses the Green--Kim--Laza--Robles character formula (Theorem 1.2):
  1 + q sum_i (x_i+x_i^-1) + sum_{n>=2} ch H(K3^[n]) q^n
    = prod_{m>=1} prod_{i=1}^{12} ((1-x_i q^m)(1-x_i^-1 q^m))^-1.
The q coefficient is exceptional; all decompositions use n >= 2.

Calculates the virtual B_12 character H(K3^[9])-Lambda^2 H(K3^[3]),
its highest-weight decomposition, Betti specialization and dimensions.
Its identification with the intersection cohomology of the dimension-18
moduli space requires hypothesis (H3) of the manuscript. Positivity of
this character does not prove (H3). The q=2 comparison uses the separate
O'Grady-ten calculation. Controls check known Hilbert-scheme characters,
Hodge specializations, and rejection of the symmetric-square alternative.

Pure Python; exact integer and rational arithmetic.
"""
from fractions import Fraction
from itertools import combinations
from math import comb, factorial
from collections import defaultdict
import sys

RANK = 12
NMAX = 9

# ---------------------------------------------------------------- one-variable series
def one_variable_series(nmax):
    """f[n][k] = coefficient of q^n x^k in prod_{m>=1} 1/((1-x q^m)(1-x^{-1} q^m))."""
    ser = {(0, 0): 1}
    for m in range(1, nmax + 1):
        for sign in (1, -1):
            new = defaultdict(int)
            for (n, k), c in ser.items():
                j = 0
                while n + j * m <= nmax:
                    new[(n + j * m, k + sign * j)] += c
                    j += 1
            ser = dict(new)
    f = [defaultdict(int) for _ in range(nmax + 1)]
    for (n, k), c in ser.items():
        f[n][k] += c
    return f

F = one_variable_series(NMAX)

def Fk(k):
    return [F[n].get(k, 0) for n in range(NMAX + 1)]

def mul_trunc(a, b, nmax=NMAX):
    out = [0] * (nmax + 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if i + j > nmax:
                    break
                out[i + j] += ai * bj
    return out

_F0_pows = {0: [1] + [0] * NMAX}
def F0_pow(e):
    if e not in _F0_pows:
        _F0_pows[e] = mul_trunc(F0_pow(e - 1), Fk(0))
    return _F0_pows[e]

def partitions_upto(total, maxlen, maxpart=None):
    """All partitions (as tuples, nonincreasing, positive parts) of size <= total."""
    out = [()]
    def rec(prefix, remaining, cap):
        for p in range(min(cap, remaining), 0, -1):
            new = prefix + (p,)
            out.append(new)
            if len(new) < maxlen:
                rec(new, remaining - p, p)
    rec((), total, total if maxpart is None else maxpart)
    return out

def pad(mu):
    return tuple(mu) + (0,) * (RANK - len(mu))

def strip(mu):
    return tuple(x for x in mu if x)

# ---------------------------------------------------------------- characters of K3^[n]
def ch_hilb(n):
    """Dominant-weight multiplicities of H^*(K3^[n]) as a T^12 character."""
    out = {}
    for mu in partitions_upto(n, RANK):
        ser = F0_pow(RANK - len(mu))
        for part in mu:
            ser = mul_trunc(ser, Fk(part))
        if ser[n]:
            out[mu] = ser[n]
    return out

# two-variable per-coordinate series for the square of a character
def Gk(k, a, b):
    """coefficient of q^a p^b x^k in F(x,q) F(x,p)."""
    return sum(F[a].get(j, 0) * F[b].get(k - j, 0) for j in range(-a, a + 1))

def ch_square(n, sign):
    """Dominant-weight multiplicities of Lambda^2 (sign=-1) or Sym^2 (sign=+1) of H(K3^[n])."""
    c = ch_hilb(n)
    out = {}
    # per coordinate polynomial in (a,b), 0<=a,b<=n
    def coord(k):
        return {(a, b): Gk(k, a, b) for a in range(n + 1) for b in range(n + 1)}
    def mul2(A, B):
        C = defaultdict(int)
        for (a, b), x in A.items():
            if x:
                for (a2, b2), y in B.items():
                    if y and a + a2 <= n and b + b2 <= n:
                        C[(a + a2, b + b2)] += x * y
        return C
    c0 = coord(0)
    pows = [{(0, 0): 1}]
    for _ in range(RANK):
        pows.append(mul2(pows[-1], c0))
    for mu in partitions_upto(2 * n, RANK):
        ser = pows[RANK - len(mu)]
        for part in mu:
            ser = mul2(ser, coord(part))
        sq = ser.get((n, n), 0)
        half = 0
        if all(p % 2 == 0 for p in mu):
            half = c.get(tuple(p // 2 for p in mu), 0)
        val = sq + sign * half
        assert val % 2 == 0
        if val:
            out[mu] = val // 2
    return out

def char_sub(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) - v
        if out[k] == 0:
            del out[k]
    return out

def char_add(*chars):
    out = defaultdict(int)
    for c in chars:
        for k, v in c.items():
            out[k] += v
    return {k: v for k, v in out.items() if v}

# ---------------------------------------------------------------- B_12 representation theory
RHO = tuple(Fraction(2 * (RANK - i) - 1, 2) for i in range(RANK))

def dominant_rep(w):
    return strip(tuple(sorted((abs(x) for x in w), reverse=True)))

def below(mu, lam):
    """mu <= lam in dominance order for B_n (partial sums)."""
    mu, lam = pad(mu), pad(lam)
    s = t = 0
    for a, b in zip(mu, lam):
        s += a
        t += b
        if s > t:
            return False
    return True

_irrep_cache = {}
def irrep(lam):
    """Dominant-weight multiplicities of the B_12 irreducible module V_lam (Freudenthal)."""
    lam = strip(lam)
    if lam in _irrep_cache:
        return _irrep_cache[lam]
    lamp = pad(lam)
    size = sum(lam)
    doms = [mu for mu in partitions_upto(size, RANK) if below(mu, lam)]
    doms.sort(key=lambda m: pad(m), reverse=True)
    norm_lam = sum((a + r) ** 2 for a, r in zip(lamp, RHO))
    mult = {lam: 1}
    roots = []
    for i in range(RANK):
        e = [0] * RANK
        e[i] = 1
        roots.append(tuple(e))
        for j in range(i + 1, RANK):
            for s in (1, -1):
                e = [0] * RANK
                e[i] = 1
                e[j] = s
                roots.append(tuple(e))
    for mu in doms:
        if mu == lam:
            continue
        mup = pad(mu)
        denom = norm_lam - sum((a + r) ** 2 for a, r in zip(mup, RHO))
        if denom == 0:
            continue
        total = 0
        for al in roots:
            k = 1
            while True:
                w = tuple(a + k * b for a, b in zip(mup, al))
                d = dominant_rep(w)
                m = mult.get(d, 0)
                if m == 0:
                    if sum(d) > size or not below(d, lam):
                        break
                    # weight inside the hull but multiplicity zero cannot occur before
                    # higher ones are exhausted; continue the string
                    k += 1
                    if k > 2 * size + 2:
                        break
                    continue
                total += m * sum(x * y for x, y in zip(w, al))
                k += 1
        val = Fraction(2 * total, 1) / denom
        assert val.denominator == 1, (lam, mu, val)
        if val:
            mult[mu] = int(val)
    _irrep_cache[lam] = mult
    return mult

def weyl_dim(lam):
    lamp = pad(lam)
    l = [a + r for a, r in zip(lamp, RHO)]
    num = Fraction(1)
    for i in range(RANK):
        num *= l[i] / RHO[i]
        for j in range(i + 1, RANK):
            num *= (l[i] - l[j]) / (RHO[i] - RHO[j])
            num *= (l[i] + l[j]) / (RHO[i] + RHO[j])
    assert num.denominator == 1
    return int(num)

def orbit_size(mu):
    mup = pad(mu)
    cnt = defaultdict(int)
    for x in mup:
        cnt[x] += 1
    size = factorial(RANK)
    for v, c in cnt.items():
        size //= factorial(c)
    size *= 2 ** sum(1 for x in mup if x)
    return size

def total_dim(char):
    return sum(v * orbit_size(mu) for mu, v in char.items())

def decompose(char):
    """Peel off irreducible B_12 characters. Returns (dict lam -> multiplicity, ok flag)."""
    work = dict(char)
    out = {}
    while work:
        lam = max(work, key=lambda m: pad(m))
        m = work[lam]
        out[lam] = m
        for mu, c in irrep(lam).items():
            work[mu] = work.get(mu, 0) - m * c
            if work[mu] == 0:
                del work[mu]
    return out

# ---------------------------------------------------------------- Hodge specialisation
def spec2(char):
    """Specialise x_2 = ... = x_11 = 1.  Returns dict (theta0, theta1) -> multiplicity."""
    out = defaultdict(int)
    for mu, mult in char.items():
        vals = list(pad(mu))
        cnt = defaultdict(int)
        for x in vals:
            cnt[x] += 1
        distinct = sorted(cnt)
        for u in distinct:
            for v in distinct:
                if u == v and cnt[u] < 2:
                    continue
                rest = dict(cnt)
                rest[u] -= 1
                rest[v] -= 1
                arr = factorial(RANK - 2)
                nonzero = 0
                for val, c in rest.items():
                    arr //= factorial(c)
                    if val:
                        nonzero += c
                arr *= 2 ** nonzero
                for su in ((1, -1) if u else (1,)):
                    for sv in ((1, -1) if v else (1,)):
                        out[(su * u, sv * v)] += mult * arr
    return dict(out)

def hodge_from_spec(spec, n):
    """h^{p,q} with p = theta0 + theta1 + n, q = theta0 - theta1 + n (GKLR (2.31))."""
    h = defaultdict(int)
    for (t0, t1), m in spec.items():
        h[(t0 + t1 + n, t0 - t1 + n)] += m
    return dict(h)

def betti_from_hodge(h, n):
    b = [0] * (4 * n + 1)
    for (p, q), m in h.items():
        b[p + q] += m
    return b

def gs_hodge(nmax):
    """Goettsche-Soergel: dict n -> {(p,q): h^{p,q}(K3^[n])}."""
    # generating series in s, t, q, truncated in q
    ser = {(0, 0, 0): 1}
    k3 = {(0, 0): 1, (2, 0): 1, (0, 2): 1, (1, 1): 20, (2, 2): 1}
    for m in range(1, nmax + 1):
        for (p, q), mult in k3.items():
            for _ in range(mult):
                new = defaultdict(int)
                for (n, a, b), c in ser.items():
                    j = 0
                    while n + j * m <= nmax:
                        new[(n + j * m, a + j * (p + m - 1), b + j * (q + m - 1))] += c
                        j += 1
                ser = dict(new)
    out = defaultdict(dict)
    for (n, a, b), c in ser.items():
        out[n][(a, b)] = c
    return out

# ---------------------------------------------------------------- so(5) = B_2 decomposition
def b2_decompose(spec):
    """Decompose a W(B_2)-invariant character in (x_0, x_1) into irreducibles.
    Uses doubled exponents.  rho = (3/2, 1/2)."""
    # Weyl denominator sum_w sgn(w) x^{w rho}, doubled exponents (3,1)
    den = {}
    for (a, b, perm_sign) in ((3, 1, 1), (1, 3, -1)):
        for sa in (1, -1):
            for sb in (1, -1):
                den[(sa * a, sb * b)] = perm_sign * sa * sb
    prod = defaultdict(int)
    for (t0, t1), m in spec.items():
        for (a, b), s in den.items():
            prod[(2 * t0 + a, 2 * t1 + b)] += m * s
    out = {}
    for (a, b), c in prod.items():
        if c and a > b > 0:
            lam = ((a - 3) // 2, (b - 1) // 2)
            out[lam] = c
    return out

# ---------------------------------------------------------------- reporting helpers
def fmt(decomp):
    items = sorted(decomp.items(), key=lambda kv: pad(kv[0]), reverse=True)
    return " + ".join(("%d*" % m if m != 1 else "") + "V" + str(tuple(lam) if lam else "(0)") for lam, m in items)

def check(name, cond):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", name))
    return cond

def main():
    ok = True
    print("== T1: control, GKLR Cor 3.2 for K3^[n] ==")
    gklr = {
        2: {(2,): 1},
        3: {(3,): 1, (1, 1): 1},
        4: {(4,): 1, (2, 1): 1, (2,): 1, (): 1},
        5: {(5,): 1, (3, 1): 1, (3,): 1, (2, 1): 1, (1, 1): 1, (1,): 1},
        6: {(6,): 1, (4, 1): 1, (4,): 1, (3, 1): 1, (3,): 1, (2, 2): 1, (2, 1): 1,
            (2,): 2, (1, 1, 1): 1, (1,): 1, (): 1},
    }
    hilb = {}
    hilb_dec = {}
    for n in range(2, NMAX + 1):
        hilb[n] = ch_hilb(n)
        hilb_dec[n] = decompose(hilb[n])
    for n in range(2, 7):
        print("  K3^[%d] = %s" % (n, fmt(hilb_dec[n])))
        ok &= check("K3^[%d] matches GKLR Cor 3.2" % n, hilb_dec[n] == gklr[n])
    euler = [1, 24, 324, 3200, 25650, 176256, 1073720, 5930496, 30178575, 143184000]
    for n in range(2, NMAX + 1):
        ok &= check("dim H(K3^[%d]) = %d (Euler number)" % (n, euler[n]),
                    total_dim(hilb[n]) == euler[n]
                    and sum(m * weyl_dim(l) for l, m in hilb_dec[n].items()) == euler[n])

    print("== T2: control, torus character vs Goettsche-Soergel Hodge numbers ==")
    gs = gs_hodge(NMAX)
    for n in (2, 3, 5, 9):
        h = hodge_from_spec(spec2(hilb[n]), n)
        ok &= check("Hodge numbers of K3^[%d] agree" % n,
                    {k: v for k, v in h.items() if v} == {k: v for k, v in gs[n].items() if v})

    print("== T3: q = 2 (O'Grady 10 singular model), unconditional control ==")
    lam2_2 = ch_square(2, -1)
    sym2_2 = ch_square(2, +1)
    chi2 = char_sub(hilb[5], lam2_2)
    dec2 = decompose(chi2)
    print("  chi_2 = ch H(K3^[5]) - Lambda^2 ch H(K3^[2]) = %s" % fmt(dec2))
    predicted2 = {(5,): 1, (3,): 1, (2, 1): 1, (1,): 1}
    ok &= check("chi_2 = V(5)+V(3)+V(2,1)+V(1) (hand prediction)", dec2 == predicted2)
    # independent route: OG10 = V_(5)+V_(2,2) of D_13, branching D_13 -> B_12 interlacing
    og10_res = char_add(*[irrep((k,)) for k in range(0, 6)],
                        irrep((2, 2)), irrep((2, 1)), irrep((2,)))
    og_route = char_sub(char_sub(og10_res, sym2_2), hilb[2])
    ok &= check("OG10 route gives the same B_12 character", og_route == chi2)
    ok &= check("dim OG10 restricted = 176904", total_dim(og10_res) == 176904)
    ok &= check("Sym^2/Lambda^2 of V(2): Sym^2 = V(4)+V(2,2)+V(2)+V(0), Lambda^2 = V(3,1)+V(1,1)",
                decompose(sym2_2) == {(4,): 1, (2, 2): 1, (2,): 1, (): 1}
                and decompose(lam2_2) == {(3, 1): 1, (1, 1): 1})
    h2 = hodge_from_spec(spec2(chi2), 5)
    b2 = betti_from_hodge(h2, 5)
    print("  ib_{2j}(q=2):", b2[0::2])
    ok &= check("ib(q=2) = 1,23,277,2346,15756,87124,... (dCRS-derived, LG2 A1')",
                b2[0::2] == [1, 23, 277, 2346, 15756, 87124, 15756, 2346, 277, 23, 1])
    # Hodge numbers of OG10 (dCRS Thm A), upper half rows k = 0..10 (p = 0..k), symmetric
    og10_h = {0: [1], 2: [1, 22, 1], 4: [1, 22, 254, 22, 1], 6: [1, 23, 276, 2299, 276, 23, 1],
              8: [1, 22, 276, 2531, 16490, 2531, 276, 22, 1],
              10: [1, 22, 254, 2299, 16490, 88024, 16490, 2299, 254, 22, 1]}
    def gsnum(n, p, q):
        return gs[n].get((p, q), 0)
    def sym2_h(p, q):
        # Hodge numbers of Sym^2 K3^[2] = invariant part of H(N x N); all classes even
        tot = 0
        diag = 0
        for a in range(0, 5):
            for b in range(0, 5):
                c, d = p - a, q - b
                if 0 <= c <= 4 and 0 <= d <= 4:
                    tot += gsnum(2, a, b) * gsnum(2, c, d)
        if p % 2 == 0 and q % 2 == 0:
            diag = gsnum(2, p // 2, q // 2)
        return (tot + diag) // 2
    good = True
    for k, row in og10_h.items():
        for p, val in enumerate(row):
            q = k - p
            ih = val - sym2_h(p - 1, q - 1) - gsnum(2, p - 3, q - 3)
            if ih != h2.get((p, q), 0):
                good = False
    ok &= check("intersection Hodge numbers (q=2) agree with OG10 - Sym^2 - N route", good)

    print("== T4: q = 4, first example (quartic K3, v_0 = (5,H,0), dim 18) ==")
    lam2_3 = ch_square(3, -1)
    sym2_3 = ch_square(3, +1)
    chi4 = char_sub(hilb[9], lam2_3)
    neg_weights = [mu for mu, v in chi4.items() if v < 0]
    ok &= check("all weight multiplicities of chi_4 are nonnegative", not neg_weights)
    dec4 = decompose(chi4)
    print("  H(K3^[9])           = %s" % fmt(hilb_dec[9]))
    print("  Lambda^2 H(K3^[3])  = %s" % fmt(decompose(lam2_3)))
    print("  chi_4               = %s" % fmt(dec4))
    ok &= check("B_12 decomposition of chi_4 has nonnegative multiplicities",
                all(m > 0 for m in dec4.values()))
    ok &= check("Verbitsky component V(9) occurs exactly once", dec4.get((9,), 0) == 1)
    ok &= check("dimension = 138065600", total_dim(chi4) == 138065600
                and sum(m * weyl_dim(l) for l, m in dec4.items()) == 138065600)
    spec4 = spec2(chi4)
    h4 = hodge_from_spec(spec4, 9)
    b4 = betti_from_hodge(h4, 9)
    print("  ib_{2j}:", b4[0:19:2])
    lg2 = [1, 23, 300, 2876, 22427, 150053, 887710, 4727859, 22432662, 81617778]
    ok &= check("ib_{2j} agree with printed Betti table", b4[0:19:2] == lg2)
    ok &= check("odd intersection Betti numbers vanish", all(x == 0 for x in b4[1::2]))
    ok &= check("ib_2 = 23", b4[2] == 23)
    sym_dims = [comb(22 + k, k) for k in range(0, 10)]
    ok &= check("ib_{2k} >= dim Sym^k H^2 for k <= 9", all(b4[2 * k] >= sym_dims[k] for k in range(10)))
    row18 = [h4.get((p, 18 - p), 0) for p in range(19)]
    print("  ih^{p,18-p}:", row18)
    ok &= check("middle Hodge row agrees with the printed Hodge row",
                row18 == [1, 21, 254, 2278, 16723, 105841, 594366, 3012475, 13605820, 46942220,
                          13605820, 3012475, 594366, 105841, 16723, 2278, 254, 21, 1])
    so5 = b2_decompose(spec4)
    ok &= check("so(5) decomposition of the Hodge character is nonnegative",
                all(m > 0 for m in so5.values()))
    print("  so(5) types occurring: %d, largest highest weight %s" % (len(so5), max(so5)))

    print("== T5: failure control, Sym^2 variant ==")
    bad = char_sub(hilb[9], sym2_3)
    hb = hodge_from_spec(spec2(bad), 9)
    bb = betti_from_hodge(hb, 9)
    decb = decompose(bad)
    fails = (bb[6] != 2876) or any(m < 0 for m in decb.values()) or any(v < 0 for v in bad.values())
    print("  Sym^2 variant: ib_6 = %d, negative B_12 multiplicities: %s" %
          (bb[6], sorted([l for l, m in decb.items() if m < 0], reverse=True)[:6]))
    ok &= check("Sym^2 variant is rejected (negative multiplicity)", any(m < 0 for m in decb.values()))
    bad2 = char_sub(hilb[5], sym2_2)
    decb2 = decompose(bad2)
    ok &= check("q = 2 Sym^2 variant is rejected", any(m < 0 for m in decb2.values()))

    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())

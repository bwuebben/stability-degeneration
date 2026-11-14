#!/usr/bin/env python3
"""Exact check of the corank-one character argument in Section 10.

The peeling algorithm sees only the restriction of
ch H(K3^[9])-Lambda^2 ch H(K3^[3]) to the torus mu_12=0.
It must recover the 23 highest weights in the manuscript and leave zero
remainder. These finite checks supplement the printed proof; they do not
establish the geometric Hodge-structure hypothesis (H3).

Failure controls: two-torus data do not determine the representation;
the highest-weight length hypothesis is necessary; length-12 ambiguity
must be detected. Uses ih_llv_test.py, pure Python and exact arithmetic.
"""
import sys, time
from collections import defaultdict
from itertools import product

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import ih_llv_test as T

RANK = T.RANK
pad, strip = T.pad, T.strip


def restrict(char):
    """Restriction of a W(B_12)-invariant character (dominant weight -> mult)
    to the torus mu_12 = 0, as a W(B_11)-invariant function:
    dominant 11-tuple -> multiplicity of that weight."""
    out = defaultdict(int)
    for nu, m in char.items():
        vals = list(pad(nu))
        for v in sorted(set(vals)):
            rest = list(vals)
            rest.remove(v)
            kappa = tuple(sorted(rest, reverse=True))
            out[kappa] += m * (2 if v else 1)
    return {k: v for k, v in out.items() if v}


_res_cache = {}
def res_irrep(lam):
    lam = strip(lam)
    if lam not in _res_cache:
        _res_cache[lam] = restrict(T.irrep(lam))
    return _res_cache[lam]


def restricted_peel(res):
    """Decompose using the restricted character only (corank-one peeling algorithm).
    Returns (decomposition, ambiguous_weights, remainder)."""
    work = dict(res)
    out = {}
    ambiguous = []
    while work:
        kappa = max(work)
        if kappa[-1] != 0:
            ambiguous.append(kappa)
            break
        lam = strip(kappa)
        m = work[kappa]
        out[lam] = m
        for mu, c in res_irrep(lam).items():
            work[mu] = work.get(mu, 0) - m * c
            if work[mu] == 0:
                del work[mu]
    return out, ambiguous, work


def check(name, cond):
    print(('PASS  ' if cond else 'FAIL  ') + name)
    return cond


def b2_restricted(lam):
    """Restricted character (to mu_2 = 0) of the B_2 module V_lam, by the Weyl
    character formula evaluated as a Laurent polynomial in t = x_1 (x_2 = 1
    by a limit: computed from the weight multiplicities via Freudenthal)."""
    # weights of B_2 irreps by brute force Freudenthal in rank 2
    from fractions import Fraction
    rho = (Fraction(3, 2), Fraction(1, 2))
    roots = [(1, 0), (0, 1), (1, 1), (1, -1)]
    def dom(w):
        return tuple(sorted((abs(x) for x in w), reverse=True))
    def below(mu, la):
        return mu[0] <= la[0] and mu[0] + mu[1] <= la[0] + la[1]
    size = sum(lam)
    doms = [(a, b) for a in range(size + 1) for b in range(a + 1) if below((a, b), lam)
            and (a + b - sum(lam)) % 1 == 0]
    doms.sort(reverse=True)
    nl = sum((a + r) ** 2 for a, r in zip(lam, rho))
    mult = {tuple(lam): 1}
    for mu in doms:
        if mu == tuple(lam):
            continue
        den = nl - sum((a + r) ** 2 for a, r in zip(mu, rho))
        if den == 0:
            continue
        tot = 0
        for al in roots:
            k = 1
            while True:
                w = (mu[0] + k * al[0], mu[1] + k * al[1])
                d = dom(w)
                if not below(d, lam):
                    break
                tot += mult.get(d, 0) * (w[0] * al[0] + w[1] * al[1])
                k += 1
        val = Fraction(2 * tot) / den
        assert val.denominator == 1
        if val:
            mult[mu] = int(val)
    out = defaultdict(int)
    for (a, b), m in mult.items():
        orbit = set()
        for sa, sb in product((1, -1), repeat=2):
            orbit.add((sa * a, sb * b))
            orbit.add((sb * b, sa * a))
        for w in orbit:
            out[w[0]] += m
    return dict(out)


def main():
    ok = True
    t0 = time.time()
    print('building characters ...', flush=True)
    hilb9 = T.ch_hilb(9)
    lam2_3 = T.ch_square(3, -1)
    chi = T.char_sub(hilb9, lam2_3)
    print('  chi built, %d dominant weights, dim %d  (%.1fs)'
          % (len(chi), T.total_dim(chi), time.time() - t0), flush=True)
    full = T.decompose(chi)
    print('  full peeling done (%.1fs)' % (time.time() - t0), flush=True)

    expected = {
        (9,): 1, (7, 1): 1, (7,): 1, (6, 1): 1, (6,): 1, (5, 2): 1, (5, 1): 1,
        (5,): 3, (4, 2): 1, (4, 1, 1): 1, (4, 1): 2, (4,): 2, (3, 2): 2,
        (3, 1): 4, (3,): 4, (2, 2, 1): 1, (2, 2): 1, (2, 1, 1): 1, (2, 1): 3,
        (2,): 3, (1, 1, 1): 1, (1, 1): 2, (1,): 4,
    }
    ok &= check('full Cartan: chi decomposes as the module (29)', full == expected)

    res = restrict(chi)
    ok &= check('restriction preserves dimension',
                sum(v * T_orbit11(k) for k, v in res.items()) == T.total_dim(chi))
    dec, amb, rem = restricted_peel(res)
    print('  restricted peeling done (%.1fs)' % (time.time() - t0), flush=True)
    ok &= check('PREDICTION: restricted peeling returns the printed decomposition', dec == expected)
    ok &= check('PREDICTION: remainder zero, no ambiguous leading weight',
                not rem and not amb)

    # C1: two-torus control
    def two_torus(mods):
        tot = defaultdict(int)
        for lam, m in mods.items():
            for k, v in T.spec2(T.irrep(lam)).items():
                tot[k] += m * v
        return {k: v for k, v in tot.items() if v}
    A = {(1, 1, 1): 1, (1,): 230}
    B = {(1, 1): 21, (): 1750}
    ok &= check('C1: V(1,1,1)+230V(1) and 21V(1,1)+1750V(0) agree on the 2-torus',
                two_torus(A) == two_torus(B))
    ra, rb = defaultdict(int), defaultdict(int)
    for lam, m in A.items():
        for k, v in res_irrep(lam).items():
            ra[k] += m * v
    for lam, m in B.items():
        for k, v in res_irrep(lam).items():
            rb[k] += m * v
    ok &= check('C1: ... but differ on the corank-one torus', dict(ra) != dict(rb))

    # C2: hypothesis needed (B_2)
    a = defaultdict(int)
    for k, v in b2_restricted((1, 1)).items():
        a[k] += v
    a[0] += 5
    b = defaultdict(int)
    for k, v in b2_restricted((1, 0)).items():
        b[k] += 3 * v
    ok &= check('C2: B_2, V(1,1)+5V(0) and 3V(1,0) agree on the torus mu_2=0',
                dict(a) == dict(b))

    # C3: ambiguity is detected
    lam12 = (1,) * 12
    dec3, amb3, rem3 = restricted_peel(res_irrep(lam12))
    ok &= check('C3: a highest weight of length 12 is reported as ambiguous',
                bool(amb3))

    print('ALL PASS' if ok else 'SOME FAILED')
    print('time %.1fs' % (time.time() - t0))
    return 0 if ok else 1


def T_orbit11(kappa):
    from math import factorial
    cnt = defaultdict(int)
    for x in kappa:
        cnt[x] += 1
    size = factorial(len(kappa))
    for v, c in cnt.items():
        size //= factorial(c)
    return size * 2 ** sum(1 for x in kappa if x)


if __name__ == '__main__':
    sys.exit(main())

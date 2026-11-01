"""Exact tensor-Cech calculations for the split lines and their twisted endpoint.

Exploratory calculation, not an independently reviewed deformation theorem.
Every coefficient is rational. Each coordinate uses the two affine charts of P1.
States 0, 1, 2 denote the first chart, second chart and their intersection.
All sections are expressed in the first-chart frame of O(d).
"""
from collections import defaultdict
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import sympy as S

LINES = ((-2, -2, 4), (1, -2, 1), (-2, 1, 1), (0, 0, 0))


def clean(a):
    return {k: v for k, v in a.items() if v}


def add(*terms):
    out = defaultdict(lambda: S.Integer(0))
    for a in terms:
        for k, v in a.items():
            out[k] += v
    return clean(out)


def scale(a, c):
    return clean({k: c*v for k, v in a.items()})


def deg(k):
    return sum(t == 2 for t in k[2])


def differential(a):
    out = defaultdict(lambda: S.Integer(0))
    for (i, j, st, ex), c in a.items():
        for r in range(3):
            if st[r] == 2:
                continue
            ss = list(st); ss[r] = 2
            sign = (-1)**(sum(t == 2 for t in st[:r]) + (st[r] == 0))
            out[i, j, tuple(ss), ex] += sign*c
    return clean(out)


def multiply(a, b):
    out = defaultdict(lambda: S.Integer(0))
    for (i, j, st, ex), c in a.items():
        for (k, m, tt, ey), v in b.items():
            if j != k:
                continue
            result = []
            for u, w in zip(st, tt):
                if u == w and u != 2:
                    result.append(u)
                elif (u, w) in ((0, 2), (2, 1)):
                    result.append(2)
                else:
                    break
            if len(result) != 3:
                continue
            sign = (-1)**sum(st[r] == 2 and tt[s] == 2
                              for r in range(3) for s in range(r))
            out[i, m, tuple(result), tuple(x+y for x, y in zip(ex, ey))] += sign*c*v
    return clean(out)


def one_ip(d, st, e):
    if st == 1 and 0 <= e <= d:
        return [(0, 1), (1, 1)]
    if st == 2 and d < e < 0:
        return [(2, 1)]
    return []


def one_h(d, st, e):
    if st != 2:
        return []
    if e >= 0:
        return [(0, -1)]
    if e <= d:
        return [(1, 1)]
    return []


def homotopy(a):
    out = defaultdict(lambda: S.Integer(0))
    for (i, j, st, ex), c in a.items():
        ds = [x-y for x, y in zip(LINES[i], LINES[j])]
        for r in range(3):
            options = [one_ip(ds[s], st[s], ex[s]) for s in range(r)]
            options += [one_h(ds[r], st[r], ex[r])]
            options += [[(st[s], 1)] for s in range(r+1, 3)]
            for entry in product(*options):
                sign = (-1)**sum(st[s] == 2 for s in range(r))
                coef = c*sign
                for _, z in entry:
                    coef *= z
                out[i, j, tuple(x[0] for x in entry), ex] += coef
    return clean(out)


def project(a):
    out = defaultdict(lambda: S.Integer(0))
    for (i, j, st, ex), c in a.items():
        ds = [x-y for x, y in zip(LINES[i], LINES[j])]
        if all((t == 1 and 0 <= e <= d) or (t == 2 and d < e < 0)
               for d, t, e in zip(ds, st, ex)):
            out[i, j, ex] += c
    return clean(out)


def include(v):
    out = {}
    for (i, j, ex), c in v.items():
        ds = [x-y for x, y in zip(LINES[i], LINES[j])]
        states = [(0, 1) if 0 <= e <= d else (2,) for d, e in zip(ds, ex)]
        for st in product(*states):
            out[i, j, st, ex] = c
    return out


def cohomology_basis():
    out = {q: [] for q in range(4)}
    for i, j in product(range(4), repeat=2):
        ds = [x-y for x, y in zip(LINES[i], LINES[j])]
        if -1 in ds:
            continue
        exps = [range(d+1) if d >= 0 else range(d+1, 0) for d in ds]
        q = sum(d <= -2 for d in ds)
        out[q].extend((i, j, ex) for ex in product(*exps))
    return out


BASIS = cohomology_basis()
A = include({(0, 1, (-2, 0, 3)): S.Integer(1)})
B = include({(1, 3, (0, -1, 1)): S.Integer(1)})
GAMMA = scale(homotopy(multiply(A, B)), -1)
N = add(A, B, GAMMA)


def perturb(a):
    out = multiply(N, a)
    for q in range(4):
        aq = {k: v for k, v in a.items() if deg(k) == q}
        out = add(out, scale(multiply(aq, N), -(-1)**q))
    return out


def twisted_d(a):
    return add(differential(a), perturb(a))


def lift(v):
    term = include(v); out = term
    for _ in range(8):
        term = scale(homotopy(perturb(term)), -1)
        if not term:
            return out
        out = add(out, term)
    raise AssertionError('Upper-triangular transfer failed to terminate')


def transferred_d(v):
    return project(perturb(lift(v)))


def transferred_project(a):
    term = a; out = project(a)
    for _ in range(8):
        term = scale(perturb(homotopy(term)), -1)
        if not term:
            return out
        out = add(out, project(term))
    raise AssertionError('Upper-triangular transfer failed to terminate')


def matrix_d(q):
    result = S.zeros(len(BASIS[q+1]), len(BASIS[q]))
    ix = {k: i for i, k in enumerate(BASIS[q+1])}
    for j, key in enumerate(BASIS[q]):
        out = transferred_d({key: S.Integer(1)})
        for k, c in out.items():
            result[ix[k], j] = c
    return result


def describe(key):
    i, j, ex = key
    return {'block': f'{i+1}{j+1}', 'exponents': list(ex)}


def run():
    assert not add(differential(N), multiply(N, N))
    tested = 0
    for q in range(4):
        for k in BASIS[q]:
            v = {k: S.Integer(1)}
            a = include(v)
            assert project(a) == v
            assert not differential(a)
            assert not add(twisted_d(lift(v)), scale(lift(transferred_d(v)), -1))
            for b in [a, perturb(a), homotopy(perturb(a))]:
                assert add(differential(homotopy(b)), homotopy(differential(b))) == add(b, scale(include(project(b)), -1))
                assert not twisted_d(twisted_d(b))
                tested += 1
    d0, d1, d2 = [matrix_d(q) for q in range(3)]
    assert d1*d0 == S.zeros(5, 4)
    ranks = [d0.rank(), d1.rank(), d2.rank()]
    nonzero = []
    for j, key in enumerate(BASIS[1]):
        if any(d1[:, j]):
            nonzero.append({'input': describe(key), 'output': [
                {'basis': describe(BASIS[2][i]), 'coefficient': str(d1[i,j])}
                for i in range(5) if d1[i,j]]})
    ans = {'status': 'exact finite checks passed; finite calculation only; analytic claims are proved in the manuscript',
           'graded_dimensions': [len(BASIS[q]) for q in range(4)],
           'transferred_differential_ranks': ranks,
           'endpoint_dimensions': [len(BASIS[q])-(ranks[q-1] if q else 0)-(ranks[q] if q<3 else 0) for q in range(4)],
           'contraction_and_chain_checks': tested,
           'nonzero_degree_one_columns': nonzero,
           'quadratic_primitive_terms': [{'block': f'{k[0]+1}{k[1]+1}', 'states':k[2], 'exponents':k[3], 'coefficient':str(v)} for k,v in GAMMA.items()]}
    dest = (Path(__file__).resolve().parents[1]/'output'/'cohomology_results.json')
    dest.write_text(json.dumps(ans, indent=2)+'\n')
    print(json.dumps(ans, indent=2))


if __name__ == '__main__':
    run()

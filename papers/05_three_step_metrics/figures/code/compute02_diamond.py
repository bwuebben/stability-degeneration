#!/usr/bin/env python3
"""Exact finite checks for R07; these do not prove its analytic estimates."""
from fractions import Fraction as F
from itertools import permutations, combinations
from pathlib import Path
import json

# Sparse rational polynomials in (lambda,k,l).
def const(x):
    return {} if x == 0 else {(0, 0, 0): F(x)}

def add(a, b):
    c = dict(a)
    for e, v in b.items():
        c[e] = c.get(e, F(0)) + v
        if not c[e]:
            del c[e]
    return c

def scale(a, s):
    return {e: v * s for e, v in a.items() if v * s}

def mul(a, b):
    c = {}
    for e, v in a.items():
        for f, w in b.items():
            g = tuple(x + y for x, y in zip(e, f))
            c[g] = c.get(g, F(0)) + v * w
    return {e: v for e, v in c.items() if v}

def det(m):
    ans = {}
    for p in permutations(range(len(m))):
        inversions = sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))
        term = const((-1) ** inversions)
        for i, j in enumerate(p):
            term = mul(term, m[i][j])
        ans = add(ans, term)
    return ans

T, K, L = ({e: F(1)} for e in [(1, 0, 0), (0, 1, 0), (0, 0, 1)])
weights = {(0, 1): scale(L, F(5, 8)), (0, 2): scale(K, F(5, 8)),
           (1, 3): scale(K, F(3, 8)), (2, 3): scale(L, F(3, 8))}
lap = [[{} for j in range(4)] for i in range(4)]
for (i, j), w in weights.items():
    for v in [i, j]:
        lap[v][v] = add(lap[v][v], w)
    lap[i][j] = lap[j][i] = scale(w, -1)
char = det([[add(T if i == j else {}, scale(lap[i][j], -1)) for j in range(4)] for i in range(4)])
expected = {(4, 0, 0): F(1), (3, 1, 0): F(-2), (3, 0, 1): F(-2),
            (2, 2, 0): F(15, 16), (2, 0, 2): F(15, 16), (2, 1, 1): F(3),
            (1, 2, 1): F(-15, 16), (1, 1, 2): F(-15, 16)}
assert char == expected, (char, expected)

factors = [(-2, -2, 4), (1, -2, 1), (-2, 1, 1), (0, 0, 0)]
def cohom(d):
    h = [0] * 4
    if -1 in d:
        return h
    q = sum(x <= -2 for x in d)
    size = 1
    for x in d:
        size *= x + 1 if x >= 0 else -x - 1
    h[q] = size
    return h

groups = []
total = [0] * 4
for i in range(4):
    for j in range(4):
        d = tuple(a - b for a, b in zip(factors[i], factors[j]))
        h = cohom(d)
        if i != j:
            assert h[0] == 0
        total = [a + b for a, b in zip(total, h)]
        groups.append({'target': i + 1, 'source': j + 1, 'degree': d, 'cohomology': h})
for i, j, h1 in [(0, 1, 8), (0, 2, 8), (1, 3, 4), (2, 3, 4)]:
    assert cohom(tuple(a - b for a, b in zip(factors[i], factors[j]))) == [0, h1, 0, 0]
assert cohom(factors[0]) == [0, 0, 5, 0]
pair_patterns = []
for j in [1, 2, 3]:
    other = [i for i in [1, 2, 3] if i != j]
    lhs = tuple(a + b for a, b in zip(factors[0], factors[j]))
    rhs = tuple(a + b for a, b in zip(factors[other[0]], factors[other[1]]))
    assert lhs != rhs
    pair_patterns.append({'pair': [1, j + 1], 'sum': lhs, 'opposite_sum': rhs})

# Slopes and moment balance, as polynomial identities.
u, v0 = scale(add(K, scale(L, 5)), F(1, 12)), scale(add(scale(K, 5), L), F(1, 12))
kappa = [add(scale(u, d[0]), scale(v0, d[1])) for d in factors]
avg = scale(add(add(kappa[0], kappa[1]), add(kappa[2], kappa[3])), F(1, 4))
assert avg == scale(add(K, L), F(-3, 8))
moment = [{} for _ in range(4)]
for (i, j), w in weights.items():
    moment[i] = add(moment[i], w)
    moment[j] = add(moment[j], scale(w, -1))
assert moment == [add(avg, scale(x, -1)) for x in kappa]
assert add(scale(add(kappa[0], kappa[1]), F(1, 2)), scale(avg, -1)) == scale(K, F(-1, 2))
assert add(scale(add(kappa[0], kappa[2]), F(1, 2)), scale(avg, -1)) == scale(L, F(-1, 2))
assert mul(weights[0, 1], weights[1, 3]) == mul(weights[0, 2], weights[2, 3])
q = {}
for w in weights.values():
    q = add(q, w)
assert q == add(K, L)

# Every cut preserves blocks of the reconstructed ambient semiconnection.
# Retained path products encode its possible (1,4) block and obstruction.
cuts = []
edges = list(weights)
for size in [1, 2, 3]:
    for subset in combinations(range(4), size):
        s = set(subset)
        keep = {e: (e[0] in s) == (e[1] in s) for e in edges}
        path1 = keep[0, 1] and keep[1, 3]
        path2 = keep[0, 2] and keep[2, 3]
        if (0 in s) != (3 in s):
            assert not path1 and not path2
        cuts.append({'subset': [i + 1 for i in subset], 'retained_paths': [int(path1), int(path2)],
                     'integrable_at_unit_coefficients': path1 == path2})
assert len(cuts) == 14
assert [1, 2, 4] in [c['subset'] for c in cuts if not c['integrable_at_unit_coefficients']]

result = {'status': 'all exact finite checks passed', 'limitations': 'No PDE, Hodge-theoretic or stability proof is automated.',
          'graded_end_cohomology': total, 'cohomology_groups': groups,
          'excluded_four_line_tensor_patterns': pair_patterns,
          'characteristic_polynomial': [{'exponents_lambda_k_l': e, 'coefficient': str(v)} for e, v in sorted(char.items())],
          'moment_balance': 'exact polynomial identity', 'subbundle_slope_differences': ['-k/2', '-l/2'],
          'all_cuts': cuts}
path = (Path(__file__).resolve().parents[1]/"output"/"compute02_diamond_results.json")
path.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'status': result['status'], 'cohomology_groups_checked': len(groups),
                  'graded_end_cohomology': total, 'cuts_checked': len(cuts),
                  'nonintegrable_cuts': [c['subset'] for c in cuts if not c['integrable_at_unit_coefficients']],
                  'output': str(path)}))

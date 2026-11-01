#!/usr/bin/env python3
"""Exact algebra for R10; no geometric classification or PDE claim is automated."""
from fractions import Fraction as F
from pathlib import Path
import json

names = ([f'{letter}{i}{j}' for letter in ['a', 'c'] for i in range(2) for j in range(4)]
         + [f'{letter}{i}{j}' for letter in ['b', 'd'] for i in range(2) for j in range(2)])
assert len(names) == 24
equations = []
for r in range(5):
    equations.append([(sign, f'{a}{i}{j}', f'{b}{i}{v}')
                      for a, b, sign in [('a', 'b', 1), ('c', 'd', -1)]
                      for i in range(2) for j in range(4) for v in range(2) if j + v == r])

def evaluate(point):
    return [sum(sign * point.get(a, 0) * point.get(b, 0) for sign, a, b in eq) for eq in equations]

def jacobian(point):
    return [[sum(sign * ((point.get(b, 0) if a == v else 0) +
                        (point.get(a, 0) if b == v else 0)) for sign, a, b in eq)
             for v in names] for eq in equations]

def rank(matrix):
    a = [[F(x) for x in row] for row in matrix]
    r = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        d = a[r][col]
        a[r] = [x / d for x in a[r]]
        for i in range(len(a)):
            if i != r:
                d = a[i][col]
                a[i] = [x - d * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r

smooth = dict.fromkeys(['a00', 'a13', 'c00', 'c13', 'b00', 'b11', 'd00', 'd11'], F(1))
old = dict.fromkeys(['a00', 'c00', 'b00', 'd00'], F(1))
assert evaluate(smooth) == evaluate(old) == [0] * 5
assert rank(jacobian(smooth)) == 5
assert rank(jacobian(old)) == 4
assert rank([[row[j] for j in range(8)] for row in jacobian(smooth)]) == 5
for m in [0, 1, 7, 1000]:
    point = {**old, 'a10': F(m)}
    assert evaluate(point) == [0] * 5

# Sparse polynomials in (k,l,a,rho).
def c(x):
    return {} if x == 0 else {(0, 0, 0, 0): F(x)}

def plus(p, q):
    r = dict(p)
    for mon, value in q.items():
        r[mon] = r.get(mon, F(0)) + value
        if not r[mon]:
            del r[mon]
    return r

def times(p, q):
    r = {}
    for x, a in p.items():
        for y, b in q.items():
            z = tuple(u + v for u, v in zip(x, y))
            r[z] = r.get(z, F(0)) + a * b
    return {z: a for z, a in r.items() if a}

def scale(p, q):
    return times(p, c(q))

def sub_a(p, value):
    out = {}
    for mon, coef in p.items():
        base = {(mon[0], mon[1], 0, mon[3]): coef}
        for _ in range(mon[2]):
            base = times(base, value)
        out = plus(out, base)
    return out

k, l, a, rho = [{tuple(int(i == j) for i in range(4)): F(1)} for j in range(4)]
sigma = plus(k, l)
p = scale(sigma, F(5, 8))
c0 = scale(plus(scale(k, 3), scale(l, -5)), F(1, 8))
weights = [a, plus(p, scale(a, -1)), plus(a, c0), plus(l, scale(a, -1))]
assert plus(plus(weights[0], weights[1]), plus(weights[2], weights[3])) == sigma
moment = [plus(weights[0], weights[1]), plus(scale(weights[0], -1), weights[2]),
          plus(scale(weights[1], -1), weights[3]), scale(plus(weights[2], weights[3]), -1)]
demand = [scale(sigma, F(5, 8)), c0, scale(plus(scale(k, -5), scale(l, 3)), F(1, 8)), scale(sigma, F(-3, 8))]
assert moment == demand
cycle = plus(times(weights[0], weights[2]), scale(times(rho, times(weights[1], weights[3])), -1))
expected = plus(plus(times(plus(c(1), scale(rho, -1)), times(a, a)),
                    times(plus(c0, times(rho, plus(p, l))), a)), scale(times(rho, times(p, l)), -1))
assert cycle == expected
special = [sub_a(w, scale(l, F(5, 8))) for w in weights]
assert special == [scale(l, F(5, 8)), scale(k, F(5, 8)), scale(k, F(3, 8)), scale(l, F(3, 8))]
assert times(special[0], special[2]) == times(special[1], special[3])
trees = [[sub_a(w, value) for w in weights] for value in [{}, p, scale(c0, -1), l]]
for i, tree in enumerate(trees):
    assert tree[i] == {}
    assert plus(plus(tree[0], tree[1]), plus(tree[2], tree[3])) == sigma

# The extra weak-mode Rayleigh quotient is an exact rational identity at samples.
weak = []
xi = [F(-1, 3), F(-1, 3), F(1), F(-1, 3)]
edges = [(0, 1), (0, 2), (1, 3), (2, 3)]
mass = sum(x * x for x in xi)
for delta in [F(1, 10), F(1, 100), F(1, 10000)]:
    w = [1 - delta, delta, F(3, 5) - delta, delta]
    invariant = w[0] * w[2] / (w[1] * w[3])
    rayleigh = sum(energy * (xi[i] - xi[j]) ** 2 for energy, (i, j) in zip(w, edges)) / mass
    assert rayleigh == F(8, 3) * delta
    assert sum(w) == F(8, 5)
    weak.append({'delta_over_l': str(delta), 'rho': str(invariant), 'trial_Rayleigh_over_l': str(rayleigh)})

result = {'status': 'all exact finite checks passed', 'ambient_complex_dimension': 24,
          'obstruction_components': [{'quartic_t_degree': r, 'terms': eq} for r, eq in enumerate(equations)],
          'rank_at_smooth_point': rank(jacobian(smooth)), 'rank_at_R09_direction': rank(jacobian(old)),
          'smooth_locus_local_dimension': 19, 'effective_torus_dimension': 3,
          'balance_and_cycle_equation': 'exact polynomial identities', 'R09_weights_recovered': True,
          'connected_tree_supports_checked': 4, 'weak_mode_exact_samples': weak,
          'limitations': 'Local quotient interpretation, bundle isomorphism classification, stable sufficiency and PDE estimates require written proofs and independent reviews.'}
out = (Path(__file__).resolve().parents[1]/"output"/"compute03_family_results.json")
out.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({key: result[key] for key in ['status', 'ambient_complex_dimension', 'rank_at_smooth_point',
                                             'rank_at_R09_direction', 'connected_tree_supports_checked']}))

"""Exact six-dimensional stress contraction and scalar curvature vertex.

These algebraic checks do not evaluate the full next heterotic coefficient.
Run from the paper directory: python figures/code/compute06_gravity.py
"""
from itertools import permutations
import json
from pathlib import Path

import sympy as sp

c = sp.symbols('c', real=True)
eta = sp.diag(-1, 1, 1, 1, 1, 1)
z = sp.sqrt(1-c*c)
p1 = sp.Matrix([1, 0, 0, 0, 0, 1])
p2 = sp.Matrix([1, 0, 0, 0, 0, -1])
p3 = sp.Matrix([1, z, 0, 0, 0, c])
p4 = sp.Matrix([1, -z, 0, 0, 0, -c])


def dot(p, q):
    return (p.T*eta*q)[0]


def stress(p, q):
    pl, ql = eta*p, eta*q
    return pl*ql.T+ql*pl.T-eta*dot(p, q)


def zero(expr):
    assert sp.simplify(expr) == 0, expr


T1, T2 = stress(p1, p3), stress(p2, p4)
contraction = sp.trace(T1*eta*T2*eta)-sp.trace(eta*T1)*sp.trace(eta*T2)/4
s0, t0, u0 = -dot(p1+p2, p1+p2), -dot(p1-p3, p1-p3), -dot(p1-p4, p1-p4)
zero(contraction+s0*u0)
for p in (p1, p2, p3, p4):
    zero(dot(p, p))
zero(s0+t0+u0)
for val in T1*(p1-p3):
    zero(val)

# All-incoming scalar momenta and their pairings. Differentiate the
# normal-coordinate quartic by enumerating external-leg assignments.
s, t, u, K = sp.symbols('s t u K', real=True)
products = {(0, 1): -s/2, (0, 2): -t/2, (0, 3): -u/2,
            (1, 2): -u/2, (1, 3): -t/2, (2, 3): -s/2}
species = ['phi', 'chi', 'phi', 'chi']
terms = [(K/6, [('chi', False), ('chi', False), ('phi', True), ('phi', True)]),
         (K/6, [('phi', False), ('phi', False), ('chi', True), ('chi', True)]),
         (-K/3, [('phi', False), ('chi', False), ('phi', True), ('chi', True)])]
curvature_vertex = 0
for coeff, factors in terms:
    for assignment in permutations(range(4)):
        if any(species[leg] != factor[0] for leg, factor in zip(assignment, factors)):
            continue
        derivative_legs = sorted(leg for leg, factor in zip(assignment, factors) if factor[1])
        curvature_vertex -= coeff*products[tuple(derivative_legs)]
zero((curvature_vertex-K*t).subs(u, -s-t))

x = sp.symbols('x', positive=True)
model = x/2/(1-x)-5*x/4/(1+3*x/4)
zero(sp.limit(model/x, x, 0)+sp.Rational(3, 4))
zero(sp.limit(model, x, sp.oo)+sp.Rational(13, 6))
record = {'status': 'PASS', 'sympy_version': sp.__version__, 'signature': 'mostly plus',
          'dimension': 6, 'stress_contraction': str(sp.simplify(contraction)),
          'stress_identity': 'T1 P T2 = -s*u',
          'normal_coordinate_vertex': str(sp.simplify(curvature_vertex.subs(u, -s-t))),
          'figure_low_energy_slope': '-3/4', 'figure_intermediate_limit': '-13/6',
          'scope': 'Exact tensor and finite vertex identities; not the full heterotic correction.'}
out = Path(__file__).resolve().parents[1]/'data'/'gravity_and_vertex.json'
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(record, indent=2)+'\n')
print(json.dumps(record, indent=2))

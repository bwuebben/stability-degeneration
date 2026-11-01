"""Exact identities for the illustrative phase system; no PDE simulation."""
import sympy as sp

t, eps, s, a, y0, v0 = sp.symbols("t eps s a y0 v0", positive=True)
v = v0 * sp.exp(-s*t)
phase = y0 + a*v0/(s-eps)
selected = phase*sp.exp(-eps*t)
y = selected-a*v0*sp.exp(-s*t)/(s-eps)
assert sp.simplify(sp.diff(y,t)+eps*y-a*v) == 0
assert sp.simplify(sp.diff(v,t)+s*v) == 0
assert sp.simplify(y.subs(t,0)-y0) == 0
assert sp.simplify(selected-y-a*v0*sp.exp(-s*t)/(s-eps)) == 0
assert sp.diff(phase,v0) == a/(s-eps)
print("PASS: exact solution, initial value, selected phase, fast error and normal response.")

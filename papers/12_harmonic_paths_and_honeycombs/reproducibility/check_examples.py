# /// script
# requires-python = ">=3.11"
# dependencies = ["sympy==1.14.0"]
# ///
"""Exact finite calculations in Section 8. Run: uv run check_examples.py.

The valued family uses b=T^s, 0<s<1. Every nonzero rational function
of b has valuation ord_b(f)*s. Affine inequalities are checked at both
ends of [0,1], thereby certifying the entire open interval, not samples.
This checks finite algebra and support geometry, not the analytic theorem.
"""
from fractions import Fraction as F
from itertools import combinations, product
import json
import sympy as sp

b = sp.symbols("b")
w = (-2, -1, 0)
bases = [sp.eye(3), sp.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]]),
         sp.Matrix([[1, 1, 0], [b, 2, 0], [1, 3, 1]])]


def order(expr):
    expr = sp.cancel(expr)
    if expr == 0:
        return None
    num, den = sp.fraction(expr)
    def least_power(p):
        return min(monomial[0] for monomial, coeff in sp.Poly(p, b).terms() if coeff)
    return least_power(num) - least_power(den)


def certify_minimum(candidates, expected):
    """Each tuple (slope,constant) represents a real affine function."""
    assert expected in candidates
    ea, ec = expected
    for a, c in candidates:
        assert c - ec >= 0 and a - ea + c - ec >= 0, (a, c, expected)


profiles = {}
minor_checks = 0
for i, j in product(range(3), repeat=2):
    if i == j:
        continue
    M = bases[j].inv() * bases[i]
    minima = [0]
    for k, target in enumerate((-2, -2, 0), start=1):
        candidates = []
        for rows in combinations(range(3), k):
            for cols in combinations(range(3), k):
                a = order(M.extract(rows, cols).det())
                if a is not None:
                    candidates.append((a, sum(w[c] for c in cols) - sum(w[r] for r in rows)))
        certify_minimum(candidates, (0, target))
        minor_checks += len(candidates)
        minima.append(target)
    logs = [minima[k-1] - minima[k] for k in range(1, 4)]
    assert logs == [2, 0, -2]
    profiles[f"{i+1}{j+1}"] = logs

triple_candidates = []
for cols in product(range(3), repeat=3):
    mat = sp.Matrix.hstack(*(bases[i][:, cols[i]] for i in range(3)))
    a = order(mat.det())
    if a is not None:
        triple_candidates.append((a, sum(w[c] for c in cols)))
certify_minimum(triple_candidates, (1, -6))

# Coordinates below are (x/p, sqrt(3)*t/p), p>0.
r = (-1, 0, 1)
vertices = {f"Y{i+1}": (F(ri), F(0)) for i, ri in enumerate(r)}
for i, j in combinations(range(3), 2):
    vertices[f"P{i+1}{j+1}"] = (F(r[i] + r[j], 2), F(r[j] - r[i], 2))
finite_edges = set()
for slope in (-1, 1):
    for ri in r:
        incident = sorted((y, name) for name, (x, y) in vertices.items()
                          if x == ri + slope*y)
        for (_, v), (_, v_next) in zip(incident, incident[1:]):
            finite_edges.add(tuple(sorted((v, v_next))))
assert len(vertices) == 6 and len(finite_edges) == 6
seen = {"Y1"}
while True:
    more = seen | {v for edge in finite_edges if set(edge) & seen for v in edge}
    if more == seen:
        break
    seen = more
assert seen == set(vertices)
assert len(finite_edges) - len(vertices) + 1 == 1
cycle = [vertices[v] for v in ("Y2", "P12", "P13", "P23")]
area = abs(sum(cycle[i][0]*cycle[(i+1)%4][1] - cycle[(i+1)%4][0]*cycle[i][1]
               for i in range(4))) / 2
assert area == F(1, 2)

print(json.dumps({
    "scope": "Exact finite algebra for every 0<s<1; normalized support for p>0",
    "pairwise_spectra": profiles,
    "nonzero_minor_candidates_checked": minor_checks,
    "triple_determinant_candidates_checked": len(triple_candidates),
    "D111": "6-s",
    "vertices": {k: [str(x), str(y)] for k, (x, y) in vertices.items()},
    "finite_edges": sorted(finite_edges),
    "exterior_rays": 9,
    "bounded_faces": 1,
    "normalized_euclidean_area": str(area),
    "status": "passed; no analytic theorem is certified by this calculation"
}, indent=2))

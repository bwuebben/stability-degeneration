#!/usr/bin/env python3
"""Exact finite checks for R02. These do not verify any analytic theorem."""
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
import json

def cohomology(degrees):
    out = [1]
    for d in degrees:
        h = [max(d + 1, 0), max(-d - 1, 0)]
        nxt = [0] * (len(out) + 1)
        for i, a in enumerate(out):
            for j, b in enumerate(h):
                nxt[i + j] += a * b
        out = nxt
    return out

factors = [(-1, -1, 2), (1, -2, 1), (0, 0, 0)]
hom = {}
for i, j in product(range(3), repeat=2):
    d = tuple(a-b for a, b in zip(factors[i], factors[j]))
    hom[f"{j+1}->{i+1}"] = {"degrees": d, "h0_h1_h2_h3": cohomology(d)}
    if i != j:
        assert cohomology(d)[0] == 0
assert hom["2->1"]["h0_h1_h2_h3"] == [0, 4, 0, 0]
assert hom["3->2"]["h0_h1_h2_h3"] == [0, 4, 0, 0]
assert hom["3->1"]["h0_h1_h2_h3"] == [0, 0, 0, 0]

# Sparse polynomials in (lambda,k1,k2), computed over the integers.
def add(*ps):
    out = {}
    for p in ps:
        for e, c in p.items():
            out[e] = out.get(e, 0) + c
    return {e:c for e,c in out.items() if c}

def scale(p, c):
    return {e:a*c for e,a in p.items() if a*c}

def mul(p, q):
    out = {}
    for e, a in p.items():
        for f, b in q.items():
            g = tuple(x+y for x,y in zip(e,f))
            out[g] = out.get(g, 0) + a*b
    return {e:c for e,c in out.items() if c}

lam, k1, k2 = ({e:1} for e in [(1,0,0),(0,1,0),(0,0,1)])
matrix = [[add(lam,scale(k1,-1)), k1, {}],
          [k1, add(lam,scale(k1,-1),scale(k2,-1)), k2],
          [{}, k2, add(lam,scale(k2,-1))]]
det = {}
for perm in permutations(range(3)):
    inversions = sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
    term = {(0,0,0):(-1)**inversions}
    for i,j in enumerate(perm):
        term = mul(term,matrix[i][j])
    det = add(det,term)
assert det == {(3,0,0):1,(2,1,0):-2,(2,0,1):-2,(1,1,1):3}

# Exact demand map using a=k1/(2pi), b=k2/(2pi).
average = [sum(F(d[i],3) for d in factors) for i in range(3)]
samples = [(F(1,10),F(1,10)),(F(1,10**4),F(1,10)),
           (F(1,10),F(1,10**20)),(F(-1,10),F(1,10))]
for a,b in samples:
    inverse_h = [1+a,1+b,F(1)]
    demands = [-sum((F(d[i])-average[i])*inverse_h[i] for i in range(3)) for d in factors]
    assert demands == [a,b-a,-b]

# Coefficient of alpha wedge beta in (dbar N + N^2)_13.
# Raw projection keeps the old induced block; recompletion recomputes it.
cuts = []
for mask in range(1,7):
    inside = [bool(mask & (1<<i)) for i in range(3)]
    keep12,keep23,keep13 = (int(inside[i]==inside[j]) for i,j in [(0,1),(1,2),(0,2)])
    raw = -keep13 + keep12*keep23
    completed = -(keep12*keep23) + keep12*keep23
    assert completed == 0
    cuts.append({"vertices":[i+1 for i in range(3) if inside[i]],
                 "raw_truncation_MC_coefficient":raw,
                 "recompleted_MC_coefficient":completed})
assert [x["vertices"] for x in cuts if x["raw_truncation_MC_coefficient"]] == [[2],[1,3]]

out = {"scope":"Exact finite algebra only; no PDE or proof validation",
       "hom_cohomology":hom,
       "characteristic_polynomial":"lambda^3 - 2*(k1+k2)*lambda^2 + 3*k1*k2*lambda",
       "exact_slope_samples":len(samples),"projection_truncation":cuts,
       "result":"PASS"}
dest = (Path(__file__).resolve().parents[1]/"output"/"compute01_chain_results.json")
dest.write_text(json.dumps(out,indent=2)+"\n")
print(f"PASS: 9 cohomology vectors, symbolic characteristic polynomial, {len(samples)} exact slope checks, 6 cuts.")
print("Raw truncation fails at cuts {2} and {1,3}; recompletion preserves the integrability identity at all six cuts.")

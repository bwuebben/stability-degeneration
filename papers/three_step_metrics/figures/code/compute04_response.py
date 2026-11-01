#!/usr/bin/env python3
"""Exact finite response identities; does not verify geometric/PDE estimates."""
from pathlib import Path
import json
import sympy as S

t=S.symbols("tau", positive=True)
D=S.Matrix([[1,1,0,0],[-1,0,1,0],[0,-1,0,1],[0,0,-1,-1]])
B=S.Matrix([[5,5],[3,-5],[-5,3],[-3,-3]])/8
L=D*S.diag(1-t,t,S.Rational(3,5)-t,t)*D.T
A=S.simplify((L+S.ones(4)/4).inv()*B)
assert S.simplify(L*A-B)==S.zeros(4,2)
assert S.ones(1,4)*A==S.zeros(1,2) or S.simplify(S.ones(1,4)*A)==S.zeros(1,2)
assert B.T*B==S.Matrix([[17,1],[1,17]])/16
G=S.simplify(A.T*A)
tr=(12800*t**4-19840*t**3+8904*t**2-816*t+459)/(512*t**2*(4*t-3)**2)
det=25*(64*t**2-48*t+27)/(512*t**2*(4*t-3)**2)
assert S.simplify(S.trace(G)-tr)==0
assert S.simplify(G.det()-det)==0
residue=S.Matrix([[5,-3],[5,-3],[-15,9],[5,-3]])/64
assert A.applyfunc(lambda x:S.limit(t*x,t,0))==residue
assert S.simplify(A*S.Matrix([3,5]))==S.Matrix([5,0,0,-5])
assert S.limit(t*t*tr,t,0)==S.Rational(51,512)
assert S.limit(det/tr,t,0)==S.Rational(25,17)
u=S.Matrix([-1,-1,3,-1])/S.sqrt(12)
v=S.Matrix([-5,3])/S.sqrt(34)
assert S.simplify(residue*v-S.sqrt(102)*u/32)==S.zeros(4,1)
assert (u.T*u)[0]==1 and (v.T*v)[0]==1
out={
 "date":"2026-10-30",
 "status":"all exact finite response identities passed",
 "sympy_version":S.__version__,
 "incidence":D.tolist(),
 "polarization_derivative":[[str(x) for x in row] for row in B.tolist()],
 "scaled_response":[[str(x) for x in row] for row in A.tolist()],
 "gram_trace":str(S.factor(tr)),
 "gram_determinant":str(S.factor(det)),
 "tau_squared_trace_limit":"51/512",
 "smaller_singular_value_squared_limit":"25/17",
 "weak_residue":[[str(x) for x in row] for row in residue.tolist()],
 "hard_direction_identity":"A0(tau)*(3,5)=(5,0,0,-5)",
 "scope":"Exact finite algebra only; no PDE or geometric-derivative verification."
}
out["incidence"]=[[int(x) for x in row] for row in out["incidence"]]
p=(Path(__file__).resolve().parents[1]/"output"/"compute04_response_results.json")
p.write_text(json.dumps(out,indent=2)+"\n")
print(out["status"])

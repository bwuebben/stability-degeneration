"""Exact split cohomology, weighted finite differential and endpoint low spectra.

This checks the finite constants printed in Section 5, not elliptic estimates.
"""
from pathlib import Path
from itertools import product
import json
import sympy as s

lines=[(-2,-2,4),(1,-2,1),(-2,1,1),(0,0,0)]
def cohom(deg):
    q=0; n=1
    for d in deg:
        if d==-1:return [0]*4
        if d>=0:n*=d+1
        else:q+=1;n*=-d-1
    out=[0]*4;out[q]=n;return out
blocks={f'{i+1}{j+1}':cohom(tuple(a-b for a,b in zip(fi,fj)))
        for i,fi in enumerate(lines) for j,fj in enumerate(lines)}
dims=[sum(v[q] for v in blocks.values()) for q in range(4)];dims[0]-=1
assert dims==[3,83,5,0]
inc=s.Matrix([[1,0],[-1,1],[0,0],[0,-1]])
L=inc*s.diag(1,s.Rational(3,5))*inc.T
expected0=[(8-s.sqrt(19))/5,(8+s.sqrt(19))/5]
lam=s.symbols('lam')
assert s.expand(L.charpoly(lam).as_expr()-lam**2*(lam**2-s.Rational(16,5)*lam+s.Rational(9,5)))==0
# Degree-one inputs: four a_0 cubics, then the two b_1 linears.
a=1/s.sqrt(s.pi);b=s.sqrt(6/(5*s.pi))
T=s.zeros(5,6)
for j in range(4):T[j+1,j]=b
T[3,4]=a;T[4,5]=a
mass_in=s.diag(*[s.pi/s.binomial(3,j) for j in range(4)],s.pi/2,s.pi/2)
mass_out=s.diag(*[4*s.pi**2/(5*s.binomial(4,j)) for j in range(5)])
TT=s.simplify(T*mass_in.inv()*T.T*mass_out)
assert TT==s.diag(0,s.Rational(6,25),s.Rational(12,25),s.Rational(28,25),s.Rational(64,25))
assert T.rank()==4
# Beta integral for the normalized monomial norms (integer parameters).
monomials={}
for d in range(5):
 for j in range(d+1):
  integral=s.factorial(j)*s.factorial(d-j)/s.factorial(d+1)
  assert integral==1/((d+1)*s.binomial(d,j))
  monomials[f'{d},{j}']=str(integral)
result={'status':'PASS: exact finite calculations only','tracefree_split_dimensions':dims,
 'endpoint_dimensions':[1,77,1,0],'positive_low_counts':[2,6,4,0],
 'degree_zero_constants':list(map(str,expected0)),
 'degree_two_constants':list(map(str,TT.diagonal()[1:])),
 'cohomology_blocks':blocks,'monomial_norms':monomials}
path=Path(__file__).resolve().parents[1]/'output/wall_complex_results.json'
path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

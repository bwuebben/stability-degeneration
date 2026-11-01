"""Exact cohomology and graph spectrum of the tensor-product diamond."""
from pathlib import Path
from itertools import product
import json
import sympy as s
lines=[(-1,1,-1,1),(1,-1,-1,1),(-1,1,1,-1),(1,-1,1,-1)]
h=[0]*5
for f,g in product(lines, repeat=2):
 d=[x-y for x,y in zip(f,g)]
 if -1 in d: continue
 q=sum(x<=-2 for x in d)
 n=1
 for x in d: n*=x+1 if x>=0 else -x-1
 h[q]+=n
assert h==[4,24,36,0,0],h
a,b,t=s.symbols('a b t')
D=s.Matrix([[1,1,0,0],[-1,0,1,0],[0,-1,0,1],[0,0,-1,-1]])
L=D*s.diag(a,b,b,a)*D.T
assert s.factor(L.charpoly(t).as_expr())==t*(t-2*a)*(t-2*b)*(t-2*a-2*b)
out={'graded_end_cohomology':h,'characteristic_polynomial':str(s.factor(L.charpoly(t).as_expr())),'status':'passed'}
(Path(__file__).resolve().parents[1]/'output/compute05_product_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)

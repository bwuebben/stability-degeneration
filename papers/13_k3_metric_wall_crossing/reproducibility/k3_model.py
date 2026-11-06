# /// script
# requires-python = ">=3.11"
# dependencies = ["sympy==1.14.0"]
# ///
"""Exact polynomial and intersection calculations for the K3 surface.

The affine Jacobian ideals verify smooth good reduction. The geometric
characteristic-zero argument appears in Section 2 of the manuscript.
"""
import json
import time
from pathlib import Path
import sympy as s

out = Path(__file__).resolve().with_name('k3-model.json')
x = s.symbols('x0:5')
x0, x1, x2, x3, x4 = x
q = x0*x3-x1*x2+x4**2
c0 = (x0**3+x0*x1**2+2*x0*x2**2+5*x1*x2*x3
      +4*x0*x1*x2+3*x1**2*x3+2*x2**2*x3)
c = c0+x4*(x0**2+x3**2+x4**2)
a,b,u,v = s.symbols('a b u v')
f = a*(u**2+v**2)+b*u*v
g = (a**2+2*b**2)*u+3*a*b*v
seg = {x0:a*u,x1:a*v,x2:b*u,x3:b*v,x4:0}
assert s.expand(c.subs(seg, simultaneous=True)-f*g) == 0
assert s.expand(q.subs(seg, simultaneous=True)) == 0
t = s.symbols('t')
p = -t**4+2*t**2+2
assert s.gcd(p,s.diff(p,t)) == 1
# The fifth intersection, t=infinity, is simple: the local leading term
# of g on f=0 is -u, after setting v=b=1.
local_g = s.factor(g.subs({b:1,v:1,a:-u/(1+u**2)}))
assert s.limit(local_g/u,u,0) == -1
# Normal derivative at the finite intersections: 4(1+t^2), modulo p.
normal = 1+(1+t**2)**2
assert s.rem(normal,p,t) == 4*(1+t**2)
assert s.gcd(normal,p) == 1

minors = [s.diff(q,x[i])*s.diff(c,x[j])-s.diff(q,x[j])*s.diff(c,x[i])
          for i in range(5) for j in range(i+1,5)]
prime = 101
charts=[]
for k in range(5):
    start=time.monotonic()
    variables=[z for i,z in enumerate(x) if i!=k]
    polynomials=[s.expand(z.subs(x[k],1)) for z in [q,c]+minors]
    basis=s.groebner(polynomials,*variables,modulus=prime,order='grevlex')
    unit = list(basis)==[s.Integer(1)]
    row={'chart':str(x[k])+'=1','unit_jacobian_ideal':unit,
         'seconds':round(time.monotonic()-start,3)}
    charts.append(row)
    print(json.dumps(row),flush=True)
    assert unit, 'Candidate surface has not passed the smoothness test'

gram=s.Matrix([[-2,5],[5,-2]])
H=s.Matrix([1,1]); D=s.Matrix([1,-1]); C=s.Matrix([2,1])
pair=lambda z,w:int((z.T*gram*w)[0])
invariants={
    'gram':[[int(gram[i,j]) for j in range(2)] for i in range(2)],
    'H_square':pair(H,H),'D_square':pair(D,D),'H_D':pair(H,D),
    'c1_square':pair(C,C),'c1_D1':pair(C,s.Matrix([1,0])),
    'rank':3,'c2':3+pair(C,C)//2,'mukai_s':0,
    'mukai_square':pair(C,C),'expected_moduli_dimension':pair(C,C)+2,
    'cross_ext1':-(2+pair(D,D)//2),
}
assert invariants['H_square']==6 and invariants['cross_ext1']==5
assert invariants['c2']==8 and invariants['expected_moduli_dimension']==12
result={'candidate_quadric':str(q),'candidate_cubic':str(c),
        'sympy_version':s.__version__,'smooth_reduction_prime':prime,
        'charts':charts,'segre_product_identity':True,
        'transverse_five_intersections':True,'invariants':invariants,
        'scope':'Finite algebra and smooth good reduction only; no moduli or metric theorem.'}
out.write_text(json.dumps(result,indent=2)+'\n')
print('Saved',out,flush=True)

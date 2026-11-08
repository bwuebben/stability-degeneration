"""Exact Pauli Gram and O'Neill coefficient with unnormalized real vectors.

Gaussian-integer matrix operations are exact in the small range used here.
Fractions handle every division. This verifies a finite calculation, not the
analytic metric comparison or physical reduction.
"""
from fractions import Fraction
import json

S = [[[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]]

def mul(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]

def comm(a,b):
    ab,ba=mul(a,b),mul(b,a)
    return [[ab[i][j]-ba[i][j] for j in range(2)] for i in range(2)]

def hs(a,b):
    value=sum(complex(a[i][j]).conjugate()*b[i][j] for i in range(2) for j in range(2))
    assert value.imag == 0 and value.real == int(value.real)
    return int(value.real)

def scale(c,a):
    return [[c*x for x in row] for row in a]

def gram(x):
    d=[comm(scale(1j,s),x) for s in S]
    return [[hs(u,v) for v in d] for u in d]

qa,qb=gram(S[2]),gram(S[0])
assert qa == [[8,0,0],[0,8,0],[0,0,0]]
assert qb == [[0,0,0],[0,8,0],[0,0,8]]
# u=sigma1, v=sigma2 in the unused loop: area squared=4, beta3=-4.
beta=[hs(comm(scale(1j,s),S[0]),S[1]) for s in S]
area=hs(S[0],S[0])*hs(S[1],S[1])-hs(S[0],S[1])**2
assert beta == [0,0,-4] and area == 4
coefficient=3*Fraction(beta[2]**2, qb[2][2]*area)
assert coefficient == Fraction(3,2)
# A wrong factor from treating the unnormalized basis as orthonormal fails.
assert 2*coefficient != Fraction(3,2)
index=4*30*24-4*496
for nh,nv in [(340,96),(338,94),(337,93)]:
    assert 4*(nh-20)-4*nv == index == 896
print(json.dumps({'Gram_a2':qa,'Gram_b2':qb,'unnormalized_beta':beta,
                  'plane_area_squared':area,'b2_times_sectional_curvature':str(coefficient),
                  'ASD_index':index,
                  'scope':'Exact finite Gram, curvature and index arithmetic'},indent=2))

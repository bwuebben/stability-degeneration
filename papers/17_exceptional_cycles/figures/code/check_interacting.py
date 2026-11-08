#!/usr/bin/env python3
"""Exact arithmetic supplement; not a substitute for any geometric/PDE proof."""
from fractions import Fraction as F
from functools import reduce
from itertools import product
from math import gcd
from pathlib import Path
import json


def dot(c, d):
    return sum(c[2*i]*d[2*i+1]+c[2*i+1]*d[2*i] for i in range(3))


def pairing(a, b):
    return dot(a[1], b[1])-a[0]*b[2]-b[0]*a[2]


def lin(*terms):
    return (sum(q*a[0] for q,a in terms),
            tuple(sum(q*a[1][i] for q,a in terms) for i in range(6)),
            sum(q*a[2] for q,a in terms))


H = (1,1,0,0,0,0)
D1 = (0,0,1,-3,0,0)
D2 = (0,0,0,-4,1,-4)
zero = (0,)*6
a0, a1, a2 = (1,zero,1), (1,D1,-2), (1,D2,-3)
avec = [a0,a1,a2]
v = lin((1,a0),(1,a1),(1,a2))
w = lin((1,a0),(1,a2))
c0 = lin((F(1,2),v),(-1,a0))
c2 = lin((F(1,2),v),(-1,a2))
f = a1
assert lin((1,c0),(1,c2)) == f
assert pairing(v,v)==2 and pairing(w,w)==0 and pairing(a1,w)==2
ns = [[dot(c,d) for d in [H,D1,D2]] for c in [H,D1,D2]]
mukai = [[pairing(c,d) for d in avec] for c in avec]
assert ns == [[2,0,0],[0,-6,-4],[0,-4,-8]]
assert mukai == [[-2,1,2],[1,-2,1],[2,1,-2]]
x=(0,(0,0,0,1,0,-1),0)
y=(0,(0,0,0,-3,0,0),-1)
assert pairing(v,x)==pairing(v,y)==0
witness_pairings=[[pairing(a,z) for a in avec] for z in [x,y]]
assert witness_pairings==[[0,1,-1],[1,-2,1]]
index_y=(0,(0,0,0,-1,0,0),0)
index_euler=[-pairing(a,index_y) for a in avec]
assert index_euler==[0,1,0] and -pairing(v,index_y)==1
restrictions=[[-pairing(a,z) for a in [a0,a2]] for z in [x,y]]
assert restrictions==[[0,1],[-1,-1]]
det=lambda B:B[0][0]*B[1][1]-B[0][1]*B[1][0]
assert abs(det(restrictions))==1
assert abs(det([[a*a for a in row] for row in restrictions]))==1
B=[[3,1],[1,3]]
Bp=[[3,-2],[-2,4]]
for matrix, functional in [(B,[1,-3]),(Bp,[6,1])]:
    assert det(matrix)==8
    assert reduce(gcd,[abs(k) for row in matrix for k in row])==1
    assert all(sum(functional[i]*matrix[i][j] for i in range(2))%8==0
               for j in range(2))
assert F(1,2)*pairing(c0,c2)**2-F(1,8)==1


def feasible(signs, levels):
    """Strict positivity of all three aggregated edge flows.

    q=x02>0, x01=(-zeta0-s02*q)/s01,
    x12=(zeta2-s02*q)/s12. Exact open interval feasibility.
    """
    s01,s12,s02=signs
    lower=F(0)
    upper=None
    for intercept,slope in [(-F(levels[0],s01),-F(s02,s01)),
                            (F(levels[2],s12),-F(s02,s12))]:
        endpoint=-intercept/slope
        if slope>0:
            lower=max(lower,endpoint)
        else:
            upper=endpoint if upper is None else min(upper,endpoint)
    return upper is None or lower<upper


acyclic=[s for s in product([-1,1],repeat=3)
         if s not in [(1,1,-1),(-1,-1,1)]]
chambers=[(1,-2,1),(-1,2,-1),(1,1,-2),(-1,-1,2),(-2,1,1),(2,-1,-1)]
orientations=[]
for levels in chambers:
    supports=[s for s in acyclic if feasible(s,levels)]
    assert len(supports)==2
    orientations.append({'levels':levels,'maximal_feasible_orientations':supports})
# The rank-one plus rank-two analytic bound forces k in {0,1}; for these
# k, the displayed quadratic formulas force j in {0,1} or {0} respectively.
Fpoly=lambda k,j:20*k+24*j-18*k*k-24*k*j-24*j*j
rank_pairs=[(k,j,Fpoly(k,j)) for k in [0,1] for j in range(-2,3)
            if Fpoly(k,j)>=0]
assert rank_pairs==[(0,0,0),(0,1,0),(1,0,2)]

result={
    'scope':'Exact finite arithmetic and oriented-support feasibility only; no geometric or asymptotic theorem is inferred from this computation.',
    'ns_gram':ns,'mukai_gram':mukai,
    'v_square':pairing(v,v),'w_square':pairing(w,w),
    'cycle_relation':'f=c0+c2',
    'curve_pairings':{k:str(val) for k,val in {
        'q(c0)':pairing(c0,c0),'q(c2)':pairing(c2,c2),
        'q(c0,c2)':pairing(c0,c2),'q(f)':pairing(f,f)}.items()},
    'witness_pairings':witness_pairings,'plane_restrictions':restrictions,
    'index_normalization_euler_pairings':index_euler,
    'intersection_matrices':[B,Bp],'smith_diagonals':[[1,8],[1,8]],
    'boundary_functionals_mod8':[[1,-3],[6,1]],
    'rank_one_candidates_after_analytic_bound':rank_pairs,
    'six_chamber_supports':orientations,'verdict':'PASS'
}
out=Path(__file__).with_name('interacting_data.json')
out.write_text(json.dumps(result,indent=2)+'\n')
print('PASS: lattice/Mukai pairings, primitive restriction witnesses, cycle relation, both boundary cokernels and all six chamber supports.')
print(out)

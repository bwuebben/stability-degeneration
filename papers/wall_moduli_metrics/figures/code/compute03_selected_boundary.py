#!/usr/bin/env python3
"""Exact finite checks for the singular endpoint; no PDE verification."""
from pathlib import Path
import json
import sympy as sp

s, t, L = sp.symbols('s t L')
x, y, u0, u1, v0, v1 = sp.symbols('x y u0 u1 v0 v1')
alpha, gamma1, gamma2, eta = sp.symbols('alpha gamma1 gamma2 eta')
P = sp.symbols('P0:6')
e0, e1 = sp.symbols('e0 e1')
a1 = alpha*L**3 + gamma1*L**2*t + gamma2*L*t**2 + x*t**3
b1 = eta*L+y*t
c0 = P[0]*L**3+P[1]*L**2*t+P[2]*L*t**2+u0*t**3
c1 = P[3]*L**3+P[4]*L**2*t+P[5]*L*t**2+u1*t**3
d0, d1 = e0*L+v0*t, e1*L+v1*t
quadratic = x*y-u0*v0-u1*v1
numerator = sp.expand(a1*b1-c0*d0-c1*d1-quadratic*t**4)
assert sp.Poly(numerator,L).coeff_monomial(1) == 0
a0 = sp.cancel(-numerator/L)
assert sp.expand(a0*L+a1*b1-c0*d0-c1*d1-quadratic*t**4) == 0
hessian = sp.hessian(quadratic,(x,y,u0,u1,v0,v1))
assert hessian.det() == -1

degrees = [3,3,3,3,1,1,1,1]
symbols, blocks = [], []
for i, deg in enumerate(degrees):
    coeff = sp.symbols(f'z{i}_0:{deg+1}')
    symbols.extend(coeff)
    blocks.append(sum(coeff[j]*s**(deg-j)*t**j for j in range(deg+1)))
A0,A1,C0,C1,B0,B1,D0,D1 = blocks
omega = sp.expand(A0*B0+A1*B1-C0*D0-C1*D1)
equations = sp.Matrix([omega.coeff(s,4-j).coeff(t,j) for j in range(5)])
jacobian = equations.jacobian(symbols)
inc = sp.Matrix([[1,1,0,0],[-1,0,1,0],[0,-1,0,1],[0,0,-1,-1]])
examples = {
    'S4': [0,s**3,0,0,s,0,0,0],
    'S5x': [0,s**3+t**3,0,0,s,0,0,0],
    'S5y': [-s**2*t,s**3,0,0,s,t,0,0],
    'S11': [s**3,s**3,s**3,0,s,0,s,0],
    'S16': [s**3,s**3+t**3,s**3,0,s,0,s,0],
}
expected = {'S4':(4,2),'S5x':(5,2),'S5y':(5,2),'S11':(4,3),'S16':(5,3)}
rank_results = {}
for name, polynomials in examples.items():
    values=[]
    for polynomial, degree in zip(polynomials,degrees):
        polynomial=sp.sympify(polynomial).expand()
        values.extend(polynomial.coeff(s,degree-j).coeff(t,j) for j in range(degree+1))
    mapping=dict(zip(symbols,values))
    assert equations.subs(mapping)==sp.zeros(5,1)
    jr=jacobian.subs(mapping).rank()
    active=[any(polynomials[i] != 0 for i in pair) for pair in [(0,1),(2,3),(4,5),(6,7)]]
    orbit_rank=inc[:,[i for i,a in enumerate(active) if a]].rank()
    assert (jr,orbit_rank)==expected[name]
    rank_results[name]={'obstruction_rank':jr,'orbit_rank':orbit_rank,'effective_stabilizer_dimension':3-orbit_rank}

l, delta = sp.symbols('l delta',positive=True)
weights=sp.Matrix([l-delta,delta,sp.Rational(3,5)*l-delta,delta])
assert inc*weights==sp.Matrix([l,-sp.Rational(2,5)*l,0,-sp.Rational(3,5)*l])
assert sum(weights)==sp.Rational(8,5)*l

# M=r^-2; all scalar functions below are independent of l.
r=sp.symbols('r',real=True)
c=sp.sqrt(sp.Rational(3,5))
h=(-sp.Rational(8,5)*r**2+sp.sqrt(sp.Rational(12,5)+sp.Rational(124,25)*r**4))/(2*(1+r**4))
assert sp.simplify((1+r**4)*h**2+sp.Rational(8,5)*r**2*h-sp.Rational(3,5))==0
hs=sp.series(h,r,0,5).removeO()
assert hs.subs(r,0)==c and hs.coeff(r,2)==-sp.Rational(4,5)
f=sp.sqrt((1-r**2*h)/(1+2*r**4))
# Series suffice for the proved endpoint speed coefficient, not just numerics.
fs=sp.series(f,r,0,5).removeO()
cs=sp.series(r*sp.sqrt(h),r,0,5).removeO()
bs=sp.series(sp.sqrt(sp.Rational(3,5)-r**2*h),r,0,5).removeO()
# Orthogonal monomial coefficients: a0 p, a1 p, a1 q, C p, B ell, D ell.
curve=sp.Matrix([fs*r**2,fs,fs*r**2,cs,bs,cs])
velocity=curve.diff(r)
speed=sp.expand(velocity.dot(velocity))
assert sp.simplify(speed.subs(r,0)-2*c)==0
assert sp.simplify(speed.coeff(r,1))==0

# Exact moment-level horizontality: each block has real coefficients.
# For arbitrary r the complex orbit pairing is one half D_inc w'(r), zero.
tau=r**2*h
w_r=sp.Matrix([1-tau,tau,sp.Rational(3,5)-tau,tau])
assert sp.simplify(inc*w_r.diff(r))==sp.zeros(4,1)
# Product equality of the two nonzero contractions, squared (positive roots).
assert sp.simplify((1-tau)*(sp.Rational(3,5)-tau)/(r**-4+2)-tau**2)==0

result={
 'kind':'exact finite algebra and series; not verification of PDE, Hodge or distance proofs',
 'slice_division':'PASS', 'six_variable_hessian_determinant':-1,
 'piece_rank_examples':rank_results,
 'moment_demand':['l','-2*l/5','0','-3*l/5'],
 'balanced_energy':'8*l/5',
 'h_at_zero':str(c), 'h_r2_coefficient':'-4/5',
 'squared_speed_divided_by_l_at_endpoint':str(2*c),
 'finite_distance_coefficient':'sqrt(2*l*sqrt(3/5))',
 'distance_rate_in_M':'M^(-1/2)',
 'moment_horizontality':'PASS', 'balanced_path_integrability':'PASS',
}
output=(Path(__file__).resolve().parents[1]/'output'/'selected_boundary_results.json')
output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

"""Rigorous interval enclosure of the real-boundary mass derivative / pi.

The exponential inputs are enclosed by exact rational Taylor sums.
All subsequent operations are outward-rounded mpmath interval arithmetic.
Only real 2x2 matrices and square roots are needed; no eigenvector solver,
matrix logarithm, numerical differentiation or quadrature enters this bound.
"""
from fractions import Fraction
from math import factorial
from pathlib import Path
import json
from mpmath import iv

iv.dps = 60


def taylor(x, n):
    return sum(((-x)**k / factorial(k) for k in range(n+1)), Fraction(0))


def exp_bounds(x, low_text, high_text, n):
    low, high = Fraction(low_text), Fraction(high_text)
    # Taylor's integral remainder has the next term's sign for exp(-x), x>=0.
    assert n % 2 == 0
    assert low < taylor(x, n+1) < taylor(x, n) < high
    return iv.mpf([low_text, high_text])


def tr(a):
    return a[0,0]+a[1,1]


def det(a):
    return a[0,0]*a[1,1]-a[0,1]*a[1,0]


def negative_projection(a):
    disc = iv.sqrt(tr(a)**2-4*det(a))
    return (iv.eye(2)-(2*a-tr(a)*iv.eye(2))/disc)/2


def display(a):
    if hasattr(a, 'rows'):
        return [[str(a[i,j]) for j in range(a.cols)] for i in range(a.rows)]
    return str(a)


b1,b2 = iv.mpf(1)/25, iv.mpf(4)
p = exp_bounds(Fraction(1,25),'0.96078943915232320','0.96078943915232322',30)
q = exp_bounds(Fraction(4),'0.01831563888873417','0.01831563888873419',80)
gap = b2-b1
rr1,rr2 = (1-p)/b1,(1-q)/b2
f1,f2 = rr1-p,rr2-q
ell = (p-q)/gap
s = iv.sqrt(b1*b2)*(rr1-rr2)/gap
c = iv.matrix([[1,2],[2,4]])
h = iv.matrix([[2*p,2*ell],[2*ell,5*q]])
g = iv.matrix([[f1-p,2*s],[2*s,4*f2-q]])
w = 4-rr1-4*rr2
assert det(h).b < 0 and det(g).b < 0 and w.a > 0
pv,pu = negative_projection(h),negative_projection(g)
b = iv.matrix([[b1,0],[0,b2]])
d = iv.matrix([[iv.mpf(1)/5,0],[0,2]])
db = (pv*b+b*pv)/2+d*pu*d
gamma = db[0,1]/gap
rotation = iv.matrix([[0,gamma],[-gamma,0]])
dc = (pv*c+c*pv)/2+c*rotation-rotation*c
db1,db2 = db[0,0],db[1,1]
dell = ((-p*db1+q*db2)*gap-(p-q)*(db2-db1))/(gap**2)
dr1,dr2 = -f1/b1*db1,-f2/b2*db2
ds = s*(db1/b1+db2/b2)/2+iv.sqrt(b1*b2)*(
       (dr1-dr2)*gap-(rr1-rr2)*(db2-db1))/(gap**2)
dh = iv.matrix([[p*(dc[0,0]-2*db1),dc[0,1]*ell+2*dell],
                [dc[1,0]*ell+2*dell,q*(dc[1,1]-5*db2)]])
dg = iv.matrix([[f1*dc[0,0]+((-f1/b1+p)+p)*db1,dc[0,1]*s+2*ds],
                [dc[1,0]*s+2*ds,f2*dc[1,1]+(4*(-f2/b2+q)+q)*db2]])
derivative_over_pi = -2*tr(pv*dh+pu*dg)
assert derivative_over_pi.a > iv.mpf('0.2395').b
assert derivative_over_pi.b < iv.mpf('0.2397').a
result = {key:display(val) for key,val in dict(p=p,q=q,H=h,G=g,w=w,detH=det(h),
    detG=det(g),Pv=pv,Pu=pu,dB=db,dC_rotated=dc,dH=dh,dG=dg,
    derivative_over_pi=derivative_over_pi).items()}
result['precision_decimal_digits'] = iv.dps
result['bound'] = '0.2395 < lim_(epsilon->0+) mass_derivative / pi < 0.2397'
output = Path(__file__).parent/'results'/'limit-sign-interval.json'
output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({key:result[key] for key in ('detH','detG','w','derivative_over_pi','bound')},indent=2))

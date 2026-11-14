"""Exact normalization checks printed in Paper 33; Python 3 + SymPy."""
import sympy as s
x, scale, radius = s.symbols('x scale radius', positive=True)
pi = s.pi
energy = s.integrate(48*scale**4/(radius**2+scale**2)**4 * 2*pi**2*radius**3, (radius,0,s.oo))
assert energy == 8*pi**2
second = s.integrate(48*scale**4/(radius**2+scale**2)**4 * 2*pi**2*radius**5, (radius,0,s.oo))
assert second == 16*pi**2*scale**2
print('BPST curvature mass:', energy)
print('BPST second moment:', second)
for n in (2,3,4,7):
 q=32*pi**4*x
 fujiki=s.factorial(2*n)*q**n/(2**n*s.factorial(n))
 volume=s.factor(fujiki/s.factorial(2*n))
 assert s.simplify(volume-(16*pi**4*x)**n/s.factorial(n))==0
 c2=6*(n+3)*s.factorial2(2*n-3)*q**(n-1)
 curvature=8*pi**2*c2/s.factorial(2*n-2)
 assert s.simplify(curvature-48*pi**2*(n+3)*(16*pi**4*x)**(n-1)/s.factorial(n-1))==0
 diameter4=(volume/(pi**(2*n)/s.factorial(2*n)))**s.Rational(1,n)
 expected=16*pi**2*x*(s.factorial(2*n)/s.factorial(n))**s.Rational(1,n)
 assert s.simplify(diameter4-expected)==0
 print('N=',n,'volume=',volume,'Chern-Weil and diameter PASS')
for r in (2,3,4,5):
 k=s.symbols('k')
 dim=s.expand(s.expand_func((2*k+r-1)*s.binomial(k+r-2,r-2)**2/(r-1)))
 poly=s.Poly(dim,k)
 assert poly.degree()==2*r-3 and poly.eval(1)==r*r-1
 leading=poly.LC()*s.factorial(2*r-3)/2**(2*r-2)
 expected=s.binomial(2*r-2,r-1)/4**(r-1)
 assert s.simplify(leading-expected)==0
 print('rank',r,'one-instanton density',expected)
for ss,n,c2,kappa in ((0,2,3,s.Rational(5,2)),(-1,4,4,s.Rational(7,2))):
 assert (2-4*ss+2)//2==n and 2+1-ss==c2 and c2-s.Rational(1,4)*2==kappa
 print('double-plane v=(2,H,%s), N=%s, c2=%s, kappa=%s PASS'%(ss,n,c2,kappa))
print('ALL CHECKS PASS')

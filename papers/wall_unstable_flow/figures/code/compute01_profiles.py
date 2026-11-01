"""Exact scalar identities used in the energy and residual estimates."""
import sympy as s

a,d,q,t,tau=s.symbols("a d q t tau",positive=True)
L=q*s.exp(-a*tau*t)/(1+d*q/(a*tau)*(1-s.exp(-a*tau*t)))
assert s.simplify(s.diff(L,t)+a*tau*L+d*L**2)==0
assert s.simplify(L.subs(t,0)-q)==0
assert s.simplify(s.limit(L,tau,0,dir="+")-q/(1+d*q*t))==0
I=s.log(1+d*q/(a*tau)*(1-s.exp(-a*tau*t)))/d
assert s.simplify(s.diff(I,t)-L)==0
assert I.subs(t,0)==0
assert s.simplify(s.limit(I,t,s.oo)-s.log(1+d*q/(a*tau))/d)==0
upper=s.simplify(I.subs({a:s.Rational(2,5),d:s.Rational(1,5)}))
expected=5*s.log(1+q/(2*tau)*(1-s.exp(-2*tau*t/5)))
assert s.simplify(upper-expected)==0
assert s.simplify(s.limit(upper,tau,0,dir="+")-5*s.log(1+q*t/5))==0
print("PASS: logistic ODE, initial value, wall limit, integral and residual bound.")

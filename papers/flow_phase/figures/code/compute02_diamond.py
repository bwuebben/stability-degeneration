"""Symbolic checks of the printed diamond weights and finite spectrum."""
import sympy as sp

k,l,s = sp.symbols("k l s", positive=True)
D = sp.Matrix([[1,1,0,0],[-1,0,1,0],[0,-1,0,1],[0,0,-1,-1]])
w = sp.Matrix([5*l,5*k,3*k,3*l])/8
b = sp.Matrix([5*(k+l),3*k-5*l,-5*k+3*l,-3*(k+l)])/8
assert D*w == b
assert sum(w) == k+l
assert sp.simplify(w[0]*w[2]-w[1]*w[3]) == 0
height = sp.Matrix([1,0,0,-1])
assert (height.T*D) == sp.ones(1,4)
L = D*sp.diag(*w)*D.T
expected = s**3-2*(k+l)*s**2+(sp.Rational(15,16)*(k*k+l*l)+3*k*l)*s-sp.Rational(15,16)*k*l*(k+l)
assert sp.simplify((s*sp.eye(4)-L).det()-s*expected) == 0
assert sp.factor(expected.subs(k,0)) == s*(-5*l+4*s)*(-3*l+4*s)/16
assert sp.simplify(sp.Rational(15,16)*l*l/(sp.Rational(3,4)*l*sp.Rational(5,4)*l)) == 1
assert sp.simplify(sum(bi**2 for bi in b)-(17*k*k+2*k*l+17*l*l)/16) == 0
print("PASS: exact balance, cycle invariant, height, characteristic polynomial, endpoint roots and weak coefficient.")

"""Exact weighted five-mode identities; no numerical or PDE approximation."""

import json
from pathlib import Path

import sympy as sp

c, s, mu = sp.symbols("c s mu", real=True)
A0, A1, A2 = sp.symbols("A0 A1 A2", positive=True)
Q = sp.Matrix([
    [2, -2, 0, 0, 0],
    [-2, 4, 0, 0, -2*c],
    [0, 0, 4, 0, -2*s],
    [0, 0, 0, 4, 0],
    [0, -2*c, -2*s, 0, 2],
])
M = sp.diag(A0, 2*A1, 2*A1, 2*A1, A2)
circle = lambda x: sp.factor(sp.expand(x).subs(s**2, 1-c**2))
characteristic = circle((Q-mu*M).det())/M.det()
determinant = sp.factor(characteristic.subs(mu, 0))
assert sp.simplify(determinant-8*(1-c*c)/(A0*A1**3*A2)) == 0
nonzero_product = sp.factor(-sp.diff(characteristic, mu).subs({mu: 0, c: 1}))
expected_product = 8*(A0+2*A1+A2)/(A0*A1**3*A2)
assert sp.simplify(nonzero_product-expected_product) == 0
small_coefficient = sp.factor((8/(A0*A1**3*A2))/nonzero_product)
assert small_coefficient == 1/(A0+2*A1+A2)
kernel = Q.subs({c: 1, s: 0}).nullspace()
assert len(kernel) == 1
assert kernel[0] == sp.Matrix([1, 1, 0, 0, 1])
assert Q.subs({c: 1, s: 0}).adjugate() == 64*kernel[0]*kernel[0].T
assert (kernel[0].T*M*kernel[0])[0] == A0+2*A1+A2
equal_mass = sp.factor(characteristic.subs({A0: 1, A1: 1, A2: 1}))
expected_equal = (2-mu)*(((2-mu)**2-2)**2-4*c*c)
assert sp.simplify(equal_mass-expected_equal) == 0
out = {
    "status": "PASS: exact weighted finite identities; actual geometric/PDE comparison is a separate proof",
    "mass": [str(x) for x in M.diagonal()],
    "characteristic": str(sp.factor(characteristic)),
    "spectral_product": str(determinant),
    "four_positive_eigenvalue_product_at_c_1": str(nonzero_product),
    "small_eigenvalue_theta_squared_coefficient": str(small_coefficient),
    "sin_squared_estimator_coefficient": str(A0*A1**3*A2/8),
    "zero_mode": [str(x) for x in kernel[0]],
}
(Path(__file__).resolve().parents[1] / "output" / "weighted-spectrum.json").write_text(json.dumps(out, indent=2)+"\n")
print(out["status"])
print("Product:", determinant)
print("Weak coefficient:", small_coefficient)

"""Exact matrix calculation for Proposition 3.5 (label prop:three-scales).

Run with uv run --with sympy python figures/code/compute_three_scales.py.
The printed matrix and analytic proof, not this check, establish the result.
"""
from pathlib import Path
import json
import sympy as sp

h, z = sp.symbols("h z", positive=True)
lam = sp.Symbol("lam")
eye = sp.eye(3)
rx = sp.Matrix([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
ry = sp.Matrix([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
edges = [(0, 1, 1, eye), (0, 2, 1, eye), (0, 3, 1, eye),
         (1, 2, 2, rx), (1, 3, 4, ry), (2, 3, 3, eye)]
q = sp.zeros(12)
for i, j, power, transport in edges:
    w = h**power
    for a in range(3):
        q[3*i+a, 3*i+a] += w
        q[3*j+a, 3*j+a] += w
        for b in range(3):
            q[3*i+a, 3*j+b] -= w*transport[a, b]
            q[3*j+b, 3*i+a] -= w*transport[a, b]
assert q == q.T
assert rx.T*rx == eye and ry.T*ry == eye

characteristic = q.charpoly(lam)
assert characteristic.gen == lam
poly = characteristic.as_expr()
determinant = sp.factor(q.det())
fast = sp.factor(sp.limit(poly.subs(lam, h*z)/h**12, h, 0, dir="+"))
middle = sp.factor(sp.limit(poly.subs(lam, h**2*z)/h**15, h, 0, dir="+"))
det_lead = sp.limit(determinant/h**17, h, 0, dir="+")
slow = sp.factor(sp.limit(poly.subs(lam, h**4*z)/h**17, h, 0, dir="+"))

assert sp.expand(fast - z**3*(z-1)**6*(z-4)**3) == 0
assert sp.expand(middle + 64*z*(z-sp.Rational(1, 2))**2) == 0
assert det_lead == 8
assert sp.expand(slow - (8-16*z)) == 0

q0 = sp.zeros(12)
for i, j, _, transport in edges[:3]:
    for a in range(3):
        q0[3*i+a, 3*i+a] += 1
        q0[3*j+a, 3*j+a] += 1
        q0[3*i+a, 3*j+a] -= 1
        q0[3*j+a, 3*i+a] -= 1
assert q0.eigenvals() == {0: 3, 1: 6, 4: 3}
diag_embedding = sp.Matrix.vstack(eye, eye, eye, eye)/2
assert diag_embedding.T*diag_embedding == eye
q2 = q.applyfunc(lambda entry: sp.expand(entry).coeff(h, 2))
compressed_middle = sp.simplify(diag_embedding.T*q2*diag_embedding)
assert compressed_middle == sp.diag(0, sp.Rational(1, 2), sp.Rational(1, 2))
q4 = q.applyfunc(lambda entry: sp.expand(entry).coeff(h, 4))
root_x = diag_embedding[:, 0]
assert (root_x.T*q4*root_x)[0] == sp.Rational(1, 2)

result = {
    "status": "PASS; exact finite algebra, not a PDE proof",
    "mass": "I_12 (four unit-area components, orthonormal Pauli basis)",
    "edge_weights": {f"{i}{j}": str(h**power) for i, j, power, _ in edges},
    "transport_12": sp.sstr(rx), "transport_13": sp.sstr(ry),
    "determinant": str(determinant),
    "rescaled_characteristic_limits": {
        "P(h*z,h)/h^12": str(fast),
        "P(h^2*z,h)/h^15": str(middle),
        "P(h^4*z,h)/h^17": str(slow)},
    "determinant_leading_coefficient_at_h17": str(det_lead),
    "middle_kernel_compression": sp.sstr(compressed_middle),
    "weakest_kernel_compression": str((root_x.T*q4*root_x)[0]),
    "predicted_positive_eigenvalue_orders": {
        "h": {"coefficient_1": 6, "coefficient_4": 3},
        "h^2": {"coefficient_1/2": 2},
        "h^4": {"coefficient_1/2": 1}},
    "geometric_path": "tau_e=exp(-2*pi/h^p_e); real path in the actual complex polydisc, not claimed holomorphic at h=0",
    "scope": "Finite algebra only. The actual spectral comparison is Theorem thm:spectrum; no nonlinear weakest-time inference."
}
out = Path(__file__).resolve().parents[1] / "output" / "three-scales-exact.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps({"status": result["status"], "determinant": str(determinant),
                  "fast": str(fast), "middle": str(middle), "slow": str(slow)}, indent=2))

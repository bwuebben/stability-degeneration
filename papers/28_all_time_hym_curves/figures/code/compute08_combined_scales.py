"""Exact finite checks for the proposed R30 combined-scale example.

This verifies finite matrices and polynomial limits, not the relative PDE
comparison or any actual HYM trajectory theorem.
"""

from pathlib import Path
import json
import sympy as sp


h, z = sp.symbols("h z", positive=True)
eye = sp.eye(3)
rotation = sp.Matrix([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
assert rotation.T * rotation == eye
rows = []
for i, j, transport in [(0, 1, eye), (1, 2, eye), (2, 3, eye), (3, 0, rotation)]:
    block = sp.zeros(3, 12)
    block[:, 3 * i : 3 * i + 3] = eye
    block[:, 3 * j : 3 * j + 3] = -transport
    rows.append(block)
path = sum((block.T * block for block in rows[:3]), sp.zeros(12))
cycle = rows[3].T * rows[3]
component = sp.zeros(12)
component[:3, :3] = sp.diag(4, 4, 0)
matrix = h * path + h**2 * cycle + h**4 * component

# The normalized constant-tuple inclusion has norm one since there are
# four vertices of area one.
constant = sp.Matrix.vstack(eye, eye, eye, eye) / 2
assert constant.T * constant == eye
cycle_compression = sp.simplify(constant.T * cycle * constant)
component_compression = sp.simplify(constant.T * component * constant)
assert cycle_compression == sp.diag(0, sp.Rational(1, 2), sp.Rational(1, 2))
assert component_compression == sp.diag(1, 1, 0)
assert path * constant == sp.zeros(12, 3)
assert cycle * constant[:, 0] == sp.zeros(12, 1)
assert (path + cycle).rank() == 11
assert (path + cycle + component).rank() == 12

x_indices = [0, 3, 6, 9]
yz_indices = [i for i in range(12) if i not in x_indices]
assert matrix.extract(x_indices, yz_indices) == sp.zeros(4, 8)
x_matrix = matrix.extract(x_indices, x_indices)
yz_matrix = matrix.extract(yz_indices, yz_indices)
char_x = x_matrix.charpoly().as_expr()
gen_x = x_matrix.charpoly().gen
char_x = sp.factor(char_x.subs(gen_x, z))
char_yz_obj = yz_matrix.charpoly()
char_yz = sp.factor(char_yz_obj.as_expr().subs(char_yz_obj.gen, z))
characteristic = sp.expand(char_x * char_yz)
determinant = sp.factor(characteristic.subs(z, 0))
assert sp.limit(determinant / h**17, h, 0, dir="+") == 16

# At the three scales these are the exact limiting characteristic
# polynomials after removing the common power of h.
high = sp.factor(sp.limit(characteristic.subs(z, h * z) / h**12, h, 0, dir="+"))
middle = sp.factor(sp.limit(characteristic.subs(z, h**2 * z) / h**15, h, 0, dir="+"))
weak = sp.factor(sp.limit(characteristic.subs(z, h**4 * z) / h**17, h, 0, dir="+"))
assert high == sp.expand(z**3 * (z - 2)**3 * (z**2 - 4*z + 2)**3).factor()
assert middle == -16 * z * (2*z - 1)**2
assert sp.expand(weak - 16 * (1 - z)) == 0

out = {
    "scope": "Exact finite identities and polynomial limits only; the R30 PDE proof is separate.",
    "graph": "Path 0-1-2-3 with identity transports and closing edge 3-0 with Ad((I+i sigma_x)/sqrt(2)).",
    "weights": ["h", "h", "h", "h^2"],
    "component_coefficient_alpha_squared_kappa": "h^4 at vertex 0, T=sigma_z",
    "mass": "Identity in Hilbert-Schmidt orthonormal Pauli coordinates; all four areas one.",
    "rank_without_component": 11,
    "rank_with_component": 12,
    "cycle_compression": [[str(v) for v in row] for row in cycle_compression.tolist()],
    "component_compression": [[str(v) for v in row] for row in component_compression.tolist()],
    "determinant": str(determinant),
    "characteristic_x_sector": str(char_x),
    "characteristic_yz_sector": str(char_yz),
    "high_scale_polynomial": str(high),
    "middle_scale_polynomial": str(middle),
    "weak_scale_polynomial": str(weak),
    "leading_modes": [
        {"order": "h", "coefficient": "2-sqrt(2)", "multiplicity": 3},
        {"order": "h", "coefficient": "2", "multiplicity": 3},
        {"order": "h", "coefficient": "2+sqrt(2)", "multiplicity": 3},
        {"order": "h^2", "coefficient": "1/2", "multiplicity": 2},
        {"order": "h^4", "coefficient": "1", "multiplicity": 1},
    ],
    "result": "PASS",
}
destination = Path(__file__).with_name("combined_scales.json")
destination.write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))

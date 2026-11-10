"""Exact finite checks; these do not prove nonlinear HYM tracking."""
import json
from pathlib import Path
import sympy as sp

I = sp.eye(2)
H2 = sp.diag(2, sp.Rational(1, 2))
H3 = sp.Matrix([[sp.Rational(5, 3), sp.Rational(4, 3)],
                [sp.Rational(4, 3), sp.Rational(5, 3)]])
commutator = H2 * H3 - H3 * H2
Q = 3 * sp.eye(3) - sp.ones(3)
assert H2.det() == H3.det() == 1
assert H2[0, 0] > 0 and H3[0, 0] > 0
assert commutator == sp.Matrix([[0, 2], [-2, 0]])
assert Q.eigenvals() == {0: 1, 3: 2}

result = {
    "status": "finite algebra passed; no analytic theorem certified",
    "H1": [[str(x) for x in I.row(i)] for i in range(2)],
    "H2": [[str(x) for x in H2.row(i)] for i in range(2)],
    "H3": [[str(x) for x in H3.row(i)] for i in range(2)],
    "commutator": [[str(x) for x in commutator.row(i)] for i in range(2)],
    "graph_eigenvalues": {str(k): v for k, v in Q.eigenvals().items()},
    "rank_two_positive_mode_count": (Q.rank()) * (I.rows**2 - 1),
    "rank_two_leading_coefficient": str(2 * sp.pi * min(k for k in Q.eigenvals() if k > 0)),
}
out = Path(__file__).with_suffix(".json")
out.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result))

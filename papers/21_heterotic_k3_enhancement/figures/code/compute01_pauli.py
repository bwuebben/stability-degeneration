#!/usr/bin/env python3
"""Finite matrix and index checks only; no HYM or string approximation is proved.

All matrix entries used below are Gaussian integers exactly representable by
Python complex numbers. The symbolic a²,b² coefficient matrices are assembled
separately, so no small-eigenvalue numerical subtraction is involved.
"""
import json

Z = [[0j, 0j], [0j, 0j]]
I = [[1, 0], [0, 1]]
S = [[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]]

def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]

def add(a, b):
    return [[a[i][j] + b[i][j] for j in range(2)] for i in range(2)]

def sub(a, b):
    return [[a[i][j] - b[i][j] for j in range(2)] for i in range(2)]

def adj(a):
    return [[complex(a[j][i]).conjugate() for j in range(2)] for i in range(2)]

def comm(a, b):
    return sub(mul(a, b), mul(b, a))

def hs(a, b):
    c = mul(adj(a), b)
    return c[0][0] + c[1][1]

def real_integer(x):
    assert x.imag == 0 and x.real == int(x.real)
    return int(x.real)

def orbit_coefficient(x):
    # The Pauli basis has HS mass 2 I_3. Divide by 2, before rank*Vol.
    return [[real_integer(hs(comm(s, x), comm(t, x)) / 2) for t in S] for s in S]

ka, kb = orbit_coefficient(S[2]), orbit_coefficient(S[0])
assert ka == [[4, 0, 0], [0, 4, 0], [0, 0, 0]]
assert kb == [[0, 0, 0], [0, 4, 0], [0, 0, 4]]
assert all(comm(s, adj(s)) == Z for s in [S[2], S[0]])
assert comm(S[2], Z) == Z and comm(S[0], Z) == Z
assert hs(comm(S[2], S[0]), comm(S[2], S[0])) == 8
# Failure controls: nonnormal X violates the real equation; noncommuting
# X,Y in the same symplectic pair violates the complex equation.
raising = [[0, 1], [0, 0]]
assert comm(raising, adj(raising)) != Z
assert comm(S[2], S[0]) != Z

# Mukai/Ext arithmetic conditional on stable, locally free K3 bundles.
def mukai(r, c2):
    s = r - c2
    square = -2 * r * s
    return {"v": [r, 0, s], "square": square, "ext1": square + 2, "chi": 2 * r - c2}

am, bm, em = mukai(3, 7), mukai(3, 10), mukai(6, 14)
assert am["ext1"] == 26 and bm["ext1"] == 44 and em["ext1"] == 98
assert 2 * 7 + 10 == 24 and 2 * (3 * 2 + 3) + 14 == 32

# Six-dimensional gravitational anomaly diagnostic for the orthogonal
# adjoint decomposition derived in the manuscript.
hypers = [
    ("adjoint U2", 13, 4), ("B moduli", 22, 1),
    ("Sym2 U2", 1, 3), ("wedge2 U2", 23, 1),
    ("B charged singlets", 4, 1), ("A B plus", 33, 2),
    ("A B minus", 33, 2), ("SO14 U2", 1, 28),
    ("SO14 B", 4, 14), ("geometric", 20, 1),
]
nh = sum(n * d for _, n, d in hypers)
nv = 14 * 13 // 2 + 4 + 1
assert nh == 340 and nv == 96 and nh - nv == 244

print(json.dumps({
    "scope": "Exact finite identities and conditional index arithmetic, not an asymptotic proof",
    "orbit_a2_coefficient_before_rank_volume": ka,
    "orbit_b2_coefficient_before_rank_volume": kb,
    "ordered_finite_eigenvalues_for_0_b_le_a": ["4*b^2/(3*Vol)", "4*a^2/(3*Vol)", "4*(a^2+b^2)/(3*Vol)"],
    "kernel_dimensions_tracefree": {"a=b=0": 3, "one_axis_nonzero": 1, "a*b_nonzero": 0},
    "A": am, "B": bm, "E": em,
    "candidate_nH": nh, "candidate_nV": nv,
    "anomaly_difference": nh - nv,
    "failure_controls": "nonnormal real-moment and noncommuting complex-moment inputs rejected"
}, indent=2))

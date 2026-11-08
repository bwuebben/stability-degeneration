#!/usr/bin/env python3
"""Reproduce exact numerical identities, not geometric existence or PDE proofs."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json


def pairing(x, y):
    r, h, d, s = x
    R, H, D, S = y
    return 2*h*H - 22*d*D - r*S - R*s


def add(x, y):
    return tuple(a+b for a, b in zip(x, y))


def mul(c, x):
    return tuple(c*a for a in x)


def c2(x):
    r, h, d, s = x
    return r + h*h - 11*d*d - s


def det3(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


a = (1, 0, 0, 1)
b = (2, 0, 1, -5)
v = add(a, b)
ell = add(a, mul(Q(-1, 2), v))
rho = mul(2, ell)
w0 = (9, 0, 3, -10)
basis = [(0, 1, 0, 0), (1, 0, 1, -6), (3, 0, 0, 4)]
gram = [[pairing(x, y) for y in basis] for x in basis]

assert (pairing(a, a), pairing(b, b), pairing(a, b), pairing(v, v)) == (-2, -2, 3, 2)
assert c2(b) == c2(v) == -4
# Rank and D coordinates give a unit minor for the factor sublattice.
# Unimodularity of the full Mukai lattice then supplies the index normalizer.
factor_minor=a[0]*b[2]-b[0]*a[2]
assert factor_minor==1
assert pairing(ell, ell) == Q(-5, 2)
assert pairing(rho, rho) == -10
assert pairing(v, w0) == 0 and pairing(ell, w0) == 1
assert pairing(rho, w0) == 2
assert all(pairing(v, x) == 0 for x in basis)
assert det3(gram) == 88
assert add(mul(3, basis[1]), mul(2, basis[2])) == w0

# Check the displayed quadratic polynomial by enough exact evaluations to
# determine its k-degree <=2 and affine n,l coefficients. The inequalities
# for all integers are proved in the notes, not by these finite evaluations.
for k, n, lost in [(0,0,0),(1,0,0),(2,0,0),(0,1,0),(0,0,1)]:
    residual = (2, 0, 1-k, -5+lost+11*k*k+n)
    assert pairing(residual, residual) == -2+44*k-66*k*k-4*(n+lost)

# Fujiki identities: integral rho^4, c2*rho^2, c2^2.
rho4, c2rho2, c2square = 300, -300, 828
plane_square = Q(rho4,64) + Q(c2rho2,96) + Q(c2square,576)
plane_c2 = Q(c2rho2,8) + Q(c2square,24)
assert plane_square == 3 and plane_c2 == -3

wall = (0, 1, 0, 0)
direction = (0, 0, 1, Q(-22,3))
assert pairing(v, wall) == pairing(v, direction) == 0
assert pairing(wall, ell) == 0 and pairing(direction, ell) == Q(22,3)
assert pairing(wall, wall) == 2 and pairing(direction, direction) == -22
area_over_pi2 = 2*pairing(direction, ell)
volume_over_pi4 = area_over_pi2**2/2
assert area_over_pi2 == Q(44,3) and volume_over_pi4 == Q(968,9)

outdir = Path(__file__).resolve().parent
notes = []
output = {
    'scope': 'Exact lattice, polynomial-coefficient, Fujiki and normalization identities only; not geometric existence, global contraction, compactness or metric proof.',
    'vectors': {'a':a, 'b':b, 'v':v, 'rho':rho, 'w0':w0},
    'orthogonal_gram': gram, 'orthogonal_determinant':det3(gram),
    'primitive_factor_unit_minor':factor_minor,
    'line_square':str(pairing(ell,ell)), 'plane_square':str(plane_square),
    'plane_c2':str(plane_c2), 'area_coefficient_over_pi2':str(area_over_pi2),
    'volume_coefficient_over_pi4':str(volume_over_pi4),
    'note_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in notes},
}
(outdir/'isolated_data.json').write_text(json.dumps(output,indent=2,default=str)+'\n')
print('Exact identities pass: lattice determinant 88; line square -5/2; plane square 3; c2 pairing -3; area 44*pi^2*|t|/3; plane volume 968*pi^4*t^2/9.')
print('Geometric existence, gluing, compactification and PDE assertions require the written proofs.')

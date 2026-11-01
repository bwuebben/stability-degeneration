#!/usr/bin/env python3
"""Exact finite checks for the balanced rank-five family; no geometric or PDE claims certified here."""
from fractions import Fraction as F
from pathlib import Path
import json


def rank(a):
    a = [list(map(F, row)) for row in a]
    i = 0
    for j in range(len(a[0])):
        piv = next((k for k in range(i, len(a)) if a[k][j]), None)
        if piv is None:
            continue
        a[i], a[piv] = a[piv], a[i]
        v = a[i][j]
        a[i] = [x / v for x in a[i]]
        for k in range(len(a)):
            if k != i:
                v = a[k][j]
                a[k] = [x - v*y for x, y in zip(a[k], a[i])]
        i += 1
        if i == len(a):
            break
    return i


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


# Independent coordinates z1,M11,M12,M21,M22,z3; z4 fixed by trace zero.
u = [[F(i == j) for j in range(6)] for i in range(6)]
z1, m11, m12, m21, m22, z3 = u
z4 = [-x for x in [1, 1, 0, 0, 1, 1]]
sub = lambda a, b: [x-y for x, y in zip(a, b)]
rows = [sub(z1, m11), m12, m21, sub(z1, m22), sub(z1, z3),
        sub(m11, z4), m21, sub(z3, z4)]
xi = list(map(F, [1, 1, 0, 0, -4, 1]))
mass = dot(xi, xi) + dot(z4, xi)**2
assert mass == 20
assert sum(dot(row, xi)**2 for row in [m21, sub(z1, m22)]) == 25


def gram(e, z):
    weights = [1+e, 1+e, e, e, (1+e)/z, 1, 1, z]
    mat = [[sum((w*row[i]*row[j] for w, row in zip(weights, rows)), F(0))
            for j in range(6)] for i in range(6)]
    return mat


cases = []
z_values = [F(3, 5), F(3, 5)+F(1, 10**40), F(61, 100), F(5, 8), F(2, 3), F(3, 4)]
for z in z_values:
    e = (5*z*z+2*z-3)/(4*z+3)
    assert 0 <= e < F(1, 4)
    assert 5*z*z+(2-4*e)*z-3-3*e == 0
    q = gram(e, z)
    qrank = rank(q)
    assert qrank == (5 if e == 0 else 6)
    assert dot(xi, [dot(row, xi) for row in q]) == 25*e
    for T in [F(1, 1000), F(1, 10**6), F(1, 10**12)]:
        A = e+1+z
        B = -2*e+7*(1+z)
        beta = 1/(1-T*B)
        p1 = 5*T*beta*A
        p2 = beta-1
        alpha = 1+p1
        assert p2 == T*beta*B and alpha > 0 and beta > 0
        D = beta*z/alpha
        # Moment of unscaled tuple, divided by beta.
        moment = [1+2*e+(1+e)/z, -e, z-(1+e)/z, -(1+z)]
        assert moment[0]+2*moment[1]+moment[2]+moment[3] == 0
        assert 3*moment[0]+2*moment[1]+5*moment[3] == 0
        # b / t^2 with T=t^2/(18*pi); pi cancels exactly.
        demands = [(F(16, 5)*p1+2*p2)/(18*T),
                   (-F(14, 5)*p1+2*p2)/(18*T),
                   (F(16, 5)*p1-4*p2)/(18*T),
                   (-F(4, 5)*p1-2*p2)/(18*T)]
        assert demands == [beta*x for x in moment]
        # Off-diagonal commutator entry squared of one row Gram with e1 e1*.
        commutator_square = (beta/alpha)**2*(1+e)*e/4
        assert (commutator_square > 0) == (e > 0)
        # All four arrow energies, before overall t^2.
        energy = beta*(2+2*e+(1+e)/z+z)
        assert energy > 0 and D > 0
        cases.append({'z': str(z), 'epsilon_squared': str(e), 'T': str(T),
                      'tracefree_rank': qrank, 'balance_exact': True,
                      'individual_noncommutation': e > 0})

out = {'method': 'exact rational arithmetic; geometric and analytic proof is separate',
       'cases': cases, 'case_count': len(cases), 'boundary_kernel_mass': str(mass),
       'first_variation_on_kernel': '25', 'weak_coefficient': str(F(25, 20)),
       'boundary_energy_coefficient': str(F(2)+F(5, 3)+F(3, 5)), 'pass': True}
(Path(__file__).resolve().parents[1] / 'output' / 'balanced_results.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps({k: v for k, v in out.items() if k != 'cases'}, indent=2))

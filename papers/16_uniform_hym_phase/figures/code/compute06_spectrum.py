#!/usr/bin/env python3
"""Exact finite check of the complex and Hermitian commutator forms. Standard library only.

The two coefficient Gram matrices are computed independently from commutators.
All complex arithmetic is Gaussian-integer arithmetic; rank computations use
Fraction. Equality of the coefficient matrices proves identities for all u,v.
"""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


def unit(i, j):
    return [[int(a == i and b == j) for b in range(3)] for a in range(3)]


def add(a, b):
    return [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]


def scale(c, a):
    return [[c * a[i][j] for j in range(3)] for i in range(3)]


def multiply(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3))
             for j in range(3)] for i in range(3)]


def commutator(a, b):
    return add(multiply(a, b), scale(-1, multiply(b, a)))


def inner(a, b):
    return sum(complex(a[i][j]).conjugate() * b[i][j]
               for i in range(3) for j in range(3))


def gram(basis, arrow=None, real=False):
    vectors = basis if arrow is None else [commutator(x, arrow) for x in basis]
    values = [[inner(x, y) for y in vectors] for x in vectors]
    if not real:
        assert all(z.imag == 0 for row in values for z in row)
    assert all(z.real == int(z.real) for row in values for z in row)
    return [[int(z.real) for z in row] for row in values]


def rank(a):
    a = [[Fraction(x) for x in row] for row in a]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        p = a[r][c]
        a[r] = [x / p for x in a[r]]
        for i in range(len(a)):
            if i != r:
                factor = a[i][c]
                a[i] = [x - factor * y for x, y in zip(a[i], a[r])]
        r += 1
    return r


def weighted(a, b, u, v):
    return [[u*x + v*y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def expected(diagonal, block):
    a = [[0] * 8 for _ in range(8)]
    for i, value in enumerate(diagonal):
        a[i][i] = value
    for i in range(2):
        for j in range(2):
            a[6+i][6+j] = block[i][j]
    return a


def main():
    e12, e23, e13 = unit(0, 1), unit(1, 2), unit(0, 2)
    d1 = [[1, 0, 0], [0, -1, 0], [0, 0, 0]]
    d2 = [[0, 0, 0], [0, 1, 0], [0, 0, -1]]
    complex_basis = [e13, e12, e23, unit(1, 0), unit(2, 1), unit(2, 0), d1, d2]
    mass_c, u_c, v_c = (gram(complex_basis), gram(complex_basis, e12),
                        gram(complex_basis, e23))
    assert mass_c == expected([1]*6, [[2, -1], [-1, 2]])
    assert u_c == expected([0, 0, 1, 2, 0, 1], [[4, -2], [-2, 1]])
    assert v_c == expected([0, 1, 0, 0, 2, 1], [[1, -2], [-2, 4]])

    real_basis = []
    for i, j in [(0, 1), (1, 2), (0, 2)]:
        e, f = unit(i, j), unit(j, i)
        real_basis += [add(e, f), scale(1j, add(e, scale(-1, f)))]
    real_basis += [d1, d2]
    mass_r, u_r, v_r = (gram(real_basis, real=True), gram(real_basis, e12, True),
                        gram(real_basis, e23, True))
    assert mass_r == expected([2]*6, [[2, -1], [-1, 2]])
    assert u_r == expected([2, 2, 1, 1, 1, 1], [[4, -2], [-2, 1]])
    assert v_r == expected([1, 1, 2, 2, 1, 1], [[1, -2], [-2, 4]])

    # These ranks confirm the already symbolic block formulas at each stratum.
    strata = [(1, 1), (1, 0), (0, 1), (0, 0)]
    complex_nullities = [8-rank(weighted(u_c, v_c, u, v)) for u, v in strata]
    real_nullities = [8-rank(weighted(u_r, v_r, u, v)) for u, v in strata]
    assert complex_nullities == [1, 4, 4, 8]
    assert real_nullities == [0, 1, 1, 8]
    assert multiply(e12, e23) == e13
    assert all(multiply(a, b) == [[0]*3 for _ in range(3)]
               for a, b in [(e23, e12), (e12, e12), (e23, e23),
                            (e12, e13), (e13, e12), (e23, e13),
                            (e13, e23), (e13, e13)])

    paper = Path(__file__).resolve().parents[2]
    sources = [paper / 'sections/06_application.tex']
    report = {
        'status': 'PASS: exact coefficient Gram identities and matrix products',
        'complex_basis': ['E13', 'E12', 'E23', 'E21', 'E32', 'E31',
                          'diag(1,-1,0)', 'diag(0,1,-1)'],
        'complex_mass': mass_c, 'complex_u_gram': u_c, 'complex_v_gram': v_c,
        'real_basis': ['E12+E21', 'i(E12-E21)', 'E23+E32', 'i(E23-E32)',
                       'E13+E31', 'i(E13-E31)', 'diag(1,-1,0)', 'diag(0,1,-1)'],
        'real_mass': mass_r, 'real_u_gram': u_r, 'real_v_gram': v_r,
        'diagonal_generalized_characteristic': 'lambda^2 - 2(u+v)lambda + 3uv',
        'strata_uv': strata,
        'complex_nullities': complex_nullities, 'real_nullities': real_nullities,
        'source_sha256': {str(p.relative_to(paper)): sha256(p.read_bytes()).hexdigest()
                          for p in sources},
        'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope': 'Finite identities only; no HYM decay, confinement or general theorem proof.'
    }
    output = paper / 'checks/finite_spectrum.json'
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(report['status'])
    print(f'Complex nullities: {complex_nullities}; Hermitian real nullities: {real_nullities}')
    print(output)


if __name__ == '__main__':
    main()

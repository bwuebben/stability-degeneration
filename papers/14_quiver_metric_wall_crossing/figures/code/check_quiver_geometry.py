#!/usr/bin/env python3
"""Exact rational checks for R01/R03; no numerical asymptotic inference.

Run from any directory with Python 3; only the standard library is used.
The output is written to checks/quiver_geometry.json relative to this project.
"""
from fractions import Fraction as F
from functools import reduce
from math import gcd
from pathlib import Path
import json


def transpose(a):
    return [list(r) for r in zip(*a)]


def mul(a, b):
    return [[sum(x*y for x, y in zip(r, c)) for c in zip(*b)] for r in a]


def zeros(n, m):
    return [[0 for _ in range(m)] for _ in range(n)]


def add_to(a, b, sign=1):
    for i, row in enumerate(b):
        for j, x in enumerate(row):
            a[i][j] += sign*x


def rank(a):
    a = [[F(x) for x in row] for row in a]
    if not a:
        return 0
    r = 0
    for c in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        d = a[r][c]
        a[r] = [x/d for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                d = a[i][c]
                a[i] = [x-d*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def paired_representation(dims):
    arrows = []
    for i in range(3):
        for j in range(i+1, 3):
            ni, nj = dims[i], dims[j]
            if not ni or not nj:
                continue
            if (ni, nj) == (2, 2):
                cs = [[[1, 0], [0, 1]], [[0, 1], [0, 0]]]
            elif (ni, nj) == (2, 1):
                cs = [[[1, 0]], [[0, 1]]]
            elif (ni, nj) == (1, 2):
                cs = [[[1], [0]], [[0], [1]]]
            else:
                cs = [[[1]], [[0]]]
            for c in cs:
                arrows.append((i, j, c, transpose(c)))
                arrows.append((i, j, c, [[-x for x in r] for r in transpose(c)]))
    return arrows


def check_representation(dims):
    pairs = paired_representation(dims)
    real = [zeros(n, n) for n in dims]
    complex_moment = [zeros(n, n) for n in dims]
    directed = []
    for i, j, a, b in pairs:
        add_to(real[j], mul(a, transpose(a)))
        add_to(real[i], mul(transpose(a), a), -1)
        add_to(real[i], mul(b, transpose(b)))
        add_to(real[j], mul(transpose(b), b), -1)
        add_to(complex_moment[j], mul(a, b))
        add_to(complex_moment[i], mul(b, a), -1)
        directed.extend([(i, j, a), (j, i, b)])
    coordinates = [(i, r, c) for i, n in enumerate(dims) for r in range(n) for c in range(n)]
    index = {key: k for k, key in enumerate(coordinates)}
    equations = []
    for i, j, a in directed:
        for r in range(dims[j]):
            for c in range(dims[i]):
                row = [0]*len(coordinates)
                for k in range(dims[j]):
                    row[index[j, r, k]] += a[k][c]
                for k in range(dims[i]):
                    row[index[i, k, c]] -= a[r][k]
                equations.append(row)
    real_zero = all(x == 0 for a in real for row in a for x in row)
    complex_zero = all(x == 0 for a in complex_moment for row in a for x in row)
    nullity = len(coordinates)-rank(equations)
    assert real_zero and complex_zero and nullity == 1, dims
    return {'dimension_vector': dims, 'real_moment_zero': real_zero,
            'complex_moment_zero': complex_zero, 'endomorphism_dimension': nullity,
            'equation_rank': rank(equations), 'endomorphism_variables': len(coordinates)}


def pairing(n, p):
    return -2*sum(x*y for x, y in zip(n, p)) + 4*sum(
        n[i]*p[j]+n[j]*p[i] for i in range(3) for j in range(i+1, 3))


def nonsimple_example():
    """The forward-only representation printed after the wall equations."""
    units = []
    for r in range(2):
        for c in range(2):
            a = zeros(2, 2)
            a[r][c] = 1
            units.append(a)
    directed = [(0, 1, a) for a in units]
    directed += [(i, 2, a) for i in (0, 1) for a in ([[1, 0]], [[0, 1]])]
    moment = [zeros(n, n) for n in (2, 2, 1)]
    for i, j, a in directed:
        add_to(moment[j], mul(a, transpose(a)))
        add_to(moment[i], mul(transpose(a), a), -1)
    level = [a[0][0] for a in moment]
    assert moment == [[[-3, 0], [0, -3]], [[1, 0], [0, 1]], [[4]]]
    assert 2*level[0]+2*level[1]+level[2] == 0
    assert all(x for x in (level[0], level[1], level[0]+level[1],
                           level[0]+2*level[1], 2*level[0]+level[1]))
    return {'real_level': level, 'complex_moment_zero': True,
            'proper_subrepresentation': [0, 0, 1], 'off_all_five_walls': True}


def main():
    gram_u = zeros(6, 6)
    for i in range(3):
        gram_u[2*i][2*i+1] = gram_u[2*i+1][2*i] = 1
    embedding = [[1, 0, 0], [1, 0, 0], [0, 1, 0],
                 [0, -2, 2], [0, 0, 1], [0, 0, -2]]
    gram = mul(mul(transpose(embedding), gram_u), embedding)
    assert gram == [[2, 0, 0], [0, -4, 2], [0, 2, -4]]
    assert [embedding[i] for i in (0, 2, 4)] == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    c = [[0, 1, 0], [0, 0, 1], [0, -1, -1]]
    square = lambda x: mul(mul([x], gram), transpose([x]))[0][0]
    differences = [square([x-y for x, y in zip(c[i], c[j])])
                   for i in range(3) for j in range(i+1, 3)]
    assert [square(x) for x in c] == [-4]*3 and differences == [-12]*3
    m = (2, 2, 1)
    first_chern = [sum(m[i]*c[i][j] for i in range(3)) for j in range(3)]
    ch2 = sum(m[i]*F(square(c[i]), 2) for i in range(3))
    c2 = F(square(first_chern), 2)-ch2
    s = sum(m)+ch2
    mukai_square = square(first_chern)-2*sum(m)*s
    assert (c2, s, mukai_square) == (8, -5, 46)
    degrees = [mul([x], gram)[0][1:] for x in c]
    mean_degree = [sum(m[i]*degrees[i][j] for i in range(3))/F(sum(m))
                   for j in range(2)]
    slopes = [[degrees[i][j]-mean_degree[j] for j in range(2)]
              for i in range(3)]
    geometric_parameter = [[-x for x in row] for row in slopes]
    assert geometric_parameter == [[F(18, 5), F(-12, 5)],
                                   [F(-12, 5), F(18, 5)],
                                   [F(-12, 5), F(-12, 5)]]
    assert rank(slopes) == 2
    assert all(sum(m[i]*slopes[i][j] for i in range(3)) == 0 for j in range(2))
    wall_normals = set()
    for i in range(3):
        for j in range(3):
            for k in range(2):
                n = (i, j, k)
                if n in ((0, 0, 0), m):
                    continue
                normal = (i-2*k, j-2*k)
                divisor = reduce(gcd, normal)
                normal = tuple(x//abs(divisor) for x in normal)
                if next(x for x in normal if x) < 0:
                    normal = tuple(-x for x in normal)
                wall_normals.add(normal)
    assert wall_normals == {(1, 0), (0, 1), (1, 1), (1, 2), (2, 1)}
    splits = [(1, 0, 0), (0, 1, 0), (1, 1, 0), (1, 2, 0), (2, 1, 0)]
    representations = {m}
    split_rows = []
    for n in splits:
        p = tuple(x-y for x, y in zip(m, n))
        representations.update((n, p))
        cross = pairing(n, p)
        tangent_dimension = pairing(n, n)+pairing(p, p)+4
        transverse_dimension = 2*(cross-1)
        assert tangent_dimension+transverse_dimension == mukai_square+2
        split_rows.append({'first': n, 'second': p, 'mukai_pairing': cross,
                           'first_square': pairing(n, n), 'second_square': pairing(p, p),
                           'tangent_dimension': tangent_dimension,
                           'transverse_dimension': transverse_dimension,
                           'projective_dimension': cross-1})
    # Single-arrow extension sign test, with coefficients in units of 2*pi*t.
    # At H+t*a: L1 has degree -4*t, L2 degree 2*t, and mean slope -t.
    # A:L1->L2 has moment (-|A|^2,+|A|^2) but L2 destabilizes its extension.
    two_degrees = [degrees[0][0], degrees[1][0]]
    two_mean = F(sum(two_degrees), 2)
    two_parameter = [two_mean-d for d in two_degrees]
    assert two_degrees == [-4, 2] and two_mean == -1
    assert two_parameter == [3, -3]
    assert two_parameter[1] < 0  # no forward-only balanced representation
    assert two_parameter[0] > 0  # reverse arrow has the correct stable-subbundle sign
    extension_dimension = 2*4*sum(m[i]*m[j] for i in range(3) for j in range(i+1, 3))
    group_dimension = sum(n*n for n in m)-1
    assert extension_dimension-2*group_dimension == mukai_square+2
    result = {'status': 'exact finite calculations only; geometric and analytic proofs separate',
              'embedding': embedding, 'gram': gram,
              'line_squares': [-4]*3, 'difference_squares': differences,
              'mukai_vector': [5, first_chern, int(s)], 'c2': int(c2),
              'mukai_square': int(mukai_square), 'expected_moduli_dimension': int(mukai_square+2),
              'extension_space_complex_dimension': extension_dimension,
              'effective_group_complex_dimension': group_dimension,
              'slope_matrix': [[str(x) for x in row] for row in slopes],
              'geometric_parameter_matrix_in_units_of_2pi':
                  [[str(x) for x in row] for row in geometric_parameter],
              'two_factor_extension_sign_test': {
                  'degrees_in_units_of_t': two_degrees,
                  'mean_slope_in_units_of_t': str(two_mean),
                  'parameter_in_units_of_2pi_t': [str(x) for x in two_parameter],
                  'destabilizing_forward_arrow_excluded': True},
              'potential_wall_normals': sorted(wall_normals),
              'number_of_chambers': 2*len(wall_normals),
              'wall_splits': split_rows,
              'stable_but_not_module_simple_example': nonsimple_example(),
              'simple_balanced_representations': [check_representation(d) for d in sorted(representations)]}
    destination = Path(__file__).resolve().parents[2]/'checks/quiver_geometry.json'
    destination.write_text(json.dumps(result, indent=2)+'\n')
    print(f'PASS: lattice, Mukai data, five potential wall lines, {len(representations)} simple balanced representations')
    print(destination)


if __name__ == '__main__':
    main()

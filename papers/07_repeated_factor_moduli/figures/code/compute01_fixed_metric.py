#!/usr/bin/env python3
"""Exact finite checks for the composable rank-five example; no PDE claim."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import itertools
import json


def transpose(a):
    return list(map(list, zip(*a)))


def mul(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def rref(a):
    a = [[F(x) for x in row] for row in a]
    pivots = []
    i = 0
    for j in range(len(a[0])):
        k = next((k for k in range(i, len(a)) if a[k][j]), None)
        if k is None:
            continue
        a[i], a[k] = a[k], a[i]
        v = a[i][j]
        a[i] = [x/v for x in a[i]]
        for k in range(len(a)):
            if k != i:
                v = a[k][j]
                a[k] = [x-v*y for x, y in zip(a[k], a[i])]
        pivots.append(j)
        i += 1
        if i == len(a):
            break
    return a, pivots


def quad(a, v):
    return sum(x*y for x, y in zip(v, [sum(aij*vj for aij, vj in zip(row, v)) for row in a]))


def gram(rows):
    return mul(transpose(rows), rows)


def line_cohom(degrees):
    if -1 in degrees:
        return [0]*4
    degree = sum(d <= -2 for d in degrees)
    dim = 1
    for d in degrees:
        dim *= d+1 if d >= 0 else -d-1
    result = [0]*4
    result[degree] = dim
    return result


types = [(-2,-2,4), (1,-2,1), (-2,1,1), (0,0,0)]
mult = [1,2,1,1]
cohom = [0]*4
blocks = {}
for i, j in itertools.product(range(4), repeat=2):
    dims = [v*mult[i]*mult[j] for v in line_cohom([a-b for a,b in zip(types[i],types[j])])]
    cohom = [a+b for a,b in zip(cohom,dims)]
    if dims[1] or dims[2]:
        blocks[f'{i+1}{j+1}'] = dims
cohom[0] -= 1
assert cohom == [6,119,5,0]

# Coordinates are z1, M11, M12, M21, M22, z3; z4 is fixed by trace zero.
embedding = [[int(i==j) for j in range(6)] for i in range(6)]
embedding.append([-1,-1,0,0,-1,-1])
mass = gram(embedding)
rows0_full = [
    [1,-1,0,0,0,0,0], [0,0,-1,0,0,0,0],  # first row arrow
    [1,0,0,0,0,-1,0],                       # scalar alternate path
    [0,1,1,0,0,0,-1], [0,0,0,1,1,0,-1],   # column (1,1)
    [0,0,0,0,0,1,-1],                       # final scalar arrow
]
rows1_full = [
    [0,0,0,-1,0,0,0], [1,0,0,0,-1,0,0],   # second row arrow
    [1,0,0,0,0,-1,0],                       # second scalar harmonic coefficient
]
rows0, rows1 = mul(rows0_full, embedding), mul(rows1_full, embedding)
rank0 = len(rref(rows0)[1])
rank_positive = len(rref(rows0+rows1)[1])
assert (rank0, rank_positive) == (5,6)
kernel = [1,1,0,5,-4,1]
assert mul(rows0, [[x] for x in kernel]) == [[0]]*6
assert quad(mass, kernel) == 45
assert quad(gram(rows1), kernel) == 50
coefficient = F(50,45)
assert coefficient == F(10,9)

# Dropping the off-diagonal matrix directions loses the boundary kernel.
scalar_columns = [0,1,4,5]
scalar_restriction = [[row[j] for j in scalar_columns] for row in rows0]
assert len(rref(scalar_restriction)[1]) == 4

# Check the full finite quadratic expression directly at rational test values.
checks = 0
for eps in [F(0), F(1,100), F(1,3)]:
    for v in [[1,2,3,4,5,6], [2,-1,F(1,2),-3,4,-2]]:
        full = [row[0] for row in mul(embedding, [[x] for x in v])]
        z,m11,m12,m21,m22,z3,z4 = full
        direct = ((z-m11)**2+m12**2+(z-z3)**2
                  +(m11+m12-z4)**2+(m21+m22-z4)**2+(z3-z4)**2
                  +eps**2*(m21**2+(z-m22)**2+(z-z3)**2))
        matrix = quad(gram(rows0),v)+eps**2*quad(gram(rows1),v)
        assert direct == matrix
        checks += 1

result = {
    'scope': 'exact finite algebra and cohomology only; not an elliptic or PDE proof',
    'cohomology': cohom,
    'blocks': blocks,
    'finite_rank_at_boundary': rank0,
    'finite_rank_for_positive_parameter': rank_positive,
    'tracefree_kernel_coordinates': kernel,
    'kernel_mass': 45,
    'first_variation_energy_in_units_of_a': 50,
    'weak_coefficient_in_units_of_a': str(coefficient),
    'scalar_only_restriction_rank': 4,
    'rational_direct_checks': checks,
    'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'pass': True,
}
output = (Path(__file__).resolve().parents[1] / 'output' / 'fixed_metric_results.json')
output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

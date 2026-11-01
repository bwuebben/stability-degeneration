"""Exact rational linear algebra for the finite geometric computations."""
from fractions import Fraction as Q

def rank(rows):
    a = [[Q(x) for x in row] for row in rows]
    if not a:
        return 0
    k = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(k, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[k], a[pivot] = a[pivot], a[k]
        p = a[k][j]
        a[k] = [x / p for x in a[k]]
        for i in range(k + 1, len(a)):
            c = a[i][j]
            a[i] = [x - c*y for x, y in zip(a[i], a[k])]
        k += 1
        if k == len(a):
            break
    return k


def determinant(rows):
    a = [[Q(x) for x in row] for row in rows]
    d = Q(1)
    for j in range(len(a)):
        k = next((i for i in range(j, len(a)) if a[i][j]), None)
        if k is None:
            return Q(0)
        if k != j:
            a[j], a[k] = a[k], a[j]
            d = -d
        p = a[j][j]
        d *= p
        for i in range(j + 1, len(a)):
            c = a[i][j] / p
            for t in range(j + 1, len(a)):
                a[i][t] -= c*a[j][t]
    return d

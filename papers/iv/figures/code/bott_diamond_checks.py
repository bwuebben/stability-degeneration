#!/usr/bin/env python3
"""Exact Chow-ring checks of the printed positive-twist obstruction identities."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import json


def check(t, u, H):
    A = ((0,0,0,0),(t,0,0,0),(t,0,0,0),(0,u,u,0))
    @lru_cache(None)
    def integral(e):
        for i in range(3, -1, -1):
            if e[i] >= 2:
                result = Q(0)
                for j in range(i):
                    if A[i][j]:
                        v = list(e)
                        v[i] -= 1
                        v[j] += 1
                        result += A[i][j]*integral(tuple(v))
                return result
        return Q(int(e == (1,1,1,1)))
    degrees = []
    for i in range(4):
        C = [-Q(x) for x in A[i]]
        C[i] = Q(2)
        degree = Q(0)
        for indices in product(range(4), repeat=3):
            coefficient = H[indices[0]]*H[indices[1]]*H[indices[2]]
            for j in range(4):
                e = [0]*4
                e[j] += 1
                for k in indices:
                    e[k] += 1
                degree += C[j]*coefficient*integral(tuple(e))
        degrees.append(degree)
    a,b,c,d = H
    expected = 6*d*(c-b)*(t*(b+c+u*d)+2*a)
    assert degrees[1]-degrees[2] == expected
    symmetric_expected = None
    if b == c:
        symmetric_expected = -4*(t-1)*d*(3*b*b+3*b*u*d+u*u*d*d)-6*a*d*(2*b+u*d)
        assert degrees[0]-degrees[1] == symmetric_expected < 0
    return {'t':t,'u':u,'H':list(map(str,H)), 'factor_degrees':list(map(str,degrees)),
            'middle_difference':str(expected),
            'symmetric_base_difference':None if symmetric_expected is None else str(symmetric_expected)}


if __name__ == '__main__':
    cases=[]
    for t,u in ((1,1),(2,1),(1,3),(4,2)):
        for H in ((Q(1),Q(1),Q(1),Q(1)),
                  (Q(2,3),Q(5,4),Q(5,4),Q(3,2)),
                  (Q(1,2),Q(2),Q(3),Q(4,3))):
            cases.append(check(t,u,H))
    out={'scope':'Twelve rational checks; uniform proof is in the manuscript.',
         'cases':cases,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dest=Path(__file__).resolve().parents[1]/'data/bott_diamond_checks.json'
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print('PASS:',len(cases),'exact Chow-ring comparisons, including nonsymmetric ample classes')

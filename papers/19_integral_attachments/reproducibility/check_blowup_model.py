"""Exact preliminary checks for the diagonal-blow-up construction.

No numerical approximation; standard library only. The geometric interpretation,
integral descent and applicability of Markman's correspondence require proofs.
"""
from fractions import Fraction as F
from pathlib import Path
import json

N = 22
PT = N + 1
Q = [[0] * N for _ in range(N)]
for k in range(3):
    Q[2*k][2*k+1] = Q[2*k+1][2*k] = 1
for off in (6, 14):
    for i in range(8): Q[off+i][off+i] = -2
    for i, j in [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(2,7)]:
        Q[off+i][off+j] = Q[off+j][off+i] = 1


def det(a):
    a = [list(map(int, row)) for row in a]
    sign, prev = 1, 1
    for k in range(len(a)-1):
        if not a[k][k]:
            j = next(j for j in range(k+1, len(a)) if a[j][k])
            a[k], a[j] = a[j], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, len(a)):
            x = a[i][k]
            for j in range(k+1, len(a)):
                numerator = pivot*a[i][j] - x*a[k][j]
                assert numerator % prev == 0
                a[i][j] = numerator // prev
            a[i][k] = 0
        prev = pivot
    return sign*a[-1][-1]


def inverse(a):
    n = len(a)
    b = [[F(x) for x in row] + [F(i==j) for j in range(n)]
         for i, row in enumerate(a)]
    for i in range(n):
        k = next(k for k in range(i,n) if b[k][i])
        b[i], b[k] = b[k], b[i]
        x = b[i][i]
        b[i] = [y/x for y in b[i]]
        for j in range(n):
            if j != i and b[j][i]:
                x = b[j][i]
                b[j] = [u-x*v for u,v in zip(b[j],b[i])]
    return [row[n:] for row in b]


QI = inverse(Q)
assert det(Q) == -1
assert all(x.denominator == 1 for row in QI for x in row)
assert all(QI[i][i] % 2 == 0 for i in range(N))


def add(*vs):
    out = {}
    for v in vs:
        for k,x in v.items(): out[k] = out.get(k,F(0))+x
    return {k:v for k,v in out.items() if v}


def scale(v, x): return {k:F(x)*a for k,a in v.items() if x*a}
def tensor(i,j): return {('t',i,j): F(1)}
def exc(i): return {('e',i): F(1)}


def surf(i,j):
    if i == 0: return j,1
    if j == 0: return i,1
    if i == PT or j == PT: return 0,0
    return PT,Q[i-1][j-1]


DIAGONAL = add(tensor(0,PT), tensor(PT,0),
               *[scale(tensor(i+1,j+1), QI[i][j])
                 for i in range(N) for j in range(N)])


def mul(a,b):
    out = {}
    for ka,va in a.items():
        for kb,vb in b.items():
            coeff = va*vb
            if ka[0] == kb[0] == 't':
                i,x = surf(ka[1],kb[1]); j,y = surf(ka[2],kb[2])
                term = scale(tensor(i,j), x*y*coeff)
            elif ka[0] == kb[0] == 'e':
                i,x = surf(ka[1],kb[1])
                term = scale(mul(DIAGONAL,tensor(i,0)), -x*coeff)
            else:
                e,t = (ka,kb) if ka[0] == 'e' else (kb,ka)
                i,x = surf(t[1],t[2]); j,y = surf(e[1],i)
                term = scale(exc(j),x*y*coeff)
            out = add(out,term)
    return out


def integral(v): return v.get(('t',PT,PT),F(0))/2
def pair(a,b): return integral(mul(a,b))
def mu(a): return add(*[scale(add(tensor(i+1,0),tensor(0,i+1)),a[i]) for i in range(N)])
DELTA = exc(0)
A = add(tensor(0,PT),tensor(PT,0))
D2 = mul(DELTA,DELTA)
C2 = add(scale(A,24),scale(D2,-3))

BASIS4 = [A]
LABELS = ['A']
for i in range(N):
    BASIS4.append(scale(exc(i+1),2)); LABELS.append(f'T{i}')
for i in range(N):
    for j in range(i+1,N):
        BASIS4.append(add(tensor(i+1,j+1),tensor(j+1,i+1)))
        LABELS.append(f'C{i}_{j}')
for i in range(N):
    BASIS4.append(add(tensor(i+1,i+1),scale(exc(i+1),-1)))
    LABELS.append(f'Q{i}')
assert len(BASIS4) == 276


def coords4(v):
    a = v.get(('t',PT,0),F(0))
    qs = [v.get(('t',i+1,i+1),F(0)) for i in range(N)]
    ts = [(v.get(('e',i+1),F(0))+qs[i])/2 for i in range(N)]
    cs = [v.get(('t',i+1,j+1),F(0)) for i in range(N) for j in range(i+1,N)]
    vals = [a]+ts+cs+qs
    assert add(*[scale(b,c) for b,c in zip(BASIS4,vals)]) == v
    return vals


def vec(*items):
    a = [F(0)]*N
    for i,x in items: a[i] = F(x)
    return a


def dot(a,b): return sum(a[i]*Q[i][j]*b[j] for i in range(N) for j in range(N))
def plus(a,b): return [x+y for x,y in zip(a,b)]
def times(a,k): return [k*x for x in a]
def mukai(a,b): return dot(a[1],b[1])-a[0]*b[2]-a[2]*b[0]

H = vec((0,1),(1,1))
D1 = vec((2,1),(3,-3))
D2K = vec((3,-4),(4,1),(5,-4))
v = (F(3),plus(D1,D2K),F(-4))
a0 = (F(1),vec(),F(1))
a1 = (F(1),D1,F(-2))
a2 = (F(1),D2K,F(-3))
ell = (F(1),H,F(2))
C = plus(plus(D1,D2K),times(H,-2))


def transform(x):
    k = mukai(x,ell)
    r,c,s = x[0]+k,plus(x[1],times(H,k)),x[2]+2*k
    return r,plus(c,times(C,-r)),s-dot(c,C)+r*dot(C,C)/2


def theta(x):
    assert x[0] == x[2]
    return add(mu(x[1]),scale(DELTA,x[0]))


def sparse(vals):
    return {k:int(v) if v.denominator == 1 else str(v)
            for k,v in zip(LABELS,vals) if v}


def main():
    gram = [[pair(a,b) for b in BASIS4] for a in BASIS4]
    assert all(x.denominator == 1 for row in gram for x in row)
    gramdet = det(gram)
    assert gramdet == 1
    assert pair(D2,D2) == 12
    assert pair(A,D2) == -1
    assert pair(C2,C2) == 828
    for i in range(N):
        for j in range(N):
            m = mul(mu(vec((i,1))),mu(vec((j,1))))
            assert pair(D2,m) == -2*Q[i][j]
            assert pair(C2,m) == 30*Q[i][j]
            assert all(x.denominator == 1 for x in coords4(m))
    assert transform(v) == (1,vec(),-1)
    cycles = []
    for ai in (a0,a2):
        curve = (v[0]/2-ai[0],plus(times(v[1],F(1,2)),times(ai[1],-1)),v[2]/2-ai[2])
        c = theta(transform(curve))
        plane = add(scale(mul(c,c),F(1,2)),scale(C2,F(1,24)))
        coeffs = coords4(plane)
        assert all(x.denominator == 1 for x in coeffs)
        cycles.append(plane)
    intersections = [[int(pair(a,b)) for b in cycles] for a in cycles]
    assert intersections == [[3,1],[1,3]]
    report = {'status':'exact algebra verified; geometric interpretation is proved in the manuscript',
              'rank_H4':len(BASIS4),'det_K3':det(Q),'det_H4':gramdet,
              'delta_fourth':int(pair(D2,D2)),'c2_squared':int(pair(C2,C2)),
              'plane_intersection_matrix':intersections,
              'plane_coordinates':[sparse(coords4(p)) for p in cycles],
              'c2_coordinates':sparse(coords4(C2)),
              'mukai_transform_v':[1,[0]*N,-1]}
    path=Path(__file__).with_name('blowup-model-check.json')
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if not k.endswith('coordinates')},indent=2))


if __name__ == '__main__': main()

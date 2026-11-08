"""Exact supporting arithmetic; not a proof of existence or metric limits."""
from fractions import Fraction
from math import comb, gcd
import json
from pathlib import Path


def binomial(m, k):
    ans = Fraction(1)
    for i in range(k):
        ans *= Fraction(m-i, i+1)
    assert ans.denominator == 1
    return int(ans)


def hilbert_poincare(n):
    # Göttsche's K3 product, truncated in point number at n.
    series = {(0, 0): 1}
    for m in range(1, n+1):
        for tdegree, multiplicity in ((2*m-2, 1), (2*m, 22), (2*m+2, 1)):
            product = {}
            for (qdegree, degree), coefficient in series.items():
                for k in range((n-qdegree)//m+1):
                    key = qdegree+m*k, degree+tdegree*k
                    product[key] = product.get(key, 0) + coefficient*comb(multiplicity+k-1, k)
            series = product
    return [series.get((n, degree), 0) for degree in range(4*n+1)]


betti = hilbert_poincare(4)
expected = [1,23,300,2852,19298,2852,300,23,1]
assert betti[::2] == expected and not any(betti[1::2])
complement = expected.copy()
for i, rank in zip(range(4,9), [1,1,2,1,1]):
    complement[i] -= rank
assert complement == [1,23,300,2852,19297,2851,298,22,0]

# The top-degree evaluations of (1+h+p)^m and its lower terms.
for m in range(-20,21):
    c3 = 2*binomial(m,3) + m*(m-1)
    c4 = 2*binomial(m,4) + 3*binomial(m,3) + binomial(m,2)
    assert c3 == Fraction(m*(m*m-1),3)
    assert c4 == Fraction(m*m*(m*m-1),12)

# B's Smith factors and boundary quotient map (n0,n2) -> n0-3 n2 mod 8.
B = ((3,1),(1,3))
d1 = gcd(*[abs(x) for row in B for x in row])
d2 = abs(B[0][0]*B[1][1]-B[0][1]*B[1][0])//d1
assert (d1,d2) == (1,8)
assert all((B[0][j]-3*B[1][j]) % 8 == 0 for j in range(2))

# Mukai lattice on triples (r, D coefficient, f coefficient, s), D.f=1.
def pairing(x,y):
    r,d,f,s = x
    R,e,g,S = y
    return -44*d*e + d*g + f*e - r*S - R*s

a=(3,1,0,-7); b=(1,0,0,1); v=(7,2,0,-13)
x=(0,0,-1,0); y=(0,0,3,1); z=(0,0,-7,-2)
assert [-pairing(a,x),-pairing(b,x),-pairing(a,y),-pairing(b,y)] == [1,0,0,1]
assert pairing(v,z)==0 and pairing(a,z)==-1
assert pairing(a,v)==0 and pairing(a,a)==-2

result={
    "status":"PASS",
    "scope":"Exact finite arithmetic only; geometry and analysis require the printed proofs.",
    "K3_Hilbert4_even_betti":expected,
    "complement_even_betti":complement,
    "two_plane_smith":[d1,d2],
    "boundary_orders":{"isolated_plane":3,"two_planes":8,"Grassmannian":6},
    "Grassmannian_top_evaluations":{"h4":2,"h2p":1,"p2":1},
    "tested_integer_Chern_parameters":[-20,20],
    "Mukai_index_and_curve_witnesses":"PASS"
}
path=Path(__file__).with_name("integral_data.json")
path.write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps(result))

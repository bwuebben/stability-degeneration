"""Exhaust the finite support patterns and verify signed characters exactly."""
from fractions import Fraction as F
from itertools import permutations

blocks=[(1,2),(3,2),(3,4),(1,3),(2,3),(2,4),(2,1),(3,1),(4,1)]
dims=[8,8,4,8,8,4,8,8,27]
assert sum(dims)==83
heights=[F(-3,2),F(-1,2),F(1,2),F(3,2)]
assert sum(heights)==0 and sum(x*x for x in heights)==5

def acyclic(edges):
    incoming={i:0 for i in range(1,5)}
    for i,j in edges: incoming[i]+=1
    ready=[i for i in incoming if incoming[i]==0]
    seen=0
    while ready:
        j=ready.pop()
        seen+=1
        for i,k in edges:
            if k==j:
                incoming[i]-=1
                if incoming[i]==0: ready.append(i)
    return seen==4

for mask in range(1<<len(blocks)):
    edges=[e for k,e in enumerate(blocks) if mask>>k&1]
    valid=[p for p in permutations(heights)
           if all(p[i-1]-p[j-1]>=1 for i,j in edges)]
    assert bool(valid)==acyclic(edges)
    for p in valid:
        assert sum(p)==0 and sum(z*z for z in p)==5
# A cycle monomial has total diagonal character zero, by telescoping.
for cycle in permutations(range(1,5),4):
    char=[0]*4
    for k,j in enumerate(cycle):
        i=cycle[(k+1)%4]
        char[i-1]+=1; char[j-1]-=1
    assert char==[0]*4
unstable=[(1,2),(1,3),(3,2),(2,4),(3,4)]
assert sum(dims[blocks.index(e)] for e in unstable)==32
b=[F(-1),F(2,5),F(0),F(3,5)]
assert [b[i-1]-b[j-1] for i,j in unstable]==[
    F(-7,5),F(-1),F(-2,5),F(-1,5),F(-3,5)]
curv=[-b[i-1] for i in [1,3,2,4]]
assert curv==[F(1),F(0),F(-2,5),F(-3,5)]
assert all(curv[i]>curv[i+1] for i in range(3))
h=[F(3,2),F(-1,2),F(1,2),F(-3,2)]
assert sum(h)==0 and sum(z*z for z in h)==5
assert all(h[i-1]-h[j-1]>=1 for i,j in unstable)
print("PASS: all 512 support patterns, height norm, signed characters and HN order.")

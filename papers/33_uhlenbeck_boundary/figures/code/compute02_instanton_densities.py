"""Exact localization sums, Nakajima--Yoshioka Theorem 2.11.
Run with Python 3 and SymPy; cancellation precedes specialization.
"""
import sympy as S
from itertools import product

def partitions(n,m=None):
 if n==0: yield (); return
 for k in range(min(n,m if m else n),0,-1):
  for q in partitions(n-k,k): yield (k,)+q

def tuples(r,n):
 if r==1:
  for q in partitions(n): yield (q,)
 else:
  for k in range(n+1):
   for q in partitions(k):
    for t in tuples(r-1,n-k): yield (q,)+t

def cells(q): return [(i,j) for i,l in enumerate(q) for j in range(l)]
def arm(q,i,j): return q[i]-j-1
def leg(q,i,j): return sum(l>j for l in q)-i-1
h=S.symbols('h'); e1=S.Integer(1);e2=1+3*h
for r,n in [(2,1),(2,2),(2,3),(3,1),(3,2),(1,3)]:
 total=0
 for ys in tuples(r,n):
  den=1
  for a,qa in enumerate(ys):
   for b,qb in enumerate(ys):
    ab=(b-a)*7*h
    for i,j in cells(qa): den*=ab-leg(qb,i,j)*e1+(arm(qa,i,j)+1)*e2
    for i,j in cells(qb): den*=ab+(leg(qa,i,j)+1)*e1-arm(qb,i,j)*e2
  if den==0: raise RuntimeError((r,n,ys))
  total+=1/den
 val=S.limit(S.cancel(total),h,0)
 expected={(2,1):S.Rational(1,2),(2,2):S.Rational(11,108),(2,3):S.Rational(61,5184),(3,1):S.Rational(3,8),(3,2):S.Rational(61,1296),(1,3):S.Rational(1,6)}
 assert val==expected[r,n], ((r,n),val)
 print((r,n),val, 'PASS', flush=True)

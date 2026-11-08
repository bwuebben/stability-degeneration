"""Check the printed integral Gram pairing and all six attachment systems.

Uses the sparse blow-up arithmetic as an interpreter and reconstructs the
pairing from the displayed formulas. Includes every boundary edge map.
"""
import sys
sys.dont_write_bytecode = True
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_blowup_model as m
import check_universal_indices as ui
import check_wall_action as wa

# The wall root and A-image have independently simplified expressions.
k = m.plus(m.times(m.H, 2), m.times(m.D2K, -1))
mu_k = m.mu(k)
r = m.add(mu_k, m.DELTA)
assert m.dot(k,k) == 0 and r == wa.ROOT
simple_A = m.add(m.A, m.scale(m.mul(mu_k, m.add(mu_k,m.DELTA)), F(1,2)))
assert simple_A == wa.IA

# Construct the degree-four pairing separately from the printed formula table.
descriptors=[('A',)] + [('T',i) for i in range(22)] + [('C',i,j) for i in range(22) for j in range(i+1,22)] + [('Q',i) for i in range(22)]
order={'A':0,'T':1,'C':2,'Q':3}
def gram_value(d,e):
    if order[d[0]] > order[e[0]]: d,e=e,d
    t,u=d[0],e[0]
    if t=='A': return int(u=='A')
    if t=='T' and u=='T': return -2*m.Q[d[1]][e[1]]
    if t=='T' and u=='C': return 0
    if t=='T' and u=='Q': return m.Q[d[1]][e[1]]
    if t=='C' and u=='C':
        i,j=d[1:]; a,b=e[1:]
        return m.Q[i][a]*m.Q[j][b]+m.Q[i][b]*m.Q[j][a]
    if t=='C' and u=='Q':
        i,j=d[1:]; a=e[1]
        return m.Q[a][i]*m.Q[a][j]
    if t=='Q' and u=='Q':
        q=m.Q[d[1]][e[1]]
        return (q*q-q)//2
    raise AssertionError((d,e))
G=[[gram_value(d,e) for e in descriptors] for d in descriptors]
cols=[]
for image in wa.IB4:
    vals=m.coords4(image)
    assert all(v.denominator==1 for v in vals)
    cols.append({i:int(v) for i,v in enumerate(vals) if v})
def apply_column(c):
    out={}
    for j,n in c.items():
        for i,v in cols[j].items(): out[i]=out.get(i,0)+n*v
    return {i:v for i,v in out.items() if v}
for j,c in enumerate(cols): assert apply_column(c)=={j:1}
for i in range(276):
    for j in range(i,276):
        computed=sum(a*b*G[s][t] for s,a in cols[i].items() for t,b in cols[j].items())
        assert computed==G[i][j], (i,j)

p0,p2=[p for _,p in ui.CURVES_PLANES]
# Explicit witnesses and all six attachment templates on the 23 divisor generators.
x=(F(0),m.vec((3,1),(5,-1)),F(0))
y0=(F(0),m.vec((3,-3)),F(-1))
a=m.theta(m.transform(x)); b=m.theta(m.transform(y0))
assert [m.mukai(x,z) for z in (m.a0,m.a1,m.a2)]==[0,1,-1]
assert [m.mukai(y0,z) for z in (m.a0,m.a1,m.a2)]==[1,-2,1]
lam=[m.theta(c) for c,p in ui.CURVES_PLANES]
def bb(t,s):
    tc=ui.coords2(t); sc=ui.coords2(s)
    return m.dot(tc[:-1],sc[:-1])-2*tc[-1]*sc[-1]
def ds(t): return [bb(t,l) for l in lam]
assert ds(a)==[0,1] and ds(b)==[-1,-1]
assert [m.pair(m.mul(a,a),p) for p in (p0,p2)]==[0,1]
assert [m.pair(m.mul(b,b),p) for p in (p0,p2)]==[1,1]

divs=wa.MU+[m.DELTA]
chambers={'A':('plane',1),'D':('plane',-1),'B':(2,0,-1),'C':(0,2,1),'E':(2,0,1),'F':(0,2,-1)}
for ch,entry in chambers.items():
    for u,v in product(divs,repeat=2):
        du,dv=ds(u),ds(v)
        prod=m.mul(u,v)
        if entry[0]=='plane':
            eps=entry[1]
            s0,s2=eps*du[0],eps*du[1]
            t0,t2=eps*dv[0],eps*dv[1]
            n=s0*t0-3*s2*t2
            rhs=m.pair(prod,m.add(p0,m.scale(p2,-3)))
        else:
            ki,ji,eps=entry
            ki=ki//2; ji=ji//2
            s,t=eps*du[ki],-eps*(du[0]+du[1])
            ss,tt=eps*dv[ki],-eps*(dv[0]+dv[1])
            n=6*s*ss+s*tt+t*ss+t*tt
            rhs=m.pair(prod,m.add([p0,p2][ji],m.scale([p0,p2][ki],5)))
            assert s*tt+t*ss+t*tt == m.pair(prod,m.add([p0,p2][ji],m.scale([p0,p2][ki],-1)))
        assert (n-rhs)%8==0,(ch,ui.coords2(u),ui.coords2(v))

# Every integral degree-four basis class, not merely divisor products.
for z in m.BASIS4:
    n0,n2=m.pair(z,p0),m.pair(z,p2)
    kap=n0-3*n2
    assert (6*n2+(n0-n2)-kap)%8==0
    assert (6*n0+(n2-n0)-5*kap)%8==0

# Corrected middle basis transformation and its discriminant functional.
B=[[3,1],[1,3]]; BP=[[3,-2],[-2,4]]; T=[[1,1],[0,1]]
for i,j in product(range(2),repeat=2):
    assert sum(T[s][i]*BP[s][t]*T[t][j] for s,t in product(range(2),repeat=2))==B[i][j]
assert [sum([6,1][i]*BP[i][j] for i in range(2))%8 for j in range(2)]==[0,0]
assert [sum([1,-3][i]*B[i][j] for i in range(2))%8 for j in range(2)]==[0,0]

# All-degree index normalization under the reflection, with the explicit line exp(r).
inputs=[m.a0]+[(F(0),m.vec((i,1)),F(0)) for i in range(22)]+[(F(0),m.vec(),F(1))]
assert ui.part(ui.normalized_ch(wa.source_reflect(ui.Y)),2)==m.scale(r,-1)
index_count=0
for z in inputs:
    old=ui.chern_classes(ui.normalized_ch(z))
    new=ui.chern_classes(m.mul(ui.normalized_ch(wa.source_reflect(z)),ui.exp(r)))
    for degree,oldc,newc in zip((2,4,6,8),old,new):
        assert wa.reflect(oldc)==newc
        ui.encode(newc,degree)
        index_count+=1

report = {'H4_integral_involutions':276,'all_H4_pairings':276*277//2,
       'six_chamber_degree2_products':6*23*23,'all_H4_boundary_functionals':276,
       'normalized_index_classes':index_count,'primitive_witnesses':'pass',
       'A_image_independent_simplification':'pass','status':'all exact checks passed'}

# Every degree-six basis vector and its full pairing with reflected H2.
for i,im in enumerate(wa.IB6):
    ui.encode(im,6)
    old=m.mul(m.mu(m.QI[i]),m.A) if i<22 else m.scale(m.mul(m.DELTA,m.A),-1)
    assert wa.reflect(im)==old
    for j,h in enumerate(wa.IMU+[wa.IDELTA]):
        assert m.integral(m.mul(h,im))==int(i==j)

def matmul(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
def matinv(a):
    d=a[0][0]*a[1][1]-a[0][1]*a[1][0]
    assert abs(d)==1
    return [[a[1][1]//d,-a[0][1]//d],[-a[1][0]//d,a[0][0]//d]]
def transpose(a): return [list(row) for row in zip(*a)]
I=[[1,0],[0,1]]; S=[[0,-1],[-1,0]]
boundary_basis={'A':I,'B':[[0,-1],[1,1]],'C':[[1,0],[-1,-1]],'D':[[-1,0],[0,-1]],'E':[[0,1],[-1,-1]],'F':[[-1,0],[1,1]]}
kappa={'A':1,'B':1,'C':5,'D':1,'E':1,'F':5}
P=I; Q=I; tors=1; rows=[]
for src,dst in zip('ABCDEF','BCDEFA'):
    edge=S if src in 'BE' else I
    action=matmul(matmul(boundary_basis[dst],edge),matinv(boundary_basis[src]))
    action5=transpose(matinv(action))
    mod8=(kappa[dst]*(5 if src in 'BE' else 1)*kappa[src])%8
    P=matmul(action,P); Q=matmul(action5,Q); tors=tors*mod8%8
    rows.append({'edge':src+dst,'H2':action,'H5':action5,'H4_multiplier_mod8':mod8})
assert P==Q==I and tors==1
report.update({'H6_integral_involutions':23,'H2_H6_pairings':529,'six_boundary_edges':rows,
               'boundary_loop':'identity in degrees 2,4,5,7'})
print(json.dumps(report,indent=2,sort_keys=True))

#!/usr/bin/env sage -python
"""Exact toric intersection certificate for the original-side nodal X'_9.

Run with SageMath: sage -python compute03_x9.py --output x9_exact.json
No numerical approximation or source-project mutation is used. Rays, coarse
fan and integral support heights are transcribed from CY VI's public x9.json;
the fine fan, intersection form and quotient spectrum are computed here.
"""
import argparse
import json
from itertools import combinations
from pathlib import Path
from sage.all import QQ, ZZ, Fan, Polyhedron, ToricVariety, matrix, vector, PolynomialRing

RAYS = [(1,0,0,0),(0,1,0,0),(-1,-1,0,0),(0,0,1,0),(0,0,0,1),
        (-4,-2,-1,0),(-4,-2,0,-1),(3,1,1,1),(2,1,1,1),(-1,0,0,0),(-2,-1,0,0)]
HEIGHT = [2,2,2,3,3,3,3,2,2,1,2]
COARSE = [[0,1,3,6],[0,1,3,7],[0,1,4,5],[0,1,4,7],[0,1,5,6],
          [0,2,3,6,7],[0,2,4,5,7],[0,2,5,6],[1,3,6,9],[1,3,7,8],
          [1,3,8,9],[1,4,5,9],[1,4,7,8],[1,4,8,9],[1,5,6,9],
          [2,3,6,7,8,10],[2,4,5,7,8,10],[2,5,6,10],
          [3,6,8,9,10],[4,5,8,9,10],[5,6,9,10]]
NODES = [(3,6,10,8),(2,6,3,7),(4,5,10,8),(2,5,4,7)]
BASIS = [2,5,6,7,8,10]

def polytope(heights):
    return Polyhedron(ieqs=[[heights[i]]+list(v) for i,v in enumerate(RAYS)], base_ring=QQ)

def active(heights, point):
    return [i for i,v in enumerate(RAYS) if heights[i]+vector(QQ,v).dot_product(point)==0]

def rows(A):
    return [[str(x) for x in row] for row in A.rows()]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args=parser.parse_args()
    P=polytope(HEIGHT)
    assert P.is_compact() and P.dim()==4
    vertices=[vector(QQ,p) for p in P.vertices_list()]
    assert sorted(active(HEIGHT,p) for p in vertices)==sorted(COARSE)
    certificate=[]
    for p in vertices:
        assert all(x in ZZ for x in p)
        cone=active(HEIGHT,p)
        margin=min(HEIGHT[i]+vector(QQ,v).dot_product(p) for i,v in enumerate(RAYS) if i not in cone)
        assert margin>=1
        certificate.append({'cone':cone,'cartier_vertex':[int(x) for x in p],'omitted_margin':str(margin)})
    fine_h=[QQ(h)-QQ(int(i in (3,4)))/100+QQ(i*i)/1000000 for i,h in enumerate(HEIGHT)]
    P_fine=polytope(fine_h)
    fine=sorted(active(fine_h,vector(QQ,p)) for p in P_fine.vertices_list())
    assert P_fine.is_compact() and len(fine)==29 and all(len(c)==4 for c in fine)
    for c in fine:
        V=matrix(ZZ,[RAYS[i] for i in c])
        assert abs(V.det())==1 and any(set(c)<=set(old) for old in COARSE)
        # Each cone is contained in a facet of the original reflexive polytope.
        ell=V.solve_right(vector(QQ,[-1]*4))
        assert all(1+vector(QQ,v).dot_product(ell)>=0 for v in RAYS)
    for a,b,c,d in NODES:
        assert vector(RAYS[a])+vector(RAYS[c])==vector(RAYS[b])+vector(RAYS[d])
        assert any({a,b,c}<=set(s) for s in fine) and any({a,c,d}<=set(s) for s in fine)
        assert not any({b,d}<=set(s) for s in fine)
    fan=Fan(fine,rays=RAYS)
    assert fan.is_complete() and fan.is_smooth()
    T=ToricVariety(fan)
    A=T.cohomology_ring(); D=A.gens(); Y=sum(D)
    H=sum(HEIGHT[i]*D[i] for i in range(11))
    assert D[9]*Y==0
    Hcoords=vector(QQ,[6,18,18,-12,-10,8])
    assert (H-sum(Hcoords[i]*D[j] for i,j in enumerate(BASIS)))*Y==0
    B=matrix(QQ,6,6,lambda i,j:T.integrate(D[BASIS[i]]*D[BASIS[j]]*H*Y))
    expected_B=matrix(QQ,[[-18,6,6,6,0,0],[6,0,14,0,0,2],[6,14,0,0,0,2],
                          [6,0,0,-22,12,0],[0,0,0,12,-16,2],[0,2,2,0,2,-6]])
    assert B==expected_B and B.det()==-5064192
    volume=T.integrate(H**3*Y)
    assert volume==8712 and Hcoords*B*Hcoords==volume
    circuits=[]
    for a,b,c,d in NODES:
        row=[0]*11
        for j in (b,d): row[j]+=1
        for j in (a,c): row[j]-=1
        circuits.append(row)
    C=matrix(QQ,[[row[j] for j in BASIS] for row in circuits])
    relation=vector(QQ,[1,-1,-1,1])
    assert C.rank()==3 and relation*C==0 and C*Hcoords==0
    K=matrix(QQ,[[1,0,0],[0,1,0],[0,1,0],[1,-1,0],[0,0,1],[0,1,1]])
    Q=matrix(QQ,6,3,lambda i,j:int(i==j))
    V=K.augment(Q)
    assert C*K==0 and V.det()!=0
    assert K*vector(QQ,[6,18,-10])==Hcoords
    Bkk=K.transpose()*B*K; Bdk=Q.transpose()*B*K
    M=-Q.transpose()*B*Q+Bdk*Bkk.inverse()*Bdk.transpose()
    E=(C*Q).transpose()*(C*Q)
    assert M==matrix(QQ,[[177,-79,-79],[-79,158,-52],[-79,-52,158]])/15
    assert E==matrix(QQ,[[2,-1,-1],[-1,2,0],[-1,0,2]])
    assert all(M[:i,:i].det()>0 and E[:i,:i].det()>0 for i in range(1,4))
    PR=PolynomialRing(QQ,'mu'); mu=PR.gen()
    characteristic=PR((E-mu*M).det()/((-1)**3*M.det()))
    assert characteristic==(mu-QQ(1)/7)*(mu**2-QQ(375)/628*mu+QQ(45)/628)
    # A unit minor detects the full primitive period image in relation^perp.
    unit_minors=[]
    for rr in combinations(range(4),3):
        for cc in combinations(range(6),3):
            if abs(C.matrix_from_rows_and_columns(rr,cc).det())==1:
                unit_minors.append([list(rr),list(cc)])
    assert unit_minors
    output={'family':'original-side nodal partial resolution Xprime9 of CY VI',
            'rays':RAYS,'coarse_support_height':HEIGHT,'coarse_fan':COARSE,
            'cartier_certificate':sorted(certificate,key=lambda x:x['cone']),
            'fine_height':[str(x) for x in fine_h],'fine_fan':fine,
            'all_fine_cones_unimodular':True,'nodes_ac_diagonal':NODES,
            'divisor_basis':BASIS,'H_coordinates':list(map(str,Hcoords)),
            'H_cube':str(volume),'B':rows(B),'B_determinant':str(B.det()),
            'periods':rows(C),'relation':list(map(str,relation)),
            'kernel_columns':rows(K),'quotient_columns':rows(Q),
            'B_kernel':rows(Bkk),'mass':rows(M),'energy':rows(E),
            'mass_determinant':str(M.det()),'characteristic':str(characteristic),
            'characteristic_factorization':str(characteristic.factor()),
            'unit_period_minors':unit_minors,
            'status':'exact finite arithmetic; geometric and analytic hypotheses require separate proof'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2)+'\n')
    print('PASS: 21 integral Cartier vertices; 29 smooth crepant fine cones; H^3=8712.')
    print('PASS: rank(periods)=3; full quotient mass; characteristic='+str(characteristic.factor()))
    print('Certificate: '+str(args.output))

if __name__=='__main__':
    main()

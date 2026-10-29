#!/usr/bin/env python3
"""Exact common toric fourfold: slope wall, anticanonical nodes and CY counts.

No source-project modules. All lattice vertices are exhaustively enumerated
in rational arithmetic. Geometric implications are proved/reviewed separately.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache, reduce
from itertools import combinations, product
from math import gcd
from pathlib import Path
import hashlib
import json
from exact_linear import rank, determinant

R=((1,0,0,0),(-1,0,0,0),(0,1,0,0),(0,-1,0,0),
   (0,0,1,0),(1,1,-1,0),(0,0,0,1),(-1,-1,0,-1))
FAN=[tuple(2*j+b for j,b in enumerate(bits)) for bits in product((0,1),repeat=4)]


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def solve(rows,rhs):
    a=[[Q(x) for x in row]+[Q(y)] for row,y in zip(rows,rhs)]
    for j in range(4):
        k=next((k for k in range(j,4) if a[k][j]),None)
        if k is None:
            return None
        a[j],a[k]=a[k],a[j]
        p=a[j][j];a[j]=[x/p for x in a[j]]
        for i in range(4):
            if i!=j:
                c=a[i][j];a[i]=[x-c*y for x,y in zip(a[i],a[j])]
    return tuple(row[-1] for row in a)


@lru_cache(None)
def integral(e):
    # Coordinate order f,g,xi1,xi2; xi_i^2=(f-g)xi_i.
    for i in (3,2):
        if e[i]>=2:
            a=list(e);b=list(e)
            a[i]-=1;b[i]-=1;a[0]+=1;b[1]+=1
            return integral(tuple(a))-integral(tuple(b))
    return Q(int(e==(1,1,1,1)))


def intersection(classes):
    ans=Q(0)
    for indices in product(range(4),repeat=4):
        c=Q(1);e=[0]*4
        for form,j in zip(classes,indices):
            c*=form[j];e[j]+=1
        if c:
            ans+=c*integral(tuple(e))
    return ans


def lattice_counts(points,normals):
    bycodim=Counter();groups=defaultdict(list)
    for v in points:
        active=tuple(i for i,h in enumerate(normals) if dot(v,h)==1)
        codim=rank([normals[i] for i in active])
        bycodim[codim]+=1;groups[active].append(v)
    return dict(sorted(bycodim.items())),groups


def build():
    delta_vertices=set()
    for ids in combinations(range(8),4):
        h=solve([R[i] for i in ids],[1]*4)
        if h is not None and all(dot(h,r)<=1 for r in R):
            assert all(x.denominator==1 for x in h)
            delta_vertices.add(tuple(map(int,h)))
    delta_vertices=sorted(delta_vertices)
    assert len(delta_vertices)==12
    for f in FAN:
        assert abs(determinant([R[i] for i in f]))==1
        h=solve([R[i] for i in f],[1]*4)
        assert h in delta_vertices
    assert all(rank([h for h in delta_vertices if dot(h,r)==1])==4 for r in R)
    delta=[v for v in product(range(-1,2),range(-1,2),range(-3,2),range(-3,2))
           if all(dot(v,r)<=1 for r in R)]
    dual=[v for v in product(range(-1,2),repeat=4)
          if all(dot(v,h)<=1 for h in delta_vertices)]
    assert len(delta)==69 and len(dual)==9
    dc,dgroups=lattice_counts(delta,R)
    qc,qgroups=lattice_counts(dual,delta_vertices)
    assert dc=={0:1,1:4,2:22,3:30,4:12} and qc=={0:1,4:8}
    correction=0
    for active,points in dgroups.items():
        if rank([R[i] for i in active])==2:
            assert len(active)==2
            length=reduce(gcd,(abs(a-b) for a,b in zip(R[active[0]],R[active[1]])))
            correction+=len(points)*(length-1)
    assert correction==0
    h11=len(dual)-5-qc.get(1,0)
    h21=len(delta)-5-dc.get(1,0)+correction
    assert (h11,h21)==(4,60)
    singular=[]
    for ids in ((0,2,4,5),(1,3,6,7)):
        assert rank([R[i] for i in ids])==3
        assert all(R[ids[0]][j]+R[ids[1]][j]==R[ids[2]][j]+R[ids[3]][j] for j in range(4))
        endpoints=[v for v in delta_vertices if all(dot(v,R[i])==1 for i in ids)]
        assert len(endpoints)==2
        length=reduce(gcd,(abs(a-b) for a,b in zip(*endpoints)))
        assert length==4
        assert any(abs(determinant([[R[i][j] for j in cols] for i in ids[:3]]))==1
                   for cols in combinations(range(4),3))
        singular.append({'circuit':ids,'dual_edge':endpoints,'anticanonical_degree':length})
    nonsimplicial=[]
    for v in delta_vertices:
        active=[i for i,r in enumerate(R) if dot(v,r)==1]
        if len(active)>4:
            assert len(active)==5
            assert sum(set(x['circuit']).issubset(active) for x in singular)==1
            nonsimplicial.append(active)
    assert len(nonsimplicial)==4
    ray_classes=(( -1,1,1,0),(0,0,1,0),(-1,1,0,1),(0,0,0,1),
                 (1,0,0,0),(1,0,0,0),(0,1,0,0),(0,1,0,0))
    slope_cases=[]
    for r in (Q(5,4),Q(3,2)):
        H=(r-1,r+1,1,1)
        coefficients=(0,1,0,1,r-1,0,r+1,0)
        vertices=[]
        for f in FAN:
            v=solve([R[i] for i in f],[-coefficients[i] for i in f])
            assert all(dot(v,R[i])==-coefficients[i] for i in f)
            assert all(dot(v,R[i])>-coefficients[i] for i in range(8) if i not in f)
            vertices.append(v)
        assert len(set(vertices))==16
        degrees=[intersection((C,H,H,H)) for C in ray_classes]
        assert degrees==[6*r*r-2]*4+[6*r]*4
        mu=sum(degrees)/4
        flats={}
        for mask in range(1,1<<8):
            ids=[i for i in range(8) if mask>>i&1]
            k=rank([R[i] for i in ids])
            if k==4:
                continue
            closed=tuple(i for i,v in enumerate(R) if rank([R[j] for j in ids]+[v])==k)
            flats[closed]=k
        assert len(flats)==32
        margins=[(sum(degrees[i] for i in ids)-k*mu,ids,k) for ids,k in flats.items()]
        maximum=max(margins)
        status='stable' if maximum[0]<0 else 'strictly semistable' if maximum[0]==0 else 'unstable'
        assert status==('stable' if r==Q(5,4) else 'unstable')
        slope_cases.append({'r':str(r),'H_in_f_g_xi1_xi2':list(map(str,H)),
                            'ray_degrees':list(map(str,degrees)),'tangent_slope':str(mu),
                            'proper_flat_count':len(flats),'maximum_margin':str(maximum[0]),
                            'flat_inequalities':[{'members':list(ids),'rank':k,
                                'degree':str(sum(degrees[i] for i in ids)),
                                'margin':str(sum(degrees[i] for i in ids)-k*mu)}
                                for ids,k in sorted(flats.items())],
                            'maximizing_flat':maximum[1],'status':status,
                            'maximal_cones':FAN,'strict_ample_support_checks':True})
    D=(0,4,2,2) # sum of the eight ray divisor classes
    basis=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))
    gram=[[intersection((a,b,D,D)) for b in basis] for a in basis]
    gram_det=determinant(gram)
    assert gram_det!=0
    curves=[(1,0,0,0)]*4+[(0,1,-1,-1)]*4
    assert rank(curves)==2
    relation=(1,1,1,-3,1,1,1,-3)
    assert all(relation) and all(sum(c*v[j] for c,v in zip(relation,curves))==0 for j in range(4))
    return {'rays':R,'smooth_refinement_fan':FAN,
            'polar_vertices':delta_vertices,'polar_lattice_points':delta,'dual_lattice_points':dual,
            'polar_lattice_codimensions':dc,'dual_lattice_codimensions':qc,
            'hodge_numbers_by_Batyrev_formula':[h11,h21],
            'singular_curve_circuits':singular,'nonsimplicial_anticanonical_maximal_cones':nonsimplicial,
            'generic_node_count_from_two_degree_four_sections':8,
            'slope_cases':slope_cases,'divisor_restriction_gram':[[str(x) for x in row] for row in gram],
            'gram_determinant':str(gram_det),'exceptional_curve_vectors':curves,
            'curve_span_rank':2,'relation_kernel_dimension':6,'all_nonzero_relation':relation,
            'smoothing_hodge_numbers':[2,66],
            'geometric_scope':'Exact lattice and intersection data; geometric conclusions use the manuscript proofs.',
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__=='__main__':
    data=build()
    dest=Path(__file__).resolve().parents[1]/'data/cy_toric_comparison.json'
    dest.write_text(json.dumps(data,indent=2)+'\n')
    print('PASS: polar lattice counts 69/9; Hodge formula (4,60); two degree-4 circuit curves; 32-flat slope comparisons; curve relation rank 2, kernel 6')
    print('Restriction Gram determinant:',data['gram_determinant'])

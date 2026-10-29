#!/usr/bin/env python3
"""Exact independent polygon and all-flat computation for the blowup wall."""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import hashlib
import json

R=((1,0),(-1,0),(0,1),(1,-1),(1,1))
CONES=((0,4),(4,2),(2,1),(1,3),(3,0))


def det(v,w):
    return v[0]*w[1]-v[1]*w[0]


def run(t):
    coefficients=(Q(4),Q(0),Q(0),Q(7)+t,Q(2))
    vertices=set()
    for i,j in combinations(range(5),2):
        d=det(R[i],R[j])
        if not d:
            continue
        a,b=-coefficients[i],-coefficients[j]
        p=((a*R[j][1]-b*R[i][1])/d,
           (R[i][0]*b-R[j][0]*a)/d)
        if all(sum(x*y for x,y in zip(p,r))>=-c for r,c in zip(R,coefficients)):
            vertices.add(p)
    lengths=[]
    active_pairs=set()
    for p in vertices:
        active=tuple(i for i,(r,c) in enumerate(zip(R,coefficients))
                     if sum(x*y for x,y in zip(p,r))==-c)
        assert len(active)==2
        active_pairs.add(frozenset(active))
    assert active_pairs=={frozenset(x) for x in CONES}
    assert all(det(R[i],R[j])==1 for i,j in CONES)
    for r,c in zip(R,coefficients):
        ends=[p for p in vertices if sum(x*y for x,y in zip(p,r))==-c]
        assert len(ends)==2
        tang=(r[1],-r[0]);k=next(i for i,x in enumerate(tang) if x)
        lengths.append(abs((ends[1][k]-ends[0][k])/tang[k]))
    assert lengths==[1+t,7+t,2,4,2] and min(lengths)>0
    mu=sum(lengths)/2
    flats={tuple(j for j,w in enumerate(R) if det(r,w)==0) for r in R}
    slopes=[(ids,sum(lengths[j] for j in ids)) for ids in sorted(flats)]
    worst=max(x[1] for x in slopes)
    status='unstable' if worst>mu else 'strictly semistable' if worst==mu else 'stable'
    expected='unstable' if t>0 else 'strictly semistable' if t==0 else 'stable'
    assert status==expected
    tight=[set(ids) for ids,slope in slopes if slope==mu]
    cone_counts=[sum(i in {0,1} for i in f) for f in CONES]
    failures=[list(f) for f,k in zip(CONES,cone_counts) if k!=1]
    if t==0:
        assert tight==[{0,1}]
        assert failures==[[4,2]]
    return {'t':str(t),'vertices':[list(map(str,p)) for p in sorted(vertices)],
            'divisor_degrees':list(map(str,lengths)),'tangent_slope':str(mu),
            'all_flats':[{'members':list(ids),'slope':str(s)} for ids,s in slopes],
            'status':status,'ray_line_cone_counts':cone_counts,
            'condition15_fails_at_cone':failures if t==0 else None}


if __name__=='__main__':
    records=[run(Q(t)) for t in ('-1/2','-1/100','0','1/100','1','3')]
    out={'rays':R,'maximal_cones':CONES,'records':records,
         'scope':'Exact finite checks; uniform parameter proof in the manuscript',
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dest=Path(__file__).resolve().parents[1]/'data/surface_wall_checks.json'
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print('PASS:',len(records),'exact fan/polygon/full-flat checks')

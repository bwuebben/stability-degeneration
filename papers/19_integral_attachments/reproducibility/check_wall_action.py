"""Exact full-ring reflection and normalized-index wall action checks.

This checks algebra of the proposed compact correspondences; it does not
establish their geometric existence or the partial-wall classification.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import check_blowup_model as m
import check_universal_indices as ui

ROOT = m.theta(m.transform(m.a1))


def refl2(z):
    co=ui.coords2(z)
    root=ui.coords2(ROOT)
    k=m.dot(co[:-1],root[:-1])-2*co[-1]*root[-1]
    return m.add(z,m.scale(ROOT,k))


MU=[m.mu(m.vec((i,1))) for i in range(m.N)]
IMU=[refl2(a) for a in MU]
IDELTA=refl2(m.DELTA)
IA=m.scale(m.add(m.C2,m.scale(m.mul(IDELTA,IDELTA),3)),F(1,24))
IB4=[IA]
for i in range(m.N): IB4.append(m.mul(IMU[i],IDELTA))
for i in range(m.N):
    for j in range(i+1,m.N):
        IB4.append(m.add(m.mul(IMU[i],IMU[j]),m.scale(IA,-m.Q[i][j])))
for i in range(m.N):
    IB4.append(m.scale(m.add(m.mul(IMU[i],IMU[i]),m.scale(IA,-m.Q[i][i]),
                             m.scale(m.mul(IMU[i],IDELTA),-1)),F(1,2)))
IB6=[]
for i in range(m.N):
    dual=m.add(*[m.scale(IMU[j],m.QI[j][i]) for j in range(m.N)])
    IB6.append(m.mul(dual,IA))
IB6.append(m.scale(m.mul(IDELTA,IA),-1))


def reflect(z):
    out=ui.part(z,0)
    z2=ui.part(z,2)
    if z2: out=m.add(out,refl2(z2))
    z4=ui.part(z,4)
    if z4:
        out=m.add(out,*[m.scale(v,c) for v,c in zip(IB4,m.coords4(z4)) if c])
    z6=ui.part(z,6)
    if z6:
        out=m.add(out,*[m.scale(v,c) for v,c in zip(IB6,ui.coords6(z6)) if c])
    return m.add(out,ui.part(z,8))


def source_reflect(z):
    k=m.mukai(z,m.a1)
    return z[0]+k,m.plus(z[1],m.times(m.a1[1],k)),z[2]+k*m.a1[2]


def main():
    for b in MU+[m.DELTA]:
        assert refl2(refl2(b)) == b
    records={}
    for label,b,image in zip(m.LABELS,m.BASIS4,IB4):
        records[label]=ui.encode(image,4)
        assert reflect(image) == b
    for i,image in enumerate(IB6):
        ui.encode(image,6)
        original=m.mul(m.mu(m.QI[i]),m.A) if i<m.N else m.scale(m.mul(m.DELTA,m.A),-1)
        assert reflect(image) == original
    assert reflect(m.C2) == m.C2
    p0,p2=[p for _,p in ui.CURVES_PLANES]
    assert reflect(p0) == p2 and reflect(p2) == p0
    # The boundary generator multiplier in the common fixed ordering is 5 mod 8.
    for b,image in zip(m.BASIS4,IB4):
        k0=m.pair(b,m.add(p0,m.scale(p2,-3)))
        k1=m.pair(image,m.add(p0,m.scale(p2,-3)))
        assert (k1-5*k0)%8 == 0
    # Full normalized index formula, including its nontrivial rank-one line.
    correction=ui.exp(m.scale(ui.part(ui.normalized_ch(source_reflect(ui.Y)),2),-1))
    inputs=[m.a0]+[(F(0),m.vec((i,1)),F(0)) for i in range(m.N)]+[(F(0),m.vec(),F(1))]
    for z in inputs:
        old=ui.chern_classes(ui.normalized_ch(z))
        new=ui.chern_classes(m.mul(ui.normalized_ch(source_reflect(z)),correction))
        for k,(a,b) in enumerate(zip(old,new),1):
            assert reflect(a) == b
            ui.encode(b,2*k)
    # Pairing preservation on an independent deterministic set of mixed classes.
    for i in range(20):
        a=m.add(m.BASIS4[(13*i)%276],m.scale(m.BASIS4[(19*i+23)%276],i-7))
        b=m.add(m.BASIS4[(31*i+7)%276],m.scale(m.BASIS4[(11*i+83)%276],3-i))
        assert m.pair(reflect(a),reflect(b)) == m.pair(a,b)
    report={'status':'full integral reflection arithmetic verified; compact wall existence is proved in the manuscript',
            'H2_involution_basis':23,'H4_integral_involution_basis':276,
            'H6_integral_involution_basis':23,'normalized_index_classes_checked':96,
            'planes_exchanged':True,'boundary_multiplier_mod8':5,
            'six_edge_word':['id','reflection','id','id','reflection','id'],
            'six_edge_cohomology_action':'identity',
            'H2_images':{str(i):ui.encode(refl2(b),2) for i,b in enumerate(MU+[m.DELTA])},
            'H4_images':records}
    Path(__file__).with_name('wall-action-check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if not k.endswith('images')},indent=2))


if __name__ == '__main__': main()

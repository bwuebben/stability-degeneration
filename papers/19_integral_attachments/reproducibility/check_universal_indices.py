"""Calculate the 24 normalized index generators in the proposed integral model.

Independent checks compare their actual plane restrictions with Paper 17's
binomial formulas. The finite calculation does not prove geometric transport.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
import check_blowup_model as m

ONE = m.tensor(0,0)


def degree(k):
    def sd(i): return 0 if i == 0 else 4 if i == m.PT else 2
    return 2+sd(k[1]) if k[0] == 'e' else sd(k[1])+sd(k[2])


def part(v,d): return {k:x for k,x in v.items() if degree(k) == d}


def exp(v):
    out,t = ONE,ONE
    for k in range(1,5):
        t = m.mul(t,v)
        out = m.add(out,m.scale(t,F(1,factorial(k))))
    return out


EXP_MINUS_E = exp(m.scale(m.DELTA,-1))


def taut_ch(z):
    r,c,s = z
    left = m.add(m.scale(ONE,r),
                 *[m.scale(m.tensor(i+1,0),c[i]) for i in range(m.N)],
                 m.scale(m.tensor(m.PT,0),s-r))
    right = m.add(m.scale(ONE,r),
                  *[m.scale(m.tensor(0,i+1),c[i]) for i in range(m.N)],
                  m.scale(m.tensor(0,m.PT),s-r))
    return m.add(m.mul(left,EXP_MINUS_E),right)


def index_ch(z):
    r,c,s=m.transform(z)
    return m.add(m.scale(ONE,r+s),m.scale(taut_ch((r,m.times(c,-1),s)),-1))


Y = (F(0),m.vec((3,-1)),F(0))
NORMALIZER = exp(m.scale(part(index_ch(Y),2),-1))


def normalized_ch(z): return m.mul(index_ch(z),NORMALIZER)


def chern_classes(ch):
    cs=[ONE]
    for n in range(1,5):
        cs.append(m.scale(m.add(*[
            m.scale(m.mul(cs[n-k],part(ch,2*k)),(-1)**(k-1)*factorial(k))
            for k in range(1,n+1)]),F(1,n)))
    return cs[1:]


def coords2(v):
    c=[v.get(('t',i+1,0),F(0)) for i in range(m.N)]
    d=v.get(('e',0),F(0))
    assert m.add(m.mu(c),m.scale(m.DELTA,d)) == v
    return c+[d]


def coords6(v):
    a=[v.get(('t',i+1,m.PT),F(0)) for i in range(m.N)]
    d=-v.get(('e',m.PT),F(0))/2
    assert m.add(m.mul(m.mu(a),m.A),m.scale(m.mul(m.DELTA,m.A),-d)) == v
    return [sum(m.Q[i][j]*a[j] for j in range(m.N)) for i in range(m.N)]+[d]


def encode(v,d):
    if d == 2:
        labels=[f'mu{i}' for i in range(m.N)]+['delta']; vals=coords2(v)
    elif d == 4:
        labels=m.LABELS; vals=m.coords4(v)
    elif d == 6:
        labels=[f'R{i}' for i in range(m.N)]+['Rdelta']; vals=coords6(v)
    elif d == 8:
        labels=['point']; vals=[m.integral(v)]
        assert m.scale(m.tensor(m.PT,m.PT),2*vals[0]) == v
    else: raise ValueError(d)
    assert all(x.denominator == 1 for x in vals), (d,vals)
    return {k:int(x) for k,x in zip(labels,vals) if x}


def curve_and_plane(ai):
    c=(m.v[0]/2-ai[0],m.plus(m.times(m.v[1],F(1,2)),m.times(ai[1],-1)),m.v[2]/2-ai[2])
    gc=m.transform(c)
    lam=m.theta(gc)
    p=m.add(m.scale(m.mul(lam,lam),F(1,2)),m.scale(m.C2,F(1,24)))
    return gc,p


CURVES_PLANES=[curve_and_plane(a) for a in (m.a0,m.a2)]


def verify_index(z):
    cs=chern_classes(normalized_ch(z))
    for k,c in enumerate(cs,1): encode(c,2*k)
    c1=coords2(cs[0])
    for ai,(gc,p) in zip((m.a0,m.a2),CURVES_PLANES):
        chi=-m.mukai(ai,z)
        assert m.dot(gc[1],c1[:-1])-2*gc[0]*c1[-1] == chi
        assert m.pair(p,cs[1]) == chi*(chi-1)/2
    return cs


def main():
    basis=[('O',m.a0)]+[(f'alpha{i}',(F(0),m.vec((i,1)),F(0))) for i in range(m.N)]+[('point',(F(0),m.vec(),F(1)))]
    records=[]
    for label,z in basis:
        cs=verify_index(z)
        records.append({'class':label,'rank':int(-m.mukai(m.v,z)),
                        'chern':{f'c{k}':encode(c,2*k) for k,c in enumerate(cs,1)},
                        'plane_c2':[int(m.pair(p,cs[1])) for _,p in CURVES_PLANES]})
    # Nonlinear/negative-rank checks, independent of merely testing a K-basis.
    extra=[Y]
    for n in range(1,8):
        extra.append((F(n-4),m.vec((0,n),(2,1-n),(7,n%3-1)),F(3-2*n)))
    for z in extra: verify_index(z)
    assert part(normalized_ch(Y),2) == {}
    report={'status':'exact index arithmetic; topological correspondence is proved in the manuscript',
            'basis_indices':len(records),'additional_virtual_indices':len(extra),
            'all_chern_coordinates_integral':True,'plane_binomial_restrictions_pass':True,
            'rank_one_normalizer_c1_zero':True,'generators':records}
    Path(__file__).with_name('universal-index-check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='generators'},indent=2))


if __name__ == '__main__': main()

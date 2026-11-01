"""Exact formal obstruction probes using the endpoint tensor-Cech contraction.

This is an exploratory probe. Nonzero higher coefficients in this particular
choice of coordinates do not by themselves disprove analytic quadraticity.
"""
from functools import lru_cache
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sympy as S
import cech_model as C


def unit(i, j, ex):
    return {(i-1, j-1, ex): S.Integer(1)}


H1 = {}
for key in C.BASIS[1]:
    i, j, ex = key
    if (i,j) == (0,1) and ex[0] == -1:
        continue
    if key in [(0,1,(-2,0,3)), (1,3,(0,-1,1))]:
        continue
    v = {key:S.Integer(1)}
    if (i,j)==(1,3) and ex[0]==1:
        v=C.add(v,unit(1,2,(-1,0,ex[2]+2)))
        # The closed lift needs the negative a0 coefficient.
        v[(0,1,(-1,0,ex[2]+2))]=-1
    name=f'{i+1}{j+1}_'+('_'.join(str(x).replace('-','m') for x in ex))
    H1[name]=v

assert len(H1)==77
assert all(not C.transferred_d(v) for v in H1.values())


def h_linear(v):
    out={}
    for (i,j,ex),coef in v.items():
        if (i,j,ex)==(0,1,(-2,0,3)):
            out=C.add(out,C.scale(unit(1,1,(0,0,0)),-coef))
        elif (i,j,ex)==(1,3,(0,-1,1)):
            out=C.add(out,C.scale(unit(4,4,(0,0,0)),coef))
        elif (i,j)==(0,3) and ex[2]>=1:
            out=C.add(out,C.scale(unit(1,2,(-1,0,ex[2]-1)),coef))
    return out


def k_twisted(a):
    out={};term=a
    for _ in range(8):
        out=C.add(out,C.homotopy(term))
        term=C.scale(C.perturb(C.homotopy(term)),-1)
        if not term:
            return out
    raise AssertionError('Nonterminating homotopy')


def k_total(a):
    return C.add(k_twisted(a),C.lift(h_linear(C.transferred_project(a))))


def obstruction(a):
    return C.transferred_project(a).get((0,3,(-1,-1,0)),S.Integer(0))


@lru_cache(None)
def representative(name):
    return C.lift(H1[name])


def qmatrix():
    names=list(H1);mat=S.zeros(77,77)
    for i,a in enumerate(names):
        for j in range(i,77):
            b=names[j]
            value=obstruction(C.add(C.multiply(representative(a),representative(b)),C.multiply(representative(b),representative(a))))
            mat[i,j]=mat[j,i]=value
    return names,mat


def series(initial,order=6):
    aa={1:initial}; obs={}
    for n in range(2,order+1):
        curvature=C.add(*(C.multiply(aa[i],aa[n-i]) for i in range(1,n)))
        obs[n]=S.expand(obstruction(curvature))
        aa[n]=C.scale(k_total(curvature),-1)
    return obs,aa


def main():
    names,mat=qmatrix()
    terms=[{'a':names[i],'b':names[j],'polar_coefficient':str(mat[i,j])}
           for i in range(77) for j in range(i,77) if mat[i,j]]
    assert mat.rank()==6
    # Exact identities for the second contraction, beyond dimension counts.
    for name in names:
        a=representative(name)
        assert not C.twisted_d(a)
        assert not k_total(a)
        assert not k_total(k_total(C.multiply(a,a)))
    x,y,r,b,c=S.symbols('x y r b c')
    selected=['12_m2_0_0','24_1_m1_0','12_m2_0_1',
              '24_0_m1_0','21_3_0_m2']
    initial=C.add(*(C.scale(representative(n),v)
                    for n,v in zip(selected,[x,y,r,b,c])))
    obs,_=series(initial,4)
    assert S.expand(obs[2]-x*y)==0
    assert S.expand(obs[3]+r*(b*y+c*x))==0
    assert S.expand(obs[4]+b**3*y-b*c*r**2)==0
    # Removing the nondegenerate pair cancels this restricted quartic term.
    xs=-S.diff(obs[3],y).subs({x:0,y:0})
    ys=-S.diff(obs[3],x).subs({x:0,y:0})
    residual=S.expand(obs[4].subs({x:0,y:0})-xs*ys)
    assert residual==0
    # The first proposed middle-arrow cubic does not give a nonzero term.
    zero_names=['13_0_m2_0','32_m1_2_0','24_0_m1_0']
    zero_obs,_=series(C.add(*(representative(n) for n in zero_names)),6)
    assert not any(zero_obs.values())
    result={'status':'formal probes; no assertion of analytic quadraticity',
            'basis_names':names,'quadratic_polar_rank':mat.rank(),
            'quadratic_polar_terms':terms,
            'mixed_probe':{'directions':selected,'variables':['x','y','r','b','c'],
                           'coefficients':{str(n):str(S.factor(v)) for n,v in obs.items()},
                           'critical_pair_second_order':[str(xs),str(ys)],
                           'restricted_morse_quartic':str(residual)},
            'unsuccessful_cubic_probe':{'directions':zero_names,
                'coefficients':{str(n):str(v) for n,v in zero_obs.items()}}}
    (Path(__file__).resolve().parents[1]/'output'/'obstruction_probe_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='basis_names'},indent=2))


if __name__=='__main__':
    main()

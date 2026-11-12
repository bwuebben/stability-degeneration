#!/usr/bin/env python3
"""Exact finite two-cycle spectrum and its plot; no PDE or flow simulation."""
from pathlib import Path
import json
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/figures'
OUT.mkdir(parents=True,exist_ok=True)

# Independently form the energy on an orthonormal real trace-free Hermitian basis.
sqrt2=sp.sqrt(2)
E11=sp.diag(1,0);E22=sp.diag(0,1);Z=sp.zeros(2)
sx=sp.Matrix([[0,1],[1,0]])/sqrt2
sy=sp.Matrix([[0,-sp.I],[sp.I,0]])/sqrt2
basis=[((E11-E22)/2,(E11-E22)/2),
       (E11/sqrt2,-E11/sqrt2),(E22/sqrt2,-E22/sqrt2),
       (sx,Z),(sy,Z),(Z,sx),(Z,sy)]
inner=lambda X,Y:sp.trace(X.conjugate().T*Y)
M=sp.Matrix(7,7,lambda i,j:sp.simplify(sum(inner(x,y) for x,y in zip(basis[i],basis[j]))))
assert M==sp.eye(7)
d,e,t=sp.symbols('d e t',real=True)
D=sp.diag(d,e)
images=[(S*D-D*T,T*D-D*S) for S,T in basis]
Q=sp.Matrix(7,7,lambda i,j:sp.simplify(2*sum(inner(x,y) for x,y in zip(images[i],images[j]))))
poly=sp.factor((t*sp.eye(7)-Q).det())
expected=t*(t-8*d*d)*(t-8*e*e)*(t-2*(d-e)**2)**2*(t-2*(d+e)**2)**2
assert sp.simplify(poly-expected)==0
qfun=sp.lambdify((d,e),Q,'numpy')
ratios=np.unique(np.r_[np.linspace(0,2,801),1.0])
maxerr=0.0
for r in ratios:
 computed=np.linalg.eigvalsh(np.asarray(qfun(1.,float(r)),dtype=float))
 exact=np.sort([0.,8.,8*r*r,2*(1-r)**2,2*(1-r)**2,2*(1+r)**2,2*(1+r)**2])
 maxerr=max(maxerr,float(np.max(np.abs(computed-exact))/max(1.,np.max(exact))))
assert maxerr<1e-12

plt.rcParams.update({'font.family':'serif','font.size':9,'axes.titlesize':10,
                     'axes.labelsize':10,'legend.fontsize':8,'pdf.fonttype':42})
fig,ax=plt.subplots(1,2,figsize=(7.2,3.05),constrained_layout=True)
s=np.linspace(0,2,801)
colors=['#555555','#997c00','#1b8a68','#345caa','#aa4b78']
ax[0].plot(s,np.zeros_like(s),color=colors[0],linestyle=':',label=r'$0$')
ax[0].plot(s,np.full_like(s,8),color=colors[1],label=r'$8$')
ax[0].plot(s,8*s*s,color=colors[2],label=r'$8s^2$')
ax[0].plot(s,2*(1-s)**2,color=colors[3],label=r'$2(1-s)^2$ (twice)')
ax[0].plot(s,2*(1+s)**2,color=colors[4],label=r'$2(1+s)^2$ (twice)')
ax[0].axvline(1,color='#bbbbbb',linewidth=.6,zorder=0)
ax[0].set(xlim=(0,2),ylim=(-.6,33),xlabel=r'$s=e/d$',ylabel=r'$V\mu/d^2$',title='(a) Complete finite spectrum')
ax[0].legend(loc='upper left',frameon=False,ncol=1,handlelength=1.8)

h=np.geomspace(1e-3,.4,501)
fast1=8*h*h;fast2=8*(h+h**3)**2;fast3=2*(2*h+h**3)**2
ax[1].loglog(h,fast1,color=colors[2],linewidth=1.2,label='Four fast modes')
ax[1].loglog(h,fast2,color=colors[2],linewidth=.8,linestyle='--')
ax[1].loglog(h,fast3,color=colors[2],linewidth=.8,linestyle=':')
# Evaluate the exact difference h^3 without floating subtraction.
ax[1].loglog(h,2*h**6,color=colors[3],linewidth=1.5,label='Two slow modes')
ax[1].text(.03,3e-4,r'$V\mu\sim8h^2$',color=colors[2])
ax[1].text(.012,2e-14,r'$V\mu=2h^6$',color=colors[3])
ax[1].text(.05,5e-18,'One exact zero',color=colors[0],fontsize=8)
ax[1].set(xlim=(1e-3,.4),ylim=(1e-18,3),xlabel=r'$h$',ylabel=r'$V\mu$',title=r'(b) $d=h,\ e=h+h^3$')
ax[1].legend(loc='upper left',frameon=False)
for a in ax:
 a.spines[['top','right']].set_visible(False)
 a.tick_params(direction='out')
 a.grid(axis='y',color='#dddddd',linewidth=.45,zorder=0)
fig.savefig(OUT/'two_cycle_spectrum.pdf',metadata={'Title':'Finite two-cycle spectrum','Author':'Bernd Johannes Wuebben','CreationDate':None,'ModDate':None})
fig.savefig(OUT/'two_cycle_spectrum.png',dpi=220)
plt.close(fig)
result={'scope':'Exact finite energy/mass calculation and illustration; no actual metric or nonlinear-flow simulation',
        'dimension':7,'basis_gram':'identity','characteristic_polynomial':str(poly),
        'normalization':'Eigenvalues shown are V times the generalized eigenvalues',
        'sampled_ratios':len(ratios),'maximum_error_scaled_by_largest_root':maxerr,
        'zero_multiplicities':{'e=0<d':2,'e=d>0':3,'generic_positive_distinct':1,'d=e=0':7},
        'result':'PASS'}
(OUT/'two_cycle_spectrum_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

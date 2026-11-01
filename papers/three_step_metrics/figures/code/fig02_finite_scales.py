#!/usr/bin/env python3
"""Finite graph and leading response only; no numerical HE solution."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
D=np.array([[1,1,0,0],[-1,0,1,0],[0,-1,0,1],[0,0,-1,-1]],dtype=float)
B=np.array([[5,5],[3,-5],[-5,3],[-3,-3]],dtype=float)/8

def compute_data():
    rho=np.logspace(0,6,241)
    # Stable positive root of (rho-1) tau^2 + 8/5 tau - 3/5=0.
    tau=(6/5)/(8/5+np.sqrt((8/5)**2+12/5*(rho-1)))
    ev=[];sv=[]
    for t in tau:
        w=np.array([1-t,t,3/5-t,t])
        L=D@np.diag(w)@D.T
        ev.append(np.linalg.eigvalsh(L)[1:])
        A=np.linalg.solve(L+np.ones((4,4))/4,B)
        assert np.max(np.abs(L@A-B))<1e-9
        sv.append(np.linalg.svd(A,compute_uv=False))
    return rho,tau,np.array(ev),np.array(sv)

def make_figure(data):
    rho,tau,ev,sv=data
    plt.rcParams.update({'font.family':'serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axs=plt.subplots(1,2,figsize=(9.0,3.4),layout='constrained')
    colors=['#2457a7','#d37420','#44815c']
    for j in range(3):axs[0].loglog(rho,ev[:,j],color=colors[j],label=rf'$\nu_{j+1}/l$')
    axs[0].set(title='Small graph spectrum',xlabel=r'Orbit invariant $\rho$',ylabel='Eigenvalue / $l$')
    axs[0].legend(frameon=False,loc='lower left')
    for j in range(2):axs[1].loglog(rho,sv[:,j],color=colors[j],label=rf'$s_{j+1}(A_0)$')
    axs[1].loglog(rho,np.sqrt(102)/(32*tau),'--',color='#555555',lw=1,label=r'$\sqrt{102}/(32\tau)$')
    axs[1].axhline(5/np.sqrt(17),ls=':',color='#555555',lw=1,label=r'$5/\sqrt{17}$')
    axs[1].set(title='Two polarization-response scales',xlabel=r'Orbit invariant $\rho$',ylabel='Leading response × $l$')
    axs[1].legend(frameon=False,loc='upper left',fontsize=8)
    for ax in axs:ax.grid(True,which='major',alpha=.18)
    out=ROOT/'figures/output';out.mkdir(exist_ok=True)
    fig.savefig(out/'fig02_finite_scales.pdf',metadata={'Title':'Finite spectral and response scales','Author':'Bernd Johannes Wuebben'})
    fig.savefig(out/'fig02_finite_scales.png',dpi=180)
    np.savetxt(out/'fig02_finite_scales.csv',np.column_stack([rho,tau,ev,sv]),delimiter=',',header='rho,tau,nu1_over_l,nu2_over_l,nu3_over_l,leading_response1,leading_response2',comments='')
    plt.close(fig)

if __name__=='__main__':
    make_figure(compute_data())
    print('Generated finite-model figure and all 241 plotted rows.')

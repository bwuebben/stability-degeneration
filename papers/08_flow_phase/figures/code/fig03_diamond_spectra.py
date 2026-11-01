from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[2]
out=root/"figures"; data=out/"data"; data.mkdir(exist_ok=True)
D=np.array([[1,1,0,0],[-1,0,1,0],[0,-1,0,1],[0,0,-1,-1]],float)
Q=np.linalg.qr(np.array([[1,1,1],[-1,1,1],[0,-2,1],[0,0,-3]],float))[0]
ratios=np.geomspace(1e-7,1,501)
eig=[]
for k in ratios:
    w=np.array([5,5*k,3*k,3])/8
    vals=np.linalg.eigvalsh(Q.T @ (D @ np.diag(w) @ D.T) @ Q)
    assert np.all(vals>0)
    assert abs(vals.sum()-2*(k+1))<1e-12
    eig.append(vals)
eig=np.asarray(eig)
norm=eig/np.column_stack([ratios,np.ones_like(ratios)*.75,np.ones_like(ratios)*1.25])
assert abs(norm[0,0]-1)<1e-6
assert np.max(abs(norm[0,1:]-1))<1e-6
np.savetxt(data/"diamond_spectra.csv",
           np.column_stack([ratios,eig,norm]),delimiter=",",
           header="k_over_l,nu1_over_l,nu2_over_l,nu3_over_l,nu1_over_k,nu2_over_075l,nu3_over_125l",comments="")
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,
                     "axes.spines.top":False,"axes.spines.right":False})
fig,ax=plt.subplots(1,2,figsize=(9,3.2),layout="constrained")
colors=["#0069aa","#c36b18","#27764e"]
for j in range(3):
    ax[0].loglog(ratios,eig[:,j],color=colors[j],label=rf"$\nu_{j+1}/l$")
    ax[1].semilogx(ratios,norm[:,j],color=colors[j],
                  label=[r"$\nu_1/k$",r"$\nu_2/(3l/4)$",r"$\nu_3/(5l/4)$"][j])
ax[0].set(xlabel=r"$k/l$",ylabel="Finite eigenvalue / l",title="One weak mode, two larger modes")
ax[1].set(xlabel=r"$k/l$",ylabel="Ratio to leading asymptotic",title="Exact asymptotic coefficients")
for aa in ax:aa.grid(alpha=.18);aa.legend(frameon=False,fontsize=8)
fig.savefig(out/"fig03_diamond_spectra.pdf")
fig.savefig(out/"fig03_diamond_spectra.png",dpi=180)
print("Saved diamond figure and 501 finite spectra; trace and endpoint checks PASS.")

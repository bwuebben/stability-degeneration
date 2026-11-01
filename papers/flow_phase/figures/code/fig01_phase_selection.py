from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[2]
out=root/"figures"; data=out/"data"; data.mkdir(exist_ok=True)
eps,s,a,y0,v0=.08,1.,.3,.4,.5
t=np.linspace(0,35,701)
phase=y0+a*v0/(s-eps)
selected=phase*np.exp(-eps*t)
fast=a*v0/(s-eps)*np.exp(-s*t)
actual=selected-fast
projected=y0*np.exp(-eps*t)
slow=a*v0/(s-eps)*(np.exp(-eps*t)-np.exp(-s*t))
np.savetxt(data/"phase_selection.csv",
           np.column_stack([t,actual,selected,projected,fast,slow]),
           delimiter=",",header="time,actual_y,selected_y,projected_y,selected_abs_error,projected_abs_error",comments="")
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,
                     "axes.spines.top":False,"axes.spines.right":False})
fig,ax=plt.subplots(1,2,figsize=(9,3.2),layout="constrained")
ax[0].plot(t,actual,color="#222222",lw=2,label="Actual finite coordinate")
ax[0].plot(t,selected,color="#0069aa",ls="--",label="Selected finite trajectory")
ax[0].plot(t,projected,color="#c36b18",ls=":",lw=2,label="Projected initial value")
ax[0].set(xlabel="Time",ylabel="Finite coordinate",title="Two choices on the invariant line")
ax[0].legend(fontsize=7.4,frameon=False)
ax[1].semilogy(t,fast,color="#0069aa",label="Selected: fast error")
ax[1].semilogy(t[1:],slow[1:],color="#c36b18",ls=":",lw=2,label="Projected: slow tail")
ax[1].set(xlabel="Time",ylabel="Absolute finite-coordinate error",title="Normal rate versus finite rate")
ax[1].legend(fontsize=8,frameon=False,loc="lower left")
for aa in ax:aa.grid(alpha=.18)
fig.savefig(out/"fig01_phase_selection.pdf")
fig.savefig(out/"fig01_phase_selection.png",dpi=180)
print("Saved phase figure and 701 exact-solution rows.")

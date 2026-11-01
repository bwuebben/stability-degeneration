"""Integral upper bounds for the finite nonlinear drift, not measured residuals."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root=Path(__file__).resolve().parents[1]
(root/"data").mkdir(exist_ok=True)
plt.rcParams.update({"font.family":"serif","font.size":10,"axes.labelsize":11,
                     "pdf.fonttype":42,"ps.fonttype":42})
u=np.geomspace(1e-3,1e5,701)
fig,ax=plt.subplots(figsize=(7.8,3.7))
rows=[]
colors=["#555555","#18577d","#238c87","#b66a21"]
for theta,color in zip([0,.001,.03,1.],colors):
    y=5*np.log1p(u/5) if theta==0 else 5*np.log1p(-np.expm1(-.4*theta*u)/(2*theta))
    ax.semilogx(u,y,color=color,ls="--" if theta==0 else "-",
                label="Wall" if theta==0 else rf"$\tau/q_0={theta:g}$")
    rows.extend(zip([theta]*len(u),u,y))
ax.set(xlabel=r"Scaled time $q_0T$",ylabel=r"Upper bound for $\int_0^T q(t)\,dt$",
       xlim=(u[0],u[-1]),ylim=(0,53))
ax.grid(True,alpha=.2);ax.legend(loc="upper left",frameon=False)
fig.tight_layout()
for ext in ["pdf","png"]:fig.savefig(root/f"fig03_residual_bound.{ext}",dpi=180)
with (root/"data/integrated_energy_bounds.csv").open("w",newline="") as f:
    w=csv.writer(f);w.writerow(["tau_over_q0","q0_time","integral_upper_bound"]);w.writerows(rows)
print("Wrote integral-bound figure and 2804 exact-profile evaluations.")

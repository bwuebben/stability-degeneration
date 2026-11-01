"""Exact comparison profiles; this is not a numerical HYM solution."""
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
u=np.geomspace(1e-4,1e3,601)
def profile(a,d,theta):
    return np.exp(-a*theta*u)/(1+d/(a*theta)*(-np.expm1(-a*theta*u)))
fig,axes=plt.subplots(1,2,figsize=(8.6,3.5),sharey=True)
rows=[]
for ax,theta in zip(axes,[.01,1.]):
    lo=profile(14/5,5,theta)
    hi=profile(2/5,1/5,theta)
    wall=1/(1+u/5)
    ax.fill_between(u,lo,hi,color="#347cad",alpha=.14)
    ax.loglog(u,hi,color="#18577d",label="Upper comparison")
    ax.loglog(u,lo,color="#b66a21",label="Lower comparison")
    ax.loglog(u,wall,ls=":",color="#777777",label="Wall upper bound")
    ax.axvline(1/theta,ls="--",lw=.8,color="#444444",alpha=.55)
    ax.set(title=rf"$\tau/q_0={theta:g}$",xlabel=r"Scaled time $q_0t$",
           xlim=(1e-4,1e3),ylim=(1e-10,1.5))
    ax.grid(True,which="major",alpha=.2)
    rows.extend(zip([theta]*len(u),u,lo,hi,wall))
axes[0].set_ylabel(r"Energy divided by $q_0$")
axes[0].legend(loc="lower left",frameon=False,fontsize=8)
fig.tight_layout()
for ext in ["pdf","png"]:fig.savefig(root/f"fig02_transition.{ext}",dpi=180)
with (root/"data/transition_profiles.csv").open("w",newline="") as f:
    w=csv.writer(f);w.writerow(["tau_over_q0","q0_time","lower_over_q0",
                                "upper_over_q0","wall_upper_over_q0"])
    w.writerows(rows)
print("Wrote comparison figure and 1202 exact-profile evaluations.")

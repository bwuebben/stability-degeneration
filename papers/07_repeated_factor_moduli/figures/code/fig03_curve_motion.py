"""Finite quotient and assigned parameter speeds of the cubic rotation."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from finite_models import curve_data

out = Path(__file__).resolve().parents[1]/"output"
out.mkdir(exist_ok=True)
t = .1
eps = np.geomspace(1e-4,.45,701)
data = np.array([curve_data(t,e)[:3] for e in eps])
bundle, parameter, lower = np.sqrt(data).T/t
assert np.all(bundle+1e-12 >= lower)
with (out/"curve_motion.csv").open("w",newline="") as stream:
    writer=csv.writer(stream)
    writer.writerow(["epsilon","t","finite_bundle_speed","parameter_speed",
                     "bundle_lower_bound"])
    writer.writerows(np.column_stack((eps,np.full(len(eps),t),
                                      t*bundle,t*parameter,t*lower)))
plt.rcParams.update({"font.family":"serif","font.size":9,
                     "axes.spines.top":False,"axes.spines.right":False,
                     "pdf.fonttype":42})
fig,ax=plt.subplots(1,2,figsize=(7.2,2.7),layout="constrained")
ax[0].plot(eps,bundle,color="#1666ac",lw=1.8,label="Bundle quotient speed")
ax[0].plot(eps,lower,color="#c76b21",lw=1.4,ls="--",
           label="Orthogonal-component bound")
ax[0].set(xlabel=r"$\varepsilon$",ylabel="Speed divided by t",
          title="Persistent bundle motion")
ax[0].legend(frameon=False,fontsize=8,loc="center left")
ax[1].loglog(eps,bundle,color="#1666ac",lw=1.8,label="Bundle")
ax[1].loglog(eps,parameter,color="#c76b21",lw=1.6,label="Polarization")
ax[1].set(xlabel=r"$\varepsilon$",ylabel="Speed divided by t",
          title="Bundle and parameter velocities")
ax[1].legend(frameon=False,loc="upper left")
for a in ax:a.grid(alpha=.15)
fig.savefig(out/"curve_motion.pdf")
fig.savefig(out/"curve_motion.png",dpi=200)
print("curve_motion: finite quotient data and PDF/PNG saved")

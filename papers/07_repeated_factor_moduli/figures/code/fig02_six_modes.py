"""Plot only the explicitly defined finite spectrum."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from finite_models import eigenvalues

out = Path(__file__).resolve().parents[1]/"output"
out.mkdir(exist_ok=True)
eps = np.geomspace(1e-5, .5, 1001)
values = np.array([eigenvalues(e) for e in eps])
linear = np.linspace(1e-5, .5, 701)
ratio = np.array([eigenvalues(e)[0]/(1.25*e*e) for e in linear])
with (out/"six_modes.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["epsilon"]+[f"theta{j}" for j in range(1,7)])
    writer.writerows(np.column_stack((eps, values)))
plt.rcParams.update({"font.family":"serif", "font.size":9,
                     "axes.spines.top":False,"axes.spines.right":False,
                     "pdf.fonttype":42})
fig, ax = plt.subplots(1,2,figsize=(7.2,2.7),layout="constrained")
for j in range(1,6):
    ax[0].loglog(eps, values[:,j],color="#687f91", lw=1.05,alpha=.8)
ax[0].loglog(eps, values[:,0],color="#1666ac",lw=1.8,label=r"$\theta_1$")
ax[0].set(xlabel=r"$\varepsilon$",ylabel="Finite eigenvalue",
          title="One vanishing mode, five positive limits")
ax[0].legend(frameon=False,loc="lower right")
ax[1].plot(linear,ratio,color="#1666ac",lw=1.8)
ax[1].axhline(1,color="#8a8a8a",lw=.7,ls=":")
ax[1].set(xlabel=r"$\varepsilon$",
          ylabel=r"$\theta_1/[(5/4)\varepsilon^2]$",
          title="Exact leading coefficient")
for a in ax:a.grid(alpha=.15)
fig.savefig(out/"six_modes.pdf")
fig.savefig(out/"six_modes.png",dpi=200)
print("six_modes: finite SVD data and PDF/PNG saved")

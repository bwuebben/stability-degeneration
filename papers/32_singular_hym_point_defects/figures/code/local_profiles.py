"""Plot the exact normalized inclusion and charge-one energy profiles."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = Path(__file__).resolve().parents[1]
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 11,
    "axes.labelsize": 12, "legend.fontsize": 10,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
r = np.linspace(0, 3.5, 701)
height = r**2 / (1 + r**2)
charge = r**4 * (r**2 + 3) / (1 + r**2)**3
fig, ax = plt.subplots(figsize=(5.9, 3.25), layout="constrained")
ax.plot(r, height, color="#195a91", linewidth=2,
        label=r"Inclusion: $R^2/(1+R^2)$")
ax.plot(r, charge, color="#a04424", linewidth=2, linestyle="--",
        label=r"Curvature energy: $R^4(R^2+3)/(1+R^2)^3$")
ax.axhline(.5, color=".6", linewidth=.8, linestyle=":")
ax.axvline(1, color=".6", linewidth=.8, linestyle=":")
ax.plot(1, .5, "o", color=".18", markersize=4)
ax.set(xlim=(0, 3.5), ylim=(0, 1.03),
       xlabel=r"Normalized radius $R=r/r_p(T)$", ylabel="Normalized value")
ax.set_xticks([0, 1, 2, 3])
ax.set_yticks([0, .25, .5, .75, 1])
ax.spines[["top", "right"]].set_visible(False)
ax.legend(loc="lower right", frameon=False)
fig.savefig(out / "local_profiles.pdf", metadata={
    "Title": "Exact local inclusion and curvature-energy profiles",
    "Author": "Bernd Johannes Wuebben", "Creator": "Matplotlib",
})
fig.savefig(out / "local_profiles.png", dpi=180)
plt.close(fig)
print("Wrote exact-model figure:", out / "local_profiles.pdf")

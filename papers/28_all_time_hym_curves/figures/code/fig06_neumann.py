"""A scalar first-angular-mode Neumann example, at rho=1.

The exact formula and integrals are independently checked by the accepted
dossier computation check_collar_angular_modes.py. This illustration is
not a numerical solution or evidence for the global matrix PDE theorem.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUTPUT = Path(__file__).resolve().parents[1] / "output"


def compute_data():
    radius = np.linspace(1., 2., 181)
    theta = np.linspace(0., 2*np.pi, 361)
    r, t = np.meshgrid(radius, theta)
    v = r - 1
    chi = 6*v**5 - 15*v**4 + 10*v**3
    chi_prime = 30*v**4 - 60*v**3 + 30*v**2
    chi_second = 120*v**3 - 180*v**2 + 60*v
    potential = (chi*r - 4*r/3 - 4/(3*r))*np.cos(t)
    source = -(chi_second*r + 3*chi_prime)*np.cos(t)
    return r*np.cos(t), r*np.sin(t), source, potential


def make_figure(data):
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                         "font.size": 11, "pdf.fonttype": 42})
    x, y, source, potential = data
    fig, axes = plt.subplots(1, 2, figsize=(6.1, 3.0))
    for ax, values, title in zip(
            axes, [source, potential],
            [r"Source $S(r,\theta)$", r"Neumann solution $u(r,\theta)$"]):
        scale = np.max(np.abs(values))
        levels = np.linspace(-scale, scale, 31)
        ax.contourf(x, y, values, levels=levels, cmap="RdBu_r")
        for radius in [1, 2]:
            ax.add_patch(plt.Circle((0, 0), radius, fill=False,
                                   ec="#222222", lw=.9))
        ax.text(0, 0, r"$|z|\leq\rho$", ha="center", va="center", fontsize=11)
        ax.set(xlim=(-2.2, 2.2), ylim=(-2.2, 2.2), aspect="equal")
        ax.axis("off")
        ax.set_title(title, fontsize=12, pad=7)
    axes[0].text(0, -2.52, r"$\int_{A_\rho}S\,dx\,dy=0$",
                 ha="center", va="center", fontsize=12)
    axes[1].text(0, -2.52, "$\\partial_nu=0$\non both circles",
                 ha="center", va="center", fontsize=11)
    fig.subplots_adjust(left=.02, right=.98, bottom=.22, top=.84, wspace=.30)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT / "fig06_neumann.pdf", bbox_inches="tight",
                metadata={"Creator": "Matplotlib", "CreationDate": None})
    fig.savefig(OUTPUT / "fig06_neumann.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    make_figure(compute_data())

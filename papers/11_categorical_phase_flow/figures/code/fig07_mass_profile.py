"""Plot exact normalized static mass and Hessian profiles and their limits.

The rationalized expression is equation (rationalprofile) in the manuscript.
It avoids subtraction of nearly equal trace norms. No trajectory simulation
or fitted asymptotic value is used. The Hessian is evaluated by analytic
differentiation of the same formula, without finite differences.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]


def compute_data():
    xi = np.linspace(0, 2, 501)
    root = np.sqrt((1 + xi**2)**2 + 4)
    limit = xi**2 + root - np.sqrt(5)
    limit_hessian = 2 + (6*xi**2+2)/root - 4*xi**2*(1+xi**2)**2/root**3
    profiles = []
    for lam in (2.0, 2.5, 3.0):
        p = np.exp(-lam**2)
        rho = -np.expm1(-lam**2) / lam**2
        eta = p / rho**2
        w = xi**2 / (1 + eta * xi**2)
        w_prime = 2*xi / (1 + eta*xi**2)**2
        w_second = 2*(1-3*eta*xi**2) / (1 + eta*xi**2)**3
        h_root = np.sqrt(5 + 2*w + w**2)
        h = w + h_root
        h_prime = 1 + (1+w)/h_root
        h_second = 4/h_root**3
        g = np.sqrt(1 + 2*p*h + 5*p**2)
        profile = (2*h + 5*p) / (np.sqrt(1 + 2*p*h + 5*p**2) + 1) - np.sqrt(5)
        hessian = (h_second*w_prime**2+h_prime*w_second)/g
        hessian -= p*(h_prime*w_prime)**2/g**3
        profiles.append(dict(lam=lam, values=profile, hessian=hessian,
                             scale=np.sqrt(2*p)/rho))
    return xi, limit, limit_hessian, profiles


def make_figure(data):
    xi, limit, limit_hessian, profiles = data
    plt.rcParams.update({
        "font.family": "serif", "font.size": 10,
        "mathtext.fontset": "cm", "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42,
    })
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), layout="constrained")
    for entry, color, style in zip(profiles,
                                  ("#4477AA", "#EE7733", "#228833"),
                                  ("--", "-.", ":")):
        for ax, key in zip(axes, ("values", "hessian")):
            ax.plot(xi, entry[key], color=color, linestyle=style, lw=1.7,
                    label=rf"$\Lambda={entry['lam']:g}$")
    axes[0].plot(xi, limit, color="black", lw=1.4, label="Limit")
    axes[1].plot(xi, limit_hessian, color="black", lw=1.4, label="Limit")
    axes[0].set(title="(a) Mass profile", ylabel=r"$P_\Lambda(\xi)$",
           xlim=(0, 2), ylim=(0, 7.5))
    axes[1].set(title="(b) Normalized Hessian", ylabel=r"$P_\Lambda''(\xi)$",
                xlim=(0, 2))
    axes[1].scatter([0, 1], [2+2/np.sqrt(5), 2+3/np.sqrt(2)],
                    s=20, color="black", zorder=5)
    for ax in axes:
        ax.set_xlabel(r"$\xi=\epsilon/s_\Lambda$")
        ax.grid(alpha=0.15)
    axes[0].legend(loc="upper left", frameon=False, fontsize=9)
    out = ROOT / "figures" / "output"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "mass-profile.pdf")
    fig.savefig(out / "mass-profile.png", dpi=180)
    plt.close(fig)
    (out / "mass-profile-data.json").write_text(json.dumps({
        "xi": xi.tolist(), "limit": limit.tolist(),
        "limit_hessian": limit_hessian.tolist(),
        "profiles": [{"lambda": row["lam"], "scale": row["scale"],
                      "values": row["values"].tolist(),
                      "hessian": row["hessian"].tolist()} for row in profiles],
        "method": "Exact mass formula and analytic derivatives; not a trajectory simulation."
    }, indent=2) + "\n")


if __name__ == "__main__":
    make_figure(compute_data())

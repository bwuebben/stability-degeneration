"""Plot the exact normalized static mass profile and its proved limit.

The rationalized expression is equation (rationalprofile) in the manuscript.
It avoids subtraction of nearly equal trace norms. No trajectory simulation
or fitted asymptotic value is used.
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
    limit = xi**2 + np.sqrt((1 + xi**2)**2 + 4) - np.sqrt(5)
    profiles = []
    for lam in (2.0, 2.5, 3.0):
        p = np.exp(-lam**2)
        r = -np.expm1(-lam**2) / lam**2
        eta = p / r**2
        w = xi**2 / (1 + eta * xi**2)
        h = w + np.sqrt(5 + 2*w + w**2)
        profile = (2*h + 5*p) / (np.sqrt(1 + 2*p*h + 5*p**2) + 1) - np.sqrt(5)
        profiles.append(dict(lam=lam, values=profile, scale=np.sqrt(2*p)/r))
    return xi, limit, profiles


def make_figure(data):
    xi, limit, profiles = data
    plt.rcParams.update({
        "font.family": "serif", "font.size": 10,
        "mathtext.fontset": "cm", "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(5.8, 3.2), layout="constrained")
    for entry, color, style in zip(profiles,
                                  ("#4477AA", "#EE7733", "#228833"),
                                  ("--", "-.", ":")):
        ax.plot(xi, entry["values"], color=color, linestyle=style, lw=1.7,
                label=rf"$\Lambda={entry['lam']:g}$")
    ax.plot(xi, limit, color="black", lw=1.4, label=r"Limit $P(\xi)$")
    ax.set(xlabel=r"Rescaled surviving parameter $\xi=\epsilon/s_\Lambda$",
           ylabel=r"Normalized mass difference $P_\Lambda(\xi)$",
           xlim=(0, 2), ylim=(0, 7.5))
    ax.legend(loc="upper left", frameon=False)
    ax.grid(alpha=0.15)
    out = ROOT / "figures" / "output"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "mass-profile.pdf")
    fig.savefig(out / "mass-profile.png", dpi=180)
    plt.close(fig)
    (out / "mass-profile-data.json").write_text(json.dumps({
        "xi": xi.tolist(), "limit": limit.tolist(),
        "profiles": [{"lambda": row["lam"], "scale": row["scale"],
                      "values": row["values"].tolist()} for row in profiles]
    }, indent=2) + "\n")


if __name__ == "__main__":
    make_figure(compute_data())

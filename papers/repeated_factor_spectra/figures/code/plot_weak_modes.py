"""Reproducible finite-model illustration; no HE differential operator is solved."""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from verify_star_model import balanced_covariances, full_spectrum

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "weak_modes.json"


def compute_data():
    rows = []
    for delta in np.geomspace(1e-10, 1.0, 121):
        covs, _, residual = balanced_covariances([delta, 1., 1.])
        values = full_spectrum(covs)
        rows.append({"delta": float(delta), "eigenvalues": values.tolist(),
                     "balance_residual": residual})
    record = {"scope": "Finite commutator spectrum at homogeneous t=1; not HE PDE eigenvalues.",
              "numpy": np.__version__, "matplotlib": matplotlib.__version__, "rows": rows}
    DATA.write_text(json.dumps(record, indent=2) + "\n")
    return record


def make_figure(record):
    delta = np.array([x["delta"] for x in record["rows"]])
    values = np.array([x["eigenvalues"] for x in record["rows"]])
    plt.rcParams.update({"font.family": "serif", "font.size": 8,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelsize": 9, "legend.fontsize": 7})
    fig, ax = plt.subplots(1, 2, figsize=(6.1, 2.65), layout="constrained")
    colors = ["#c45e00", "#0072b2", "#009e73", "#56b4e9", "#6e519a", "#333333"]
    for j in range(6):
        ax[0].loglog(delta, values[:,j], color=colors[j], lw=1.7 if j==0 else 1.0,
                     label=rf"$\nu_{j+1}$")
    ax[0].set(xlabel=r"Weak demand $\delta$", ylabel=r"Finite eigenvalue $\nu_j$",
              title="Six complex modes", xlim=(1e-10, 1), ylim=(5e-11, 5))
    ax[0].legend(ncol=3, loc="lower right", frameon=False, handlelength=1.4,
                 columnspacing=0.8)
    ax[1].semilogx(delta, values[:,0]/delta, color=colors[0], lw=1.8)
    ax[1].axhline(1.25, color="#444444", lw=0.9, ls="--", label=r"$5/4$")
    ax[1].set(xlabel=r"Weak demand $\delta$", ylabel=r"Weak-mode coefficient $\nu_1/\delta$",
              title="Leading coefficient", xlim=(1e-10, 1), ylim=(0.5, 1.34))
    ax[1].legend(loc="lower left", frameon=False)
    for axis in ax:
        axis.set_xticks([1e-10, 1e-5, 1])
        axis.grid(which="major", color="#dddddd", lw=0.5, zorder=0)
    for suffix in ["pdf", "png"]:
        fig.savefig(ROOT / "output" / f"weak_modes.{suffix}", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args()
    record = compute_data() if args.recompute or not DATA.exists() else json.loads(DATA.read_text())
    make_figure(record)
    print(json.dumps({"data": str(DATA), "rows": len(record["rows"]),
                      "matplotlib": matplotlib.__version__,
                      "max_balance_residual": max(r["balance_residual"] for r in record["rows"])}))

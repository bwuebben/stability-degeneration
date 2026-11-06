# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy==2.3.2", "matplotlib==3.10.5"]
# ///
"""Plot the finite-model distance law of Section 3."""
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)


def distance(z, order):
    """F_1(b(z)), with a cancellation-free expression for b squared."""
    b = np.sqrt(2 * z * z / (np.sqrt(1 + 4 * z * z) + 1))
    nodes, weights = leggauss(order)
    u = b[:, None] * (nodes[None, :] + 1) / 2
    return b * (np.sqrt((1 + 2 * u * u) / (1 + u * u)) @ weights) / 2


z = np.logspace(-4, 4, 401)
profile = distance(z, 256)
comparison = distance(z, 512)
relative_check = float(np.max(np.abs(profile - comparison) / comparison))
if relative_check > 1e-9:
    raise RuntimeError(f"Quadrature refinement discrepancy: {relative_check}")

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "font.size": 10,
        "mathtext.fontset": "cm",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)
fig, ax = plt.subplots(figsize=(6.1, 3.45), layout="constrained")
ax.loglog(z, profile, color="#174f80", lw=2.0, label="Exact model profile")
small = z <= 1
large = z >= 1
ax.loglog(z[small], z[small], "--", color="#555555", lw=1.25,
          label=r"$s_2/\sqrt{\kappa}$ regime")
ax.loglog(z[large], np.sqrt(2 * z[large]), ":", color="#555555", lw=1.6,
          label=r"$\sqrt{2s_2}$ regime")
ax.axvline(1, color="#bcbcbc", lw=0.8, zorder=0)
ax.set_xlabel(r"$s_2/\kappa$, where $\kappa=|\tau|$")
ax.set_ylabel(r"$d_\tau(\,\cdot\,,D_\tau)/\sqrt{\kappa}$")
ax.set_title("Distance to the rank-one locus", fontsize=11)
ax.set_xlim(1e-4, 1e4)
ax.set_ylim(7e-5, 2e2)
ax.legend(frameon=False, loc="upper left", fontsize=9)
ax.grid(which="major", color="#e7e7e7", lw=0.6)
ax.tick_params(which="minor", length=2)
for ext in ("pdf", "svg", "png"):
    fig.savefig(OUT / f"distance-profile.{ext}", dpi=220)
plt.close(fig)

report = {
    "mathematical_source": "Section 3, distance to the rank-one locus, with kappa=1",
    "quantity": "F_1(sqrt((sqrt(1+4*z^2)-1)/2))",
    "sampling_interval": [1e-4, 1e4],
    "sample_count": len(z),
    "quadrature_orders": [256, 512],
    "maximum_relative_refinement_difference": relative_check,
    "numpy": np.__version__,
    "matplotlib": matplotlib.__version__,
    "status": "illustration of an analytic theorem, not proof by numerical fit",
}
(ROOT / "reproducibility/distance-profile.json").write_text(
    json.dumps(report, indent=2) + "\n"
)
print(json.dumps(report, indent=2))

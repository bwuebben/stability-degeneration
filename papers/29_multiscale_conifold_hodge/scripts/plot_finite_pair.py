"""Plot the exact two-dimensional finite pair; no Ricci-flat simulation."""
from fractions import Fraction
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Exact coefficients of 5 det(diag(w_f,w_g) - mu M).
m11, m12 = Fraction(36, 5), Fraction(-4, 5)
assert 5 * (m11**2 - m12**2) == 256
assert 5 * m11 == 36

M = np.array([[36.0, -4.0], [-4.0, 36.0]]) / 5
mass_values, mass_vectors = np.linalg.eigh(M)
inverse_root = (
    mass_vectors @ np.diag(mass_values**-0.5) @ mass_vectors.T
)
c = np.logspace(-4, 4, 801)
roots = np.empty((len(c), 2))
angles = np.empty_like(roots)
max_formula_error = 0.0
max_residual = 0.0
for index, ratio in enumerate(c):
    E = np.diag([ratio, 1.0])
    eigenvalues, orthogonal_vectors = np.linalg.eigh(inverse_root @ E @ inverse_root)
    vectors = inverse_root @ orthogonal_vectors
    for j in range(2):
        if vectors[0, j] < 0:
            vectors[:, j] *= -1
    # Evaluate the smaller root by the product to avoid subtraction loss.
    larger = (9 * (ratio + 1) + np.sqrt(81 * (ratio - 1)**2 + 4 * ratio)) / 128
    smaller = 5 * ratio / (256 * larger)
    formula = np.array([smaller, larger])
    max_formula_error = max(max_formula_error, float(np.max(abs(eigenvalues / formula - 1))))
    residual = E @ vectors - (M @ vectors) * eigenvalues
    denominator = np.linalg.norm(E) * np.linalg.norm(vectors)
    max_residual = max(max_residual, float(np.linalg.norm(residual) / denominator))
    assert np.allclose(vectors.T @ M @ vectors, np.eye(2), rtol=0, atol=2e-14)
    roots[index] = formula / (1 + ratio)
    angles[index] = np.rad2deg(np.arctan2(vectors[1], vectors[0]))

assert max_formula_error < 2e-11
assert max_residual < 2e-14
assert np.max(roots[:, 0] / roots[:, 1]) <= 0.8 + 1e-14

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "mathtext.fontset": "cm",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.65), layout="constrained")
colors = ["#a33c28", "#235a83"]
for j, label in enumerate([r"$\mu_-$", r"$\mu_+$"]):
    axes[0].loglog(c, roots[:, j], color=colors[j], label=label, linewidth=1.6)
    axes[1].semilogx(c, angles[:, j], color=colors[j], label=label, linewidth=1.6)
for ax in axes:
    ax.axvline(1, color="0.55", linestyle=":", linewidth=0.8)
    ax.set_xlim(c[0], c[-1])
    ax.set_xlabel(r"Capacity ratio $c=w_f/w_g$")
    ax.grid(axis="y", color="0.9", linewidth=0.6)
    ax.legend(frameon=False, fontsize=9)
axes[0].set_title("(a) Normalized finite eigenvalues", fontsize=9)
axes[0].set_ylabel(r"$\mu_\pm/(w_f+w_g)$")
axes[0].set_ylim(8e-6, 0.3)
axes[1].set_title("(b) Finite eigenline angles", fontsize=9)
axes[1].set_ylabel(r"Angle in $(f,g)$ coordinates (degrees)")
axes[1].set_yticks([-90, -45, 0, 45, 90])
axes[1].set_ylim(-95, 95)
fig.savefig(OUT / "finite_pair.pdf", metadata={"Title": "Finite conifold eigenvalues and eigenlines", "Author": "Bernd Johannes Wuebben", "Creator": "Matplotlib", "CreationDate": None, "ModDate": None})
fig.savefig(OUT / "finite_pair.png", dpi=220)
plt.close(fig)

evidence = {
    "scope": "Exact rational characteristic coefficients and numerical checks of the finite pair only",
    "samples": len(c),
    "capacity_ratio_range": [float(c[0]), float(c[-1])],
    "max_relative_formula_error": max_formula_error,
    "max_normalized_generalized_residual": max_residual,
    "max_eigenvalue_ratio": float(np.max(roots[:, 0] / roots[:, 1])),
    "result": "PASS",
}
(OUT / "finite_pair_check.json").write_text(json.dumps(evidence, indent=2) + "\n")
print(json.dumps(evidence, indent=2))

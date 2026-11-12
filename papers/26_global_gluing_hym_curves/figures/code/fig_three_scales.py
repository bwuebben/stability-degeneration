"""Draw the actual K4 attachment data; no numerical eigenvalues are plotted."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / "output"
OUT.mkdir(parents=True, exist_ok=True)

def make_figure():
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                         "font.size": 12, "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(6.4, 4.0), layout="constrained")
    pos = {0: (0, 0.1), 1: (0, 1.8), 2: (-1.8, -1.0), 3: (1.8, -1.0)}
    edges = [(0, 1, "#525252"), (0, 2, "#525252"), (0, 3, "#525252"),
             (1, 2, "#0072B2"), (1, 3, "#D55E00"), (2, 3, "#009E73")]
    for i, j, color in edges:
        ax.plot([pos[i][0], pos[j][0]], [pos[i][1], pos[j][1]],
                color=color, lw=1.8, zorder=1)
    for i, (x, y) in pos.items():
        ax.scatter([x], [y], s=390, facecolor="white", edgecolor="black",
                   linewidth=1.2, zorder=3)
        ax.text(x, y, str(i), ha="center", va="center", zorder=4)
    labels = [
        (0.14, .95, r"$h,\ I$", "left"),
        (-.85, -.29, r"$h,\ I$", "right"),
        (.85, -.29, r"$h,\ I$", "left"),
        (-1.18, .57, r"$h^2,\ R_x$", "right"),
        (1.18, .57, r"$h^4,\ R_y$", "left"),
        (0, -1.24, r"$h^3,\ I$", "center"),
    ]
    for x, y, label, align in labels:
        ax.text(x, y, label, ha=align, va="center",
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5})
    ax.set_aspect("equal")
    ax.set_xlim(-2.6, 2.6)
    ax.set_ylim(-1.55, 2.05)
    ax.axis("off")
    fig.savefig(OUT / "three-scales.pdf", bbox_inches="tight")
    fig.savefig(OUT / "three-scales.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    make_figure()

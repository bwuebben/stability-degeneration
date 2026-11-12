"""Draw the exact prism and K3,3 attachments; no analytic data are simulated."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / "output"


def make_figure():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                         "font.size": 12, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.4), layout="constrained")
    prism = [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5),
             (0, 3), (1, 4), (2, 5)]
    bipartite = [(i, j) for i in range(3) for j in range(3, 6)]
    ppos = {0: (-1, .7), 1: (1, .7), 2: (0, -1.05),
            3: (-.44, .34), 4: (.44, .34), 5: (0, -.42)}
    bpos = {0: (-.85, .8), 1: (-.85, 0), 2: (-.85, -.8),
            3: (.85, .8), 4: (.85, 0), 5: (.85, -.8)}
    for ax, edges, pos, title, group in [
        (axes[0], prism, ppos, "Triangular prism",
         r"$\mathbb{Z}/5\oplus\mathbb{Z}/15$"),
        (axes[1], bipartite, bpos, r"$K_{3,3}$",
         r"$(\mathbb{Z}/3)^2\oplus\mathbb{Z}/9$"),
    ]:
        for i, j in edges:
            ax.plot([pos[i][0], pos[j][0]], [pos[i][1], pos[j][1]],
                    color="#3b5266", lw=1.5, zorder=1)
        for i, (x, y) in pos.items():
            ax.scatter([x], [y], s=280, facecolor="white",
                       edgecolor="black", linewidth=1.0, zorder=3)
            ax.text(x, y, str(i), ha="center", va="center", zorder=4)
        ax.set_title(title, pad=10)
        ax.text(0, -1.43, group, ha="center", va="center", fontsize=14)
        ax.set_aspect("equal")
        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.6, 1.05)
        ax.axis("off")
    fig.savefig(OUT / "monodromy-graphs.pdf", bbox_inches="tight")
    fig.savefig(OUT / "monodromy-graphs.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    make_figure()

"""Draw the actual K_{2,2} attachment graph; no numerical data."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

OUTPUT = Path(__file__).resolve().parents[1] / "output"


def make_figure():
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                         "font.size": 11, "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(6.1, 3.7))
    positions = [(0, 2), (4, 2), (4, 0), (0, 0)]
    names = [r"$E_+$", r"$R_p$", r"$E_-$", r"$R_q$"]
    for a, b in [(0, 1), (1, 2), (2, 3), (3, 0)]:
        x0, y0 = positions[a]
        x1, y1 = positions[b]
        ax.plot([x0, x1], [y0, y1], color="#39434f", lw=1.4, zorder=1)
    for k, (x, y) in enumerate(positions):
        ax.add_patch(Circle((x, y), 0.29, fc="white",
                           ec="#273a4b", lw=1.5, zorder=3))
        ax.text(x, y, names[k], ha="center", va="center", zorder=4, fontsize=14)
    bbox = dict(fc="white", ec="none", pad=3)
    for x, y, label in [
        (2, 2, r"$w_1=2\pi/L_1$"),
        (4, 1, r"$w_2=2\pi/L_2$"),
        (2, 0, r"$w_3=2\pi/L_3$"),
        (0, 1, r"$w_4=2\pi/L_4$")]:
        ax.text(x, y, label, ha="center", va="center", bbox=bbox, zorder=5)
    ax.text(0, 2.51, r"$\beta_+=\alpha_+^2\kappa_+$",
            ha="center", color="#245b7a", fontsize=13)
    ax.text(4, -.53, r"$\beta_-=\alpha_-^2\kappa_-$",
            ha="center", color="#245b7a", fontsize=13)
    ax.text(2, 1.28, r"Each edge carries $U_e\in\mathrm{SU}(r)$",
            ha="center", va="center", fontsize=10)
    ax.text(2, .72, r"$u_ev_e=\tau_e,\qquad A_i=1$",
            ha="center", va="center", fontsize=12)
    ax.set(xlim=(-.95, 4.95), ylim=(-.85, 2.85), aspect="equal")
    ax.axis("off")
    fig.tight_layout(pad=.25)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT / "fig02_family.pdf", bbox_inches="tight",
                metadata={"Creator": "Matplotlib", "CreationDate": None})
    fig.savefig(OUTPUT / "fig02_family.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    make_figure()

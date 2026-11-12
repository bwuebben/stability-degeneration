"""Exact construction diagram; not a numerical experiment or PDE simulation."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

plt.rcParams.update({'font.family': 'serif', 'font.size': 10,
                     'mathtext.fontset': 'cm', 'pdf.fonttype': 42})
out = Path(__file__).resolve().parents[1] / 'output'
out.mkdir(parents=True, exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.45),
                         gridspec_kw={'width_ratios': [1, 1.32]})
blue, ink = '#235787', '#20262b'
for ax in axes:
    ax.set(xlim=(-.12, 1.12), ylim=(-.3, 1.2))
    ax.axis('off')
ax = axes[0]
ax.set_title('Four-cycle on the cover', fontsize=11, pad=8)
pts = [(0, .45), (.5, .95), (1, .45), (.5, -.05)]
for j, p in enumerate(pts):
    q = pts[(j+1) % 4]
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=10,
                                shrinkA=10, shrinkB=10, color=blue, lw=1.2))
    ax.plot(*p, 'o', color=ink, ms=5)
    dx, dy = [(-.08, 0), (0, .09), (.08, 0), (0, -.09)][j]
    ax.text(p[0]+dx, p[1]+dy, rf'$D_{j}$', ha='center', va='center')
    mx, my = (p[0]+q[0])/2, (p[1]+q[1])/2
    dx, dy = [(-.06, .05), (.06, .05), (.06, -.05), (-.06, -.05)][j]
    ax.text(mx+dx, my+dy, '$R$', color=blue, ha='center', va='center')
ax.plot([-.02, 1.02], [.45, .45], ls='--', lw=.8, color='#969da3')
ax.text(.5, .59, r'$U_\theta=R^4$', ha='center', fontsize=11)
ax.text(.5, .30, r'$D_j\longmapsto D_{-j}$', ha='center')
ax.text(.5, -.29, r'$A_3=A_1$', ha='center', va='bottom')
ax = axes[1]
ax.set_title('Invariant weak endomorphisms', fontsize=11, pad=8)
for x, label in [(.04, r'$a u$'), (.5, r'$X\in\mathbb{R}^3$'), (.96, r'$b v_\theta$')]:
    ax.plot(x, .68, 'o', color=ink, ms=5)
    ax.text(x, .82, label, ha='center', va='center', fontsize=11)
for lo, hi in [(.04, .5), (.5, .96)]:
    ax.plot([lo, hi], [.68, .68], lw=1.5, color=blue, zorder=0)
    ax.text((lo+hi)/2, .55, 'weight 2', ha='center', color=blue, fontsize=9)
for x, label in [(.04, '$A_0$'), (.5, '$2A_1$'), (.96, '$A_2$')]:
    ax.text(x, .35, label, ha='center', fontsize=11)
ax.text(.5, .15, r'$u=(1,0,0),\quad v_\theta=(\cos\theta,\sin\theta,0)$',
        ha='center', fontsize=10)
ax.text(.5, -.08, r'$q_\theta=2\left(\|au-X\|^2+\|X-bv_\theta\|^2\right)$',
        ha='center', fontsize=11)
ax.text(.5, -.28, 'Mass weights shown below the vertices', ha='center', fontsize=9)
fig.subplots_adjust(left=.04, right=.98, bottom=.12, top=.88, wspace=.22)
for suffix in ['pdf', 'png']:
    fig.savefig(out / f'gluing-constraints.{suffix}', dpi=180,
                metadata={'Title': 'Invariant gluing constraints'} if suffix=='pdf' else None)
print(out / 'gluing-constraints.pdf')

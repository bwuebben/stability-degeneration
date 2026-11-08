"""Exact finite-model curves, not numerical HE solutions or string data."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / 'output'

def compute_data():
    ratio = np.geomspace(1e-3, 1, 501)
    return ratio, (ratio**2, np.ones_like(ratio), 1+ratio**2), 1.5/ratio**2

def make_figure(data):
    ratio, masses, curvature = data
    plt.rcParams.update({'font.family':'serif', 'font.size':10,
                         'axes.spines.top':False, 'axes.spines.right':False,
                         'pdf.fonttype':42})
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), layout='constrained')
    for curve, label, color, style in zip(masses,
            [r'$b^2/a^2$', r'$1$', r'$1+b^2/a^2$'],
            ['#2166ac','#555555','#b2182b'], ['-','--','-']):
        ax[0].loglog(ratio, curve, label=label, color=color, ls=style, lw=1.8)
    ax[0].set(ylabel=r'$(3 V_0 L^2/8a^2)\,\widehat m_j^2$',
              xlabel=r'Angular ratio $b/a$', title='Leading vector spectrum',
              xlim=(1e-3,1), ylim=(5e-7,5))
    ax[0].legend(loc='lower right', frameon=False)
    ax[1].loglog(ratio, curvature, color='#2166ac', lw=1.8)
    ax[1].set(ylabel=r'$a^2 K_{g_{\mathcal{C}}}=3/[2(b/a)^2]$',
              xlabel=r'Angular ratio $b/a$', title='Exact cone curvature',
              xlim=(1e-3,1))
    for axis in ax:
        axis.grid(which='major', color='#dddddd', lw=.6)
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT/'finite_scales.pdf', bbox_inches='tight')
    fig.savefig(OUT/'finite_scales.png', dpi=160, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    make_figure(compute_data())

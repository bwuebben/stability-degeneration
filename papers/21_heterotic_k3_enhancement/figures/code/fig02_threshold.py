"""Light-vector exchange across its mass scale; exact finite-model formula.

No numerical HE solution, fitted error constant, or string remainder is used.
The plot uses t=-s/4, u=-3s/4, b/a=0.1, theta_H=0.8 and delta=0.1.
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / 'output'


def compute_data():
    delta = 0.1
    lower = (1-delta)/(1+delta)
    upper = (1+delta)/(1-delta)
    cutoff = 0.8/0.1**2
    branches = [np.geomspace(0.01, lower, 501),
                np.geomspace(upper, cutoff, 501)]
    values = [x/2/(1-x)-5*x/4/(1+3*x/4) for x in branches]
    return branches, values, lower, upper, cutoff


def make_figure(data):
    branches, values, lower, upper, cutoff = data
    plt.rcParams.update({'font.family': 'serif', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(7.5, 3.7), layout='constrained')
    ax.axvspan(lower, upper, color='#eeeeee', zorder=0)
    ax.axvline(1, color='#999999', ls=':', lw=1)
    for i, (x, value) in enumerate(zip(branches, values)):
        ax.semilogx(x, value, color='#2166ac', lw=2,
                    label='Light-vector exchange' if i == 0 else None)
    low = np.geomspace(0.01, 0.65, 301)
    ax.semilogx(low, -0.75*low, color='#b2182b', ls='--', lw=1.5,
                label=r'Low-energy term $3t/m_0^2$')
    ax.hlines(-13/6, 4, cutoff, color='#555555', ls='-.', lw=1.4,
              label=r'Intermediate limit $-13/6$')
    ax.axhline(0, color='#bbbbbb', lw=0.6)
    ax.text(1, 1.7, 'Pole exclusion', rotation=90, ha='center', va='center',
            color='#555555', fontsize=8.5)
    ax.set(xlabel=r'$s/m_0^2$ at $t=-s/4$, $u=-3s/4$',
           ylabel=r'$\mathcal{M}_0/e^2$', xlim=(0.01, cutoff), ylim=(-4.8, 3),
           title='Scalar scattering with the light vector retained')
    ax.grid(axis='x', which='major', color='#dddddd', lw=0.6)
    ax.legend(loc='lower left', frameon=False, fontsize=9)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT/'threshold_amplitude.pdf', bbox_inches='tight')
    fig.savefig(OUT/'threshold_amplitude.png', dpi=160, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    make_figure(compute_data())

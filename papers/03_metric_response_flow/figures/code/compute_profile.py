#!/usr/bin/env python3
"""Sample the proved signed scalar profile and write its native TikZ figure.

Standard library only. The points are values of the printed closed formula;
they are not numerical solutions of a bundle PDE.
"""
from pathlib import Path
import json
import math

ROOT = Path(__file__).resolve().parents[1]


def profile(delta, time):
    if delta == 0:
        return 1/(2*time)
    return delta/(-math.expm1(-2*delta*time))


def main():
    curves = []
    for delta, color in ((.5, 'red!75!black'), (0, 'black'), (-.5, 'blue!75!black')):
        points = [(1+4*j/120, profile(delta, 1+4*j/120)) for j in range(121)]
        assert all(y > 0 for x,y in points)
        assert all(a[1] > b[1] for a,b in zip(points,points[1:]))
        curves.append({'delta': delta, 'color': color, 'points': points})
    for time in (1,2,4):
        assert abs(profile(.5,time)-profile(-.5,time)-.5) < 1e-14
        assert abs(profile(1e-9,time)-profile(0,time)) < 1e-8
    lines = [r'\begin{figure}[htbp]',r'\centering',
             r'\begin{tikzpicture}[x=1.6cm,y=4.4cm,>=Stealth,every node/.style={font=\small}]',
             r'\draw[->] (.8,0)--(5.35,0) node[right] {\(\tau\)};',
             r'\draw[->] (1,0)--(1,.91) node[above] {\(P_\delta(\tau)\)};']
    for x in range(1,6):
        lines.append(f'\\draw ({x},.012)--({x},-.012) node[below] {{{x}}};')
    for y,label in ((.25,'1/4'),(.5,'1/2'),(.75,'3/4')):
        lines.append(f'\\draw (.97,{y})--(1.03,{y});\\node[left] at (.97,{y}) {{$\\tfrac{label[0]}{label[-1]}$}};')
    for curve in curves:
        points=' '.join(f'({x:.10f},{y:.10f})' for x,y in curve['points'])
        lines.append(f'\\draw[thick,{curve["color"]}] plot coordinates {{{points}}};')
    for y,curve,label in ((.84,curves[0],r'\delta=1/2'),(.73,curves[1],r'\delta=0'),(.62,curves[2],r'\delta=-1/2')):
        lines.append(f'\\draw[thick,{curve["color"]}] (3.7,{y})--(4.02,{y});\\node[right] at (4.08,{y}) {{$'+label+'$};')
    lines += [r'\end{tikzpicture}',
              r'\caption{The proved two-factor profile on the stable side, at the wall, and on the unstable side. These are evaluations of formula~\eqref{eq05:profile}, not measured HYM trajectories. Uniform relative comparison in Theorem~\ref{thm05:joint} includes the entire future interval after its late restart.}',
              r'\label{fig05:profile}',r'\end{figure}']
    (ROOT/'code/fig05_profile.tex').write_text('\n'.join(lines)+'\n')
    (ROOT/'data/profile_values.json').write_text(json.dumps({'formula':'delta/(1-exp(-2 delta time)), wall 1/(2 time)', 'curves':curves},indent=2)+'\n')
    print('PASS: scalar formula controls and native figure generated.')


if __name__ == '__main__':
    main()

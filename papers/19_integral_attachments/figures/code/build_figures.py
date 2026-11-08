"""Compile the manuscript's TikZ sources and render previews; no data generation."""
from pathlib import Path
import subprocess

here = Path(__file__).resolve().parent
out = here.parent / 'output'
out.mkdir(exist_ok=True)
preamble = r'''\documentclass[border=6pt]{standalone}
\usepackage{amsmath,amssymb,tikz,adjustbox}
\usetikzlibrary{arrows.meta,calc,positioning}
\newcommand{\Z}{\mathbb Z}
\newcommand{\PP}{\mathbb P}
\newcommand{\FF}{\mathbb F}
\newcommand{\id}{\operatorname{id}}
\begin{document}
\begin{adjustbox}{max width=440pt,center}
'''
for name in ('chambers', 'configurations', 'excision'):
    tex = out / (name + '.tex')
    tex.write_text(preamble + (here / (name + '.tex')).read_text()
                   + '\n\\end{adjustbox}\n\\end{document}\n')
    with (out / (name + '-build.log')).open('w') as log:
        subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
                        tex.name], cwd=out, stdout=log, stderr=subprocess.STDOUT,
                       check=True)
    subprocess.run(['pdftoppm', '-singlefile', '-r', '145', '-png',
                    str(out / (name + '.pdf')), str(out / name)], check=True)
    print(out / (name + '.png'))

#!/usr/bin/env python3
"""Reproduce the finite calculations accompanying the manuscripts."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
programs = sorted(p for p in (root / 'papers').glob('*/figures/code/*.py') if p.name != 'exact_linear.py')
for program in programs:
    print(str(program.relative_to(root)), flush=True)
    if program.parents[2].name == 'repeated_factor_spectra':
        command = ['uv', 'run', '--no-project', '--with-requirements', str(program.parents[2] / 'requirements.txt'), 'python', '-B', str(program)]
    else:
        command = [sys.executable, '-B', str(program)]
    subprocess.run(command, cwd=program.parents[2], check=True)
print(f'Completed {len(programs)} finite computations.')

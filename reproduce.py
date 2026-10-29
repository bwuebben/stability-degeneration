#!/usr/bin/env python3
"""Reproduce the finite calculations accompanying the three manuscripts."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
programs = sorted(p for p in (root / 'papers').glob('*/figures/code/*.py') if p.name != 'exact_linear.py')
for program in programs:
    print(str(program.relative_to(root)), flush=True)
    subprocess.run([sys.executable, '-B', str(program)], cwd=program.parents[2], check=True)
print(f'Completed {len(programs)} finite computations.')

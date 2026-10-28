#!/usr/bin/env python3
"""Reproduce the finite calculations accompanying the three manuscripts."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
programs = sorted((root / 'papers').glob('*/figures/code/*.py'))
for program in programs:
    print(str(program.relative_to(root)), flush=True)
    subprocess.run([sys.executable, '-B', str(program)], cwd=program.parents[2], check=True)
print(f'Completed {len(programs)} finite computations.')

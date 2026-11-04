#!/usr/bin/env python3
"""Reproduce the finite calculations accompanying the manuscripts."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
programs = sorted(p for p in (root / 'papers').glob('*/figures/code/*.py') if p.name not in {'exact_linear.py', 'cech_model.py', 'finite_models.py'})
for program in programs:
    print(str(program.relative_to(root)), flush=True)
    if (program.parents[2] / 'requirements.txt').exists():
        command = ['uv', 'run', '--no-project', '--with-requirements', str(program.parents[2] / 'requirements.txt'), 'python', '-B', str(program)]
    else:
        command = [sys.executable, '-B', str(program)]
    subprocess.run(command, cwd=program.parents[2], check=True)
examples = sorted((root / 'papers').glob('*/reproducibility/check_examples.py'))
for program in examples:
    print(str(program.relative_to(root)), flush=True)
    subprocess.run(['uv', 'run', '--no-project', str(program)], cwd=program.parents[1], check=True)
print(f'Completed {len(programs) + len(examples)} finite computations.')

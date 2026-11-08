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
k3_programs = sorted((root / 'papers').glob('*/reproducibility/k3_model.py')) + sorted((root / 'papers').glob('*/reproducibility/plot_distance_profile.py'))
for program in k3_programs:
    print(str(program.relative_to(root)), flush=True)
    subprocess.run(['uv', 'run', '--no-project', str(program)], cwd=program.parents[1], check=True)
wrappers = sorted((root / 'papers').glob('*/reproducibility/run_all.py'))
for program in wrappers:
    print(str(program.relative_to(root)), flush=True)
    subprocess.run([sys.executable, '-B', str(program)], cwd=program.parents[1], check=True)
print(f'Completed {len(programs) + len(examples) + len(k3_programs) + len(wrappers)} finite computations or suites.')

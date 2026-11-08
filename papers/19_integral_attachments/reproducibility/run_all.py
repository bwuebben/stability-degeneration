"""Run the finite exact checks from any working directory."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
for name, output in (
    ('check_blowup_model.py', 'blowup-run.log'),
    ('check_universal_indices.py', 'indices-run.log'),
    ('check_wall_action.py', 'wall-run.log'),
    ('check_all_attachments.py', 'all-attachments-check.json'),
):
    with (root / output).open('w') as stream:
        subprocess.run([sys.executable, str(root / name)], stdout=stream,
                       stderr=subprocess.STDOUT, check=True)
    print(name + ': PASS', flush=True)

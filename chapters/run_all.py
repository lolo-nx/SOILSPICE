"""The four chapters in the order of the paper.

    python chapters/run_all.py           # from the saved LTspice results (a few minutes)
    python chapters/run_all.py --full    # also rerun LTspice for chapters 1 and 2 (hours)

Chapters 3 and 4 use the equal-height sweep of chapter 2; chapter 1 uses it for the three model levels.
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for ch in ("1_validation.py", "2_simulation.py", "3_field_test.py", "4_economics.py"):
    print(f"\n========== {ch} ==========", flush=True)
    args = sys.argv[1:] if ch[0] in "12" else []
    if subprocess.run([sys.executable, str(HERE / ch), *args]).returncode:
        sys.exit(f"{ch} failed")

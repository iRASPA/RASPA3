#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
raspa = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("raspa3")
subprocess.check_call([str(raspa), "simulation.json"], cwd=here)

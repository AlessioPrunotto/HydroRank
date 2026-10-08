"""Seed an isolated cache-reuse trial before replacing this process with the adapter."""

import os
import shutil
import sys
from pathlib import Path

source = Path(sys.argv[1])
folder = Path.cwd() / "results"
folder.mkdir()
shutil.copy2(source, folder / "observations.npz")
os.execv(sys.argv[2], sys.argv[2:])

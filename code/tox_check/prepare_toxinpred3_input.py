#!/usr/bin/env python3
from pathlib import Path
import runpy


if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).parents[1] / "prepare_indexed_fasta.py"), run_name="__main__")

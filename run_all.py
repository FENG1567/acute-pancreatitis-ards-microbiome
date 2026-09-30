#!/usr/bin/env python3
"""Run the public analysis-ready reproduction and verification workflow."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def run(command: list[str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-figures", action="store_true", help="skip R figure regeneration")
    parser.add_argument("--require-figures", action="store_true", help="fail if Rscript is unavailable")
    args = parser.parse_args()
    if args.no_figures and args.require_figures:
        parser.error("--no-figures and --require-figures cannot be used together")

    run([sys.executable, "scripts/analysis/reproduce_primary.py"])
    run([sys.executable, "scripts/analysis/reproduce_secondary.py"])
    run([sys.executable, "scripts/analysis/verify_reproduction.py"])

    if not args.no_figures:
        rscript = shutil.which("Rscript")
        if rscript:
            run([rscript, "scripts/figures/make_figures.R", str(ROOT)])
        elif args.require_figures:
            raise SystemExit("Rscript not found; install the environment in environment.yml")
        else:
            print("Rscript not found: statistical reproduction passed; figure regeneration skipped.")

    run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"])
    print("Reproduction workflow completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python
"""Rebuild the workbook and figures from the distributed data and run all quality control."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))


def main():
    from analysis import tables

    status = 0
    print("1. supplementary datasets")
    loaded = tables.load_all()
    print(f"   {len(loaded)} datasets, {sum(len(v) for v in loaded.values()):,} rows")

    print("2. data dictionary")
    print(f"   {len(tables.data_dictionary())} fields described")

    print("3. consolidated workbook")
    tables.build_workbook(ROOT / "Supplementary_Tables.xlsx")
    verification = tables.verify_workbook(ROOT / "Supplementary_Tables.xlsx")
    mismatches = int(verification.mismatches.sum())
    print(f"   {int(verification.cells.sum()):,} cells written, {mismatches} mismatches")
    status |= 1 if mismatches else 0

    print("4. figures")
    status |= subprocess.call([sys.executable, str(ROOT / "reproduce_figures.py")])

    print("5. quality control")
    status |= subprocess.call([sys.executable, str(ROOT / "run_qc.py")])

    print("PASS" if status == 0 else "FAIL")
    return status


if __name__ == "__main__":
    sys.exit(main())

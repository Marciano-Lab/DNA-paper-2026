#!/usr/bin/env python
"""Run every quality-control check and write the report."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))
REPORT = ROOT / "audit" / "QC_REPORT.txt"


def main():
    from analysis import claims, qc, tables

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    sections = [("structural quality control", qc.structural_qc(), "structural_qc.tsv"),
                ("data dictionary", qc.dictionary_qc(), "dictionary_qc.tsv"),
                ("semantic quality control", qc.semantic_qc(), "semantic_qc.tsv"),
                ("retired references", qc.reference_qc(ROOT), "reference_qc.tsv"),
                ("manuscript claims", claims.recompute(), "claim_audit.tsv")]
    geometry = ROOT / "audit" / "figure_geometry.csv"
    if geometry.exists():
        sections.append(("figure geometry", qc.figure_qc(geometry), "figure_qc.tsv"))
    workbook = ROOT / "Supplementary_Tables.xlsx"
    if workbook.exists():
        verification = tables.verify_workbook(workbook)
        verification.to_csv(ROOT / "audit" / "workbook_verification.csv", index=False)
        sections.append(("workbook verification", pd.DataFrame([dict(
            test="every workbook cell matches its canonical table",
            status="PASS" if int(verification.mismatches.sum()) == 0 else "FAIL",
            expected=f"{int(verification.cells.sum())} cells",
            observed=f"{int(verification.mismatches.sum())} mismatches")]),
            "workbook_qc.tsv"))
    manifest = ROOT / "release_manifest.tsv"
    if manifest.exists():
        sections.append(("release checksums", qc.checksum_qc(manifest), "checksum_qc.tsv"))

    lines, failures = [], 0
    for title, frame, filename in sections:
        frame.to_csv(ROOT / "audit" / filename, sep="\t", index=False)
        bad = int((frame.status == "FAIL").sum())
        failures += bad
        lines.append(f"{title}: {int((frame.status == 'PASS').sum())} passed, {bad} failed")
        for record in frame.itertuples():
            status = record.status
            if status == "PASS":
                continue
            label = getattr(record, "test", None) or getattr(record, "claim", "")
            detail = (f"expected {record.expected_value}, recomputed {record.recomputed_value}"
                      if hasattr(record, "expected_value")
                      else f"expected {record.expected}, observed {record.observed}")
            lines.append(f"  [{status}] {label} - {detail}")
        lines.append("")
    lines.append("PASS" if failures == 0 else f"FAIL: {failures} failing checks")
    REPORT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines[-1:]))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

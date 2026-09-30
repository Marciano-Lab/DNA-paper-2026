"""Load, verify and compile the distributed supplementary datasets."""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "supplementary_data"
INPUT_DIR = ROOT / "data" / "analysis_inputs"
TABLES = ["S1", "S2", "S3a", "S3b", "S4", "S5", "S6"]
BOOLEAN_WORDS = {"true": True, "false": False, "yes": True, "no": False, "1": True, "0": False}


def sha256(path) -> str:
    """SHA-256 of a file, read in chunks."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def is_boolean_column(series) -> bool:
    """True when every non-missing value is a recognized boolean word."""
    values = {str(v).strip().lower() for v in series.dropna().unique()}
    return bool(values) and values <= set(BOOLEAN_WORDS)


def as_bool(series):
    """Coerce a recognized boolean column to real booleans."""
    return series.astype(str).str.strip().str.lower().map(BOOLEAN_WORDS).fillna(False)


def load(name: str, typed: bool = True) -> pd.DataFrame:
    """Load one supplementary dataset by its table label."""
    frame = pd.read_csv(DATA_DIR / f"{name}.csv")
    if typed:
        for column in frame.columns:
            if is_boolean_column(frame[column]):
                frame[column] = as_bool(frame[column])
    return frame


def load_all(typed: bool = True) -> dict:
    """Load every distributed dataset."""
    return {name: load(name, typed=typed) for name in TABLES}


def index() -> pd.DataFrame:
    """The supplementary table index as distributed."""
    return pd.read_csv(DATA_DIR / "Supplementary_Table_Index.csv")


def data_dictionary() -> pd.DataFrame:
    """The data dictionary as distributed."""
    return pd.read_csv(DATA_DIR / "data_dictionary.csv")


def analysis_input(name: str) -> pd.DataFrame:
    """Load one frozen analysis input by file name."""
    return pd.read_csv(INPUT_DIR / name, sep="\t")


def build_workbook(path) -> Path:
    """Write the publication workbook from the canonical delimited tables."""
    import openpyxl
    from openpyxl.utils.dataframe import dataframe_to_rows

    book = openpyxl.Workbook()
    book.remove(book.active)
    for label, frame in [("Index", index())] + [(n, load(n, typed=False)) for n in TABLES] + [
            ("Data_dictionary", data_dictionary())]:
        sheet = book.create_sheet(label)
        for row in dataframe_to_rows(frame, index=False, header=True):
            sheet.append(row)
    book.save(path)
    return Path(path)


def verify_workbook(path) -> pd.DataFrame:
    """Compare every workbook cell against the canonical table it came from."""
    import openpyxl

    book = openpyxl.load_workbook(path, read_only=True)
    rows = []
    for name in TABLES:
        frame = load(name, typed=False)
        sheet = book[name]
        values = list(sheet.iter_rows(values_only=True))
        header, body = values[0], values[1:]
        mismatches = 0
        if list(header) != list(frame.columns) or len(body) != len(frame):
            mismatches = max(len(body), len(frame))
        else:
            for written, (_, expected) in zip(body, frame.iterrows()):
                for a, b in zip(written, expected.tolist()):
                    if (a is None and pd.isna(b)) or str(a) == str(b):
                        continue
                    try:
                        if abs(float(a) - float(b)) < 1e-9:
                            continue
                    except (TypeError, ValueError):
                        pass
                    mismatches += 1
        rows.append(dict(table=name, rows=len(frame), columns=frame.shape[1],
                         cells=len(frame) * frame.shape[1], mismatches=mismatches))
    return pd.DataFrame(rows)

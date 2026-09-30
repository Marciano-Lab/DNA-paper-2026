"""Structural, semantic and release quality control for the distributed archive."""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from . import tables
from .tables import as_bool, is_boolean_column, load

RESIDUE_FIELDS = r"reference_amino_acid|alternate_amino_acid|product_amino_acid"
RESIDUE_PATTERN = r"[A-Z*]"  # the stop symbol is a legitimate product residue
CODON_FIELDS = r"reference_codon|product_codon"
RETIRED_TERMS = ["strict_neurological_convergence_subset", "displayed_in_Figure3d",
                 "passes_correction", "exact_codon_reachable_by_c_to_a",
                 "transcript_variant_exclusion", "genomic_variant_exclusion"]
RETIRED_REFERENCES = ["Supplementary Figure S3", "Figure 3d", "Figure 4d",
                      "Supplementary Table S7", "Supplementary Table S8",
                      "Supplementary Table S9"]


def _row(test, ok, expected="", observed=""):
    return dict(test=test, status="PASS" if ok else "FAIL",
                expected=str(expected), observed=str(observed))


def structural_qc() -> pd.DataFrame:
    """File-level and field-level validity of every distributed dataset."""
    ix = tables.index()
    typed = {n: load(n) for n in tables.TABLES}
    raw = {n: load(n, typed=False) for n in tables.TABLES}
    rows = []
    missing = [n for n in tables.TABLES if not (tables.DATA_DIR / f"{n}.csv").exists()]
    rows.append(_row("every indexed dataset file exists", not missing, len(tables.TABLES),
                     missing or "all present"))
    rows.append(_row("the index lists exactly the distributed datasets",
                     list(ix.worksheet) == tables.TABLES, tables.TABLES, list(ix.worksheet)))
    duplicates = [n for n in tables.TABLES
                  if len(set(typed[n].columns)) != len(typed[n].columns)]
    rows.append(_row("no duplicate column names", not duplicates, "none", duplicates or "none"))
    leaked = [n for n in tables.TABLES
              if any(str(c).startswith("Unnamed") for c in typed[n].columns)]
    rows.append(_row("no row-index column", not leaked, "none", leaked or "none"))
    shape = [r.worksheet for r in ix.itertuples()
             if len(typed[r.worksheet]) != int(r.rows)
             or len(typed[r.worksheet].columns) != int(r.columns)]
    rows.append(_row("index row and column counts match the files", not shape, "none",
                     shape or "none"))
    empty = [f"{n}.{c}" for n in tables.TABLES for c in typed[n].columns
             if typed[n][c].isna().all()]
    rows.append(_row("no entirely empty column", not empty, "none", empty or "none"))
    retired = [f"{n}.{c}" for n in tables.TABLES for c in typed[n].columns
               if str(c) in RETIRED_TERMS]
    rows.append(_row("no retired field is distributed", not retired, "none", retired or "none"))
    bad_codon = [f"{n}.{c}" for n in tables.TABLES for c in typed[n].columns
                 if re.fullmatch(CODON_FIELDS, str(c))
                 and not raw[n][c].dropna().astype(str).str.fullmatch("[ACGT]{3}").all()]
    rows.append(_row("codon fields contain three unambiguous bases", not bad_codon, "none",
                     bad_codon or "none"))
    bad_residue = [f"{n}.{c}" for n in tables.TABLES for c in typed[n].columns
                   if re.fullmatch(RESIDUE_FIELDS, str(c))
                   and not raw[n][c].dropna().astype(str).str.fullmatch(RESIDUE_PATTERN).all()]
    rows.append(_row("residue fields are single-letter codes or the stop symbol", not bad_residue, "none",
                     bad_residue or "none"))
    bad_bool = [f"{n}.{c}" for n in tables.TABLES for c in typed[n].columns
                if is_boolean_column(raw[n][c]) and typed[n][c].dtype != bool]
    rows.append(_row("boolean fields coerce cleanly", not bad_bool, "none", bad_bool or "none"))
    bad_q = [f"{n}.{c}" for n in tables.TABLES for c in typed[n].columns
             if str(c) in ("p_value", "q_value")
             and not typed[n][c].dropna().between(0, 1).all()]
    rows.append(_row("probability fields lie between zero and one", not bad_q, "none",
                     bad_q or "none"))
    return pd.DataFrame(rows)


def dictionary_qc() -> pd.DataFrame:
    """Every distributed field is described once, and every description has a field."""
    dictionary = tables.data_dictionary()
    described = set(zip(dictionary.dataset, dictionary.column))
    distributed = {(n, c) for n in tables.TABLES for c in load(n).columns}
    rows = [_row("every distributed field has a dictionary entry",
                 not distributed - described, "none",
                 sorted(distributed - described) or "none"),
            _row("every dictionary entry describes a distributed field",
                 not described - distributed, "none",
                 sorted(described - distributed) or "none"),
            _row("no dictionary entry is duplicated",
                 len(described) == len(dictionary), len(described), len(dictionary)),
            _row("no dictionary field is blank",
                 not dictionary.isna().any().any()
                 and not (dictionary.astype(str).apply(lambda s: s.str.strip() == "")).any().any(),
                 "no blanks", "checked every cell")]
    return pd.DataFrame(rows)


def semantic_qc() -> pd.DataFrame:
    """Biological consistency within and between the distributed datasets."""
    s1, s2, s4, s5, s6 = load("S1"), load("S2"), load("S4"), load("S5"), load("S6")
    model = tables.analysis_input("gene_context_model_input.tsv")
    rows = []
    missense = s1[s1.consequence_class.str.lower() == "missense"]
    pairs = set(zip(missense.reference_amino_acid, missense.product_amino_acid))
    observed = set(zip(s2.reference_amino_acid, s2.alternate_amino_acid))
    rows.append(_row("every disease state uses an accessible substitution direction",
                     observed <= pairs, "subset of the 18 directions",
                     sorted(observed - pairs) or "all within"))
    rows.append(_row("route classes partition the disease states",
                     int(s2.route_class.notna().sum()) == len(s2), len(s2),
                     int(s2.route_class.notna().sum())))
    higher = as_bool(s2.higher_confidence_subset)
    rows.append(_row("the higher-confidence subset requires at least two review stars",
                     bool((s2[higher].clinvar_review_stars >= 2).all()), "all at least two",
                     int(s2[higher].clinvar_review_stars.min())))
    rows.append(_row("no higher-confidence state is splice-proximal",
                     not bool(as_bool(s2[higher].splice_proximal).any()), "none",
                     int(as_bool(s2[higher].splice_proximal).sum())))
    rows.append(_row("accessible states never exceed eligible states per gene",
                     bool((model.accessible_states <= model.eligible_states).all()), "always",
                     int((model.accessible_states > model.eligible_states).sum())))
    direction = s4.direction.where(s4.odds_ratio.isna(),
                                   s4.odds_ratio.map(lambda v: "enriched" if v > 1 else "depleted"))
    rows.append(_row("the pathway direction field matches the odds ratio",
                     bool((s4.direction == direction).all()), "always",
                     int((s4.direction != direction).sum())))
    rows.append(_row("the pathway correction flag matches the q value",
                     bool((as_bool(s4.passes_fdr) == (s4.q_value < 0.05)).all()), "always",
                     int((as_bool(s4.passes_fdr) != (s4.q_value < 0.05)).sum())))
    rows.append(_row("every displayed gene set passes correction",
                     bool((s4[as_bool(s4.displayed_in_Figure3c)].q_value < 0.05).all()),
                     "always", int((s4[as_bool(s4.displayed_in_Figure3c)].q_value >= 0.05).sum())))
    rows.append(_row("every displayed gene set carries a selection reason",
                     bool(s4[as_bool(s4.displayed_in_Figure3c)]
                          .display_selection_reason.notna().all()), "always",
                     int(s4[as_bool(s4.displayed_in_Figure3c)]
                         .display_selection_reason.isna().sum())))
    categories = set(s5.oxidized_guanine_support_category)
    rows.append(_row("lesion-support categories are the four defined values",
                     categories == {"endogenous", "induced-only", "not detected", "unresolved"},
                     "four defined categories", sorted(categories)))
    rows.append(_row("lesion-support categories partition the observed states",
                     int(s5.oxidized_guanine_support_category.notna().sum()) == len(s5),
                     len(s5), int(s5.oxidized_guanine_support_category.notna().sum())))
    rows.append(_row("observed states flagged as matching ClinVar are present in the state table",
                     all(((s2.gene == r.gene) & (s2.protein_position == r.protein_position)
                          & (s2.alternate_amino_acid == r.alternate_amino_acid)).any()
                         for r in s5[as_bool(s5.clinvar_match)].itertuples()),
                     "all matched states found", "checked against S2"))
    supported = s6[s6.analysis.str.contains("interaction", na=False)]
    rows.append(_row("every reported interaction estimate carries an interval",
                     bool(supported[["ci_lower", "ci_upper"]].notna().all().all()), "always",
                     int(supported[["ci_lower", "ci_upper"]].isna().sum().sum())))
    return pd.DataFrame(rows)


def reference_qc(root) -> pd.DataFrame:
    """No retired figure, table or field name survives in the release.

    The analysis package is excluded because the retired terms are defined there, and
    generated output directories are excluded because this check writes the terms it
    finds into its own report.
    """
    root = Path(root)
    offenders = {}
    skip_directories = {"__pycache__", "analysis", "audit", "derived"}
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix.lower() in {".png", ".pdf", ".svg", ".xlsx", ".gz"}:
            continue
        if skip_directories & set(path.relative_to(root).parts):
            continue
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        for term in RETIRED_REFERENCES + RETIRED_TERMS:
            if term in text:
                offenders.setdefault(term, []).append(str(path.relative_to(root)))
    return pd.DataFrame([_row("no retired reference survives in the release", not offenders,
                              "none", {k: v[:3] for k, v in offenders.items()} or "none")])


def checksum_qc(manifest_path) -> pd.DataFrame:
    """Recompute the SHA-256 of every byte-stable file named in the release manifest."""
    manifest = pd.read_csv(manifest_path, sep="\t")
    root = Path(manifest_path).resolve().parent
    if "byte_stable" in manifest.columns:
        manifest = manifest[manifest.byte_stable == "yes"]
    bad = []
    for record in manifest.itertuples():
        path = root / record.relative_path
        if not path.exists():
            bad.append(f"{record.relative_path} missing")
        elif tables.sha256(path) != record.sha256:
            bad.append(f"{record.relative_path} changed")
    return pd.DataFrame([_row("every byte-stable file matches its recorded checksum", not bad,
                              f"{len(manifest)} files", bad or "all match")])


def figure_qc(path) -> pd.DataFrame:
    """Read the geometry report written by the figure driver."""
    frame = pd.read_csv(path)
    return pd.DataFrame([
        _row("every figure is at most 175 mm wide", bool((frame.width_mm <= 175).all()),
             "at most 175 mm", f"maximum {frame.width_mm.max():.1f} mm"),
        _row("no label below 7.5 pt", int(frame.labels_under_7_5pt.sum()) == 0, 0,
             int(frame.labels_under_7_5pt.sum())),
        _row("no overlapping labels", int(frame.overlapping_pairs.sum()) == 0, 0,
             int(frame.overlapping_pairs.sum()))])

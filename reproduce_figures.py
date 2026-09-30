#!/usr/bin/env python
"""Regenerate every figure from the archived figure inputs and check their geometry."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code" / "figures"))
INPUT_DIR = ROOT / "data" / "figure_inputs"
FIGURE_DIR = ROOT / "derived" / "figures"
REPORT = ROOT / "audit" / "figure_geometry.csv"

FIGURES = [("figure_1", "build_F1", "figure1_input.json", "figure1_genetic_code"),
           ("figure_2", "build_F2", "figure2_input.json", "figure2_disease_convergence"),
           ("figure_3", "build_F3", "figure3_input.json", "figure3_gene_context"),
           ("figure_4", "build_F4", "figure4_input.json", "figure4_proteomic_observation"),
           ("supplementary_figure_s1", "build_S1", "supplementary_figure_s1_input.json",
            "figureS1_specification_sensitivity"),
           ("supplementary_figure_s2", "build_S2", "supplementary_figure_s2_input.json",
            "figureS2_lesion_evidence")]
MAXIMUM_WIDTH_MM = 175.0
MINIMUM_FONT_PT = 7.5


def restore(stem, payload):
    """Restore the container types each builder expects, preserving every stored key."""
    payload = dict(payload)
    if stem == "figure1_genetic_code":
        payload["PAIRS"] = {tuple(p) for p in payload["PAIRS"]}
        payload["ANYP"] = {tuple(p) for p in payload["ANYP"]}
        for key in ("REPL", "INTRO", "DIRN"):
            payload[key] = set(payload[key])
    if stem == "figure2_disease_convergence":
        payload["SUB"] = [(k, tuple(v)) for k, v in payload["SUB"]]
        payload["CVC"] = [tuple(v) for v in payload["CVC"]]
        payload["FAM"] = [tuple(v) for v in payload["FAM"]]
    return payload


def measure(figure, axes):
    """Report undersized text and overlapping labels at final rendered size."""
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    boxes, small = [], []
    for ax in axes:
        items = list(ax.texts)
        if ax.axison:
            items += [t for t in ax.get_xticklabels() + ax.get_yticklabels() if t.get_text()]
        legend = ax.get_legend()
        if legend:
            items += legend.get_texts()
        items += [ax.xaxis.label, ax.yaxis.label]
        for text in items:
            if not text.get_text().strip():
                continue
            if text.get_fontsize() < MINIMUM_FONT_PT:
                small.append(text.get_text()[:20])
            try:
                boxes.append((text.get_text()[:22], text.get_window_extent(renderer=renderer)))
            except Exception:
                pass
    overlaps = []
    for i in range(len(boxes)):
        for k in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[k][1]
            wide = min(a.x1, b.x1) - max(a.x0, b.x0)
            high = min(a.y1, b.y1) - max(a.y0, b.y0)
            if wide > 0 and high > 0 and wide * high > 12:
                overlaps.append(f"{boxes[i][0]} x {boxes[k][0]}")
    return small, overlaps


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    import figure_style

    rows, failures = [], 0
    for module, builder, payload_name, stem in FIGURES:
        payload = restore(stem, json.loads((INPUT_DIR / payload_name).read_text()))
        namespace = {}
        exec((ROOT / "code" / "figures" / f"{module}.py").read_text(), namespace)
        figure, axes = namespace[builder](payload)
        small, overlaps = measure(figure, list(axes))
        for extension in ("png", "pdf"):
            figure.savefig(FIGURE_DIR / f"{stem}.{extension}",
                           dpi=600 if extension == "png" else None,
                           bbox_inches="tight", pad_inches=0.02)
        plt.close(figure)
        image = Image.open(FIGURE_DIR / f"{stem}.png")
        width_mm = image.width / 600 * 25.4
        rows.append(dict(figure=stem, font=figure_style.FONT_FAMILY,
                         width_mm=round(width_mm, 1),
                         height_mm=round(image.height / 600 * 25.4, 1),
                         labels_under_7_5pt=len(small), overlapping_pairs=len(overlaps),
                         detail="; ".join(overlaps[:3])))
        failures += int(width_mm > MAXIMUM_WIDTH_MM) + len(small) + len(overlaps)
    pd.DataFrame(rows).to_csv(REPORT, index=False)
    print(f"   {len(rows)} figures regenerated, {failures} geometry failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

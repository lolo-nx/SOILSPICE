"""Compose rectangle I-V curves, curve errors and the model comparison.

Run: python Nature-builder-validation.py
Writes a 183 × 170 mm fig_validation_panel.pdf (vector) and .png in results/1_validation/.
"""
import argparse
import importlib.util
from pathlib import Path

import pymupdf
import pandas as pd

from soilspice import config, fig_validation


def without_panel_letters(source, output):
    """Remove standalone figure letters from a copy, preserving plots and text."""
    with pymupdf.open(source) as doc:
        page = doc[0]
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                for span in line["spans"]:
                    if (span["text"] in tuple("abcdefgh")
                            and "bold" in span["font"].lower()):
                        page.add_redact_annot(span["bbox"], fill=None, cross_out=False)
        page.apply_redactions(images=0, graphics=0)
        doc.save(output, garbage=3, deflate=True)


def match_panel_frame(source, reference, output):
    """Fit a vector PDF into the same page dimensions as the reference panel."""
    with pymupdf.open(source) as src, pymupdf.open(reference) as ref, pymupdf.open() as doc:
        rect = ref[0].rect
        page = doc.new_page(width=rect.width, height=rect.height)
        # Keep the common frame when the compositor trims empty margins.
        page.draw_rect(page.rect, color=None, fill=(1, 1, 1))
        page.show_pdf_page(page.rect, src, 0, keep_proportion=True)
        doc.save(output, garbage=3, deflate=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=config.CH1_DIR / "fig_validation_panel")
    args = parser.parse_args()
    data = config.CH1_DIR / "data"
    inputs = {"a": data / "iv_grid_rectangle.csv",
              "b": config.VALIDATION_DIR / "curve_errors_summary.csv",
              "c": data / "three_models_measured.csv"}
    for source in inputs.values():
        if not source.is_file():
            parser.error(f"Source data not found: {source}")

    data.mkdir(parents=True, exist_ok=True)
    # Redraw compact versions from saved results to fit the fixed page size.
    # Standalone figures retain their original dimensions.
    sources = {"a": data / "iv_grid_rectangle_compact.pdf",
               "b": data / "errors_grid_compact.pdf"}
    fig_validation.iv_grid(pd.read_csv(inputs["a"]), sources["a"].with_suffix(""),
                           figsize=(7.2, 3.15), compact_layout=True)
    fig_validation.errors_grid(pd.read_csv(inputs["b"]), sources["b"].with_suffix(""),
                               figsize=(3.4, 2.8))
    panels = {}
    for letter in ("a", "b"):
        panels[letter] = data / f"{sources[letter].stem}_panel.pdf"
        without_panel_letters(sources[letter], panels[letter])

    model_source = data / "fig_model_hierarchy_compact"
    fig_validation.model_hierarchy(pd.read_csv(inputs["c"]), model_source, figsize=(3.5, 2.8))
    panels["c"] = data / "fig_model_hierarchy_panel.pdf"
    match_panel_frame(model_source.with_suffix(".pdf"), panels["b"], panels["c"])

    spec = importlib.util.spec_from_file_location(
        "arrange_figure", config.REPO_DIR / "soilspice" / "arrange-fig.py")
    composer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(composer)
    composer.PANELS = panels
    composer.LAYOUT = [["a", "a"], ["b", "c"]]
    composer.COL_FRACTIONS = [0.5, 0.5]
    composer.PANEL_SIZE = {}
    composer.FIG_WIDTH_MM = 183.0
    composer.GAP_X_MM = 5.0
    composer.GAP_Y_MM = 3.0
    composer.LABEL_SPACE_MM = 4.0
    # Total height: 90 mm + 3 mm gap + 77 mm = 170 mm.
    composer.ROW_HEIGHTS_MM = [90.0, 77.0]
    composer.LABEL_SIZE_PT = 11
    composer.LABEL_FONT = "hebo"
    composer.OUTPUT = str(args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    composer.build()
    # Keep an existing viewer-reset copy synchronized with the main figure.
    upright = args.output.with_name(args.output.stem + "_upright.pdf")
    if upright.exists():
        upright.write_bytes(args.output.with_suffix(".pdf").read_bytes())


if __name__ == "__main__":
    main()

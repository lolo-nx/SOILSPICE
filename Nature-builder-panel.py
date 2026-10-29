"""Compose soiling photographs, electrical layouts and the simulation framework.

Run from the repository: python Nature-builder-panel.py
Writes a 183 × 148 mm PDF (retaining source vectors) and a 600 dpi PNG.
The legacy --without-photos layout retains its 183 × 170 mm page.
"""
import argparse
import importlib.util
from pathlib import Path

import pymupdf
from matplotlib import font_manager

from soilspice import config, fig_soiling, schematics


def framework_without_step_labels(source, output):
    """Centre the framework headings without step prefixes, retaining vectors."""
    # Horizontal centres of the four illustrations in the source PDF (points).
    illustration_centres = {
        "Band on the module": 52.685,
        "Cell coverage": 189.85,
        "LTspice netlist": 322.045,
        "I–V curve and losses": 454.235,
    }
    font_path = font_manager.findfont(font_manager.FontProperties(
        family="DejaVu Sans", weight="bold"))
    font = pymupdf.Font(fontfile=font_path)
    with pymupdf.open(source) as doc:
        page = doc[0]
        headings = []
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                for span in line["spans"]:
                    title = span["text"]
                    if title.startswith(("(a) ", "(b) ", "(c) ", "(d) ")):
                        title = title[4:]
                    if title in illustration_centres:
                        headings.append((title, span["size"], span["origin"][1]))
                        page.add_redact_annot(span["bbox"], fill=(1, 1, 1), cross_out=False)
        page.apply_redactions(images=0, graphics=0)
        for title, size, baseline in headings:
            centre = illustration_centres[title] * page.rect.width / 514.572937
            x = centre - font.text_length(title, fontsize=size) / 2
            page.insert_text((x, baseline), title, fontsize=size,
                             fontname="FrameworkTitle", fontfile=font_path)
        doc.save(output, garbage=3, deflate=True)


def main():
    assets = config.INPUT_DIR / "figure_assets"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--framework", type=Path, default=assets / "simulation_framework.pdf")
    parser.add_argument("--output", type=Path,
                        default=config.CH1_DIR / "fig_soiling_layout_framework")
    parser.add_argument("--photos", type=Path, default=assets / "Soiling-images.png")
    parser.add_argument("--without-photos", action="store_true",
                        help="Reproduce the previous two-panel layout without photographs")
    args = parser.parse_args()
    if not args.framework.is_file():
        parser.error(f"Source figure not found: {args.framework}")

    electrical = config.CH1_DIR / "pv_module_schematics"
    schematics.draw(electrical)
    # Architecture identifiers sit beside the module names at the original text size.
    electrical_name = "pv_module_schematics_panel" if args.without_photos else "pv_module_schematics_soiling_panel"
    electrical_panel = config.CH1_DIR / "data" / electrical_name
    schematics.draw(electrical_panel, text_scale=1.0 if args.without_photos else 1.25,
                    panel_labels=False, title_labels=True)
    framework_panel = config.CH1_DIR / "data" / "simulation_framework_panel.pdf"
    framework_without_step_labels(args.framework, framework_panel)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if not args.without_photos:
        fig_soiling.draw(args.photos, electrical_panel.with_suffix(".pdf"), framework_panel, args.output)
        return

    # Reuse the repository's PDF/raster compositor, with a layout specific to
    # this figure. Both source PDFs remain vector.
    spec = importlib.util.spec_from_file_location(
        "arrange_figure", config.REPO_DIR / "soilspice" / "arrange-fig.py")
    composer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(composer)
    composer.PANELS = {"a": electrical_panel.with_suffix(".pdf"), "b": framework_panel}
    composer.LAYOUT = [["a"], ["b"]]
    composer.COL_FRACTIONS = [1.0]
    composer.PANEL_SIZE = {}
    composer.FIG_WIDTH_MM = 183.0
    composer.GAP_Y_MM = 3.0
    composer.LABEL_SPACE_MM = 4.0
    # Total height: 100 mm + 3 mm gap + 67 mm = 170 mm.
    composer.ROW_HEIGHTS_MM = [100.0, 67.0]
    composer.LABEL_SIZE_PT = 11
    composer.LABEL_FONT = "hebo"  # Helvetica Bold
    composer.OUTPUT = str(args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    composer.build()


if __name__ == "__main__":
    main()

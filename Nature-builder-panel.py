"""Compose soiling photographs, electrical layouts and the simulation framework.

Run from the repository: python Nature-builder-panel.py
Writes a 183 × 148 mm PDF (retaining source vectors) and a 600 dpi PNG.
The legacy --without-photos layout retains its 183 × 170 mm page.
Use --framework-only for a standalone framework PDF and 600 dpi PNG.
"""
import argparse
import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory

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


def standalone_framework(source, output):
    """Export only the vector framework at 183 mm wide, with 3 mm margins."""
    mm = fig_soiling.MM
    with pymupdf.open(source) as workflow, pymupdf.open() as doc:
        clip = fig_soiling._bounds(workflow[0])
        height_mm = 177.0 * clip.height / clip.width + 6.0
        page = doc.new_page(width=183.0 * mm, height=height_mm * mm)
        page.show_pdf_page(fig_soiling._rect(3, 3, 177, height_mm - 6),
                           workflow, 0, clip=clip)
        doc.save(output.with_suffix(".pdf"), garbage=3, deflate=True)
        page.get_pixmap(dpi=600).save(output.with_suffix(".png"))
    print(f"Saved {output.with_suffix('.pdf')} and .png (183 × {height_mm:.1f} mm)")


def main():
    assets = config.INPUT_DIR / "figure_assets"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--framework", type=Path, default=assets / "simulation_framework.pdf")
    parser.add_argument("--output", type=Path,
                        help="Output file stem (default depends on the selected layout)")
    parser.add_argument("--photos", type=Path, default=assets / "Soiling-images.png")
    layout = parser.add_mutually_exclusive_group()
    layout.add_argument("--without-photos", action="store_true",
                        help="Reproduce the previous two-panel layout without photographs")
    layout.add_argument("--framework-only", action="store_true",
                        help="Export only the simulation framework as PDF and PNG")
    args = parser.parse_args()
    if not args.framework.is_file():
        parser.error(f"Source figure not found: {args.framework}")
    if args.output is None:
        name = "fig_simulation_framework" if args.framework_only else "fig_soiling_layout_framework"
        args.output = config.CH1_DIR / name
    args.output.parent.mkdir(parents=True, exist_ok=True)

    if args.framework_only:
        with TemporaryDirectory(prefix="soilspice-framework-") as temporary:
            framework_panel = Path(temporary) / "framework.pdf"
            framework_without_step_labels(args.framework, framework_panel)
            standalone_framework(framework_panel, args.output)
        return

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

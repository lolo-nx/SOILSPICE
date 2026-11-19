"""Compose soiling photographs, electrical layouts and the simulation framework.

Run from the repository: python Nature-builder-panel.py
Writes a 183 × 148 mm PDF (retaining source vectors) and a 600 dpi PNG.
The legacy --without-photos layout retains its 183 × 170 mm page.
Use --framework-only for a standalone framework PDF and 600 dpi PNG.
"""
import argparse
import importlib.util
import re
from pathlib import Path
from tempfile import TemporaryDirectory

import pymupdf
from matplotlib import font_manager

from soilspice import config, fig_soiling, schematics


def annotate_module_geometry(page):
    """Redraw the schematic with a narrow, identifiable glass margin."""
    sx, sy = page.rect.width / 514.572937, page.rect.height / 175.956406
    def point(x, y):
        return pymupdf.Point(x * sx, y * sy)
    def rect(x0, y0, x1, y1):
        return pymupdf.Rect(x0*sx, y0*sy, x1*sx, y1*sy)
    # Replace only the first illustration, preserving headings, workflow and
    # explanatory text. Keep the lower glass margin comparable to the top.
    page.add_redact_annot(rect(10, 20, 114, 123),
                          fill=(1, 1, 1), cross_out=False)
    page.apply_redactions(images=0, graphics=2)
    font_path = font_manager.findfont(font_manager.FontProperties(family="DejaVu Sans"))
    font = pymupdf.Font(fontfile=font_path)
    page.insert_font(fontname="GeometryAnnotation", fontbuffer=font.buffer)
    blue, amber, ink = (0.16, 0.36, 0.58), (0.60, 0.35, 0.08), (0.18, 0.18, 0.18)
    # Warm glass surrounds blue active cells; the lower glass strip is now
    # visible independently of the dust band, rather than a one-pixel gap.
    page.draw_rect(rect(16, 25, 72, 117), color=ink,
                   fill=(1.0, 0.94, 0.78), width=0.9*sx)
    for row in range(10):
        for col in range(6):
            x, y = 19 + col*8.5, 29 + row*8.5
            page.draw_rect(rect(x, y, x+7.5, y+7.5), color=blue,
                           fill=(0.87, 0.92, 0.97), width=0.3*sx)
    # The band starts at the bottom edge and covers both inactive glass and
    # part of the bottom cells. Its h is not the height of the inactive strip.
    page.draw_rect(rect(16, 105, 72, 117), color=(0.35, 0.35, 0.35),
                   fill=(0.35, 0.35, 0.35), fill_opacity=0.32, width=0.5*sx)
    for x, y, text, color in ((84, 37, "Active cell", blue), (84, 44, "area", blue),
                              (84, 86, "Inactive glass", amber), (84, 93, "margin", amber),
                              (85, 111, "Band height h", ink)):
        page.insert_text(point(x, y), text, fontname="GeometryAnnotation",
                         fontsize=6.1*sx, color=color)
    # Short, colour-matched elbow leaders terminate in the actual target area.
    for vertices, color in ((((82, 42), (76, 42), (65, 48)), blue),
                             (((82, 94), (74, 94), (74, 115), (65, 115)), amber)):
        for start, end in zip(vertices, vertices[1:]):
            page.draw_line(point(*start), point(*end), color=color, width=0.55*sx)
        page.draw_circle(point(*vertices[-1]), radius=0.7*sx, color=color, fill=color, width=0)
    # A separate dimension arrow makes the full band height unambiguous.
    page.draw_line(point(79, 105), point(79, 117), color=ink, width=0.6*sx)
    for y, direction in ((105, 1), (117, -1)):
        for x in (77.8, 80.2):
            page.draw_line(point(79, y), point(x, y+direction*2), color=ink, width=0.6*sx)
    page.insert_text(point(36.0, 112.1), "τ", fontname="GeometryAnnotation",
                     fontsize=7.5*sx, color=(0, 0, 0))
    page.insert_text(point(40.6, 113.4), "I", fontname="GeometryAnnotation",
                     fontsize=5.0*sx, color=(0, 0, 0))
    # Replace the original three-line paragraph with a concise input summary.
    # Stop above the bracket stroke at y=156.096 so it survives redaction.
    page.add_redact_annot(rect(17, 125, 99, 155.1), fill=(1, 1, 1), cross_out=False)
    page.apply_redactions(images=0, graphics=2)
    for y, text in ((134, "Module geometry + band"),
                    (143, "Height h, transmittance τ")):
        size = 6.4*sx
        x = 53*sx - font.text_length(text, fontsize=size)/2
        page.insert_text((x, y*sy), text, fontname="GeometryAnnotation",
                         fontsize=size, color=ink)


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
        # Translate only the final workflow arrow, preserving the nearby
        # rotated Current label and all other vector artwork.
        for xref in page.get_contents():
            stream = doc.xref_stream(xref)
            pattern = rb"381\.329091 105\.3064 m.*?\nB"
            shifted, count = re.subn(pattern,
                lambda m: b"q\n1 0 0 1 -4 0 cm\n" + m.group() + b"\nQ",
                stream, flags=re.DOTALL)
            if count:
                doc.update_stream(xref, shifted)
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
        page.insert_font(fontname="FrameworkTitle", fontbuffer=font.buffer)
        for title, size, baseline in headings:
            centre = illustration_centres[title] * page.rect.width / 514.572937
            x = centre - font.text_length(title, fontsize=size) / 2
            page.insert_text((x, baseline), title, fontsize=size,
                             fontname="FrameworkTitle")
        annotate_module_geometry(page)
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

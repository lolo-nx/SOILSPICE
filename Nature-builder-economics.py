"""Stack the reference design-lever and cleaning figures as a vector Nature panel.

Run: python Nature-builder-economics.py
Only RH=15 mm/event: 183 × 170 mm PDF and 600 dpi PNG, panels a--g.
Existing standalone figures, numerical results and physical models are untouched.
"""
import argparse
import hashlib
import json
from pathlib import Path

import pymupdf

from soilspice import config

WIDTH_MM, HEIGHT_MM = 183., 170.
MARGIN_MM, GAP_MM = 2.5, 4.
LABEL_PT, DPI = 9., 600
MM = 72 / 25.4


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_without_letters(path, expected):
    """Redact only standalone bold subplot letters in an in-memory PDF copy."""
    doc = pymupdf.open(path)
    if len(doc) != 1 or doc[0].get_images():
        raise ValueError(f"Expected one vector source page: {path}")
    page = doc[0]
    letters = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line["spans"]:
                if span["text"] in tuple(expected) and "bold" in span["font"].lower():
                    letters.append(span)
    letters.sort(key=lambda span: span["origin"][0])
    if "".join(span["text"] for span in letters) != expected:
        raise ValueError(f"Expected source labels {expected!r}: {path}")
    for span in letters:
        page.add_redact_annot(span["bbox"], fill=None, cross_out=False)
    page.apply_redactions(images=0, graphics=0)
    return doc, letters


def main():
    scenario = config.CH4_DIR / "total_scenarios/PH15mm"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=scenario / "fig_design_cleaning_panel")
    args = parser.parse_args()
    paths = [scenario / "fig_design_levers.pdf", scenario / "fig_cleaning.pdf"]
    for path in paths:
        if not path.is_file():
            parser.error(f"Missing source: {path}")
    protected = {path: digest(path) for path in config.CH4_DIR.rglob("*")
                 if path.is_file() and path.suffix in (".csv", ".tex", ".pdf", ".png")
                 and path not in {args.output.with_suffix(".pdf"), args.output.with_suffix(".png")}}
    docs, source_letters = zip(*(source_without_letters(path, expected)
                                for path, expected in zip(paths, ("abcd", "abc"))))
    available_w = (WIDTH_MM - 2 * MARGIN_MM) * MM
    available_h = (HEIGHT_MM - 2 * MARGIN_MM - GAP_MM) * MM
    natural_heights = [available_w * doc[0].rect.height / doc[0].rect.width for doc in docs]
    fit = min(1., available_h / sum(natural_heights))
    widths = [available_w * fit] * 2
    heights = [height * fit for height in natural_heights]
    y = (HEIGHT_MM * MM - sum(heights) - GAP_MM * MM) / 2
    result = pymupdf.open()
    page = result.new_page(width=WIDTH_MM * MM, height=HEIGHT_MM * MM)
    placements, letter_index = [], 0
    for doc, letters, width, height, path in zip(docs, source_letters, widths, heights, paths):
        x = (WIDTH_MM * MM - width) / 2
        box = pymupdf.Rect(x, y, x + width, y + height)
        source_rect = doc[0].rect
        scale = width / source_rect.width
        page.show_pdf_page(box, doc, 0, keep_proportion=True)
        for span in letters:
            position = (x + span["origin"][0] * scale, y + span["origin"][1] * scale)
            page.insert_text(position, "abcdefg"[letter_index], fontname="hebo", fontsize=LABEL_PT)
            letter_index += 1
        placements.append(dict(source=str(path.relative_to(config.REPO_DIR)), sha256=digest(path),
                               scale=scale, box_pdf_points=list(box)))
        y += height + GAP_MM * MM
    result.set_metadata({"title": "Design levers and cleaning — RH = 15 mm/event",
                         "subject": "Reference economics scenario; existing vector figures, unchanged calculations"})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.save(args.output.with_suffix(".pdf"), garbage=4, deflate=True)
    page.get_pixmap(dpi=DPI, alpha=False).save(args.output.with_suffix(".png"))
    assert not page.get_images(), "The composite PDF must remain entirely vector"
    labels = [span["text"] for block in page.get_text("dict")["blocks"]
              for line in block.get("lines", []) for span in line["spans"]
              if span["text"] in tuple("abcdefg") and "bold" in span["font"].lower()]
    assert labels == list("abcdefg"), labels
    for path, sha in protected.items():
        assert digest(path) == sha, f"Unexpected change to an existing output: {path}"
    metadata = dict(R_H_mm_event=15, width_mm=WIDTH_MM, height_mm=HEIGHT_MM, png_dpi=DPI,
                    panels="abcdefg", row_order=["design levers", "cleaning"], vector_pdf=True,
                    source_figures=placements, label_size_pt=LABEL_PT,
                    preexisting_exports_unchanged=len(protected),
                    model_or_simulation_changes=False,
                    pdf_sha256=digest(args.output.with_suffix(".pdf")),
                    png_sha256=digest(args.output.with_suffix(".png")))
    args.output.with_suffix(".json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    for doc in docs:
        doc.close()
    result.close()
    print(f"Saved {args.output}.pdf and .png ({WIDTH_MM:g} × {HEIGHT_MM:g} mm, {DPI} dpi).")
    print(f"Preserved {len(protected)} existing figures and data files; labels a–g verified.")


if __name__ == "__main__":
    main()

"""Compose the existing soiling photographs, electrical layouts and framework.

Photographs are clipped from the supplied annotated montage, without retouching
or synthesising image content. Electrical layouts and framework retain vectors.
"""
import hashlib
import json
from pathlib import Path

import pymupdf

MM = 72 / 25.4
WIDTH_MM, HEIGHT_MM = 183.0, 148.0
TOP_TRIM_MM = 10.0


def _rect(x, y, w, h):
    return pymupdf.Rect(x * MM, y * MM, (x + w) * MM, (y + h) * MM)


def _bounds(page):
    bounds = pymupdf.Rect()
    for _, box, *_ in page.get_bboxlog():
        bounds |= pymupdf.Rect(box)
    return bounds & page.rect


def _fit(page, source, box):
    clip = _bounds(source[0])
    scale = min(box.width / clip.width, box.height / clip.height)
    target = pymupdf.Rect(0, 0, clip.width * scale, clip.height * scale)
    target += (box.x0 + (box.width - target.width) / 2,
               box.y0 + (box.height - target.height) / 2) * 2
    page.show_pdf_page(target, source, 0, clip=clip)
    return scale


def draw(photo, electrical, framework, output):
    """Write a 183 × 148 mm, three-panel PDF and 600 dpi PNG."""
    output = Path(output)
    paths = [Path(photo), Path(electrical), Path(framework)]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    with pymupdf.open(photo) as image_doc:
        photo_pdf = image_doc.convert_to_pdf()
    with pymupdf.open("pdf", photo_pdf) as photographs, \
            pymupdf.open(electrical) as layouts, pymupdf.open(framework) as workflow, \
            pymupdf.open() as doc:
        page = doc.new_page(width=WIDTH_MM * MM, height=HEIGHT_MM * MM)
        # Coordinates in the original 1308 × 522 montage. Only its outer
        # borders and separator are cropped; annotations and insets are kept.
        clips_px = [(20, 27, 645, 494), (650, 27, 1278, 494)]
        photo_width = 43.0
        photo_info = []
        for crop, top in zip(clips_px, (16.0, 54.0)):
            x0, y0, x1, y1 = crop
            source = photographs[0].rect
            clip = pymupdf.Rect(x0 / 1308 * source.width, y0 / 522 * source.height,
                                x1 / 1308 * source.width, y1 / 522 * source.height)
            height = photo_width * (y1 - y0) / (x1 - x0)
            page.show_pdf_page(_rect(5, top - TOP_TRIM_MM, photo_width, height), photographs, 0, clip=clip)
            photo_info.append(dict(source_crop_pixels=crop, width_mm=photo_width,
                                   effective_dpi=(x1 - x0) / (photo_width / 25.4)))
        layout_box = _rect(51, 9 - TOP_TRIM_MM, 129, 87)
        scale_layout = _fit(page, layouts, layout_box)
        layout_top_mm = (layout_box.y0 + (layout_box.height - _bounds(layouts[0]).height * scale_layout) / 2) / MM
        scale_framework = _fit(page, workflow, _rect(3, 96 - TOP_TRIM_MM, 177, 59))
        label_positions = (("a", 2, 14 - TOP_TRIM_MM), ("b", 50, layout_top_mm - 2),
                           ("c", 2, 92 - TOP_TRIM_MM))
        for label, x, y in label_positions:
            page.insert_text((x * MM, y * MM), label, fontname="hebo", fontsize=9)
        doc.save(output.with_suffix(".pdf"), garbage=3, deflate=True)
        page.get_pixmap(dpi=600).save(output.with_suffix(".png"))
    for path in paths:
        if hashlib.sha256(path.read_bytes()).hexdigest() != hashes[str(path)]:
            raise ValueError(f"Source changed during composition: {path}")
    with pymupdf.open(output.with_suffix(".pdf")) as check:
        page = check[0]
        if (len(check) != 1 or abs(page.rect.width / MM - WIDTH_MM) > .001
                or abs(page.rect.height / MM - HEIGHT_MM) > .001):
            raise ValueError(f"The final figure must be one {WIDTH_MM:g} × {HEIGHT_MM:g} mm page")
        text = page.get_text()
        for required in ("Full-cell module", "Half-cell module", "Band on the module", "LTspice netlist"):
            if required not in text:
                raise ValueError(f"Missing editable panel text: {required}")
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                for span in line["spans"]:
                    if not page.rect.contains(pymupdf.Rect(span["bbox"])):
                        raise ValueError(f"Text outside the final page: {span['text']}")
    info = dict(width_mm=WIDTH_MM, height_mm=HEIGHT_MM, png_dpi=600,
                panels=dict(a="existing soiling photographs", b="electrical layouts", c="simulation framework"),
                photo_source=str(photo), photo_source_is_annotated_montage=True,
                original_photographs_available=False, photo_content_retouched=False,
                photos=photo_info, electrical_scale=scale_layout, framework_scale=scale_framework,
                panel_label_positions_mm={label: [x, y] for label, x, y in label_positions},
                electrical_and_framework_vectors_preserved=True, source_sha256=hashes)
    output.with_suffix(".json").write_text(json.dumps(info, indent=2) + "\n")
    output.with_name(output.name + "_caption.tex").write_text(
        r"Observed bottom-edge soiling, module architectures and simulation workflow. "
        r"\textbf{a}, Photographs of arched and rectangular dust bands, including a close-up of the latter. "
        r"\textbf{b}, Electrical layouts of the full-cell standard (A), full-cell unbalanced-substring (B) "
        r"and half-cell butterfly (C) architectures. "
        r"\textbf{c}, Module geometry and band properties determine cell coverage and the cell-level "
        r"SPICE netlist, from which I--V curves and electrical losses are calculated." + "\n")
    print(f"Saved {output.with_suffix('.pdf')} and .png ({WIDTH_MM:g} × {HEIGHT_MM:g} mm)")

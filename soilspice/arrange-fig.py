"""
Assemble several panels into one Nature-style figure (double column, 183 mm).

Each panel can be a PDF (kept as VECTOR: text and lines stay sharp and
editable) or a raster image (PNG, JPG, TIF). Panel letters (a, b, c, ...)
are added in 8 pt bold Helvetica, as required by Nature journals.

Usage
-----
1. Fill in the paths in PANELS below (one file per panel).
2. Choose the arrangement in LAYOUT (rows of panel letters).
3. Run:   python assemble_figure.py
   -> figure_1.pdf (for submission / LaTeX) and figure_1.png (preview)

Requirements:  pip install pymupdf pillow

Notes
-----
* EPS is not supported as input: export your panels as PDF instead.
* For the text to respect Nature's 5-7 pt rule, each panel should be made
  close to its final size. The script prints, for each panel, the scale
  factor applied: a 10 pt font in a panel scaled by 0.5 ends up at 5 pt.
"""

from pathlib import Path

import pymupdf  # PyMuPDF
from PIL import Image, ImageChops

# =====================================================================
# 1. INPUT FILES -- one per panel (PDF = vector, or PNG/JPG = image)
# =====================================================================
PANELS = {
    "a": "/Users/lorrainenaux/Desktop/Valorisation-master/soilspice_paper/results/2_simulation/regime_map.pdf",
    "b": "/Users/lorrainenaux/Desktop/Valorisation-master/soilspice_paper/results/2_simulation/regime_atlas.pdf",
    # "c": "/Users/lorrainenaux/Desktop/Valorisation-master/soilspice_paper/results/2_simulation/fig-modules.png",
    # "d": "/Users/lorrainenaux/Desktop/Valorisation-master/soilspice_paper/results/2_simulation/Soiling-images.png",
}
# =====================================================================
# 2. LAYOUT
# =====================================================================
# Rows of panels, top to bottom. A panel can span several columns by
# repeating its letter, e.g. [["a", "a"], ["b", "c"]].
LAYOUT = [
    ["a", "b"],
    # ["c", "d"],
]
COL_FRACTIONS = [0.5, 0.5]   # relative widths of the columns
ROW_HEIGHTS_MM = None        # None = automatic (from the panel aspect ratios)
                             # or e.g. [70, 60] to force the row heights

# Optional: shrink individual panels, as a fraction of their cell width
# (1.0 = fill the cell). Example: {"b": 0.75} makes panel b 25% smaller,
# which also lowers the height of its row.
PANEL_SIZE = {"b": 0.75}
PANEL_ALIGN = "center"       # "left", "center" or "right" inside the cell

OUTPUT = "figure_1"          # -> figure_1.pdf and figure_1.png

# =====================================================================
# 3. FORMAT (Nature: 89 mm single column, 183 mm double column)
# =====================================================================
FIG_WIDTH_MM = 183.0
GAP_X_MM = 4.0               # horizontal gap between panels
GAP_Y_MM = 4.0               # vertical gap between rows
LABEL_SPACE_MM = 4.0         # space reserved above each panel for its letter
LABEL_SIZE_PT = 8            # panel letter: 8 pt bold (Nature)
LABEL_FONT = "hebo"          # Helvetica-Bold (built-in PDF font)
TRIM_WHITESPACE = True       # crop empty margins around each panel
PREVIEW_DPI = 450            # resolution of the PNG preview

MM = 72 / 25.4               # 1 mm in PDF points


# ---------------------------------------------------------------- helpers
def content_rect_pdf(page):
    """Bounding box of everything drawn on a PDF page."""
    r = pymupdf.Rect()
    for item in page.get_bboxlog():
        r |= pymupdf.Rect(item[1])
    r &= page.rect
    return r if not r.is_empty else page.rect


def trim_image(path):
    """Crop uniform (white) margins of a raster image; return a PIL image."""
    im = Image.open(path).convert("RGB")
    if not TRIM_WHITESPACE:
        return im
    bg = Image.new("RGB", im.size, im.getpixel((0, 0)))
    box = ImageChops.difference(im, bg).getbbox()
    return im.crop(box) if box else im


class Panel:
    """One input file: its usable area and aspect ratio."""

    def __init__(self, letter, path):
        self.letter = letter
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"panel {letter}: file not found: {path}")
        self.is_pdf = self.path.suffix.lower() == ".pdf"
        if self.is_pdf:
            self.doc = pymupdf.open(self.path)
            page = self.doc[0]
            self.clip = content_rect_pdf(page) if TRIM_WHITESPACE else page.rect
            self.w, self.h = self.clip.width, self.clip.height
        else:
            self.image = trim_image(self.path)
            self.w, self.h = self.image.size

    @property
    def aspect(self):
        return self.h / self.w


# ------------------------------------------------------------ geometry
def column_edges(width_pt):
    """x positions (pt) of the left edge of each column, and column widths."""
    ncol = len(COL_FRACTIONS)
    usable = width_pt - (ncol - 1) * GAP_X_MM * MM
    total = sum(COL_FRACTIONS)
    widths = [usable * f / total for f in COL_FRACTIONS]
    lefts, x = [], 0.0
    for w in widths:
        lefts.append(x)
        x += w + GAP_X_MM * MM
    return lefts, widths


def panel_spans(row):
    """[(letter, first_col, last_col)] for one row of LAYOUT."""
    spans, c = [], 0
    while c < len(row):
        letter, start = row[c], c
        while c + 1 < len(row) and row[c + 1] == letter:
            c += 1
        spans.append((letter, start, c))
        c += 1
    return spans


def build():
    panels = {k: Panel(k, v) for k, v in PANELS.items()}
    width_pt = FIG_WIDTH_MM * MM
    lefts, widths = column_edges(width_pt)
    label_h = LABEL_SPACE_MM * MM

    # cell geometry for each panel
    rows = []
    for r, row in enumerate(LAYOUT):
        cells = []
        for letter, c0, c1 in panel_spans(row):
            x0 = lefts[c0]
            w = lefts[c1] + widths[c1] - x0
            cells.append((letter, x0, w))
        if ROW_HEIGHTS_MM:
            h = ROW_HEIGHTS_MM[r] * MM
        else:  # tallest panel of the row, when fitted to its (shrunk) width
            h = max(w * PANEL_SIZE.get(l, 1.0) * panels[l].aspect
                    for l, _, w in cells) + label_h
        rows.append((cells, h))

    height_pt = sum(h for _, h in rows) + (len(rows) - 1) * GAP_Y_MM * MM
    out = pymupdf.open()
    page = out.new_page(width=width_pt, height=height_pt)

    y = 0.0
    print(f"Figure: {FIG_WIDTH_MM:.0f} x {height_pt / MM:.0f} mm")
    for cells, row_h in rows:
        for letter, x0, cell_w in cells:
            p = panels[letter]
            box_h = row_h - label_h
            max_w = cell_w * PANEL_SIZE.get(letter, 1.0)
            scale = min(max_w / p.w, box_h / p.h)
            w, h = p.w * scale, p.h * scale
            # horizontal alignment inside the cell, top-aligned
            free = cell_w - w
            dx = {"left": 0.0, "center": free / 2, "right": free}[PANEL_ALIGN]
            rect = pymupdf.Rect(x0 + dx, y + label_h,
                                x0 + dx + w, y + label_h + h)
            if p.is_pdf:
                page.show_pdf_page(rect, p.doc, 0, clip=p.clip)
                print(f"  {letter}: vector PDF, scale x{scale:.2f} "
                      f"(a 10 pt font becomes {10 * scale:.1f} pt)")
            else:
                dpi = p.image.width / (w / 72)
                tmp = Path(f"_panel_{letter}.png")
                p.image.save(tmp)
                page.insert_image(rect, filename=str(tmp))
                tmp.unlink()
                flag = "" if dpi >= 450 else "   <-- below 450 dpi"
                print(f"  {letter}: image, {dpi:.0f} dpi at final size{flag}")
            page.insert_text((x0, y + LABEL_SIZE_PT * 0.85), letter,
                             fontname=LABEL_FONT, fontsize=LABEL_SIZE_PT)
        y += row_h + GAP_Y_MM * MM

    out.save(f"{OUTPUT}.pdf", garbage=3, deflate=True)
    page.get_pixmap(dpi=PREVIEW_DPI).save(f"{OUTPUT}.png")
    print(f"saved: {OUTPUT}.pdf and {OUTPUT}.png")


if __name__ == "__main__":
    build()
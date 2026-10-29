"""
Electrical layout of three PV module architectures (portrait view, front side)

  A    Full-cell module, standard           : 6 x 10 cells, 3 substrings (2/2/2 columns)
  B    Full-cell module, unbalanced         : 8 x 10 cells, 3 substrings (2/4/2 columns)
  C/D  Half-cell module, butterfly design   : 6 x 20 half-cells, top and bottom halves
                                              mirrored and connected in parallel,
                                              bypass diodes on the centre line
  E    Matrix shingle module (installed, landscape): 2 blocks x 13 rows x 9
                                              shingled half-cells; the half-cells of a
                                              row are in parallel, the rows in series;
                                              one bypass diode per block (j-boxes on
                                              the short edges)

Wiring follows the usual representation in the literature:
  * cells of a column are connected in series by ribbons (vertical lines);
  * adjacent columns are joined by cross-connector ribbons at the module edge
    (U-turns), so the whole module forms one serpentine string;
  * taps are taken at the cross-connectors located on the junction-box side
    and each bypass diode is placed between two consecutive taps
    (anti-parallel to its substring: cathode towards the + terminal);
  * half-cell (butterfly) module: each half is the mirror image of the other,
    U-turns at the outer edges, taps and diodes on the centre line, and the
    two halves are in parallel through the shared taps.
Arrows give the conventional current direction in operation (from - to +).

Output (draw): <path>.pdf (vector, for LaTeX) and <path>.png.
"""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

# PDF export: embed real (Type 42) fonts so the text stays editable
# (Illustrator, Inkscape) and renders correctly in LaTeX.
matplotlib.rcParams["pdf.fonttype"] = 42

# ----------------------------------------------------------------- style
BLUE = "#1f4e9c"
RED = "#c0392b"
BLACK = "#1a1a1a"
GAP_FRAC = 0.13          # gap between cells (fraction of pitch)
WIRE_LW = 1.3
EDGE = 0.25              # distance of the edge ribbons from the cell grid

MODULE_W = 6.0           # identical footprint for the three panels
MODULE_H = 10.0
XLIM = (-1.3, MODULE_W + 1.3)
YLIM = (-2.1, MODULE_H + 1.7)


# --------------------------------------------------------------- helpers
def wire(ax, xs, ys, lw=WIRE_LW):
    ax.plot(xs, ys, color=RED, lw=lw, solid_capstyle="butt",
            solid_joinstyle="miter", zorder=3)


def dot(ax, x, y):
    ax.plot(x, y, "o", color=RED, ms=3.6, zorder=4)


def arrow(ax, x, y, up, size=0.22):
    """Current-direction arrowhead centred at (x, y)."""
    s = size if up else -size
    ax.add_patch(Polygon([(x - abs(s) * 0.55, y - s / 2),
                          (x + abs(s) * 0.55, y - s / 2),
                          (x, y + s / 2)],
                         closed=True, facecolor=RED, edgecolor=RED, zorder=4))


def diode(ax, x, y, s=0.42):
    """Horizontal diode, cathode (bar) on the left, i.e. towards '+'."""
    ax.add_patch(Polygon([(x + s / 2, y - s / 2.3), (x + s / 2, y + s / 2.3),
                          (x - s / 2, y)],
                         closed=True, facecolor="white", edgecolor=BLACK,
                         lw=1.1, zorder=5))
    ax.plot([x - s / 2, x - s / 2], [y - s / 2.3, y + s / 2.3],
            color=BLACK, lw=1.6, zorder=5)


def cells(ax, ncols, nrows, dx, dy, y0):
    w, h = dx * (1 - GAP_FRAC), dy * (1 - GAP_FRAC)
    for c in range(ncols):
        for r in range(nrows):
            ax.add_patch(FancyBboxPatch(
                (c * dx + (dx - w) / 2, y0 + r * dy + (dy - h) / 2), w, h,
                boxstyle=f"round,pad=0,rounding_size={0.12 * min(w, h):.3f}",
                lw=0, facecolor=BLUE, zorder=2))


def tap_positions(col_x, group_sizes):
    """x of the taps: first column, junctions between substrings, last column."""
    taps, c = [col_x[0]], 0
    for g in group_sizes[:-1]:
        c += g
        taps.append((col_x[c - 1] + col_x[c]) / 2)
    taps.append(col_x[-1])
    return taps


def diode_row(ax, taps, y):
    """Diode line with + / - leads at both ends."""
    wire(ax, [taps[0] - 0.75, taps[-1] + 0.75], [y, y])
    for a, b in zip(taps, taps[1:]):
        diode(ax, (a + b) / 2, y)
    for t in taps:
        dot(ax, t, y)
    ax.text(taps[0] - 0.9, y, "+", ha="right", va="center",
            fontsize=17, fontweight="bold", color=BLACK)
    ax.text(taps[-1] + 0.9, y, "−", ha="left", va="center",
            fontsize=17, fontweight="bold", color=BLACK)


def finish(ax, label, name):
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.text(MODULE_W / 2 + 0.3, -0.85, name, ha="center", va="top",
            multialignment="center", fontsize=12.5, fontweight="bold",
            linespacing=1.3)
    ax.text(0.1, -0.75, label, ha="center", va="top", fontsize=17,
            fontweight="bold")


def diode_v(ax, x, y, s=0.42, cathode_up=True):
    """Vertical diode, cathode (bar) up or down."""
    sgn = 1 if cathode_up else -1
    ax.add_patch(Polygon([(x - s / 2.3, y - sgn * s / 2), (x + s / 2.3, y - sgn * s / 2),
                          (x, y + sgn * s / 2)],
                         closed=True, facecolor="white", edgecolor=BLACK,
                         lw=1.1, zorder=5))
    ax.plot([x - s / 2.3, x + s / 2.3], [y + sgn * s / 2] * 2, color=BLACK, lw=1.6, zorder=5)


# ---------------------------------------------------- full-cell module
def full_cell(ax, group_sizes, label, name):
    n, rows = sum(group_sizes), 10
    dx, dy = MODULE_W / n, MODULE_H / rows
    cells(ax, n, rows, dx, dy, 0.0)
    col_x = [dx * (i + 0.5) for i in range(n)]
    y_lo, y_hi = -EDGE, MODULE_H + EDGE
    y_d = MODULE_H + 0.95                       # diode row (junction box)

    for i, x in enumerate(col_x):
        wire(ax, [x, x], [y_lo, y_hi])          # series ribbons of the column
        # current flows from the last column (-) to the first one (+)
        arrow(ax, x, MODULE_H / 2, up=((n - 1 - i) % 2 == 1))
    for i in range(0, n - 1, 2):                # bottom U-turns (1-2, 3-4, ...)
        wire(ax, [col_x[i], col_x[i + 1]], [y_lo, y_lo])
    for i in range(1, n - 1, 2):                # top cross-connectors (2-3, ...)
        wire(ax, [col_x[i], col_x[i + 1]], [y_hi, y_hi])

    taps = tap_positions(col_x, group_sizes)
    for t in taps:                              # taps up to the junction box
        wire(ax, [t, t], [y_hi, y_d])
    for t in taps[1:-1]:
        dot(ax, t, y_hi)
    diode_row(ax, taps, y_d)
    finish(ax, label, name)


# ---------------------------------------------------- half-cell module
def half_cell(ax, group_sizes, label, name):
    n, rows_half = sum(group_sizes), 10
    gap = 1.4                                   # central band (junction boxes)
    dx = MODULE_W / n
    dy = (MODULE_H - gap) / (2 * rows_half)
    y_c = MODULE_H / 2
    top0, bot1 = y_c + gap / 2, y_c - gap / 2   # inner edges of the two halves
    cells(ax, n, rows_half, dx, dy, top0)
    cells(ax, n, rows_half, dx, dy, 0.0)
    col_x = [dx * (i + 0.5) for i in range(n)]

    y_hi, y_lo = MODULE_H + EDGE, -EDGE         # outer U-turns
    y_ti, y_bi = top0 - 0.2, bot1 + 0.2         # inner cross-connectors

    for i, x in enumerate(col_x):
        wire(ax, [x, x], [y_ti, y_hi])          # top half column
        wire(ax, [x, x], [y_lo, y_bi])          # bottom half column
        away = (n - 1 - i) % 2 == 0             # current leaves the centre line
        arrow(ax, x, (top0 + MODULE_H) / 2, up=away)
        arrow(ax, x, bot1 / 2, up=not away)
    for i in range(0, n - 1, 2):                # outer U-turns (mirrored)
        wire(ax, [col_x[i], col_x[i + 1]], [y_hi, y_hi])
        wire(ax, [col_x[i], col_x[i + 1]], [y_lo, y_lo])
    for i in range(1, n - 1, 2):                # inner cross-connectors
        wire(ax, [col_x[i], col_x[i + 1]], [y_ti, y_ti])
        wire(ax, [col_x[i], col_x[i + 1]], [y_bi, y_bi])

    taps = tap_positions(col_x, group_sizes)
    for t in taps:                              # shared taps: halves in parallel
        wire(ax, [t, t], [y_bi, y_ti])
    for t in taps[1:-1]:
        dot(ax, t, y_ti)
        dot(ax, t, y_bi)
    diode_row(ax, taps, y_c)
    finish(ax, label, name)


# ---------------------------------------------------- matrix shingle module
def matrix_shingle(ax, label, name, rows=13, per_row=9, width=18.0):
    """Installed (landscape) view, glass 1980 x 386 mm drawn to scale.
    Each block: `rows` rows in series, each row = `per_row` half-cells in
    parallel (shingled, brick pattern). Left block: + at the top; right block:
    + at the bottom; the two blocks are joined along the bottom edge."""
    height = width * 386 / 1980
    margin_x, margin_y, gap = width * 33 / 1980, height * 16.95 / 386, width * 20 / 1980
    blk_w = (width - 2 * margin_x - gap) / 2
    dy = (height - 2 * margin_y) / rows
    dx = blk_w / per_row
    pad = 0.06 * dy
    x_blocks = (margin_x, margin_x + blk_w + gap)

    ax.add_patch(plt.Rectangle((0, 0), width, height, fill=False, ec=BLACK, lw=1.0, zorder=1))
    for x0 in x_blocks:
        for r in range(rows):
            y = margin_y + r * dy
            shift = dx / 2 if r % 2 else 0.0               # brick pattern
            edges = sorted({0.0, blk_w, *[min(max(k * dx - shift, 0.0), blk_w) for k in range(per_row + 2)]})
            for a, b in zip(edges, edges[1:]):
                if b - a > 1e-6:
                    ax.add_patch(plt.Rectangle((x0 + a + pad, y + pad), b - a - 2 * pad, dy - 2 * pad,
                                               lw=0, facecolor=BLUE, zorder=2))

    y_top, y_bot = height - margin_y / 2, margin_y / 2      # busbars, just outside the cells
    left_x, right_x = x_blocks[0], x_blocks[1] + blk_w
    for x0 in x_blocks:
        wire(ax, [x0, x0 + blk_w], [y_top, y_top], lw=2.2)
    wire(ax, [left_x, right_x], [y_bot, y_bot], lw=2.2)     # bottom busbars + link between the blocks

    # j-boxes on the short edges: bypass diode across each block, + left, - right
    xd_l, xd_r = -0.55, width + 0.55                         # j-boxes outside the glass
    wire(ax, [xd_l, left_x], [y_top, y_top]); wire(ax, [xd_l, left_x], [y_bot, y_bot])
    wire(ax, [xd_l, xd_l], [y_bot, y_top])
    wire(ax, [right_x, xd_r], [y_top, y_top]); wire(ax, [right_x, xd_r], [y_bot, y_bot])
    wire(ax, [xd_r, xd_r], [y_bot, y_top])
    diode_v(ax, xd_l, height / 2, s=0.42, cathode_up=True)      # cathode towards +
    diode_v(ax, xd_r, height / 2, s=0.42, cathode_up=False)     # cathode towards the link node
    for x, y in ((xd_l, y_top), (xd_l, y_bot), (xd_r, y_top), (xd_r, y_bot)):
        dot(ax, x, y)
    wire(ax, [xd_l - 0.7, xd_l], [y_top, y_top]); wire(ax, [xd_r, xd_r + 0.7], [y_top, y_top])
    ax.text(xd_l - 0.85, y_top, "+", ha="right", va="center", fontsize=17, fontweight="bold", color=BLACK)
    ax.text(xd_r + 0.85, y_top, "−", ha="left", va="center", fontsize=17, fontweight="bold", color=BLACK)

    ax.set_xlim(-2.2, width + 2.2)
    ax.set_ylim(-1.5, height + 0.3)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.text(width / 2, -0.35, name, ha="center", va="top", multialignment="center",
            fontsize=12.5, fontweight="bold", linespacing=1.3)
    ax.text(-1.5, -0.3, label, ha="center", va="top", fontsize=17, fontweight="bold")


# ------------------------------------------------------------ assemble
def draw(path, include_matrix=False, *, text_scale=1.0, panel_labels=True, title_labels=False):
    """Save the measured-module layouts as PDF/PNG; optionally add matrix panel E."""
    if include_matrix:
        fig = plt.figure(figsize=(14, 11.2))
        gs = fig.add_gridspec(2, 3, height_ratios=[7.4, 3.4], hspace=0.02, wspace=0.05)
        axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
        matrix_shingle(fig.add_subplot(gs[1, :]), "E",
                       "Matrix shingle module: 2 blocks × 13 rows in series, "
                       "9 half-cells in parallel per row")
    else:
        fig, axes = plt.subplots(1, 3, figsize=(14, 7.6), gridspec_kw={"wspace": 0.05})
    full_cell(axes[0], [2, 2, 2], "A", "Full-cell module\n(standard)")
    full_cell(axes[1], [2, 4, 2], "B", "Full-cell module\n(unbalanced substrings)")
    half_cell(axes[2], [2, 2, 2], "C/D", "Half-cell module\n(butterfly design)")
    for label, ax in zip("ABC", axes):
        for text in ax.texts:
            if title_labels and text.get_text().startswith(("Full-cell module", "Half-cell module")):
                text.set_text(f"{label}  {text.get_text()}")
            if not panel_labels and text.get_text() in ("A", "B", "C/D"):
                text.set_visible(False)
            else:
                text.set_fontsize(text.get_fontsize() * text_scale)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"), format="pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".png"), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved {path.with_suffix('.pdf')}")

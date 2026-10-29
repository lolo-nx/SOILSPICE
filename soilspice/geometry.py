"""From a soiling band drawn on the glass to a covered fraction per cell.
Calculation of the covered fraction per cell from a given geometry of band (rectangular; triangular)

Coordinates are in mm, origin at the bottom-left corner of the module as it is
installed (portrait or landscape), y pointing up. The band always starts at
the bottom edge.

    band (height h, shape)  ->  f_cov[i, j] = covered area / cell area
                            ->  A_cov       = mean of f_cov over all cells

Inactive cells (matrix shingle contact rows: covered by the band but producing
no current) get f_cov = NaN: they are skipped by A_cov and by the netlist.
"""
from dataclasses import dataclass

import numpy as np

_trapz = getattr(np, "trapezoid", None) or np.trapz

SHAPES = ("rectangle", "triangle", "triangle_left")   # triangle: apex (highest point) on the right
ORIENTATIONS = ("portrait", "landscape")


@dataclass(frozen=True)
class Layout:
    """Cell rectangles of a module in a given orientation (row 0 = bottom)."""
    width: float
    height: float
    x0: np.ndarray
    x1: np.ndarray
    y0: np.ndarray
    y1: np.ndarray

    @property
    def n_rows(self):
        return len(self.y0)

    @property
    def n_cols(self):
        return len(self.x0)


def layout(module, orientation):
    """Positions of all cells. Landscape is the portrait module turned so that
    its bottom edge ends up on the left: portrait column j becomes landscape
    row j, portrait row i becomes landscape column i."""
    left, right, bottom, top = module.offsets
    if orientation == "portrait":
        W, H = module.width, module.height
        n_cols, n_rows = module.n_cols, module.n_rows
        cw, ch = module.cell_width, module.cell_height
        gap_col, gap_row = 0.0, module.center_gap
        off_left, off_bottom = left, bottom
    elif orientation == "landscape":
        W, H = module.height, module.width
        n_cols, n_rows = module.n_rows, module.n_cols
        cw, ch = module.cell_height, module.cell_width
        gap_col, gap_row = module.center_gap, 0.0
        off_left, off_bottom = bottom, left
    else:
        raise ValueError(f"orientation must be one of {ORIENTATIONS}, got {orientation!r}")

    g = module.gap
    x0 = off_left + np.arange(n_cols) * (cw + g)
    y0 = off_bottom + np.arange(n_rows) * (ch + g)
    # the central gap (between the two half-modules) shifts the second half
    x0 = x0 + np.where(np.arange(n_cols) >= n_cols // 2, gap_col, 0.0)
    y0 = y0 + np.where(np.arange(n_rows) >= n_rows // 2, gap_row, 0.0)
    return Layout(W, H, x0, x0 + cw, y0, y0 + ch)


def check_dimensions(module, tol_mm=5.0):
    """Declared module size vs cells + gaps + margins + frame (portrait)."""
    left, right, bottom, top = module.offsets
    cells_x = module.n_cols * module.cell_width + (module.n_cols - 1) * module.gap
    cells_y = module.n_rows * module.cell_height + (module.n_rows - 1) * module.gap + module.center_gap
    dW = module.width - (left + cells_x + right)
    dH = module.height - (bottom + cells_y + top)
    return {"module": module.name, "dW_mm": dW, "dH_mm": dH,
            "ok": abs(dW) <= tol_mm and abs(dH) <= tol_mm}


def band_height(x, shape, height_mm, width, module_height, apex="right"):
    """Upper edge of the band at abscissa x.
    rectangle: constant height. triangle: 0 on one side, height_mm on the apex side."""
    h = float(np.clip(height_mm, 0.0, module_height))
    if shape == "triangle_left":               # same triangle, mirrored: highest on the left
        shape, apex = "triangle", ("left" if apex == "right" else "right")
    if shape == "rectangle":
        return np.full_like(x, h, dtype=float)
    if shape == "triangle":
        frac = x / width if apex == "right" else 1.0 - x / width
        return h * np.clip(frac, 0.0, 1.0)
    raise ValueError(f"shape must be one of {SHAPES}, got {shape!r}")


def fcov_grid(module, shape, orientation, height_mm, apex="right", n_sub=400):
    """Covered fraction of every cell, shape (n_rows, n_cols), row 0 = bottom.

    Each cell is sampled on n_sub vertical slices; in each slice the covered
    height is clip(h(x), y0, y1) - y0, integrated with the trapezoidal rule."""
    lay = layout(module, orientation)
    xs = np.linspace(lay.x0, lay.x1, n_sub, axis=1)                 # (n_cols, n_sub)
    h = band_height(xs, shape, height_mm, lay.width, lay.height, apex)
    y0, y1 = lay.y0[:, None, None], lay.y1[:, None, None]
    covered = np.clip(h[None], y0, y1) - y0                          # (n_rows, n_cols, n_sub)
    area = (lay.x1 - lay.x0)[None, :] * (lay.y1 - lay.y0)[:, None]
    fcov = _trapz(covered, xs[None], axis=2) / area
    return fcov * inactive_mask(module, orientation, nan=True)


def inactive_mask(module, orientation, nan=False):
    """Grid (installed orientation) of 1 for active cells and 0 (or NaN) for the
    inactive ones."""
    lay = layout(module, orientation)
    mask = np.ones((lay.n_rows, lay.n_cols))
    for i in range(lay.n_rows):
        for j in range(lay.n_cols):
            row, col = (i, j) if orientation == "portrait" else (j, i)
            if module.is_inactive(row, col):
                mask[i, j] = np.nan if nan else 0.0
    return mask


def acov(fcov):
    """Covered fraction of the module, as the mean over the active cells (not area-weighted)."""
    return float(np.nanmean(fcov))


def height_for_acov(module, shape, orientation, target, apex="right", tol=1e-4, max_iter=60):
    """Band height [mm] giving A_cov = target, found by bisection.
    A_cov grows monotonically with the height for both shapes.
    Returns (height_mm, achieved A_cov)."""
    top = layout(module, orientation).height

    def a(h):
        return acov(fcov_grid(module, shape, orientation, h, apex))

    a_max = a(top)
    if target >= a_max:
        print(f"  [!] A_cov={target:.1%} out of reach for {module.name} ({shape}, {orientation}),"
              f" max {a_max:.1%} -> band set to the full height")
        return top, a_max

    lo, hi = 0.0, top
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        a_mid = a(mid)
        if abs(a_mid - target) <= tol:
            return mid, a_mid
        lo, hi = (mid, hi) if a_mid < target else (lo, mid)
    return mid, a(mid)


def cell_fcov(module, fcov, orientation):
    """{SPICE instance: f_cov} from a grid in the installed orientation."""
    n_rows, n_cols = fcov.shape
    values = {}
    for i in range(n_rows):
        for j in range(n_cols):
            row, col = (i, j) if orientation == "portrait" else (j, i)
            if module.is_inactive(row, col):
                continue
            name = module.instance(row, col)
            if name is not None:
                values[name] = float(fcov[i, j])
    n_active = fcov.size - len(module.inactive_cells)
    if len(values) < n_active:
        raise RuntimeError(f"{module.name}: only {len(values)} of {n_active} active cells are mapped "
                           "to a SPICE instance, check the cell map")
    return values


# --- severity indicators used to put different modules on the same x axis ---------

def bypass_group_coverage(module, fcov, orientation):
    """Mean f_cov of each bypass sub-string (all its cells, shaded or not).
    bypass_groups are ranges of portrait columns, or of portrait rows when
    module.bypass_axis == "rows" (the matrix shingle: one diode per half)."""
    out = []
    for first, last in module.bypass_groups:
        idx = list(range(first, last + 1))
        along_cols = (module.bypass_axis == "cols") == (orientation == "portrait")
        sub = fcov[:, idx] if along_cols else fcov[idx, :]
        out.append(float(np.nanmean(sub)))
    return out


def band_severity(module, fcov, shape, orientation):
    """Single covered fraction that drives the electrical response.

    rectangle: every cell of the bottom row has the same f_cov -> f_cov of that row.
    triangle : f_cov varies along the row -> coverage of the most shaded bypass
               sub-string.
    Returns (severity, index of that sub-string or None)."""
    if shape == "rectangle":
        return float(np.nanmean(fcov[0, :])), None
    per_group = bypass_group_coverage(module, fcov, orientation)
    k = int(np.argmax(per_group))
    return per_group[k], k


def tau_eff(severity, tau_i):
    """Effective transmittance seen by the limiting cell or sub-string."""
    return 1.0 - severity * (1.0 - tau_i)

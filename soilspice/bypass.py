"""Simplified electrical model of a band: ideal cells, ideal bypass diodes, no reverse-bias physics.

Each cell gives a photocurrent g = 1 - f_cov (1 - tau_I) (clean = 1). A bypass sub-string carries at most the
current of its most shaded cell; in half-cell (butterfly) modules a sub-string is two half-strings in parallel,
each carrying half the current. At a string current I every sub-string able to carry I contributes its share of
the clean voltage (cells in it / all cells), the others are bypassed (no voltage, diode drop neglected).
P(I) = I x sum of those voltages; the loss is 1 - max_I P(I) / P_clean.
Used to compare three model levels: uniform (loss = S), this ideal-bypass model, and SPICE.
"""
import numpy as np

from . import geometry, modules

HALF_CELL = {"JASOLAR", "LONGI"}


def portrait_grid(mod, fcov, orientation):
    """f_cov in portrait indices (row 0 = bottom)."""
    return fcov if orientation == "portrait" else fcov.T


def bypass_loss(module_name, shape, orientation, h, tau, groups=None, half_cell=None):
    """groups (column ranges per bypass diode) and half_cell override the module's own topology."""
    mod = modules.get(module_name)
    groups = mod.bypass_groups if groups is None else groups
    half_cell = module_name in HALF_CELL if half_cell is None else half_cell
    f = portrait_grid(mod, geometry.fcov_grid(mod, shape, orientation, h), orientation)
    g = 1.0 - np.nan_to_num(f) * (1.0 - tau)
    n_cells = np.isfinite(f).sum()
    subs = []                                                     # (current capacity, voltage share)
    for first, last in groups:
        block = g[:, first:last + 1]
        cells = np.isfinite(f[:, first:last + 1]).sum()
        if half_cell:
            half = mod.n_rows // 2
            cap = 0.5 * block[:half].min() + 0.5 * block[half:].min()
            subs.append((cap, cells / 2 / n_cells * 2))            # two parallel halves share the voltage
        else:
            subs.append((block.min(), cells / n_cells))
    caps = np.array([c for c, _ in subs]); v = np.array([s for _, s in subs]); v = v / v.sum()
    p = max(I * v[caps >= I - 1e-12].sum() for I in caps)
    return 1.0 - p

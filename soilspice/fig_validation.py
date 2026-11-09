"""Figures of chapter 1 (validation): measured and simulated I-V curves of the four modules, whole-curve errors,
and the measured loss against uniform soiling and SPICE (model hierarchy)."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle
from matplotlib.ticker import MaxNLocator, MultipleLocator
from matplotlib.transforms import Bbox

from . import experiments, geometry, modules, plots, schematics, validation
from .style import LETTER_KW, MODULE_COLOR, ORDER, ORIENTATIONS, STYLE, letter_inside, module_label, save

SHAPES = ("rectangle", "triangle")
METRICS = (("nRMSE_PV_%", r"nRMSE, $P$–$V$ [%]"), ("nMBE_IV_%", r"nMBE, $I$–$V$ [%]"))
# Warm yellow-to-brown progression with wider lightness spacing than the original.
# Keep this palette local to the validation grid and preserve the shared markers.
IV_TAU_COLORS = {1.00: "#E1B83D", 0.71: "#CF8425", 0.54: "#AA501A",
                 0.44: "#6F300F", 0.00: "#26130B"}
# Architecture references in pv_module_schematics: standard full-cell,
# unequal-substring full-cell, and butterfly half-cell, respectively.
ELECTRICAL_LAYOUT = {"YINGLI": "A", "SUNPOWER": "B", "JASOLAR": "C", "LONGI": "C"}
CELL_FILL, CELL_EDGE = "#F5F5F5", "#808080"


def iv_tau_style(tau):
    color, marker = plots.tau_style(tau)
    anchor = min(IV_TAU_COLORS, key=lambda value: abs(value - tau))
    if abs(anchor - tau) <= plots.TAU_TOL:
        color = IV_TAU_COLORS[anchor]
    return color, marker


def collect(shape):
    """Runs the validation configurations and returns every curve as one long table."""
    rows = []
    for name in ORDER:
        module = modules.get(name)
        for orientation in ORIENTATIONS:
            r = validation.run_configuration(module, shape, orientation)
            for m in r.measured:
                rows.append(pd.DataFrame({"module": name, "orientation": orientation, "tau_I": m.tau_i,
                                          "source": "EXP", "V": m.v, "I": m.i}))
                sim = r.simulated.get(m.tau_i)
                if sim is not None:
                    rows.append(pd.DataFrame({"module": name, "orientation": orientation, "tau_I": m.tau_i,
                                              "source": "SPICE", "V": sim.v, "I": sim.i}))
    return pd.concat(rows, ignore_index=True)


def mini_schematic(ax, module, shape, orientation, fontsize=None):
    """Small schematic of the module in a corner of a panel, in the style of the slides:
    cells (drawn slightly smaller so that the gaps show), series ribbons, U-turns between
    the two columns of each pair, bypass diodes, the band and the f_cov of the covered cells.

    Everything is drawn in portrait coordinates and swapped for landscape, which is the
    convention of geometry.layout (portrait column j -> landscape row j, row i -> column i)."""
    lay = geometry.layout(module, "portrait")
    fcov = validation.band_fcov(module, shape, orientation)
    land = orientation == "landscape"
    T = (lambda x, y: (y, x)) if land else (lambda x, y: (x, y))
    cw, ch = lay.x1[0] - lay.x0[0], lay.y1[0] - lay.y0[0]
    cx = (lay.x0 + lay.x1) / 2
    k = 0.09 * min(cw, ch)                               # visual shrink of each cell
    half = module.junction_box != "top"                  # butterfly (half-cell) layout
    nr, nc = lay.n_rows, lay.n_cols
    wire_color = schematics.RED
    wire = dict(color=wire_color, lw=0.65, solid_capstyle="butt", zorder=4)

    def line(xs, ys, **kw):
        pts = [T(x, y) for x, y in zip(xs, ys)]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], **{**wire, **kw})

    def diode(x, y, size):
        s = size / 2
        # Full-cell cathodes face left, as in the architecture overview.
        # Apply T afterwards so the symbol rotates with its module in landscape.
        direction = 1 if land else -1
        tri = [T(x - direction * s, y - s), T(x - direction * s, y + s), T(x + direction * s, y)]
        ax.add_patch(Polygon(tri, facecolor="white", edgecolor=wire_color, lw=0.65, zorder=5))
        line([x + direction * s] * 2, [y - s, y + s], lw=0.9, zorder=5)

    # cells and the f_cov of the covered ones
    for r in range(nr):
        for c in range(nc):
            x0, y0 = T(lay.x0[c] + k, lay.y0[r] + k)
            x1, y1 = T(lay.x1[c] - k, lay.y1[r] - k)
            ax.add_patch(Rectangle((min(x0, x1), min(y0, y1)), abs(x1 - x0), abs(y1 - y0),
                                   facecolor=CELL_FILL, edgecolor=CELL_EDGE, lw=0.30, zorder=3))
            f = fcov[c, r] if land else fcov[r, c]
            if f > 5e-3 and fontsize:
                ax.text((x0 + x1) / 2, (y0 + y1) / 2, f"{f:.1f}" if f >= 0.095 else f"{f:.2f}",
                        ha="center", va="center", fontsize=fontsize, zorder=6,
                        rotation=90 if land else 0)

    # series ribbons between consecutive cells of a column (not across the centre gap)
    for c in range(nc):
        for r in range(nr - 1):
            if half and r == nr // 2 - 1:
                continue
            line([cx[c], cx[c]], [lay.y1[r] - k, lay.y0[r + 1] + k])

    # U-turns joining the two columns of each pair, at the outer edge(s)
    pairs = [(c, c + 1) for c in range(0, nc - 1, 2)]
    u = 0.45 * min(cw, ch)
    edges = [(lay.y0[0] + k, -1)] + ([(lay.y1[-1] - k, +1)] if half else [])
    for y_edge, sgn in edges:
        for a, b in pairs:
            line([cx[a], cx[a], cx[b], cx[b]], [y_edge, y_edge + sgn * u, y_edge + sgn * u, y_edge])

    # bypass diodes: above the top row (full-cell, j-box on top) or on the centre line (half-cell)
    if half:
        yb = (lay.y1[nr // 2 - 1] + lay.y0[nr // 2]) / 2
        for c in range(nc):
            line([cx[c], cx[c]], [lay.y1[nr // 2 - 1] - k, lay.y0[nr // 2] + k])
    else:
        yb = lay.y1[-1] + 0.9 * min(cw, ch)
        for c in range(nc):
            line([cx[c], cx[c]], [lay.y1[-1] - k, yb])
    line([cx[0], cx[-1]], [yb, yb])
    for first, last in module.bypass_groups:
        diode((cx[first] + cx[last]) / 2, yb, 0.70 * min(cw, ch))   # centred on its group

    # band: grey, translucent, from the lower edge of the module (display coordinates)
    lay_d = geometry.layout(module, orientation)
    W, H = lay_d.width, lay_d.height
    h = min(experiments.BAND_HEIGHT_MM, H)
    band = ([(0, 0), (W, 0), (W, h), (0, h)] if shape == "rectangle"
            else [(0, 0), (W, 0), (W, h)] if experiments.TRIANGLE_APEX == "right" else [(0, 0), (W, 0), (0, h)])
    ax.add_patch(Polygon(band, facecolor=to_rgba("0.35", 0.80), edgecolor="0.25", lw=0.5, zorder=7))

    pad = 1.0 * min(cw, ch)
    xs, ys = zip(*[T(x, y) for x, y in ((0, -pad), (lay.width, lay.height + pad))])
    ax.set_xlim(min(xs + (0,)), max(xs + (W,)))
    ax.set_ylim(min(ys + (0,)), max(ys + (H,)))
    ax.set_aspect("equal")
    ax.set_anchor("W")
    ax.axis("off")


def tau_groups(curves):
    """One legend entry per measured transmittance, common to all modules."""
    # Keep the exact transmittances as keys; round only the displayed labels.
    return {
        float(tau): rf"$\tau_I$ = {tau:.2f}"
        for tau in sorted(curves["tau_I"].unique(), reverse=True)
    }


def iv_grid(curves, path, shape="rectangle", *, figsize=(7.2, 3.9), compact_layout=False):
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, len(ORDER), figsize=figsize, sharex="col", sharey="col")
        for col, name in enumerate(ORDER):
            for row, orientation in enumerate(ORIENTATIONS):
                ax = axes[row, col]
                sub = curves[(curves["module"] == name) & (curves["orientation"] == orientation)]
                for tau, g in sub.groupby("tau_I", sort=False):
                    color, marker = iv_tau_style(tau)
                    exp, sim = g[g["source"] == "EXP"], g[g["source"] == "SPICE"]
                    ax.plot(exp["V"], exp["I"], ls="none", marker=marker, ms=2.6, mew=0,
                            markevery=0.04, color=color)
                    ax.plot(sim["V"], sim["I"], ls="--", lw=1.1, color=color)
                ax.set_ylim(bottom=0)
                ax.set_xlim(left=0)
                # integer ticks and small padding: narrower tick labels, so wider panels
                ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
                ax.xaxis.set_major_locator(MultipleLocator(20))
                ax.tick_params(pad=1.5, labelsize=7)
                letter = "abcdefgh"[row * len(ORDER) + col]
                ax.text(0.04, 0.05, letter, transform=ax.transAxes, **LETTER_KW)
                # module drawing in the empty lower-left corner (below the plateau, left of the knee)
                if orientation == "portrait":
                    # Use the bottom padding to lower the drawing without clipping
                    # the band or the outer U-turns against the current axis.
                    bottom = -0.025 if name in ("YINGLI", "SUNPOWER") else -0.005
                    height = 0.61 if name == "YINGLI" else 0.64
                    box = (0.08, bottom, 0.43, height)
                else:
                    box = (0.08, 0.015, 0.64, 0.50)
                inset = ax.inset_axes(box)
                mini_schematic(inset, modules.get(name), shape, orientation)
                if orientation == "landscape":
                    # Swapping the physical dimensions of the portrait frame keeps
                    # cell sizes equal after rotation, at any figure/export size.
                    portrait = axes[0, col].child_axes[0]

                    def same_scale_position(_, renderer, reference=portrait, parent=ax):
                        ref = reference.get_window_extent(renderer)
                        plot = parent.get_window_extent(renderer)
                        frame = Bbox.from_bounds(plot.x0 + 0.08 * plot.width,
                                                 plot.y0 + 0.015 * plot.height,
                                                 ref.height, ref.width)
                        return frame.transformed(parent.figure.transFigure.inverted())

                    inset.set_axes_locator(same_scale_position)
                inset.text(1.06, 0.5, ELECTRICAL_LAYOUT[name], transform=inset.transAxes,
                           ha="left", va="center", fontsize=7, fontweight="bold",
                           color="0.25", clip_on=False)
                if row == 0:
                    ax.set_title(plots.module_label(name))
                else:
                    ax.set_xlabel("Voltage [V]")
                if col == 0:
                    ax.set_ylabel(r"$\bf{" + orientation.capitalize() + r"}$" + "\nCurrent [A]", labelpad=1)

        labels = tau_groups(curves)
        handles = [
            Line2D([], [], color=iv_tau_style(tau)[0],
                   marker=iv_tau_style(tau)[1], ms=4,
                   ls="none", label=label)
            for tau, label in labels.items()
        ]
        handles += [Line2D([], [], color="0.3", marker="o", ms=3, ls="none", label="measured"),
                    Line2D([], [], color="0.3", ls="--", lw=1.1, label="SPICE")]
        fig.legend(handles=handles, loc="lower center", ncol=len(handles), bbox_to_anchor=(0.5, -0.04),
                   frameon=False, handletextpad=0.3, columnspacing=1.0)
        fig.tight_layout(rect=(0, 0.05, 1, 1), h_pad=0.6, w_pad=0.0)
        if compact_layout:
            # Spend less vertical space on margins so the curves and their
            # embedded schematics stay readable in the fixed article panel.
            fig.subplots_adjust(top=0.92, bottom=0.09, hspace=0.06)
            # Keep the anchor in figure coordinates when tight PNG/PDF exports
            # temporarily change the canvas bounds.
            fig.legends[0].set_bbox_to_anchor((0.5, -0.085), transform=fig.transFigure)
        plots.save(fig, path)


def errors_grid(df, path, *, figsize=(3.4, 3.6)):
    groups = [(o, s) for o in ORIENTATIONS for s in SHAPES]
    x, w = np.arange(len(groups)), 0.8 / len(ORDER)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
        for ax, (col, ylabel), letter in zip(axes, METRICS, "ab"):
            for k, name in enumerate(ORDER):
                vals = [df.loc[(df.module == name) & (df.orientation == o) & (df["shape"] == s), col].mean()
                        for o, s in groups]
                ax.bar(x + (k - (len(ORDER) - 1) / 2) * w, vals, w, color=plots.MODULE_COLOR[name],
                       edgecolor="black", lw=0.4, zorder=3)
            ax.axhline(0, color="black", lw=0.6, zorder=4)
            ax.axvline(len(SHAPES) - 0.5, color="0.6", lw=0.6, ls="--")
            ax.set_xlim(-0.5, len(groups) - 0.5)
            ax.set_ylabel(ylabel, labelpad=2)
            ax.tick_params(pad=1.5)
            # Leave the tall first Yingli bar clear in the lower panel.
            label_x = 0.14 if letter == "b" else 0.015
            ax.text(label_x, 0.97, letter, transform=ax.transAxes, va="top", **LETTER_KW)
        axes[0].set_ylim(0, 3.2)
        ax = axes[-1]
        ax.set_xticks(x)
        ax.set_xticklabels([s.capitalize() for _, s in groups])
        ax.tick_params(axis="x", length=0)
        for i, o in enumerate(ORIENTATIONS):
            ax.text(i * len(SHAPES) + (len(SHAPES) - 1) / 2, -0.13, o.capitalize(), transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontweight="bold")
        handles = [Rectangle((0, 0), 1, 1, facecolor=plots.MODULE_COLOR[n], edgecolor="black", lw=0.4,
                             label=plots.module_label(n)) for n in ORDER]
        fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.55, 1.0), frameon=False,
                   handlelength=1.2, handletextpad=0.4, columnspacing=1.2)
        fig.tight_layout(rect=(0, 0, 1, 0.92), h_pad=0.4)
        plots.save(fig, path)




def model_hierarchy(me, path, *, figsize=(3.5, 3.0)):
    """Predicted against measured loss of P_mp for uniform soiling (L = S, open squares) and SPICE (filled circles),
    one colour per module; RMSE of each model in the legend. One column wide.
    me: one row per shaded flash (columns module, L_meas, S, L_spice)."""
    models = [("S", "Uniform soiling approximation", "s", "none"), ("L_spice", "SPICE", "o", None)]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot([0, 55], [0, 55], color="black", lw=0.6, ls="--", zorder=1)
        for name in ORDER:
            g = me[me.module == name]
            col = MODULE_COLOR[name]
            ax.scatter(g.L_meas, g.S, s=12, marker="s", facecolors="none", edgecolors=col, lw=0.7, zorder=2)
            ax.scatter(g.L_meas, g.L_spice, s=12, marker="o", color=col, lw=0, alpha=0.9, zorder=3)
        hd = [Line2D([], [], ls="none", marker="o", ms=4, color=MODULE_COLOR[n], label=module_label(n)) for n in ORDER]
        for c, lab, mk, fc in models:
            e = me[c] - me.L_meas
            hd.append(Line2D([], [], ls="none", marker=mk, ms=4, color="0.3", mfc=fc or "0.3",
                             label=f"{lab}: RMSE {np.sqrt((e ** 2).mean()):.1f} pp"))
        ax.set_xlim(0, 50); ax.set_ylim(0, 55)
        ax.set_xlabel("Measured loss of $P_{mp}$ (%)"); ax.set_ylabel("Predicted loss of $P_{mp}$ (%)")
        ax.legend(handles=hd, loc="upper left", handletextpad=0.2, fontsize=6)
        save(fig, path)
